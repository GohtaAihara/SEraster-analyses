## This file integrates rasterization with nnSVG to analyze 10X Genomics Xenium mouse whole pup datasets. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(R.utils)
library(ggplot2)
library(gridExtra)
library(here)
library(tidyr)
library(tibble)
library(dplyr)

par(mfrow=c(1,1))

dataset_name = "xenium_mousePup"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = "outputs/xenium_mousePup_preprocessed.RDS")

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  geom_point(size = 0.1) +
  theme_classic()

ct_labels <- colData(spe)$cluster

# Run methods -------------------------------------------------------------

res_list <- list(50, 100, 200, 400)

## Global nnSVG
nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
  print(paste0("Resolution: ", res))
  if (res == "singlecell") {
    num_points = dim(spe)[2]
    
    ## nnSVG
    spe <- nnSVG::nnSVG(
      spe,
      assay_name = "lognorm",
      BPPARAM = BiocParallel::MulticoreParam()
    )
    df <- tibble::rownames_to_column(as.data.frame(rowData(spe)), var = "gene")
  } else {
    ## rasterization
    spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
    num_points = dim(spe_rast)[2]
    
    ## nnSVG
    spe_rast <- nnSVG::nnSVG(
      spe_rast,
      assay_name = "pixelval",
      BPPARAM = BiocParallel::MulticoreParam()
    )
    df <- tibble::rownames_to_column(as.data.frame(rowData(spe_rast)), var = "gene")
  }
  return(data.frame(dataset = dataset_name, resolution = res, num_cells = dim(spe)[2], num_pixels = dim(spe_rast)[2], df))
}))
## save results
saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global.RDS")))

## cluster-specific nnSVG
nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
  print(paste0("Resolution = ", res))
  out <- do.call(rbind, lapply(levels(ct_labels), function(ct_label) {
    print(paste0("Cluster = ", ct_label))
    ## subset cluster of interest
    spe_sub <- spe[,spe$cluster == ct_label]
    print(paste0("# of cells in cluster ", ct_label, " : ", dim(spe_sub)[2]))
    
    ## rasterization
    spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
    
    ## nnSVG
    ## using try() to handle error
    spe_sub_rast_nnsvg <- try({
      withTimeout(
        {nnSVG::nnSVG(
          spe_sub_rast,
          assay_name = "pixelval",
          BPPARAM = BiocParallel::MulticoreParam()
        )},
        timeout = 1800,
        onTimeout = c("error")
      )
    })
    
    if (class(spe_sub_rast_nnsvg) == "try-error") {
      ## do not save anything for clusters that caused error in nnSVG
      return(NULL)
    } else {
      df <- tibble::rownames_to_column(as.data.frame(rowData(spe_sub_rast_nnsvg)), var = "gene")
      return(data.frame(
        cluster = ct_label, 
        num_cells = dim(spe_sub)[2], 
        num_pixels = dim(spe_sub_rast)[2], 
        df))
    }
  }))
  return(cbind(dataset = dataset_name, resolution = res, out))
}))
## save results
saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_ct_specific.RDS")))


# Plot --------------------------------------------------------------------

## Figure (visual inspection of top SVGs (based on rank from nnSVG))
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global.RDS")))

## count the number of s.s. SVGs
alpha <- 0.05
sum(df$padj <= 0.05) ## 379 at resolution = 100 (all genes are SVGs)

res <- 100
spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

num_rank <- 25
df_ranked <- df[order(df$rank, decreasing = FALSE),]
top_svgs <- df_ranked[1:num_rank,]$gene
df_plt <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], as.matrix(t(assay(spe_rast)[top_svgs,]))) %>%
  pivot_longer(!c(x, y), names_to = "gene", values_to = "exp")
ggplot(df_plt, aes(x = x, y = y, fill = exp)) +
  facet_wrap(~ gene) +
  geom_tile() +
  scale_fill_viridis_c(name = "gene expression") +
  labs(title = paste0("Top ", num_rank, " SVGs based on nnSVG rank")) +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_nnsvg_global_top_", num_rank, "_svgs.pdf")), width = 8, heigh = 8, dpi = 300)

bottom_svgs <- df_ranked[(dim(df_ranked)[1]-num_rank+1):dim(df_ranked)[1],]$gene
df_plt <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], as.matrix(t(assay(spe_rast)[bottom_svgs,]))) %>%
  pivot_longer(!c(x, y), names_to = "gene", values_to = "exp")
ggplot(df_plt, aes(x = x, y = y, fill = exp)) +
  facet_wrap(~ gene) +
  geom_tile() +
  scale_fill_viridis_c(name = "gene expression") +
  labs(title = paste0("Bottom ", num_rank, " SVGs based on nnSVG rank")) +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_nnsvg_global_bottom_", num_rank, "_svgs.pdf")), width = 8, heigh = 8, dpi = 300)

## concordance of cluster-specific DEGs and nnSVGs
## load DEG analysis
df_deg <- read.csv("~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_mousePup/analysis/diffexp/gene_expression_graphclust/differential_expression.csv")

## Figure (number of s.s. SVG)
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_ct_specific_v1.RDS")))
alpha <- 0.05
df <- df %>%
  mutate(
    resolution = factor(resolution, levels = c(50, 100, 200, 400)),
    cluster = factor(cluster, levels = levels(ct_labels))
  )
df$svg_boolean <- df$padj <= alpha

# check which clusters were analyzed in each resolution
for (res in unique(df$resolution)) {
  ## print out analyzed clusters
  print(paste0("Resolution: ", res))
  print(paste0(unique(df[df$resolution == res,]$cluster)))
}

## number of s.s. SVGs identified in each cluster for each resolution
df_num_svg <- df %>%
  group_by(resolution, cluster) %>%
  summarise(num_svg = sum(svg_boolean), num_pixels = mean(num_pixels))
ggplot(df_num_svg, aes(x = cluster, y = resolution, size = num_svg, col = log10(num_pixels))) +
  geom_point() +
  scale_color_viridis_c() +
  labs(title = "Number of SVGs",
       x = "Cluster",
       y = "Resolution",
       size = "Number of SVGs",
       col = "log10(Number of pixels)") +
  theme_bw()
## save plot
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_num_svg.pdf")), width = 10, height = 4, dpi = 300)

## Figure (visual inspection of top 25 s.s. SVGs for each cluster)
## load relevant analysis results
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_ct_specific_v1.RDS")))
alpha <- 0.05
df <- df %>%
  mutate(
    resolution = factor(resolution, levels = c(50, 100, 200, 400)),
    cluster = factor(cluster, levels = levels(ct_labels))
  )
df$svg_boolean <- df$padj <= alpha

res <- 200
df_sub <- df[df$resolution == res,]

## load relevant rasterized data

for (ct_label in unique(df_sub$cluster)) {
  temp <- df_sub[df_sub$cluster == ct_label,]
  
  num_rank <- 25
  df_ranked <- temp[order(temp$rank, decreasing = FALSE),]
  top_svgs <- df_ranked[1:num_rank,]$gene
  df_plt <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], as.matrix(t(assay(spe_rast)[top_svgs,]))) %>%
    pivot_longer(!c(x, y), names_to = "gene", values_to = "exp")
  ggplot(df_plt, aes(x = x, y = y, fill = exp)) +
    facet_wrap(~ gene) +
    geom_tile() +
    scale_fill_viridis_c(name = "gene expression") +
    labs(title = paste0("Top ", num_rank, " SVGs based on nnSVG rank")) +
    theme_bw()
  ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_nnsvg_global_top_", num_rank, "_svgs.pdf")), width = 8, heigh = 8, dpi = 300)
  
}

# Further exploration -----------------------------------------------------

## plot the number of cells in each cluster
barplot(table(colData(spe)$cluster),
        xlab = "Cluster",
        ylab = "Number of cells")

## plot the proportion of each cluster
barplot(prop.table(table(colData(spe)$cluster)),
        xlab = "Cluster",
        ylab = "% of cells")

## run nnSVG on specific cluster
ct_label <- 1
spe_sub <- spe[,spe$cluster == ct_label]

df <- data.frame(spatialCoords(spe_sub), colData(spe_sub))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  geom_point(size = 0.1) +
  theme_classic()

## rasterization
res <- 100
spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

df <- data.frame(spatialCoords(spe_sub_rast), colData(spe_sub_rast), transcripts = colMeans(assay(spe_sub_rast)))
ggplot(df, aes(x = x, y = y, fill = transcripts)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_classic()

## nnSVG
## using try() to handle error and withTimeout() to handle when nnSVG takes forever to run
nnsvg_output <- try({
  withTimeout(
    {nnSVG::nnSVG(
      spe_sub_rast,
      assay_name = "pixelval",
      BPPARAM = BiocParallel::MulticoreParam()
    )},
    timeout = 1800,
    onTimeout = c("error")
  )
})
class(nnsvg_output)
