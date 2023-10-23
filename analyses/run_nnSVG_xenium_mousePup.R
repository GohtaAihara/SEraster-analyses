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

dataset_name <- "xenium_mousePup"
method <- "nnSVG"

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

## Figure (deg analysis)
## load DEG analysis
df_deg <- read.csv("~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_mousePup/analysis/diffexp/gene_expression_graphclust/differential_expression.csv")
rownames(df_deg) <- df_deg$Feature.Name

## input organ-specific genes suggested by 10X Genomics (https://www.10xgenomics.com/resources/datasets/mouse-pup-preview-data-xenium-mouse-tissue-atlassing-panel-1-standard)
marker_genes <- list(
  Brain = c("Gfap"),
  Kidney = c("Lrp2"),
  Pancreas = c("Tm4sf4", "Dcdc2a"),
  Liver = c("Clec4f", "Folr2"),
  Lung = c("Mfap2", "Sod3"),
  Heart = c("Pln", "Trdnepi"),
  Skin = c("Dsc3", "Sbsn"),
  Colon = c("Cdhr5", "Rbp2"),
  Thymus = c("Arpp21", "Cd8a")
)

df_marker_genes <- marker_genes %>%
  enframe(name = "organ", value = "gene") %>%
  unnest(gene) %>%
  mutate(gene = factor(gene, levels = unlist(marker_genes))) %>%
  as.data.frame()
## some marker genes don't exist in the dataset
df_marker_genes <- df_marker_genes[df_marker_genes$gene %in% rownames(spe),]

df_plt <- do.call(rbind, lapply(levels(ct_labels), function(ct_label) {
  ## subset by cluster
  spe_sub <- spe[,spe$cluster == ct_label]
  df_deg_sub <- df_deg[,grepl(paste0("Cluster.", ct_label, "."), colnames(df_deg), fixed = TRUE)]
  
  ## compute metrics
  out <- do.call(rbind, lapply(seq(dim(df_marker_genes)[1]), function(idx) {
    gene <- df_marker_genes[idx,]$gene
    return(data.frame(
      cluster = factor(ct_label, levels = levels(ct_labels)),
      organ = df_marker_genes[idx,]$organ,
      gene = factor(gene, levels = levels(df_marker_genes$gene)),
      mean_counts = mean(assay(spe_sub, "counts")[gene,]),
      mean_lognorm = mean(assay(spe_sub, "lognorm")[gene,]),
      fraction = mean(assay(spe_sub, "counts")[gene,] > 0)*100,
      log2fc = df_deg_sub[gene,grepl("Log2.fold.change", colnames(df_deg_sub))],
      padj = df_deg_sub[gene,grepl("Adjusted.p.value", colnames(df_deg_sub))]
    ))
  }))
}))

ggplot(df_plt, aes(x = gene, y = cluster, col = padj, size = fraction)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

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

## plot top n s.s. SVGs for each resolution and cluster
alpha <- 0.05
num_rank <- 25
for (res in unique(df$resolution)) {
  df_sub <- df[df$resolution == res,]
  for (ct_label in unique(df_sub$cluster)) {
    ## subset cluster of interest
    spe_sub <- spe[,spe$cluster == ct_label]
    
    ## rasterization
    spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
    
    top_svgs <- df_sub %>%
      filter(cluster == ct_label) %>%
      arrange(rank) %>%
      filter(rank <= num_rank, padj <= alpha)
  
    df_plt <- data.frame(x = spatialCoords(spe_sub_rast)[,1], y = spatialCoords(spe_sub_rast)[,2], as.matrix(t(assay(spe_sub_rast, "pixelval")[top_svgs$gene,]))) %>%
      pivot_longer(!c(x, y), names_to = "gene", values_to = "exp")
    ggplot(df_plt, aes(x = x, y = y, fill = exp)) +
      facet_wrap(~ gene) +
      coord_fixed() +
      geom_tile() +
      scale_fill_viridis_c(name = "gene expression") +
      labs(title = paste0("Cluster ", ct_label, ", top ", num_rank, " SVGs (padj <= ", alpha, ") based on nnSVG rank"),
           x = "x (um)",
           y = "y (um)") +
      theme_bw()
    ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_nnsvg_ct_specific_resolution_", res, "_cluster_", ct_label, "_top_", num_rank, "_svgs.pdf")), width = 8, heigh = 8, dpi = 300)
  }
}

## Figure (global nnSVG vs. ct-specific nnSVG)
## switch organ
## cluster 39 (kidney)
ct_label <- 39
genes <- c("Calb1", "Clcnka", "Cryab", "Ptn", "Scin", "Tfcp2l1")
# Calb1 = kidney collecting duct epithelial cell
# Clcnka = kidney loop of Henle ascending limb epithelial cell
# Cryab = podocyte (sensu Diptera) (cells in bowman's capsule)
# Ptn = transition zone lymphatic endothelial cell
# Scin = kidney cortex artery cell
# Tfcp2l1 = kidney collecting duct principal cell
res <- 100

## cluster 23 (brain) ## not doing it because nnSVG failed at resolution = 100, and many SVGs are pancreas related
ct_label <- 23
# genes <- c("Hap1", "Lhx2", "Nap1l5", "Neurod1", "Nnat", "Pnmal2", "Scg2", "Scg5", "Vsnl1")
## weird SVGs (Hap1 = pancreatic D cell, Lhx2 = keratinocyte stem cell, Nap1l5 = pancreatic A cell, Neurod1 = pancreatic A cell, Pnmal2 = pancreatic D cell, Scg2 = type B pancreatic cell, Scg5 = pancreatic PP cell, Vsnl1 = pericyte)
res <- 200

## cluster 11 (liver)
ct_label <- 11
genes <- c("Hrg", "Tat")
# Hrg = hepatocytes
# Pck1 = kidney proximal convoluted tubule epithelial cell
# Rbp2 = enterocyte of epithelium of large intestine
# Serpina3n = fibroblast of cardiac tissue
# Tat = hepatocyte
genes <- c("Hpx", "Gc", "Tat", "Slc10a1", "Abcb11", "Hrg") ## hepatocyte markers
res <- 100

## cluster 26 (intestines)
ct_label <- 26
# genes <- c("Apoc2", "Krt19", "Rbp2", "X2610528A11Rik")
# Apoc2 = intermediate monocyte
# Krt19 = club cell (usually found in lung)
# Rbp2 = enterocyte of epithelium of large intestine
# X2610528A11Rik = enterocyte of epithelium of large intestine
genes <- c("Oit1", "Krt8", "Gpx2", "Cldn7", "Cdx1", "Cdhr5", "Hmgcs2", "B4galnt2", "Mgat4c")
res <- 100

## cluster 22 (skin)
ct_label <- 22
genes <- c("Krt16", "Cnfn", "S100a14", "Krt13", "Serpinb5", "Them5")
res <- 100

num_genes <- length(genes)

## load nnSVG results
df_nnsvg_global <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global.RDS")))
df_nnsvg_ct_specific <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_ct_specific_v1.RDS")))

## subset ct-specific nnSVG results
df_nnsvg_ct_specific <- df_nnsvg_ct_specific[df_nnsvg_ct_specific$resolution == res & df_nnsvg_ct_specific$cluster == ct_label,]

## subset cluster of interest
spe_sub <- spe[,spe$cluster == ct_label]
## rasterization
spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

## extract data
df_global <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], as.matrix(t(assay(spe_rast, "pixelval")[genes,])))
rank_global <- df_nnsvg_global[df_nnsvg_global$gene %in% genes,c("gene", "rank", "padj")]
rank_global$method <- "global"

df_ct_specific <- data.frame(x = spatialCoords(spe_sub_rast)[,1], y = spatialCoords(spe_sub_rast)[,2], as.matrix(t(assay(spe_sub_rast, "pixelval")[genes,])))
rank_ct_specific <- df_nnsvg_ct_specific[df_nnsvg_ct_specific$gene %in% genes,c("gene", "rank", "padj")]
rank_ct_specific$method <- "cluster specific"

## plot global SVGs
df_global <- df_global %>%
  pivot_longer(-c("x","y"), names_to = "gene", values_to = "exp") %>%
  mutate(gene = factor(gene, levels = genes[order(rank_ct_specific$rank)]))

ggplot(df_global, aes(x = x, y = y, fill =　exp)) +
  facet_wrap(~ gene) +
  coord_fixed() +
  geom_tile() +
  scale_fill_viridis_c(name = "Log-normalized\nexpression") +
  labs(title = paste0("Global SVG (Resolution = ", res, ")"),
       x = "x (um)",
       y = "y (um)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_nnsvg_global_resolution_", res, "_cluster_", ct_label, "_visualization.pdf")), width = 6, height = 6, dpi = 300)

## plot cluster-specific SVGs
df_ct_specific <- df_ct_specific %>%
  pivot_longer(-c("x","y"), names_to = "gene", values_to = "exp") %>%
  mutate(gene = factor(gene, levels = genes[order(rank_ct_specific$rank)]))

ggplot(df_ct_specific, aes(x = x, y = y, fill =　exp)) +
  facet_wrap(~ gene) +
  coord_fixed() +
  geom_tile() +
  scale_fill_viridis_c(name = "Log-normalized\nexpression") +
  labs(title = paste0("Cluster-specific SVG (Resolution = ", res, ")"),
       x = "x (um)",
       y = "y (um)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_nnsvg_ct_specific_resolution_", res, "_cluster_", ct_label, "_visualization.pdf")), width = 6, heigh = 3, dpi = 300)

## plot rank comparison
df_rank <- rbind(rank_global, rank_ct_specific)

ggplot(df_rank, aes(x = gene, y = rank, col = method, shape = method)) +
  geom_point(size = 3) +
  labs(title = paste0("Comparison of nnSVG rank (Resolution = ", res, ")"),
       x = "Gene",
       y = "Rank") +
  theme_bw()
ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_nnsvg_global_ct_specific_comparison_resolution_", res, "_cluster_", ct_label, ".pdf")), width = 6, heigh = 3, dpi = 300)

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
