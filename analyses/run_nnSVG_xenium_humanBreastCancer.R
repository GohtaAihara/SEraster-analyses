## This file integrates rasterization with nnSVG to analyze 10X Genomics Xenium human breast cancer datasets. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(rhdf5)
library(here)
library(tibble)
library(dplyr)
library(reshape2)
library(pheatmap)

par(mfrow=c(1,1))

dataset_name = "xenium_humanBreastCancer"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = "outputs/xenium_humanBreastCancer_preprocessed.RDS")

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = celltype)) +
  geom_point(size = 0.5, stroke = 0) +
  theme_classic()

ct_labels <- colData(spe)$celltype

# Run methods -------------------------------------------------------------

res_list <- c(50, 100, 200, 400)

# ct_labels <- factor(c(1,16), levels = c(1, 16))

nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
  out <- do.call(rbind, lapply(levels(ct_labels), function(ct_label) {
    # ## subset cluster of interest
    # spe_sub <- spe[,spe$cluster == ct_label]
    # print(paste0("# of cells in cluster ", ct_label, " : ", dim(spe_sub)[2]))
    
    ## subset celltype of interest
    spe_sub <- spe[,spe$celltype == ct_label]
    print(paste0("# of cells in celltype ", ct_label, " : ", dim(spe_sub)[2]))
    
    ## rasterization
    spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
    
    ## nnSVG
    ## using try() to handle error
    spe_sub_rast_nnsvg <- try({nnSVG::nnSVG(
      spe_sub_rast,
      assay_name = "pixelval",
      BPPARAM = BiocParallel::MulticoreParam()
    )})
    
    if (class(spe_sub_rast_nnsvg) == "try-error") {
      ## do not save anything for celltypes that caused error in nnSVG
      return(NULL)
    } else {
      df <- tibble::rownames_to_column(as.data.frame(rowData(spe_sub_rast_nnsvg)), var = "gene")
      return(data.frame(
        celltype = ct_label, 
        num_cells = dim(spe_sub)[2], 
        num_pixels = dim(spe_sub_rast)[2], 
        df))
    }
  }))
  return(cbind(dataset = dataset_name, resolution = res, out))
}))

## save results
saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_ct_specific_supervised_celltypes.RDS")))


# Plot --------------------------------------------------------------------

## Figure (single cell visualizations)
df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = celltype)) +
  coord_fixed() +
  geom_point(size = 0.1) +
  labs(title = "Single cell",
       col = "Cell-type") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
## save plot
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_celltypes_singlecell.pdf")))

## Figure (cell-type specific SVGs)
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_ct_specific_v1.RDS")))
alpha <- 0.05
df <- df %>%
  mutate(
    resolution = factor(resolution, levels = c(50, 100, 200, 400)),
    celltype = factor(celltype, levels = levels(ct_labels))
    )
df$svg_boolean <- df$padj <= alpha

# check which celltypes were analyzed in each resolution
for (res in unique(df$resolution)) {
  ## print out analyzed celltypes
  print(paste0("Resolution: ", res))
  print(paste0(unique(df[df$resolution == res,]$celltype)))
}

## number of s.s. SVGs identified in each celltype for each resolution
df_num_svg <- df %>%
  group_by(resolution, celltype) %>%
  summarise(num_svg = sum(svg_boolean), num_pixels = mean(num_pixels))
ggplot(df_num_svg, aes(x = celltype, y = resolution, size = num_svg, col = log10(num_pixels))) +
  geom_point() +
  scale_color_viridis_c() +
  labs(title = "Number of SVGs",
       x = "Cell-type",
       y = "Resolution",
       size = "Number of SVGs",
       col = "log10(Number of pixels)") +
  theme_bw()
## save plot
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_num_svg.pdf")), width = 6, height = 4, dpi = 300)

## heatmap + hiearchical clustering
## test one resolution
res <- 50
df_sub <- df[df$resolution == res,]
df_sub_heatmap <- dcast(df_sub, gene ~ cluster, value.var = "padj")
df_plt <- df_sub_heatmap[,-1]
rownames(df_plt) <- df_sub_heatmap$gene
pheatmap(
  df_plt,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  fontsize_row = 6,
  angle_col = 0,
  filename = here("plots", dataset_name, paste0(dataset_name, "_resolution_", res, "_heatmap.pdf")),
  width = 6,
  height = 20
  )

## see unique SVGs for each cluster
res <- 50
df_sub <- df[df$resolution == res,]
unique_genes <- df_sub %>%
  group_by(gene) %>%
  filter(sum(svg_boolean) == 1, svg_boolean == TRUE) %>%
  ungroup()

gene_num_cluster <- df_sub %>%
  group_by(gene) %>%
  summarise(sum(svg_boolean))

## visualize unique genes
for (i in 1:dim(unique_genes)[1]) {
  print(paste0("Visualizing cluster = ", unique_genes[i,]$cluster, " , gene = ", unique_genes[i,]$gene))
  spe_sub <- spe[,spe$cluster == unique_genes[i,]$cluster]
  spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
  gene <- unique_genes[i,]$gene
  df <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], gene = assay(spe_sub, "counts")[gene,])
  p1 <- ggplot(df, aes(x = x, y = y, col = gene)) +
    geom_point(size = 0.1) +
    scale_color_distiller(name = gene, palette = "OrRd", direction = 1) +
    labs(title = "Single cell (raw counts)") +
    theme_classic()
  df <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], gene = assay(spe_sub, "lognorm")[gene,])
  p2 <- ggplot(df, aes(x = x, y = y, col = gene)) +
    geom_point(size = 0.1) +
    scale_color_distiller(name = gene, palette = "OrRd", direction = 1) +
    labs(title = "Single cell (pseudolog normalized)") +
    theme_classic()
  df <- data.frame(x = spatialCoords(spe_sub_rast)[,1], y = spatialCoords(spe_sub_rast)[,2], gene = assay(spe_sub_rast, "pixelval")[gene,])
  p3 <- ggplot(df, aes(x = x, y = y, fill = gene)) +
    geom_tile() +
    scale_fill_distiller(name = gene, palette = "OrRd", direction = 1) +
    labs(title = paste0("Resolution = ", res)) +
    theme_classic()
  
  plt <- grid.arrange(p1, p2, p3, ncol = 3, top = paste0("Cluster = ", ct_label, ", Gene = ", gene))
  ## save plot
  ggsave(plt, filename = here("plots", dataset_name, paste0(dataset_name, "_cluster_", ct_label, "_gene_", gene, ".pdf")), width = 18, height = 5, dpi = 300)
}

ct_label <- 2
spe_sub <- spe[,spe$cluster == ct_label]
spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
gene <- "CLIC6"
df <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], gene = assay(spe_sub, "counts")[gene,])
p1 <- ggplot(df, aes(x = x, y = y, col = gene)) +
  geom_point(size = 0.1) +
  scale_color_distiller(name = gene, palette = "OrRd", direction = 1) +
  labs(title = "Single cell (raw counts)") +
  theme_classic()
df <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], gene = assay(spe_sub, "lognorm")[gene,])
p2 <- ggplot(df, aes(x = x, y = y, col = gene)) +
  geom_point(size = 0.1) +
  scale_color_distiller(name = gene, palette = "OrRd", direction = 1) +
  labs(title = "Single cell (pseudolog normalized)") +
  theme_classic()
df <- data.frame(x = spatialCoords(spe_sub_rast)[,1], y = spatialCoords(spe_sub_rast)[,2], gene = assay(spe_sub_rast, "pixelval")[gene,])
p3 <- ggplot(df, aes(x = x, y = y, fill = gene)) +
  geom_tile() +
  scale_fill_distiller(name = gene, palette = "OrRd", direction = 1) +
  labs(title = paste0("Resolution = ", res)) +
  theme_classic()

plt <- grid.arrange(p1, p2, p3, ncol = 3, top = paste0("Cluster = ", ct_label, ", Gene = ", gene))
## save plot
ggsave(plt, filename = here("plots", dataset_name, paste0(dataset_name, "_cluster_", ct_label, "_gene_", gene, ".pdf")), width = 18, height = 5, dpi = 300)

# Further exploration -----------------------------------------------------

## filter by cluster --> rasterize --> nnSVG
ct_label <- 1

## subset cluster of interest
spe_sub <- spe[,spe$cluster == ct_label]

df <- data.frame(spatialCoords(spe_sub), colData(spe_sub))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  geom_point(size = 0.1) +
  theme_classic()

## nnSVG at single cell
# spe_sub <- nnSVG::nnSVG(
#   spe_sub,
#   assay_name = "lognorm",
#   BPPARAM = BiocParallel::MulticoreParam(),
#   verbose = TRUE
# )
# View(as.data.frame(rowData(spe_sub)))

## rasterization
res <- 25
spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
dim(spe_sub_rast)

df <- data.frame(x = spatialCoords(spe_sub_rast)[,1], y = spatialCoords(spe_sub_rast)[,2], transcripts = colMeans(assay(spe_sub_rast, "pixelval")))
ggplot(df, aes(x = x, y = y, fill = transcripts)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean transcripts/bin") +
  theme_classic()

## nnSVG
spe_sub_rast <- try({nnSVG::nnSVG(
  spe_sub_rast,
  assay_name = "pixelval",
  BPPARAM = BiocParallel::MulticoreParam()
)})
df <- tibble::rownames_to_column(as.data.frame(rowData(spe_sub_rast)), var = "gene")
View(df)

## visual inspection
gene <- "SERPINA3"
df <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], gene = assay(spe_sub, "counts")[gene,])
p1 <- ggplot(df, aes(x = x, y = y, col = gene)) +
  geom_point(size = 0.1) +
  scale_color_distiller(name = gene, palette = "OrRd", direction = 1) +
  labs(title = "Single cell (raw counts)") +
  theme_classic()
df <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], gene = assay(spe_sub, "lognorm")[gene,])
p2 <- ggplot(df, aes(x = x, y = y, col = gene)) +
  geom_point(size = 0.1) +
  scale_color_distiller(name = gene, palette = "OrRd", direction = 1) +
  labs(title = "Single cell (pseudolog normalized)") +
  theme_classic()
df <- data.frame(x = spatialCoords(spe_sub_rast)[,1], y = spatialCoords(spe_sub_rast)[,2], gene = assay(spe_sub_rast, "pixelval")[gene,])
p3 <- ggplot(df, aes(x = x, y = y, fill = gene)) +
  geom_tile() +
  scale_fill_distiller(name = gene, palette = "OrRd", direction = 1) +
  labs(title = paste0("Resolution = ", res)) +
  theme_classic()

grid.arrange(p1, p2, p3, ncol = 3, top = paste0(gene))
