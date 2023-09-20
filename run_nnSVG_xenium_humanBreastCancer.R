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

par(mfrow=c(1,1))

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = "outputs/xenium_humanBreastCancer_preprocessed.RDS")

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  geom_point(size = 0.1) +
  theme_classic()

ct_labels <- colData(spe)$cluster


# Run methods -------------------------------------------------------------

## filter by cluster --> rasterize --> nnSVG
ct_label <- 1

## subset cluster of interest
spe_sub <- spe[,spe$cluster == ct_label]

df <- data.frame(spatialCoords(spe_sub), colData(spe_sub))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  geom_point(size = 0.1) +
  theme_classic()

## rasterization
res <- 100
spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

df <- data.frame(x = spatialCoords(spe_sub_rast)[,1], y = spatialCoords(spe_sub_rast)[,2], transcripts = colMeans(assay(spe_sub_rast, "lognorm")))
ggplot(df, aes(x = x, y = y, fill = transcripts)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean transcripts/bin") +
  theme_classic()

## nnSVG
spe_sub_rast <- nnSVG::nnSVG(
  spe_sub_rast,
  assay_name = "pixelval",
  BPPARAM = BiocParallel::MulticoreParam()
)
View(as.data.frame(rowData(spe_sub_rast)))

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
