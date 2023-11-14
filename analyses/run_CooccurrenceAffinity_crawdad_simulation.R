## This file integrates rasterization with CooccurrenceAffinity to analyze CRAWDAD simulated dataset. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::document()
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(ggrastr)
library(gridExtra)
library(here)
library(CooccurrenceAffinity)
library(tidyr)
library(tibble)
library(dplyr)
library(reshape)

par(mfrow=c(1,1))

dataset_name = "crawdad_simulation"
method = "CooccurrenceAffinity"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_preprocessed.RDS")))

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = celltypes)) +
  geom_point(size = 1, stroke = 0) +
  theme_classic()

ct_labels <- as.factor(colData(spe)$celltypes)

# Run methods -------------------------------------------------------------

## compute relative enrichment and binarize for various resolutions
res_list <- seq(20, 1000, by = 20)

## parameters for CooccurrenceAffinity
## set confidence interval method ("CP", "Blaker", "midQ", or "midP")
CI_method <- "Blaker"
## create a dictionary and index of CI method
CI_method_dict <- c("CP" = 6, "Blaker" = 7, "midQ" = 8, "midP" = 9)
## set confidence interval level (default = 0.95)
CI_lev <- 0.95
## set pval method ("Blaker" or "midP", default = "Blaker")
pval_method <- "Blaker"

## permutate by rotations
n_rotation <- 10
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)

affinity_results <- do.call(rbind, lapply(res_list, function(res) {
  deg_results <- do.call(rbind, lapply(angle_deg_list, function(deg) {
    print(paste0("Resolution = ", res))
    print(paste0("Rotation (degrees): ", deg))
    
    ## rotate xy coordinates
    spe_rotated <- SpatialExperiment::SpatialExperiment(
      spatialCoords = rotateAroundCenter(spatialCoords(spe), deg),
      colData = colData(spe)
    )
    
    ## rasterize
    spe_rast <- SEraster::rasterizeCellType(spe_rotated, "celltypes", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
    
    ## compute relative enrichment (RE) metric (the output is dense matrix)
    mat <- assay(spe_rast, "pixelval")
    mat_re <- do.call(rbind, lapply(rownames(spe_rast), function(ct_label) {
      ## relative enrichment = celltype observed / celltype expected = celltype observed / (celltype frequency * total # of cells in the pixel)
      mat[ct_label,] / (sum(mat[ct_label,]) / sum(mat) * colSums(mat))
    }))
    rownames(mat_re) <- rownames(mat)
    
    ## binarize (1 if RE >= 1, 0 if RE < 1)
    mat_bin <- ifelse(mat_re >= 1, 1, 0)
    
    # ## add RE and binary layers to SpatialExperiment object
    # assays(spe_rast) <- list(pixelval = assay(spe_rast, "pixelval"), re = mat_re, bin = mat_bin)
    # 
    # ## save updated SpatialExperiment object
    # saveRDS(spe_rast, file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))
    
    ## compute affinity MLE using CooccurrenceAffinity (store the log-affinity metric, CI, and pvalue for each pair)
    ## get pair combinations
    non_self <- combn(levels(ct_labels), 2, simplify = FALSE)
    self <- lapply(levels(ct_labels), function(ct_label) c(ct_label, ct_label))
    pairs <- c(non_self, self)
    
    ## multiply pixel values for each pair of cell types
    affinity_single <- do.call(rbind, lapply(pairs, function(pair) {
      ## create a 2x2 contingency table of counts
      cont_tab <- table(factor(mat_bin[pair[1],], levels = c(0,1)), factor(mat_bin[pair[2],], levels = c(0,1)))
      X <- cont_tab[2,2]
      mA <- sum(cont_tab[2,])
      mB <- sum(cont_tab[,2])
      N <- sum(cont_tab)
      
      out <- CooccurrenceAffinity::ML.Alpha(X,c(mA,mB,N), lev = CI_lev, pvalType = pval_method)
      
      ## set index for the chosen CI method
      CI_idx <- CI_method_dict[CI_method][[1]]
      
      return(data.frame(
        resolution = res,
        rotation_deg = deg,
        pair = paste(pair, collapse = " & "),
        celltypeA = factor(pair[1], levels = levels(ct_labels)),
        celltypeB = factor(pair[2], levels = levels(ct_labels)),
        X = X,
        mA = mA,
        mB = mB,
        N = N,
        alpha = out$est, ci.min = out[CI_idx][[1]][1], 
        ci.max = out[CI_idx][[1]][2], 
        pval = out$pval)
      )
    }))
  }))
}))
## save results
saveRDS(affinity_results, file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_n_rotation_", n_rotation, ".RDS")))


# Plots -------------------------------------------------------------------

col_clu <- gg_color_hue(length(levels(ct_labels)))
names(col_clu) <- levels(ct_labels)

## Figure 3a (cell type single-cell resolution)
df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], celltypes = colData(spe)$celltypes)
ggplot(df, aes(x = x, y = y, col = celltypes)) +
  coord_fixed() +
  rasterise(geom_point(size = 2, stroke = 0), dpi = 300) +
  scale_color_manual(name = "Cell type", values = col_clu) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Single-cell Resolution") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_ct_labels_sc.pdf")))

## Figure (resolution vs. alpha)
n_rotation <- 10
df <- readRDS(file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_n_rotation_", n_rotation, ".RDS")))

p1 <- ggplot(df, aes(x = resolution, y = alpha, col = pair)) +
  geom_line() +
  geom_errorbar(aes(ymin = ci.min, ymax = ci.max)) +
  geom_point() +
  # scale_x_continuous(breaks = unique(df$resolution)) +
  theme_bw()
show(p1)

plotly::ggplotly(p1)

alpha <- 0.05
p2 <- ggplot(df, aes(x = resolution, y = pval, col = pair)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  # scale_x_continuous(breaks = unique(df$resolution)) +
  theme_bw()
show(p2)

plotly::ggplotly(p2)


# Further exploration -----------------------------------------------------

resolution <- 1000

## create bbox
pos <- spatialCoords(spe)
bbox <- sf::st_bbox(c(
  xmin = floor(min(pos[,1])-resolution/2), 
  xmax = ceiling(max(pos[,1])+resolution/2), 
  ymin = floor(min(pos[,2])-resolution/2), 
  ymax = ceiling(max(pos[,2])+resolution/2)
))

## create grid for rasterization
grid <- sf::st_make_grid(bbox, cellsize = resolution)
grid_coord <- st_coordinates(grid)

## plot single cell (visualize how it would be split)
df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = celltypes)) +
  coord_fixed() +
  geom_point(size = 1, stroke = 0) +
  # scale_color_viridis_c() +
  geom_hline(yintercept = grid_coord[,2], linetype = "solid", color = "black") +
  geom_vline(xintercept = grid_coord[,1], linetype = "solid",color = "black") +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )

## rasterize
spe_rast <- SEraster::rasterizeCellType(spe, "celltypes", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())

## compute relative enrichment (RE) metric (the output is dense matrix)
mat <- assay(spe_rast, "pixelval")
mat_re <- do.call(rbind, lapply(rownames(spe_rast), function(ct_label) {
  ## relative enrichment = celltype observed / celltype expected = celltype observed / (celltype frequency * total # of cells in the pixel)
  mat[ct_label,] / (sum(mat[ct_label,]) / sum(mat) * colSums(mat))
}))
rownames(mat_re) <- rownames(mat)

## binarize (1 if RE >= 1, 0 if RE < 1)
mat_bin <- ifelse(mat_re >= 1, 1, 0)

## add RE and binary layers to SpatialExperiment object
assays(spe_rast) <- list(pixelval = assay(spe_rast, "pixelval"), re = mat_re, bin = mat_bin)

ct_label <- "B"

df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], pixelval = assay(spe_rast, "pixelval")[ct_label,], re = assay(spe_rast, "re")[ct_label,], bin = factor(assay(spe_rast, "bin")[ct_label,], levels = c(0,1)))

p1 <- ggplot(df, aes(x = x, y = y, fill = pixelval)) +
  coord_fixed() +
  geom_tile() +
  scale_fill_viridis_c(name = "cells/pixel") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p2 <- ggplot(df, aes(x = x, y = y, fill = re)) +
  coord_fixed() +
  geom_tile() +
  scale_fill_viridis_c(name = "RE") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p3 <- ggplot(df, aes(x = x, y = y, fill = bin)) +
  coord_fixed() +
  geom_tile() +
  scale_fill_viridis_d(name = "Binarized") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
grid.arrange(p1, p2, p3, ncol = 3, top = paste0("Cluster ", ct_label))
# ggsave(plt_comb, filename = here("plots", dataset_name, paste0(dataset_name, "_rasterized_data.pdf")), width = 8, height = 3, dpi = 300)
