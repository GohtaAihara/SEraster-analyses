## This file integrates rasterization with CooccurrenceAffinity to analyze 10X Genomics Xenium mouse whole pup datasets. Produce figures for 

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
library(ggstar)

par(mfrow=c(1,1))

dataset_name = "xenium_mousePup"
method = "CooccurrenceAffinity"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = "outputs/xenium_mousePup_preprocessed.RDS")

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  geom_point(size = 0.1) +
  theme_classic()

ct_labels <- colData(spe)$cluster

# Run methods -------------------------------------------------------------

## compute relative enrichment and binarize for various resolutions
res_list <- list(50, 100, 200, 400)

## parameters for CooccurrenceAffinity
## set confidence interval method ("CP", "Blaker", "midQ", or "midP")
CI_method <- "Blaker"
## create a dictionary and index of CI method
CI_method_dict <- c("CP" = 6, "Blaker" = 7, "midQ" = 8, "midP" = 9)
## set confidence interval level (default = 0.95)
CI_lev <- 0.95
## set pval method ("Blaker" or "midP", default = "Blaker")
pval_method <- "Blaker"

for (res in res_list) {
  print(paste0("Resolution = ", res))
  ## rasterize
  spe_rast <- SEraster::rasterizeCellType(spe, "cluster", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
  
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
  
  ## save updated SpatialExperiment object
  saveRDS(spe_rast, file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))
  
  ## compute affinity MLE using CooccurrenceAffinity (store the log-affinity metric, CI, and pvalue for each pair)
  ## get pair combinations
  non_self <- combn(levels(ct_labels), 2, simplify = FALSE)
  self <- lapply(levels(ct_labels), function(ct_label) c(ct_label, ct_label))
  pairs <- c(non_self, self)
  
  ## multiply pixel values for each pair of cell types
  affinity_results <- do.call(rbind, lapply(pairs, function(pair) {
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
  ## save results
  saveRDS(affinity_results, file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
}


# ## rasterize
# res <- 100
# spe_rast <- SEraster::rasterizeCellType(spe, "cluster", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
# 
# ## compute relative enrichment (RE) metric (the output is dense matrix)
# mat <- assay(spe_rast, "pixelval")
# mat_re <- do.call(rbind, lapply(rownames(spe_rast), function(ct_label) {
#   ## relative enrichment = celltype observed / celltype expected = celltype observed / (celltype frequency * total # of cells in the pixel)
#   mat[ct_label,] / (sum(mat[ct_label,]) / sum(mat) * colSums(mat))
# }))
# rownames(mat_re) <- rownames(mat)
# 
# ## binarize (1 if RE >= 1, 0 if RE < 1)
# mat_bin <- ifelse(mat_re >= 1, 1, 0)
# 
# ## add RE and binary layers to SpatialExperiment object
# assays(spe_rast) <- list(pixelval = assay(spe_rast, "pixelval"), re = mat_re, bin = mat_bin)
# ## save updated SpatialExperiment object
# saveRDS(spe_rast, file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))
# 
# ## compute affinity MLE using CooccurrenceAffinity
# mat_bin <- assay(spe_rast, "bin")
# 
# ## set confidence interval method ("CP", "Blaker", "midQ", or "midP")
# CI_method <- "Blaker"
# ## create a dictionary and index of CI method
# CI_method_dict <- c("CP" = 6, "Blaker" = 7, "midQ" = 8, "midP" = 9)
# ## set confidence interval level (default = 0.95)
# CI_lev <- 0.95
# ## set pval method ("Blaker" or "midP", default = "Blaker")
# pval_method <- "Blaker"
# 
# ## store the log-affinity metric, CI, and pvalue for each pair
# ## get pair combinations
# non_self <- combn(levels(ct_labels), 2, simplify = FALSE)
# self <- lapply(levels(ct_labels), function(ct_label) c(ct_label, ct_label))
# pairs <- c(non_self, self)
# 
# ## multiply pixel values for each pair of cell types
# affinity_results <- do.call(rbind, lapply(pairs, function(pair) {
#   ## create a 2x2 contingency table of counts
#   cont_tab <- table(factor(mat_bin[pair[1],], levels = c(0,1)), factor(mat_bin[pair[2],], levels = c(0,1)))
#   X <- cont_tab[2,2]
#   mA <- sum(cont_tab[2,])
#   mB <- sum(cont_tab[,2])
#   N <- sum(cont_tab)
#   
#   out <- CooccurrenceAffinity::ML.Alpha(X,c(mA,mB,N), lev = CI_lev, pvalType = pval_method)
#   
#   ## set index for the chosen CI method
#   CI_idx <- CI_method_dict[CI_method][[1]]
#   
#   return(data.frame(
#     pair = paste(pair, collapse = " & "),
#     celltypeA = factor(pair[1], levels = levels(ct_labels)),
#     celltypeB = factor(pair[2], levels = levels(ct_labels)),
#     alpha = out$est, ci.min = out[CI_idx][[1]][1], 
#     ci.max = out[CI_idx][[1]][2], 
#     pval = out$pval)
#   )
# }))
# ## save results
# saveRDS(affinity_results, file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res, ".RDS")))

# Plot --------------------------------------------------------------------

col_clu <- gg_color_hue(length(levels(ct_labels)))

## Figure (single cell spatial visualization)
df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  coord_fixed() +
  rasterize(geom_point(size = 0.1, stroke = 0), dpi = 300) +
  scale_color_manual(values = col_clu) +
  labs(title = "Single cell",
       col = "Cluster") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
## save plot
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_clusters_singlecell.pdf")), width = 6, height = 7, dpi = 300)

## Figure (rasterized spatial visualization)
res <- 100
spe_rast <- readRDS(file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))
ct_label <- 1

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
plt_comb <- grid.arrange(p1, p2, p3, ncol = 3, top = paste0("Cluster ", ct_label))
ggsave(plt_comb, filename = here("plots", dataset_name, paste0(dataset_name, "_rasterized_data.pdf")), width = 8, height = 3, dpi = 300)

 ## Figure (alpha MLE vs. -log10(p value))
res <- 100
df <- readRDS(file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
cluster <- 26
## subset by cluster
df_sub <- df[(df$celltypeA == cluster | df$celltypeB == cluster) & (df$celltypeA != df$celltypeB),]
## subset by alpha MLE and pval
alpha_mle <- 0
pval <- 0.05
df_sub <- df_sub[df_sub$alpha >= alpha_mle & df_sub$pval <= pval,]

## plot alpha MLE vs. -log10(p value)
ggplot(df_sub, aes(x = alpha, y = -log10(pval), col = pair)) +
  geom_point() +
  geom_errorbar(xmin = df_sub$ci.min, xmax = df_sub$ci.max) +
  theme_bw()

## plot single-cell
for (i in seq(dim(df_sub)[1])) {
  pair <- factor(unlist(strsplit(df_sub[i,"pair"], " & ")), levels = levels(ct_labels))
  ## subset by clusters
  spe_sub <- spe[,spe$cluster %in% pair]
  
  ## plot
  df_plt <- data.frame(spatialCoords(spe))
  df_plt_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], cluster = colData(spe_sub)$cluster)
  ggplot(df_plt, aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 0.1), dpi = 300) +
    rasterise(geom_point(data = df_plt_sub, aes(x = x, y = y, col = cluster), size = 0.1), dpi = 300) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
    labs(title = paste0("Alpha MLE = ", sprintf("%.2f", df_sub[i,"alpha"]), ", p value = ", sprintf("%.2f", df_sub[i,"pval"])),
         x = "x (um)",
         y = "y (um)",
         col = "Clusters") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(filename = here("plots", dataset_name, method, paste0("singlecell_pairs_", paste(pair, collapse = "_"), ".pdf")), width = 4, height = 5, dpi = 300)
}


## Figure (number of s.s. pair-wise celltype colocalization)
res_list <- list(50, 100, 200, 400)
alpha <- 0.05
df_plt <- do.call(rbind, lapply(res_list, function(res) {
  df <- readRDS(file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
  return(data.frame(resolution = res, num_localized_pairs = sum(df$pval <= alpha), prop_localized_pairs = sum(df$pval <= alpha)/length(unique(df$pair))))
}))
ggplot(df_plt, aes(x = resolution, y = num_localized_pairs)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = unique(df_plt$resolution)) +
  labs(title = paste0("Number of localized pairs (p value <= ", res, ")"),
       x = "Resolution (um)",
       y = "Number of localized pairs") +
  theme_bw()

ggplot(df_plt, aes(x = resolution, y = prop_localized_pairs)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = unique(df_plt$resolution)) +
  labs(title = paste0("Proportion of localized pairs (p value <= ", res, ")"),
       x = "Resolution (um)",
       y = "Proportion of localized pairs") +
  theme_bw()

## Figure (CooccurrenceAffinity heatmap)
res_list <- list(50, 100, 200, 400)
for (res in res_list) {
  df <- readRDS(file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
  ggplot(df, aes(x = celltypeA, y = celltypeB, fill = alpha)) +
    geom_tile() +
    scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red") +
    labs(title = paste0("Pair-wise cell type colocalization (Resolution = ", res, ")"),
         x = "Cluster A",
         y = "Cluster B") +
    theme_bw()
  ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_heatmap_alpha_without_clustering_resolution_", res, ".pdf")), width = 12, height = 10, dpi = 300)
  
  ## create symmetric data
  df_flipped <- df[df$celltypeA != df$celltypeB,]
  df_flipped[,c("celltypeA", "celltypeB")] <- df_flipped[,c("celltypeB", "celltypeA")]
  df_sym <- rbind(df, df_flipped)
  
  ## use non-symmetric (non-redundant) data
  ## reset label order
  df_sym <- df_sym %>%
    mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
           celltypeB = factor(celltypeB, levels(ct_labels)))
  ## reorganize into matrix
  df_heatmap_non_sym <- cast(df, celltypeA ~ celltypeB, value = "alpha")
  df_heatmap_non_sym <- df_heatmap_non_sym[,-1]
  ## cluster
  hc_non_sym <- hclust(dist(df_heatmap_non_sym))
  ## reorder labels
  df_sym$celltypeA <- factor(df_sym$celltypeA, levels = rownames(df_heatmap_non_sym)[hc_non_sym$order])
  df_sym$celltypeB <- factor(df_sym$celltypeB, levels = colnames(df_heatmap_non_sym)[hc_non_sym$order])
  ## plot
  ggplot(df_sym, aes(x = celltypeA, y = celltypeB, fill = alpha, col = pval <= alpha)) +
    coord_fixed() +
    geom_tile(linewidth = 0.5) +
    scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red") +
    scale_color_manual(name = paste0("p value <= ", alpha), values = c("grey", "black")) +
    labs(title = paste0("Pair-wise cell type colocalization (Resolution = ", res, ")"),
         x = "Cluster A",
         y = "Cluster B") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_heatmap_alpha_with_nonsym_clustering_resolution_", res, ".pdf")), width = 12, height = 10, dpi = 300)
  
  ## use symmetric (redundant) data
  ## reset label order
  df_sym <- df_sym %>%
    mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
           celltypeB = factor(celltypeB, levels(ct_labels)))
  ## reorganize into matrix
  df_heatmap_sym <- cast(df_sym, celltypeA ~ celltypeB, value = "alpha")
  df_heatmap_sym <- df_heatmap_sym[,-1]
  isSymmetric.matrix(as.matrix(df_heatmap_sym))
  ## cluster
  hc_sym <- hclust(dist(df_heatmap_sym))
  ## reorder labels
  df_sym$celltypeA <- factor(df_sym$celltypeA, levels = rownames(df_heatmap_sym)[hc_sym$order])
  df_sym$celltypeB <- factor(df_sym$celltypeB, levels = colnames(df_heatmap_sym)[hc_sym$order])
  ## plot
  ggplot(df_sym, aes(x = celltypeA, y = celltypeB, fill = alpha, col = pval <= alpha)) +
    coord_fixed() +
    geom_tile(linewidth = 0.5) +
    scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red") +
    scale_color_manual(name = paste0("p value <= ", alpha), values = c("grey", "black")) +
    labs(title = paste0("Pair-wise cell type colocalization (Resolution = ", res, ")"),
         x = "Cluster A",
         y = "Cluster B") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_heatmap_alpha_with_sym_clustering_resolution_", res, ".pdf")), width = 12, height = 10, dpi = 300)
}

## Figure (CooccurrenceAffinity heatmap, pyramid, vertical)
res <- 100
alpha <- 0.05
## load data
df <- readRDS(file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
## create symmetric data
df_flipped <- df[df$celltypeA != df$celltypeB,]
df_flipped[,c("celltypeA", "celltypeB")] <- df_flipped[,c("celltypeB", "celltypeA")]
df_sym <- rbind(df, df_flipped)
## use symmetric (redundant) data
## reset label order
df_sym <- df_sym %>%
  mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
         celltypeB = factor(celltypeB, levels(ct_labels)))
## reorganize into matrix
df_heatmap_sym <- cast(df_sym, celltypeA ~ celltypeB, value = "alpha")
df_heatmap_sym <- df_heatmap_sym[,-1]
isSymmetric.matrix(as.matrix(df_heatmap_sym))
## cluster
hc_sym <- hclust(dist(df_heatmap_sym))
## reorder labels
df_sym$celltypeA <- factor(df_sym$celltypeA, levels = rownames(df_heatmap_sym)[hc_sym$order])
df_sym$celltypeB <- factor(df_sym$celltypeB, levels = colnames(df_heatmap_sym)[hc_sym$order])
## plot
ggplot(df_sym, aes(x = celltypeA, y = celltypeB, fill = alpha, col = pval <= alpha)) +
  coord_fixed() +
  geom_tile(linewidth = 0.5) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red") +
  scale_color_manual(name = paste0("p value <= ", alpha), values = c("grey", "black")) +
  labs(x = "Cluster A",
       y = "Cluster B") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0),
        axis.text.y = element_text(angle = 45, vjust = 0, hjust=1))
ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_heatmap_alpha_with_sym_clustering_resolution_", res, "_pyramid.pdf")), width = 12, height = 10, dpi = 300)

# ## Figure (CooccurrenceAffinity heatmap, pyramid, horizontal)
# res <- 100
# alpha <- 0.05
# ## load data
# df <- readRDS(file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
# ## create symmetric data
# df_flipped <- df[df$celltypeA != df$celltypeB,]
# df_flipped[,c("celltypeA", "celltypeB")] <- df_flipped[,c("celltypeB", "celltypeA")]
# df_sym <- rbind(df, df_flipped)
# ## use symmetric (redundant) data
# ## reset label order
# df_sym <- df_sym %>%
#   mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
#          celltypeB = factor(celltypeB, levels(ct_labels)))
# ## reorganize into matrix
# df_heatmap_sym <- cast(df_sym, celltypeA ~ celltypeB, value = "alpha")
# df_heatmap_sym <- df_heatmap_sym[,-1]
# isSymmetric.matrix(as.matrix(df_heatmap_sym))
# ## cluster
# hc_sym <- hclust(dist(df_heatmap_sym))
# ## reorder labels
# df_sym$celltypeA <- factor(df_sym$celltypeA, levels = rownames(df_heatmap_sym)[hc_sym$order])
# df_sym$celltypeB <- factor(df_sym$celltypeB, levels = colnames(df_heatmap_sym)[hc_sym$order])
# ## plot
# ggplot(df_sym, aes(x = celltypeA, y = celltypeB, fill = alpha, col = pval <= alpha)) +
#   coord_fixed() +
#   geom_tile(linewidth = 0.5) +
#   scale_x_discrete(position = "bottom") +
#   scale_y_discrete(position = "right") +
#   scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red") +
#   scale_color_manual(name = paste0("p value <= ", alpha), values = c("grey", "black")) +
#   labs(x = "Cluster A",
#        y = "Cluster B") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 135, vjust = 0.5, hjust=0.5),
#         axis.text.y = element_text(angle = 135, vjust = 0, hjust=0))
# ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_heatmap_alpha_with_sym_clustering_resolution_", res, "_pyramid_horizontal.pdf")), width = 12, height = 10, dpi = 300)

## Figure (visual inspection of spatial niches)
niches <- list(
  factor(c(6, 8, 12, 15, 20, 56), levels = levels(ct_labels)),
  factor(c(10, 22, 46, 42, 49), levels = levels(ct_labels)),
  factor(c(13, 30, 51, 54, 57, 64), levels = levels(ct_labels)),
  factor(c(11, 25, 28, 31, 38, 45, 60), levels = levels(ct_labels)),
  factor(c(24, 26, 33, 35, 36, 55), levels = levels(ct_labels)),
  factor(c(39, 41, 47, 66), levels = levels(ct_labels)),
  factor(c(29, 50, 58, 59), levels = levels(ct_labels)),
  factor(c(9, 23), levels = levels(ct_labels))
)

for (clusters in niches) {
  ## subset by clusters
  spe_sub <- spe[,spe$cluster %in% clusters]
  
  ## plot
  df <- data.frame(spatialCoords(spe))
  df_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], cluster = colData(spe_sub)$cluster)
  ggplot(df, aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 0.01), dpi = 300) +
    rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = cluster), size = 0.01), dpi = 300) +
    scale_color_manual(values = col_clu[clusters]) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
    labs(title = paste0("Clusters = ", paste(clusters, collapse = ", ")),
         x = "x (um)",
         y = "y (um)",
         col = "Cluster") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(filename = here("plots", dataset_name, method, paste0("singlecell_niche_clusters_", paste(clusters, collapse = "_"), ".pdf")), width = 4, height = 5, dpi = 300)
}

## for poster
for (clusters in niches) {
  ## subset by clusters
  spe_sub <- spe[,spe$cluster %in% clusters]
  
  ## plot
  df <- data.frame(spatialCoords(spe))
  df_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], cluster = colData(spe_sub)$cluster)
  ggplot(df, aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 0.1, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = cluster), size = 0.1, stroke = 0), dpi = 300) +
    scale_color_manual(values = col_clu[clusters]) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
    labs(title = paste0("Clusters = ", paste(clusters, collapse = ", ")),
         x = "x (um)",
         y = "y (um)") +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(filename = here("plots", dataset_name, method, paste0("singlecell_niche_clusters_", paste(clusters, collapse = "_"), "_v2.pdf")), width = 6, height = 12, dpi = 300)
}

## Supplementary figure 2
## heatmap clustered at resolution 100 --> apply the same orders to other resolutions
res_interest <- 100
df <- readRDS(file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res_interest, ".RDS")))

## create symmetric data
df_flipped <- df[df$celltypeA != df$celltypeB,]
df_flipped[,c("celltypeA", "celltypeB")] <- df_flipped[,c("celltypeB", "celltypeA")]
df_sym <- rbind(df, df_flipped)

## use symmetric (redundant) data
## reset label order
df_sym <- df_sym %>%
  mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
         celltypeB = factor(celltypeB, levels(ct_labels)))
## reorganize into matrix
df_heatmap_sym <- cast(df_sym, celltypeA ~ celltypeB, value = "alpha")
df_heatmap_sym <- df_heatmap_sym[,-1]
rownames(df_heatmap_sym) <- colnames(df_heatmap_sym)
isSymmetric.matrix(as.matrix(df_heatmap_sym))
## cluster
hc_sym_interest <- hclust(dist(df_heatmap_sym))
## reorder labels (use hc_sym_interest)
df_sym$celltypeA <- factor(df_sym$celltypeA, levels = rownames(df_heatmap_sym)[hc_sym_interest$order])
df_sym$celltypeB <- factor(df_sym$celltypeB, levels = colnames(df_heatmap_sym)[hc_sym_interest$order])
## plot
cutoff <- min(abs(range(df_sym$alpha)))
lim <- c(-cutoff,cutoff)
# df_plt <- df_sym %>%
#   mutate(alpha = Winsorize(alpha, min(lim), max(lim)))
# ggplot(df_plt, aes(x = celltypeA, y = celltypeB, fill = alpha, col = pval <= alpha)) +
#   coord_fixed() +
#   geom_tile(linewidth = 0.5) +
#   scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red", limits = lim) +
#   scale_color_manual(name = paste0("p value ≤ ", alpha), values = c("grey", "black")) +
#   labs(title = paste0("Pair-wise cell type colocalization (Resolution = ", res_interest, ")"),
#        x = "Cluster A",
#        y = "Cluster B") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df_plt <- df_sym %>%
  mutate(
    alpha = Winsorize(alpha, min(lim), max(lim)),
    significance = case_when(
      pval <= 0.001 ~ "***",
      pval <= 0.01 ~ "**",
      pval <= 0.05 ~ "*"
    )
  )
ggplot(df_plt, aes(x = celltypeA, y = celltypeB, fill = alpha, label = significance)) +
  coord_fixed() +
  geom_tile(linewidth = 0.5) +
  geom_text(angle = 45) +
  scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red", limits = lim) +
  labs(title = paste0("Pair-wise cell type colocalization (Resolution = ", res_interest, ")"),
       x = "Cluster A",
       y = "Cluster B") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

res_list <- list(50, 100, 200, 400)
for (i in seq_along(res_list)) {
  res <- res_list[[i]]
  df <- readRDS(file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
  
  ## create symmetric data
  df_flipped <- df[df$celltypeA != df$celltypeB,]
  df_flipped[,c("celltypeA", "celltypeB")] <- df_flipped[,c("celltypeB", "celltypeA")]
  df_sym <- rbind(df, df_flipped)
  
  ## use symmetric (redundant) data
  ## reset label order
  df_sym <- df_sym %>%
    mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
           celltypeB = factor(celltypeB, levels(ct_labels)))
  ## reorganize into matrix
  df_heatmap_sym <- cast(df_sym, celltypeA ~ celltypeB, value = "alpha")
  df_heatmap_sym <- df_heatmap_sym[,-1]
  rownames(df_heatmap_sym) <- colnames(df_heatmap_sym)
  isSymmetric.matrix(as.matrix(df_heatmap_sym))
  
  ## reorder labels (use hc_sym_interest)
  df_sym$celltypeA <- factor(df_sym$celltypeA, levels = rownames(df_heatmap_sym)[hc_sym_interest$order])
  df_sym$celltypeB <- factor(df_sym$celltypeB, levels = colnames(df_heatmap_sym)[hc_sym_interest$order])
  ## plot
  cutoff <- min(abs(range(df_sym$alpha)))
  lim <- c(-cutoff,cutoff)
  df_plt <- df_sym %>%
    mutate(
      alpha = Winsorize(alpha, min(lim), max(lim)),
      significance = case_when(
        pval <= 0.05 ~ "*"
      )
    )
  # df_plt <- df_sym %>%
  #   mutate(alpha = Winsorize(alpha, min(lim), max(lim)))
  # df_plt_sig <- df_plt[df_plt$pval <= 0.05,]
  ggplot(df_plt, aes(x = celltypeA, y = celltypeB, fill = alpha, label = significance)) +
    coord_fixed() +
    geom_tile(color = "gray") +
    geom_text(size = 5, vjust = 0.8, hjust = 0.5) +
    # geom_star(data = df_plt_sig, aes(x = celltypeA, y = celltypeB), fill = "black", size = 1.25) +
    scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red", limits = lim) +
    labs(title = paste0("Pair-wise cell type colocalization (Resolution = ", res, ")"),
         x = "Cluster A",
         y = "Cluster B") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_heatmap_alpha_with_sym_clustering_order_", res_interest, "_resolution_", res, ".pdf")), width = 12, height = 10, dpi = 300)
}
