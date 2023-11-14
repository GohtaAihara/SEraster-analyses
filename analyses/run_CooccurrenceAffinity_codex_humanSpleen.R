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

par(mfrow=c(1,1))

dataset_name = "codex_humanSpleen"
method = "CooccurrenceAffinity"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_preprocessed.RDS")))

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = celltypes)) +
  geom_point(size = 0.5, stroke = 0) +
  theme_classic()

ct_labels <- as.factor(colData(spe)$celltypes)

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


# Plot --------------------------------------------------------------------

col_clu <- gg_color_hue(length(levels(ct_labels)))
names(col_clu) <- levels(ct_labels)

## Figure 3a (cell type single-cell resolution)
df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], celltypes = colData(spe)$celltypes)
ggplot(df, aes(x = x, y = y, col = celltypes)) +
  coord_fixed() +
  rasterise(geom_point(size = 0.5, stroke = 0), dpi = 300) +
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

## Figure 3 (spatial plots)
res_list <- list("singlecell", 50, 100, 200, 400)

for (res in res_list) {
  if (res == "singlecell") {
    ## plot
    df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], transcripts = colSums(assay(spe, "lognorm")))
    plt <- ggplot(df, aes(x = x, y = y, col = transcripts)) +
      coord_fixed() +
      geom_point(size = 0.1) +
      scale_color_viridis_c() +
      labs(title = "Single cell") +
      theme_bw() +
      theme(
        legend.position="none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
      )
  } else {
    ## rasterize
    spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
    
    ## plot
    df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], transcripts = colSums(assay(spe_rast)))
    plt <- ggplot(df, aes(x = x, y = y, fill = transcripts)) +
      coord_fixed() +
      geom_tile() +
      scale_fill_viridis_c() +
      labs(title = paste0("Resolution = ", res, " um")) +
      theme_bw() +
      theme(
        legend.position="none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
      )
  }
  ## save plot
  ggsave(plot = plt, filename = here("plots", dataset_name, paste0(dataset_name, "_tot_lognorm_", res, ".pdf")))
}

## Figure 3B (rasterized --> re --> bin)
res <- 100
spe_rast <- readRDS(file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))
ct_label <- "Fol B cells"

df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], pixelval = assay(spe_rast, "pixelval")[ct_label,], re = assay(spe_rast, "re")[ct_label,], bin = factor(assay(spe_rast, "bin")[ct_label,], levels = c(0,1)))

p1 <- ggplot(df, aes(x = x, y = y, fill = pixelval)) +
  coord_fixed() +
  geom_tile() +
  scale_fill_gradient(name = "cells/pixel", low = "gray", high = col_clu[["Fol B cells"]]) +
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
  scale_fill_gradient(name = "RE", low = "gray", high = col_clu[["Fol B cells"]]) +
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
  scale_fill_manual(name = "Binarized", values = c("gray",col_clu[["Fol B cells"]])) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
plt_comb <- grid.arrange(p1, p2, p3, ncol = 3, top = paste0("Cell type: ", ct_label))
ggsave(plt_comb, filename = here("plots", dataset_name, paste0(dataset_name, "_rasterized_data_ct_col.pdf")), width = 8, height = 3, dpi = 300)

p1 <- ggplot(df, aes(x = x, y = y, fill = pixelval)) +
  coord_fixed() +
  geom_tile() +
  scale_fill_viridis_c(name = "cells/pixel", option = "inferno") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p2 <- ggplot(df, aes(x = x, y = y, fill = re)) +
  coord_fixed() +
  geom_tile() +
  scale_fill_viridis_c(name = "RE", option = "inferno") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p3 <- ggplot(df, aes(x = x, y = y, fill = bin)) +
  coord_fixed() +
  geom_tile() +
  scale_fill_viridis_d(name = "Binarized", option = "inferno") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
plt_comb <- grid.arrange(p1, p2, p3, ncol = 3, top = paste0("Cell type: ", ct_label))
ggsave(plt_comb, filename = here("plots", dataset_name, paste0(dataset_name, "_rasterized_data_inferno.pdf")), width = 8, height = 3, dpi = 300)


## Figure (CooccurrenceAffinity heatmap)
res_list <- list(50, 100, 200, 400)
alpha <- 0.05
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
  rownames(df_heatmap_non_sym) <- colnames(df_heatmap_non_sym)
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
  rownames(df_heatmap_sym) <- colnames(df_heatmap_sym)
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

## order based on 1 resolution
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
ggplot(df_sym, aes(x = celltypeA, y = celltypeB, fill = alpha, col = pval <= alpha)) +
  coord_fixed() +
  geom_tile(linewidth = 0.5) +
  scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red") +
  scale_color_manual(name = paste0("p value <= ", alpha), values = c("grey", "black")) +
  labs(title = paste0("Pair-wise cell type colocalization (Resolution = ", res_interest, ")"),
       x = "Cluster A",
       y = "Cluster B") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


res_list_sub <- res_list[res_list != res_interest]
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
  ggplot(df_sym, aes(x = celltypeA, y = celltypeB, fill = alpha, col = pval <= alpha)) +
    coord_fixed() +
    geom_tile(linewidth = 0.5) +
    scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red", limits = c(-10,10)) +
    scale_color_manual(name = paste0("p value <= ", alpha), values = c("grey", "black")) +
    labs(title = paste0("Pair-wise cell type colocalization (Resolution = ", res, ")"),
         x = "Cluster A",
         y = "Cluster B") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_heatmap_alpha_with_sym_clustering_order_", res_interest, "_resolution_", res, ".pdf")), width = 12, height = 10, dpi = 300)
}

## Figure 3c (heatmap 100 um)
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
ggplot(df_sym, aes(x = celltypeA, y = celltypeB, fill = alpha, col = pval <= alpha)) +
  coord_fixed() +
  geom_tile(linewidth = 0.5) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red", limits = c(-10,10)) +
  scale_color_manual(name = paste0("p value <= ", alpha), values = c("grey", "black")) +
  labs(x = "Cluster A",
       y = "Cluster B") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0),
        axis.text.y = element_text(angle = 45, vjust = 0, hjust=1))
ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_heatmap_alpha_with_sym_clustering_resolution_", res, "_pyramid.pdf")), width = 12, height = 10, dpi = 300)


## Figure (visual inspection of spatial niches)
niches <- list(
  factor(c("Fol B cells", "CD4 Memory T cells", "Podoplanin", "indistinct"), levels = levels(ct_labels)),
  factor(c("Blood endothelial", "Macrophages", "CD8 Memory T cells", "Ki67 proliferating", "Myeloid cells", "B cells, red pulp", "Sinusoidal cells", "Neutrophils/Monocytes"), levels = levels(ct_labels))
)

for (i in seq_along(niches)) {
  ## subset by clusters
  spe_sub <- spe[,spe$celltypes %in% niches[[i]]]
  
  ## plot
  df <- data.frame(spatialCoords(spe))
  df_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], celltype = colData(spe_sub)$celltypes)
  ggplot(df, aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 0.2, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = celltype), size = 0.2, stroke = 0), dpi = 300) +
    scale_color_manual(values = col_clu[niches[[i]]]) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
    labs(x = "x (um)",
         y = "y (um)",
         col = "") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(filename = here("plots", dataset_name, method, paste0("singlecell_niche_", i, ".pdf")), width = 6, height = 5, dpi = 300)
}
