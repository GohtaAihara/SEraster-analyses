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
library(gridExtra)
library(here)
library(CooccurrenceAffinity)

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

## Figure (single cell spatial visualization)
df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  coord_fixed() +
  geom_point(size = 0.1) +
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
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_clusters_singlecell.pdf")))

## Figure (rasterized spatial visualization)
res <- 100
spe_rast <- readRDS(file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))

## Figure (CooccurrenceAffinity heatmap)
res <- 100
df <- readRDS(file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
ggplot(df, aes(x = celltypeA, y = celltypeB, fill = alpha)) +
  geom_tile() +
  scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red") +
  labs(title = paste0("Pair-wise cell type colocalization (Resolution = ", res, ")"),
       x = "Cluster A",
       y = "Cluster B") +
  theme_bw()

