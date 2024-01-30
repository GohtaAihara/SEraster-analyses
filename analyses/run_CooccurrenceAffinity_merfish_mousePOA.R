## This file integrates rasterization with CooccurrenceAffinity to analyze MERFISH mouse preoptic area (mousePOA) datasets. Produce figures for 

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
library(DescTools)

par(mfrow=c(1,1))

dataset_name <- "merfish_mousePOA"
method <- "CooccurrenceAffinity"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_preprocessed.RDS")))

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = celltype)) +
  geom_point(size = 0.5, stroke = 0) +
  theme_classic()

ct_labels <- as.factor(colData(spe)$celltype)

# specify
animal <- 1
sex <- "Female"
behavior <- "Naive"
bregma <- "-0.29"


plot(spatialCoords(spe), pch=".", asp=1)

ct_labels <- as.factor(colData(spe)$celltype)

# Run method --------------------------------------------------------------

## modify code from "run_CooccurrenceAffinity_codex_humanSpleen.R"

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
  # col_name error
  spe_rast <- SEraster::rasterizeCellType(spe, "celltype", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
  
  ##spe_rast <- SEraster::rasterizeCellType(spe, col_name = spe$celltype, resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
  
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
  ## saveRDS(spe_rast, file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))
  
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
    
    ## out <- CooccurrenceAffinity::ML.Alpha(X,c(mA,mB,N), lev = CI_lev, pvalType = pval_method)
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

## runtime analysis
n_itr <- 5

runtime_results <- do.call(rbind, lapply(res_list, function(res) {
  out <- do.call(rbind, lapply(seq(n_itr), function(i) {
    print(paste0("Resolution: ", res, ", trial: ", i))
    start1 <- Sys.time()
    ## rasterize
    spe_rast <- SEraster::rasterizeCellType(spe, "celltype", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
    runtime_rast <- difftime(Sys.time(), start1, units = "secs")
    
    start2 <- Sys.time()
    ## compute relative enrichment (RE) metric (the output is dense matrix)
    mat <- assay(spe_rast, "pixelval")
    mat_re <- do.call(rbind, lapply(rownames(spe_rast), function(ct_label) {
      ## relative enrichment = celltype observed / celltype expected = celltype observed / (celltype frequency * total # of cells in the pixel)
      mat[ct_label,] / (sum(mat[ct_label,]) / sum(mat) * colSums(mat))
    }))
    rownames(mat_re) <- rownames(mat)
    
    ## binarize (1 if RE >= 1, 0 if RE < 1)
    mat_bin <- ifelse(mat_re >= 1, 1, 0)
    
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
    end <- Sys.time()
    runtime_coloc <- difftime(end, start2, units = "secs")
    runtime_total <- difftime(end, start1, units = "secs")
    
    return(data.frame(trial = i, num_pixels = dim(spe_rast)[2], runtime_rast = runtime_rast, runtime_coloc = runtime_coloc, runtime_total = runtime_total))
    
  }))
  return(data.frame(dataset = dataset_name, resolution = res, out))
}))


saveRDS(runtime_results, file = here("outputs", paste0(dataset_name, "_CooccurrenceAffinity_runtime.RDS")))

# Plot --------------------------------------------------------------

col_clu <- gg_color_hue(length(levels(ct_labels)))
names(col_clu) <- levels(ct_labels)

#########
df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = celltype)) +
  geom_point(size = 0.5, stroke = 0) +
  theme_classic()


## Figure 3a (cell type single-cell resolution)
df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], celltype = colData(spe)$celltype)
ggplot(df, aes(x = x, y = y, col = celltype)) +
  rasterise(geom_point(size = 0.5, stroke = 0), dpi = 300) +
  # scale_color_manual(name = "Cell type", values = col_clu) +
  # guides(col = guide_legend(override.aes = list(size = 3))) +
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


## modify code from "run_CooccurrenceAffinity_codex_humanSpleen.R" (hiearchical clustering of heatmap)


# Questions ---------------------------------------------------------------

## can we identify cell-type cooccurrence patterns in mPOA? if so, can they be validated by visualizing them at single-cell resolution?
## are cell-type cooccurrence patterns consistent across mice for each bregma?
