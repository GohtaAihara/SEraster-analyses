## This file integrates rasterization with CooccurrenceAffinity to analyze CODEX human intestine datasets. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::document()
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(BiocParallel)
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
library(ComplexHeatmap)

par(mfrow=c(1,1))

dataset_name = "codex_humanIntestine"
method = "CooccurrenceAffinity"

ct_labels <- factor(c("NK", "Enterocyte", "MUC1+ Enterocyte", "TA", "CD66+ Enterocyte", "Paneth", "Smooth muscle", "M1 Macrophage", "Goblet", "Neuroendocrine", "CD57+ Enterocyte", "Lymphatic", "CD8+ T", "DC", "M2 Macrophage", "B", "Neutrophil", "Endothelial", "Cycling TA", "Plasma", "CD4+ T cell", "Stroma", "Nerve", "ICC", "CD7+ Immune"))
donors <- c("B004", "B005", "B006", "B008", "B009", "B010", "B011", "B012")
tissue_locations <- c("Transverse", "Proximal Jejunum", "Duodenum", "Ascending", "Ileum", "Mid-jejunum", "Descending", "Descending - Sigmoid")
# celltype <- factor(c("NK", "Enterocyte", "MUC1+ Enterocyte", "TA", "CD66+ Enterocyte", "Paneth", "Smooth muscle", "M1 Macrophage", "Goblet", "Neuroendocrine", "CD57+ Enterocyte", "Lymphatic", "CD8+ T", "DC", "M2 Macrophage", "B", "Neutrophil", "Endothelial", "Cycling TA", "Plasma", "CD4+ T cell", "Stroma", "Nerve", "ICC", "CD7+ Immune"))

# Load dataset ------------------------------------------------------------
# 
# donors <- c("B009", "B010", "B011", "B012")
# tissue_locations <- c("Transverse", "Proximal Jejunum", "Duodenum", "Ascending", "Ileum", "Mid-jejunum", "Descending", "Descending - Sigmoid")
# 
# donor <- donors[[1]]
# tissue_location <- tissue_locations[[1]]

## "CL" dataset in Fig. 3 of the original paper
donor <- "B006"
tissue_location <- "Ascending"
spe <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_loc_", tissue_location, "_preprocessed.RDS")))

# ## "SB" dataset in Fig. 3 of the original paper
# donor <- "B005"
# tissue_location <- "Proximal Jejunum"
# spe <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_loc_", tissue_location, "_preprocessed.RDS")))

# Run methods -------------------------------------------------------------

## set resolution parameters
res_list <- c(10, 25, 50, 75, 100, 150, 200, 300, 400)
# res_list <- c(75, 150, 300)

## set CooccurrenceAffinity parameters
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
  spe_rast <- SEraster::rasterizeCellType(spe, "Cell.Type", resolution = res, fun = "sum", BPPARAM = MulticoreParam())
  
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
  saveRDS(spe_rast, file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_rasterized_resolution_", res, "_binarized.RDS")))
  
  ## compute affinity MLE using CooccurrenceAffinity (store the log-affinity metric, CI, and pvalue for each pair)
  ## get pair combinations
  non_self <- combn(levels(ct_labels), 2, simplify = FALSE)
  self <- lapply(levels(ct_labels), function(ct_label) c(ct_label, ct_label))
  pairs <- c(non_self, self)
  
  ## multiply pixel values for each pair of cell types
  affinity_results <- do.call(rbind, lapply(pairs, function(pair) {
    ## check if cell type labels in pair exist in the dataset
    if (pair[1] %in% rownames(mat_bin) & pair[2] %in% rownames(mat_bin)) {
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
        alpha = out$est,
        ci.min = out[CI_idx][[1]][1], 
        ci.max = out[CI_idx][[1]][2], 
        pval = out$pval)
      )
    }
  }))
  ## save results
  saveRDS(affinity_results, file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
}


# Plot --------------------------------------------------------------------

## Figure 3a (single cell visualizations)
col_clu <- gg_color_hue(length(levels(ct_labels)))
names(col_clu) <- levels(ct_labels)

for (donor in donors) {
  for (tissue_location in tissue_locations) {
    spe <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_loc_", tissue_location, "_preprocessed.RDS")))
    df <- data.frame(spatialCoords(spe), celltype = colData(spe)$Cell.Type)
    ggplot(df, aes(x = x, y = y, col = celltype)) +
      coord_fixed() +
      rasterise(geom_point(size = 0.5, stroke = 0), dpi = 300) +
      scale_color_manual(name = "Cell type", values = col_clu) +
      guides(col = guide_legend(override.aes = list(size = 3))) +
      labs(title = paste0("Donor: ", donor, ", Tissue: ", tissue_location)) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
      )
    ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_sc_ct_donor_", donor, "_tissue_", tissue_location, ".pdf")))
  }
}

## Figure 3c (rasterized visualizations)
## "CL" dataset in Fig. 3 of the original paper
donor <- "B006"
tissue_location <- "Ascending"
resolution <- 400
spe_rast <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_rasterized_resolution_", resolution, "_binarized.RDS")))
ct_interest <- "Enterocyte"
plotRaster(spe_rast, assay_name = "pixelval", feature_name = ct_interest, showLegend = TRUE , name = "cells/pixel", option = "inferno", breaks = scales::pretty_breaks()) +
  theme(legend.position = "bottom")
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_rasterized_ct_donor_", donor, "_tissue_", tissue_location, "_resolution_", resolution, "_pixelval_ct_", ct_interest, ".pdf")))
plotRaster(spe_rast, assay_name = "re", feature_name = ct_interest, showLegend = TRUE , name = "RE", option = "inferno") +
  theme(legend.position = "bottom")
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_rasterized_ct_donor_", donor, "_tissue_", tissue_location, "_resolution_", resolution, "_re_ct_", ct_interest, ".pdf")))
plotRaster(spe_rast, assay_name = "bin", feature_name = ct_interest, factor_levels = c(0,1), showLegend = TRUE , name = "Binarized", option = "inferno") +
  theme(legend.position = "bottom")
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_rasterized_ct_donor_", donor, "_tissue_", tissue_location, "_resolution_", resolution, "_bin_ct_", ct_interest, ".pdf")))


## Figure 3d-x
df_res_k <- data.frame(resolution = c(50, 300, 400), k = c(9, 6, 6))
# hclust_methods <- c("complete", "ward.D2")
hclust_methods <- c("complete")

for (i in seq(nrow(df_res_k))) {
  res <- df_res_k[i,]$resolution
  k <- df_res_k[i,]$k
  print(res)
  print(k)
  
  df <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
  
  ## create symmetric data
  df_flipped <- df[df$celltypeA != df$celltypeB,]
  df_flipped[,c("celltypeA", "celltypeB")] <- df_flipped[,c("celltypeB", "celltypeA")]
  df_sym <- rbind(df, df_flipped)
  
  ## reset label order
  df_sym <- df_sym %>%
    mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
           celltypeB = factor(celltypeB, levels(ct_labels)),
           significance = case_when(pval <= 0.05 ~ "*", pval > 0.05 ~ ""))
  
  for (hclust_method in hclust_methods) {
    df_plt <- as.matrix(cast(df_sym, celltypeA ~ celltypeB, value = "alpha"))
    df_plt_sig <- as.matrix(cast(df_sym, celltypeA ~ celltypeB, value = "significance"))
    
    # cluster with hclust
    hc_sym <- hclust(dist(df_plt), method = hclust_method)
    # plot dendrogram
    plot(hc_sym, labels = rownames(df_plt))
    rect.hclust(hc_sym, k = k, border = "red")
    # get niches
    niches <- cutree(hc_sym, k = k)
    names(niches) <- rownames(df_plt)
    
    # pyramid shape (without split and dendrogram)
    set.seed(0)
    heatmap <- ComplexHeatmap::Heatmap(
      df_plt,
      row_title = "Cell-type A",
      column_title = "Cell-type B",
      clustering_distance_rows = "euclidean",
      clustering_method_rows = hclust_method,
      clustering_distance_columns = "euclidean",
      clustering_method_columns = hclust_method,
      row_names_side = "left",
      column_names_side = "top",
      heatmap_legend_param = list(title = "Alpha MLE"),
      column_split = k,
      row_order = hc_sym$order,
      show_row_dend = FALSE,
      column_gap = unit(2.5, "mm"),
      bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = gg_color_hue(k)),
                                                             labels = unique(niches[hc_sym$order]),
                                                             labels_gp = gpar(col = "black", fontsize = 12))),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(df_plt_sig[i, j], x = x, y = y, just = "center", gp = grid::gpar(col = "black", cex = 2.5))
      }
    )
    tidyHeatmap::save_pdf(heatmap, filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_heatmap_resolution_", res, "hclust_method_", hclust_method, "_pyramid.pdf")), width = 12, height = 12)
    
    
    # square shape (with split and dendrogram)
    set.seed(0)
    heatmap <- ComplexHeatmap::Heatmap(
      df_plt,
      row_title = "Cell-type A",
      column_title = "Cell-type B",
      clustering_distance_rows = "euclidean",
      clustering_method_rows = hclust_method,
      clustering_distance_columns = "euclidean",
      clustering_method_columns = hclust_method,
      row_names_side = "left",
      column_names_side = "top",
      heatmap_legend_param = list(title = "Alpha MLE"),
      column_split = k,
      row_order = hc_sym$order,
      show_row_dend = FALSE,
      column_gap = unit(2.5, "mm"),
      bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = gg_color_hue(k)),
                                                             labels = unique(niches[hc_sym$order]),
                                                             labels_gp = gpar(col = "black", fontsize = 12))),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(df_plt_sig[i, j], x = x, y = y, just = "center", gp = grid::gpar(col = "black", cex = 2.5))
      }
    )
    tidyHeatmap::save_pdf(heatmap, filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_heatmap_resolution_", res, "hclust_method_", hclust_method, ".pdf")), width = 12, height = 12)
    
    # plot single-cell visualizations of cell-types in each niche
    df <- data.frame(spatialCoords(spe))
    for (i in unique(niches)) {
      ## subset by clusters
      spe_sub <- spe[,spe$Cell.Type %in% names(niches[niches == i])]

      ## plot
      df_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], celltype = colData(spe_sub)$Cell.Type)
      ggplot(df, aes(x = x, y = y)) +
        coord_fixed() +
        rasterise(geom_point(color = "lightgray", size = 0.5, stroke = 0), dpi = 300) +
        rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = celltype), size = 0.5, stroke = 0), dpi = 300) +
        geom_rect(aes(xmin = 11700, xmax = 13800, ymin = 15600, ymax = 16400), color = "black", fill = NA, linetype = "solid") +
        scale_color_manual(values = col_clu[names(niches[niches == i])]) +
        guides(col = guide_legend(override.aes = list(size = 3))) +
        # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
        labs(title = paste0("Donor: ", donor, ", Tissue: ", tissue_location, "\nNiche ", i),
             x = "x (um)",
             y = "y (um)",
             col = "Cell type") +
        theme_bw() +
        theme(
          legend.position = "bottom",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
        )
      ggsave(filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_niche_", i, ".pdf")), dpi = 300)
    }
    
    # plot single-cell visualizations of niches
    df_niches <- data.frame(Cell.Type = names(niches), niche = niches)
    rownames(df_niches) <- NULL
    df <- data.frame(spatialCoords(spe), colData(spe)) %>%
      left_join(df_niches, by = "Cell.Type") %>%
      mutate(niche = factor(niche))
    
    ggplot(df, aes(x = x, y = y)) +
      coord_fixed() +
      rasterise(geom_point(color = "lightgray", size = 0.5, stroke = 0), dpi = 300) +
      rasterise(geom_point(data = df, aes(x = x, y = y, col = niche), size = 0.5, stroke = 0), dpi = 300) +
      guides(col = guide_legend(override.aes = list(size = 3))) +
      # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
      labs(title = paste0("Donor: ", donor, ", Tissue: ", tissue_location),
           x = "x (um)",
           y = "y (um)",
           col = "Niche") +
      theme_bw() +
      theme(
        # legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
      )
    ggsave(filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_niches_summary.pdf")), dpi = 300)
    
    # plot niches vs. neighborhood
    # compute counts/proportions of neighborhoods in the entire dataset
    count_neighborhood_in_all <- df %>%
      group_by(Neighborhood) %>%
      summarise(count_neighborhood_in_all = n())
    
    prop_neighborhood_in_all <- df %>%
      group_by(Neighborhood) %>%
      summarise(prop_neighborhood_in_all = n()/dim(df)[1])
    
    # compute counts/proportions of niches in the entire dataset
    count_cells_in_niche <- df %>%
      group_by(niche) %>%
      summarise(count_cells_in_niche = n())
    
    # compute total cells in each neighborhood for a given niche
    count_neighborhood_in_niche <- df %>%
      group_by(niche, Neighborhood) %>%
      summarise(count_neighborhood_in_niche = n())
    
    # compute fold change
    df_plt <- count_neighborhood_in_niche %>%
      left_join(count_cells_in_niche, by = "niche") %>%
      left_join(prop_neighborhood_in_all, by = "Neighborhood") %>%
      mutate(prop_neighborhood_in_niche = count_neighborhood_in_niche / count_cells_in_niche,
             fold_change = (count_neighborhood_in_niche / count_cells_in_niche) / prop_neighborhood_in_all) %>%
      mutate(log2fold_change = log2(fold_change)) %>%
      select(niche, Neighborhood, log2fold_change) %>%
      spread(Neighborhood, log2fold_change) %>%
      column_to_rownames(var = "niche")
    
    heatmap <- ComplexHeatmap::Heatmap(
      t(as.matrix(df_plt)),
      name = "log2(fold enrichment)",
      row_title = "Neighborhoods",
      row_names_side = "left",
      column_title = "Niches",
      column_title_side = "bottom",
      column_names_rot = 0,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      column_order = rownames(df_plt)
    )
    tidyHeatmap::save_pdf(heatmap, filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_niches_vs_neighborhoods_resolution_", res, "hclust_method_", hclust_method, ".pdf")))
    
    # plot niches vs. communities
    # compute counts/proportions of communities in the entire dataset
    count_community_in_all <- df %>%
      group_by(Community) %>%
      summarise(count_community_in_all = n())
    
    prop_community_in_all <- df %>%
      group_by(Community) %>%
      summarise(prop_community_in_all = n()/dim(df)[1])
    
    # compute counts/proportions of niches in the entire dataset
    count_cells_in_niche <- df %>%
      group_by(niche) %>%
      summarise(count_cells_in_niche = n())
    
    # compute total cells in each neighborhood for a given niche
    count_community_in_niche <- df %>%
      group_by(niche, Community) %>%
      summarise(count_community_in_niche = n())
    
    # compute fold change
    df_plt <- count_community_in_niche %>%
      left_join(count_cells_in_niche, by = "niche") %>%
      left_join(prop_community_in_all, by = "Community") %>%
      mutate(prop_community_in_niche = count_community_in_niche / count_cells_in_niche,
             fold_change = (count_community_in_niche / count_cells_in_niche) / prop_community_in_all) %>%
      mutate(log2fold_change = log2(fold_change)) %>%
      select(niche, Community, log2fold_change) %>%
      spread(Community, log2fold_change) %>%
      column_to_rownames(var = "niche")
    
    heatmap <- ComplexHeatmap::Heatmap(
      t(as.matrix(df_plt)),
      name = "log2(fold enrichment)",
      row_title = "Communities",
      row_names_side = "left",
      column_title = "Niches",
      column_title_side = "bottom",
      column_names_rot = 0,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      column_order = rownames(df_plt)
    )
    tidyHeatmap::save_pdf(heatmap, filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_niches_vs_communities_resolution_", res, "hclust_method_", hclust_method, ".pdf")))
    
    # plot niches vs. tissue units
    # compute counts/proportions of communities in the entire dataset
    count_tissueunit_in_all <- df %>%
      group_by(Tissue.Unit) %>%
      summarise(count_tissueunit_in_all = n())
    
    prop_tissueunit_in_all <- df %>%
      group_by(Tissue.Unit) %>%
      summarise(prop_tissueunit_in_all = n()/dim(df)[1])
    
    # compute counts/proportions of niches in the entire dataset
    count_cells_in_niche <- df %>%
      group_by(niche) %>%
      summarise(count_cells_in_niche = n())
    
    # compute total cells in each neighborhood for a given niche
    count_tissueunit_in_niche <- df %>%
      group_by(niche, Tissue.Unit) %>%
      summarise(count_tissueunit_in_niche = n())
    
    # compute fold change
    df_plt <- count_tissueunit_in_niche %>%
      left_join(count_cells_in_niche, by = "niche") %>%
      left_join(prop_tissueunit_in_all, by = "Tissue.Unit") %>%
      mutate(prop_tissueunit_in_niche = count_tissueunit_in_niche / count_cells_in_niche,
             fold_change = (count_tissueunit_in_niche / count_cells_in_niche) / prop_tissueunit_in_all) %>%
      mutate(log2fold_change = log2(fold_change)) %>%
      select(niche, Tissue.Unit, log2fold_change) %>%
      spread(Tissue.Unit, log2fold_change) %>%
      column_to_rownames(var = "niche")
    
    heatmap <- ComplexHeatmap::Heatmap(
      t(as.matrix(df_plt)),
      name = "log2(fold enrichment)",
      row_title = "Tissue Units",
      row_names_side = "left",
      column_title = "Niches",
      column_title_side = "bottom",
      column_names_rot = 0,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      column_order = rownames(df_plt)
    )
    tidyHeatmap::save_pdf(heatmap, filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_niches_vs_tissueunits_resolution_", res, "hclust_method_", hclust_method, ".pdf")))
  }
}

## Figure 3D-
# res <- 50
res <- 400
df <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_resolution_", res, ".RDS")))

# create symmetric data
df_flipped <- df[df$celltypeA != df$celltypeB,]
df_flipped[,c("celltypeA", "celltypeB")] <- df_flipped[,c("celltypeB", "celltypeA")]
df_sym <- rbind(df, df_flipped)

# reset label order
df_sym <- df_sym %>%
  mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
         celltypeB = factor(celltypeB, levels(ct_labels)),
         significance = case_when(pval <= 0.05 ~ "*", pval > 0.05 ~ ""))

df_plt <- as.matrix(cast(df_sym, celltypeA ~ celltypeB, value = "alpha"))
df_plt_sig <- as.matrix(cast(df_sym, celltypeA ~ celltypeB, value = "significance"))

# cluster with hclust
hc_sym <- hclust(dist(df_plt), method = "complete")

# heatmap, pyramid shape (without split and dendrogram)
set.seed(0)
heatmap <- ComplexHeatmap::Heatmap(
  df_plt,
  row_title = "Cell-type A",
  column_title = "Cell-type B",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  row_names_side = "left",
  column_names_side = "top",
  heatmap_legend_param = list(title = "Alpha MLE"),
  row_order = rev(hc_sym$order),
  column_order = hc_sym$order,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(df_plt_sig[i, j], x = x, y = y, just = "center", rot = 45, gp = grid::gpar(col = "black", cex = 2.5))
  },
  row_names_rot = 45,
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 2),
  show_heatmap_legend = FALSE
)
tidyHeatmap::save_pdf(heatmap, filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_heatmap_resolution_", res, "hclust_method_complete_pyramid.pdf")), width = 12, height = 12)

# plot single-cell visualizations of cell-types in each niche
# # resolution = 50
# niches <- list(
#   factor(c("CD4+ T cell", "DC", "CD8+ T", "Neutrophil"), levels = levels(ct_labels)),
#   factor(c("Enterocyte", "Neuroendocrine", "TA"), levels = levels(ct_labels))
# )
# resolution = 400
niches <- list(
  factor(c("Stroma", "Nerve", "Smooth muscle", "Endothelial", "M2 Macrophage"), levels = levels(ct_labels)),
  factor(c("CD7+ Immune", "CD66+ Enterocyte", "B", "Plasma", "NK", "ICC", "DC", "Lymphatic", "CD8+ T", "CD4+ T cell", "M1 Macrophage", "CD57+ Enterocyte", "Enterocyte", "Cycling TA", "TA", "Neuroendocrine", "Neutrophil", "Goblet", "MUC1+ Enterocyte"), levels = levels(ct_labels)),
  factor(c("CD8+ T", "CD4+ T cell", "M1 Macrophage", "CD57+ Enterocyte", "Enterocyte", "Cycling TA", "TA", "Neuroendocrine", "Neutrophil", "Goblet", "MUC1+ Enterocyte"), levels = levels(ct_labels))
)

for (i in seq_along(niches)) {
  df <- data.frame(spatialCoords(spe))
  ## subset by clusters
  spe_sub <- spe[,spe$Cell.Type %in% niches[[i]]]
  
  ## plot whole tissue
  df_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], celltype = colData(spe_sub)$Cell.Type)
  ggplot(df, aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 1, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = celltype), size = 1, stroke = 0), dpi = 300) +
    geom_rect(aes(xmin = 11700, xmax = 13800, ymin = 15600, ymax = 16400), color = "black", fill = NA, linetype = "solid") +
    scale_color_manual(values = col_clu[niches[[i]]]) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
    labs(title = paste0("Donor: ", donor, ", Tissue: ", tissue_location),
         x = "x (um)",
         y = "y (um)",
         col = "Cell type") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  # ggsave(filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_celltypes_", paste(as.character(niches[[i]]), collapse = "_"), ".pdf")), dpi = 300)
  ggsave(filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_niche_", i, ".pdf")), dpi = 300)
  
  ## plot cloes-up
  df <- df %>%
    filter(x >= 11700 & x <= 13800, y >= 15600 & y <= 16400)
  df_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], celltype = colData(spe_sub)$Cell.Type) %>%
    filter(x >= 11700 & x <= 13800, y >= 15600 & y <= 16400)
  ggplot(df, aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 4, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = celltype), size = 4, stroke = 0), dpi = 300) +
    # geom_rect(aes(xmin = 11700, xmax = 13800, ymin = 15500, ymax = 16500), color = "black", fill = NA, linetype = "solid") +
    scale_color_manual(values = col_clu[niches[[i]]]) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
    labs(title = paste0("Donor: ", donor, ", Tissue: ", tissue_location),
         x = "x (um)",
         y = "y (um)",
         col = "Cell type") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  # ggsave(filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_celltypes_", paste(as.character(niches[[i]]), collapse = "_"), "_closeup.pdf")), dpi = 300)
  ggsave(filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_niche_", i, "_closeup.pdf")), dpi = 300)
}

## Supplementary Figure 3




## Figure x (plot close-up single cell visualizations)
# res <- 50
# k <- 9
# niches_interest <- c(2,4,6)

res <- 400
k <- 6
niches_interest <- c(2)

df <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_resolution_", res, ".RDS")))

# create symmetric data
df_flipped <- df[df$celltypeA != df$celltypeB,]
df_flipped[,c("celltypeA", "celltypeB")] <- df_flipped[,c("celltypeB", "celltypeA")]
df_sym <- rbind(df, df_flipped)

# reset label order
df_sym <- df_sym %>%
  mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
         celltypeB = factor(celltypeB, levels(ct_labels)),
         significance = case_when(pval <= 0.05 ~ "*", pval > 0.05 ~ ""))

df_plt <- as.matrix(cast(df_sym, celltypeA ~ celltypeB, value = "alpha"))
df_plt_sig <- as.matrix(cast(df_sym, celltypeA ~ celltypeB, value = "significance"))

# cluster with hclust
hc_sym <- hclust(dist(df_plt), method = "complete")
# plot dendrogram
plot(hc_sym, labels = rownames(df_plt))
rect.hclust(hc_sym, k = k, border = "red")
# get niches
niches <- cutree(hc_sym, k = k)
names(niches) <- rownames(df_plt)

# plot single-cell visualizations of cell-types in each niche
df <- data.frame(spatialCoords(spe))
for (i in niches_interest) {
  ## subset by clusters
  spe_sub <- spe[,spe$Cell.Type %in% names(niches[niches == i])]
  
  ## plot
  df <- df %>%
    filter(x >= 11700 & x <= 13800, y >= 15600 & y <= 16400)
  df_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], celltype = colData(spe_sub)$Cell.Type) %>%
    filter(x >= 11700 & x <= 13800, y >= 15600 & y <= 16400)
  ggplot(df, aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 2.5, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = celltype), size = 2.5, stroke = 0), dpi = 300) +
    # geom_rect(aes(xmin = 11700, xmax = 13800, ymin = 15500, ymax = 16500), color = "black", fill = NA, linetype = "solid") +
    scale_color_manual(values = col_clu[names(niches[niches == i])]) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
    labs(title = paste0("Donor: ", donor, ", Tissue: ", tissue_location, "\nNiche ", i),
         x = "x (um)",
         y = "y (um)",
         col = "Cell type") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_niche_", i, "_closeup.pdf")), dpi = 300)
}



## Figure x (plot dendograms and compare complete vs. ward.D2 as hclust methods)
df_res_k <- data.frame(resolution = c(50, 300, 400), k = c(7, 6, 6))
hclust_methods <- c("complete", "ward.D2")

for (i in seq(nrow(df_res_k))) {
  res <- df_res_k[i,]$resolution
  k <- df_res_k[i,]$k
  
  df <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
  
  ## create symmetric data
  df_flipped <- df[df$celltypeA != df$celltypeB,]
  df_flipped[,c("celltypeA", "celltypeB")] <- df_flipped[,c("celltypeB", "celltypeA")]
  df_sym <- rbind(df, df_flipped)
  
  ## reset label order
  df_sym <- df_sym %>%
    mutate(celltypeA = factor(celltypeA, levels(ct_labels)),
           celltypeB = factor(celltypeB, levels(ct_labels)))
  
  niches_compare <- do.call(cbind, lapply(hclust_methods, function(hclust_method) {
    df_plt <- as.matrix(cast(df_sym, celltypeA ~ celltypeB, value = "alpha"))
    
    # cluster with hclust
    hc_sym <- hclust(dist(df_plt), method = hclust_method)
    
    # plot dendrogram
    pdf(file = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_dendogram_resolution_", res, "hclust_method_", hclust_method, ".pdf")))
    plot(hc_sym, labels = rownames(df_plt))
    dev.off()
    
    # get niches
    niches <- cutree(hc_sym, k = k)
    names(niches) <- rownames(df_plt)
    
    return(niches)
  }))
  colnames(niches_compare) <- hclust_methods
  
  print(paste0("Resolution: ", res, ", k = ", k))
  print(niches_compare)
}


## Figure 3c (heatmaps)
res_list <- c(10, 25, 75, 50, 100, 150, 200, 300, 400)
for (res in res_list) {
  df <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
  
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
  df_plt <- df_sym %>%
    mutate(
      alpha = Winsorize(alpha, min(lim), max(lim)),
      significance = case_when(pval <= 0.05 ~ "*")
    )
  ggplot(df_plt, aes(x = celltypeA, y = celltypeB, fill = alpha, label = significance)) +
    coord_fixed() +
    geom_tile(color = "gray") +
    geom_text(size = 10, angle = 45) +
    scale_x_discrete(position = "top") +
    scale_fill_gradient2(name = "Alpha MLE", low = "blue", mid = "white", high = "red", limits = lim) +
    labs(x = "Cluster B",
         y = "Cluster A") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0),
          axis.text.y = element_text(angle = 45, vjust = 0, hjust=1))
  ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_heatmap_alpha_with_sym_clustering_resolution_", res, "_pyramid.pdf")), width = 12, height = 10, dpi = 300)
}

## Figure 3x
# identify niches manually
niches <- list(
  factor(c("NK", "CD57+ Enterocyte", "Paneth"), levels = levels(ct_labels)),
  factor(c("Enterocyte", "Neuroendocrine", "TA"), levels = levels(ct_labels)),
  factor(c("B", "M1 Macrophage", "Plasma"), levels = levels(ct_labels)),
  factor(c("CD4+ T cell", "DC", "CD8+ T", "Neutrophil"), levels = levels(ct_labels)),
  factor(c("Cycling TA", "MUC1+ Enterocyte", "CD66+ Enterocyte", "Goblet"), levels = levels(ct_labels))
)
names(niches) <- paste0("niche ", seq_along(niches))

# identify niches using dendrogram
res <- 300
df <- readRDS(file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_CooccurrenceAffinity_resolution_", res, ".RDS")))
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
# plot dendrogram
plot(hc_sym_interest, labels = rownames(df_heatmap_sym))
k_select <- 5
rect.hclust(hc_sym_interest, k = k_select, border = "red")
# get niches
niches <- cutree(hc_sym_interest, k = k_select)
names(niches) <- rownames(df_heatmap_sym)

# for (i in seq_along(niches)) {
#   ## subset by clusters
#   spe_sub <- spe[,spe$Cell.Type %in% niches[[i]]]
#   
#   ## plot
#   df <- data.frame(spatialCoords(spe))
#   df_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], celltype = colData(spe_sub)$Cell.Type)
#   ggplot(df, aes(x = x, y = y)) +
#     coord_fixed() +
#     rasterise(geom_point(color = "lightgray", size = 0.5, stroke = 0), dpi = 300) +
#     rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = celltype), size = 0.5, stroke = 0), dpi = 300) +
#     scale_color_manual(values = col_clu[niches[[i]]]) +
#     guides(col = guide_legend(override.aes = list(size = 3))) +
#     # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
#     labs(title = paste0("Donor: ", donor, ", Tissue: ", tissue_location, "\nNiche ", i),
#          x = "x (um)",
#          y = "y (um)",
#          col = "Cell type") +
#     theme_bw() +
#     theme(
#       # legend.position = "bottom",
#       panel.grid = element_blank(),
#       axis.title = element_blank(),
#       axis.text = element_blank(),
#       axis.ticks = element_blank(),
#     )
#   ggsave(filename = here("plots", dataset_name, method, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_niche_", i, ".pdf")), dpi = 300)
# }

for (i in unique(niches)) {
  ## subset by clusters
  spe_sub <- spe[,spe$Cell.Type %in% names(niches[niches == i])]

  ## plot
  df <- data.frame(spatialCoords(spe))
  df_sub <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], celltype = colData(spe_sub)$Cell.Type)
  ggplot(df, aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 0.5, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = celltype), size = 0.5, stroke = 0), dpi = 300) +
    scale_color_manual(values = col_clu[names(niches[niches == i])]) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
    labs(title = paste0("Donor: ", donor, ", Tissue: ", tissue_location, "\nNiche ", i),
         x = "x (um)",
         y = "y (um)",
         col = "Cell type") +
    theme_bw() +
    theme(
      # legend.position = "bottom",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_niche_", i, ".pdf")), dpi = 300)
}

## plot summary of niches
# make it adaptable later
# df <- data.frame(spatialCoords(spe), colData(spe)) %>%
#   mutate(niche = case_when(
#     Cell.Type %in% niches[["niche 1"]] ~ "niche 1",
#     Cell.Type %in% niches[["niche 2"]] ~ "niche 2",
#     Cell.Type %in% niches[["niche 3"]] ~ "niche 3",
#     Cell.Type %in% niches[["niche 4"]] ~ "niche 4",
#     Cell.Type %in% niches[["niche 5"]] ~ "niche 5",
#     TRUE ~ "unknown" # Default case for cell types not found in any niche
#   )) %>%
#   mutate(niche = factor(niche, levels = c(names(niches), "unknown")))

df_niches <- data.frame(Cell.Type = names(niches), niche = niches)
rownames(df_niches) <- NULL
df <- data.frame(spatialCoords(spe), colData(spe)) %>%
  left_join(df_niches, by = "Cell.Type") %>%
  mutate(niche = factor(niche))

ggplot(df, aes(x = x, y = y)) +
  coord_fixed() +
  rasterise(geom_point(color = "lightgray", size = 0.5, stroke = 0), dpi = 300) +
  rasterise(geom_point(data = df, aes(x = x, y = y, col = niche), size = 0.5, stroke = 0), dpi = 300) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  # scale_color_manual(name = "Clusters", values = gg_color_hue(67)) +
  labs(title = paste0("Donor: ", donor, ", Tissue: ", tissue_location),
       x = "x (um)",
       y = "y (um)",
       col = "Niche") +
  theme_bw() +
  theme(
    # legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
ggsave(filename = here("plots", dataset_name, method, res, paste0(dataset_name, "_donor_", donor, "_tissue_location_", tissue_location, "_singlecell_niches_summary.pdf")), dpi = 300)

## compare niches vs. neighborhoods in the original paper
# compute counts/proportions of neighborhoods in the entire dataset
count_neighborhood_in_all <- df %>%
  group_by(Neighborhood) %>%
  summarise(count_neighborhood_in_all = n())

prop_neighborhood_in_all <- df %>%
  group_by(Neighborhood) %>%
  summarise(prop_neighborhood_in_all = n()/dim(df)[1])

# compute counts/proportions of niches in the entire dataset
count_cells_in_niche <- df %>%
  group_by(niche) %>%
  summarise(count_cells_in_niche = n())

# compute total cells in each neighborhood for a given niche
count_neighborhood_in_niche <- df %>%
  group_by(niche, Neighborhood) %>%
  summarise(count_neighborhood_in_niche = n())

# compute fold change
df_plt <- count_neighborhood_in_niche %>%
  left_join(count_cells_in_niche, by = "niche") %>%
  left_join(prop_neighborhood_in_all, by = "Neighborhood") %>%
  mutate(prop_neighborhood_in_niche = count_neighborhood_in_niche / count_cells_in_niche,
         fold_change = (count_neighborhood_in_niche / count_cells_in_niche) / prop_neighborhood_in_all)

# plot with regular ggplot
set.seed(0)
ggplot(df_plt, aes(x = log2(fold_change), y = niche, col = Neighborhood, shape = niche)) +
  geom_violin(color = "black") +
  geom_boxplot(color = "black", width = 0.1) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(x = "log2(fold enrichment)",
       y = "Niche",
       col = "Neighborhood",
       shape = "Niche") +
  theme_bw()
ggsave(file = here("plots", dataset_name, method, res, paste0(dataset_name, "_niche_neighborhood_fold_change_resolution_50.pdf")), width = 8, height = 8, dpi = 300)

# ggplot(df_plt, aes(x = Neighborhood, y = niche, fill = log2(df_plt))) +
#   geom_tile() +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   labs(x = "Neighborhoods",
#        y = "Niches") +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )

# plot with ComplexHeatmap
df_plt2 <- df_plt %>%
  mutate(log2fold_change = log2(fold_change)) %>%
  select(niche, Neighborhood, log2fold_change) %>%
  spread(Neighborhood, log2fold_change) %>%
  column_to_rownames(var = "niche")

ComplexHeatmap::Heatmap(
  as.matrix(df_plt2),
  name = "log2(fold enrichment)",
  row_title = "Niches",
  column_title = "Neighborhoods"
)

## compare niches vs. communities in the original paper
# compute counts/proportions of communities in the entire dataset
count_community_in_all <- df %>%
  group_by(Community) %>%
  summarise(count_community_in_all = n())

prop_community_in_all <- df %>%
  group_by(Community) %>%
  summarise(prop_community_in_all = n()/dim(df)[1])

# compute counts/proportions of niches in the entire dataset
count_cells_in_niche <- df %>%
  group_by(niche) %>%
  summarise(count_cells_in_niche = n())

# compute total cells in each neighborhood for a given niche
count_community_in_niche <- df %>%
  group_by(niche, Community) %>%
  summarise(count_community_in_niche = n())

# compute fold change
df_plt <- count_community_in_niche %>%
  left_join(count_cells_in_niche, by = "niche") %>%
  left_join(prop_community_in_all, by = "Community") %>%
  mutate(prop_community_in_niche = count_community_in_niche / count_cells_in_niche,
         fold_change = (count_community_in_niche / count_cells_in_niche) / prop_community_in_all)

# plot with regular ggplot
set.seed(0)
ggplot(df_plt, aes(x = log2(fold_change), y = niche, col = Community, shape = niche)) +
  geom_violin(color = "black") +
  geom_boxplot(color = "black", width = 0.1) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(x = "log2(fold enrichment)",
       y = "Niche",
       col = "Community",
       shape = "Niche") +
  theme_bw()
ggsave(file = here("plots", dataset_name, method, res, paste0(dataset_name, "_niche_neighborhood_fold_change_resolution_50.pdf")), width = 8, height = 8, dpi = 300)

# plot with ComplexHeatmap
df_plt2 <- df_plt %>%
  mutate(log2fold_change = log2(fold_change)) %>%
  select(niche, Community, log2fold_change) %>%
  spread(Community, log2fold_change) %>%
  column_to_rownames(var = "niche")

ComplexHeatmap::Heatmap(
  as.matrix(df_plt2),
  name = "log2(fold enrichment)",
  row_title = "Niches",
  column_title = "Communities"
)

## compare niches vs. tissue units in the original paper
# compute counts/proportions of communities in the entire dataset
count_tissueunit_in_all <- df %>%
  group_by(Tissue.Unit) %>%
  summarise(count_tissueunit_in_all = n())

prop_tissueunit_in_all <- df %>%
  group_by(Tissue.Unit) %>%
  summarise(prop_tissueunit_in_all = n()/dim(df)[1])

# compute counts/proportions of niches in the entire dataset
count_cells_in_niche <- df %>%
  group_by(niche) %>%
  summarise(count_cells_in_niche = n())

# compute total cells in each neighborhood for a given niche
count_tissueunit_in_niche <- df %>%
  group_by(niche, Tissue.Unit) %>%
  summarise(count_tissueunit_in_niche = n())

# compute fold change
df_plt <- count_tissueunit_in_niche %>%
  left_join(count_cells_in_niche, by = "niche") %>%
  left_join(prop_tissueunit_in_all, by = "Tissue.Unit") %>%
  mutate(prop_tissueunit_in_niche = count_tissueunit_in_niche / count_cells_in_niche,
         fold_change = (count_tissueunit_in_niche / count_cells_in_niche) / prop_tissueunit_in_all)

# plot with regular ggplot
set.seed(0)
ggplot(df_plt, aes(x = log2(fold_change), y = niche, col = Tissue.Unit, shape = niche)) +
  geom_violin(color = "black") +
  geom_boxplot(color = "black", width = 0.1) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(x = "log2(fold enrichment)",
       y = "Niche",
       col = "Tissue Unit",
       shape = "Niche") +
  theme_bw()
ggsave(file = here("plots", dataset_name, method, paste0(dataset_name, "_niche_neighborhood_fold_change_resolution_50.pdf")), width = 8, height = 8, dpi = 300)

# plot with ComplexHeatmap
df_plt2 <- df_plt %>%
  mutate(log2fold_change = log2(fold_change)) %>%
  select(niche, Tissue.Unit, log2fold_change) %>%
  spread(Tissue.Unit, log2fold_change) %>%
  column_to_rownames(var = "niche")

ComplexHeatmap::Heatmap(
  as.matrix(df_plt2),
  name = "log2(fold enrichment)",
  row_title = "Niches",
  column_title = "Tissue Unit"
)