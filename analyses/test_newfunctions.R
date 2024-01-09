setwd("~/Desktop/SEraster")
devtools::document()
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

library(here)

dataset_name = "codex_humanSpleen"
method = "CooccurrenceAffinity"

res <- 100
spe_rast <- readRDS(file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))
ct_label <- "Fol B cells"

plotRaster(spe_rast, assay_name = "pixelval", feature_name = "Fol B cells", name = "cell-type counts", option = "inferno")
plotRaster(spe_rast, assay_name = "re", feature_name = "Fol B cells", name = "RE", option = "inferno")
plotRaster(spe_rast, assay_name = "bin", feature_name = "Fol B cells", factor_levels = c(0,1), name = "binarized", option = "inferno")

data("merfish_mousePOA")
class(merfish_mousePOA)

## for tutorials
rastCt <- SEraster::rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = 50)

# extract cell-type labels
ct_labels <- as.factor(colData(merfish_mousePOA)$celltype)

# compute relative enrichment (RE) metric
mat <- assay(rastCt, "pixelval")
mat_re <- do.call(rbind, lapply(rownames(rastCt), function(ct_label) {
  mat[ct_label,] / (sum(mat[ct_label,]) / sum(mat) * colSums(mat))
}))
rownames(mat_re) <- rownames(mat)

# binarize
mat_bin <- ifelse(mat_re >= 1, 1, 0)

# add RE and binarized layers to SpatialExperiment object
assays(rastCt) <- list(pixelval = assay(rastCt, "pixelval"), re = mat_re, bin = mat_bin)

ct_interest <- "Ependymal"

plotRaster(rastCt, assay_name = "pixelval", feature_name = ct_interest, name = "cell-type counts", option = "inferno")
plotRaster(rastCt, assay_name = "re", feature_name = ct_interest, name = "RE", option = "inferno")
plotRaster(rastCt, assay_name = "bin", feature_name = ct_interest, factor_levels = c(0,1), name = "binarized", option = "inferno")

ct_coocc <- CooccurrenceAffinity::affinity(data = mat_bin, row.or.col = "row", squarematrix = c("all"))
CooccurrenceAffinity::plotgg(data = ct_coocc, variable = "alpha_mle", legendlimit = "datarange")

# visual inspection
plotRaster(rastCt, assay_name = "bin", feature_name = "Inhibitory", factor_levels = c(0,1), name = "binarized", option = "inferno")
plotRaster(rastCt, assay_name = "bin", feature_name = "Excitatory", factor_levels = c(0,1), name = "binarized", option = "inferno")
