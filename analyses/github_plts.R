## This file produces necessary plots for the README.md file on https://github.com/JEFworks/SEraster

## Load packages
setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

# library(SEraster)
library(SpatialExperiment)
library(here)
library(ggplot2)
library(nnSVG)
library(CooccurrenceAffinity)

dir_plt <- here(file.path("plots", "SEraster_github"))

## load example dataset
data("merfish_mousePOA")
dim(merfish_mousePOA)
length(unique(colData(merfish_mousePOA)$celltype))

## visualize at single-cell resolution
## cell-type
df <- data.frame(spatialCoords(merfish_mousePOA), celltype = colData(merfish_mousePOA)$celltype, neurontype = colData(merfish_mousePOA)$neurontype)
ggplot(df, aes(x = x, y = y, col = celltype)) +
  geom_point(size = 1) +
  labs(x = "x (μm)",
       y = "y (μm)",
       col = "Cell-types") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(filename = here(dir_plt, "singlecell_celltypes.png"), dpi = 300)

## rasterize gene expression
rastGexp <- SEraster::rasterizeGeneExpression(merfish_mousePOA, assay_name="volnorm", resolution = 50)
dim(rastGexp)

# plot total rasterized gene expression
SEraster::plotRaster(rastGexp, name = "Total rasterized gene expression")
ggsave(filename = here(dir_plt, "rasterized_gexp_total.png"), dpi = 300)

# plot specific gene
SEraster::plotRaster(rastGexp, feature_name = "Esr1", name = "Esr1")
ggsave(filename = here(dir_plt, "rasterized_gexp_esr1.png"), dpi = 300)

## rasterize cell-type labels
rastCt <- SEraster::rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = 50)
dim(rastCt)

# plot total cell counts
SEraster::plotRaster(rastCt, name = "cell counts", option = "inferno")
ggsave(filename = here(dir_plt, "rasterized_ct_total.png"), dpi = 300)

# plot specific cell-type
SEraster::plotRaster(rastCt, feature_name = "Inhibitory", name = "Inhibitory neuron counts", option = "inferno")
ggsave(filename = here(dir_plt, "rasterized_ct_inhibitory.png"), dpi = 300)

## SVG analysis
# run nnSVG
set.seed(0)
rastGexp <- nnSVG(rastGexp, assay_name = "pixelval")

# number of significant SVGs
table(rowData(rastGexp)$padj <= 0.05)

# plot rasterized gene expression of top-ranked SVG
top_svg <- which(rowData(rastGexp)$rank == 1)
top_svg_name <- rownames(rowData(rastGexp))[top_svg]
SEraster::plotRaster(rastGexp, feature_name = top_svg_name, name = top_svg_name)
ggsave(filename = here(dir_plt, "rasterized_gexp_top_svg.png"), dpi = 300)

# subset data
ct_interest <- "Excitatory"
spe_sub <- merfish_mousePOA[,merfish_mousePOA$celltype == ct_interest]

# run SEraster
rastGexp_sub <- SEraster::rasterizeGeneExpression(spe_sub, assay_name="volnorm", resolution = 50)

# run nnSVG
set.seed(0)
rastGexp_sub <- nnSVG(rastGexp_sub, assay_name = "pixelval")

# number of significant SVGs
table(rowData(rastGexp_sub)$padj <= 0.05)

# plot rasterized gene expression of top-ranked SVG
top_svg <- which(rowData(rastGexp_sub)$rank == 1)
top_svg_name <- rownames(rowData(rastGexp_sub))[top_svg]
SEraster::plotRaster(rastGexp_sub, feature_name = top_svg_name, name = top_svg_name)
ggsave(filename = here(dir_plt, "rasterized_gexp_sub_top_svg.png"), dpi = 300)

## cell-type cooccurrence
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

# plot pixel value for a cell-type of interest
plotRaster(rastCt, assay_name = "pixelval", feature_name = ct_interest, name = "cell-type counts", option = "inferno")
ggsave(filename = here(dir_plt, "rasterized_ct_sub_total.png"), dpi = 300)

# plot RE value for a cell-type of interest
plotRaster(rastCt, assay_name = "re", feature_name = ct_interest, name = "RE", option = "inferno")
ggsave(filename = here(dir_plt, "rasterized_ct_sub_re.png"), dpi = 300)

# plot binarized value for a cell-type of interest
plotRaster(rastCt, assay_name = "bin", feature_name = ct_interest, factor_levels = c(0,1), name = "binarized", option = "inferno")
ggsave(filename = here(dir_plt, "rasterized_ct_sub_bin.png"), dpi = 300)

# run CooccurrenceAffinity
ct_coocc <- CooccurrenceAffinity::affinity(data = mat_bin, row.or.col = "row", squarematrix = c("all"))

# plot maximum likelihood estimates of affinity metric (alpha MLE)
CooccurrenceAffinity::plotgg(data = ct_coocc, variable = "alpha_mle", legendlimit = "datarange")
ggsave(filename = here(dir_plt, "coocc_heatmap.png"), dpi = 300)
