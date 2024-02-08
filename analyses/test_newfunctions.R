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

## load example dataset
data("merfish_mousePOA")
dim(merfish_mousePOA)

## rasterize gene expression
rastGexp <- SEraster::rasterizeGeneExpression(merfish_mousePOA, assay_name="volnorm", resolution = 50)

resolution <- 100
data <- assay(merfish_mousePOA)
pos <- spatialCoords(merfish_mousePOA)
bbox <- sf::st_bbox(c(
  xmin = floor(min(pos[,1])-resolution/2), 
  xmax = ceiling(max(pos[,1])+resolution/2), 
  ymin = floor(min(pos[,2])-resolution/2), 
  ymax = ceiling(max(pos[,2])+resolution/2)
))
test <- SEraster::rasterizeMatrix(data, pos, bbox, resolution = resolution)
data_sub <- data[,1:100]
test <- SEraster::rasterizeMatrix(data_sub, pos, bbox, resolution = resolution)


# Test hexagonal pixels ---------------------------------------------------

data("merfish_mousePOA")
dim(merfish_mousePOA)

resolution = 50
pos <- spatialCoords(merfish_mousePOA)
bbox <- sf::st_bbox(c(
  xmin = floor(min(pos[,1])-resolution/2), 
  xmax = ceiling(max(pos[,1])+resolution/2), 
  ymin = floor(min(pos[,2])-resolution/2), 
  ymax = ceiling(max(pos[,2])+resolution/2)
))
grid_sq <- sf::st_make_grid(bbox, cellsize = resolution)
grid_hx <- sf::st_make_grid(bbox, cellsize = resolution, square = FALSE)

## rasterize gene expression
# square
rastGexp_sq <- SEraster::rasterizeGeneExpression(merfish_mousePOA, assay_name="volnorm", resolution = resolution)
dim(rastGexp_sq)
df_sq <- data.frame(spatialCoords(rastGexp_sq), tot_gexp = colSums(assay(rastGexp_sq)))
ggplot() +
  coord_fixed() +
  geom_sf(data = grid_sq) +
  geom_point(data = df_sq, aes(x = x, y = y, col = tot_gexp)) +
  scale_color_viridis_c() +
  theme_bw()

# hexagon
rastGexp_hx <- SEraster::rasterizeGeneExpression(merfish_mousePOA, assay_name="volnorm", resolution = resolution, square = FALSE)
dim(rastGexp_hx)
df_hx <- data.frame(spatialCoords(rastGexp_hx), tot_gexp = colSums(assay(rastGexp_hx)))
ggplot() +
  coord_fixed() +
  geom_sf(data = grid_hx) +
  geom_point(data = df_hx, aes(x = x, y = y, col = tot_gexp)) +
  scale_color_viridis_c() +
  theme_bw()

## rasterize celltypes
# square
rastCt_sq <- SEraster::rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = resolution)
dim(rastCt_sq)
df_sq <- data.frame(spatialCoords(rastCt_sq), Inhibitory = assay(rastCt_sq)["Inhibitory",])
ggplot() +
  coord_fixed() +
  geom_sf(data = grid_sq) +
  geom_point(data = df_sq, aes(x = x, y = y, col = Inhibitory)) +
  scale_color_viridis_c(option = "inferno") +
  theme_bw()

# hexagon
rastCt_hx <- SEraster::rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = resolution, square = FALSE)
dim(rastCt_hx)
df_hx <- data.frame(spatialCoords(rastCt_hx), Inhibitory = assay(rastCt_hx)["Inhibitory",])
ggplot() +
  coord_fixed() +
  geom_sf(data = grid_hx) +
  geom_point(data = df_hx, aes(x = x, y = y, col = Inhibitory)) +
  scale_color_viridis_c(option = "inferno") +
  theme_bw()