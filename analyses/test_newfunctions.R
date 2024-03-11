
# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::document()
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

library(here)
library(ggplot2)


# Test plotting function in general --------------------------------------------------

dataset_name = "codex_humanSpleen"
method = "CooccurrenceAffinity"

res <- 100
spe_rast <- readRDS(file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))
ct_label <- "Fol B cells"

plotRaster(spe_rast, assay_name = "pixelval", feature_name = "Fol B cells", name = "cell-type counts", option = "inferno")
plotRaster(spe_rast, assay_name = "re", feature_name = "Fol B cells", name = "RE", option = "inferno")
plotRaster(spe_rast, assay_name = "bin", feature_name = "Fol B cells", factor_levels = c(0,1), name = "binarized", option = "inferno")


# Test tutorial codes -----------------------------------------------------

data("merfish_mousePOA")
class(merfish_mousePOA)

## for tutorials
rastCt <- SEraster::rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = 50, square = FALSE)

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


# Test when bigger than dataset resolution --------------------------------

data("merfish_mousePOA")
class(merfish_mousePOA)

## get range
diff(range(spatialCoords(merfish_mousePOA)[,1]))
diff(range(spatialCoords(merfish_mousePOA)[,2]))

## for tutorials
rastCt <- SEraster::rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = 5000, verbose = TRUE)
plotRaster(rastCt, assay_name = "pixelval", feature_name = "Inhibitory", name = "cell-type counts", option = "inferno")


# Test memory on sfc_POLYGON ----------------------------------------------

## merfish mousePOA
# data("merfish_mousePOA")
# spe <- merfish_mousePOA
## merfish mouseBrain
# spe <- readRDS(file = here("outputs", "merfish_mouseBrain_preprocessed.RDS"))
## xenium mousePup
spe <- readRDS(file = here("outputs", "xenium_mousePup_preprocessed.RDS"))

# set parallelization backend
bpparam <- BiocParallel::MulticoreParam()

## run SEraster to rasterize
resolution <- 100
# square pixel
rastGexp_sq <- SEraster::rasterizeGeneExpression(spe, resolution = resolution, BPPARAM = bpparam)
dim(rastGexp_sq)
# hexagonal pixel
rastGexp_hx <- SEraster::rasterizeGeneExpression(spe, resolution = resolution, square = FALSE, BPPARAM = bpparam)
dim(rastGexp_hx)

## compute size of spe
print(object.size(rastGexp_sq), units = "Mb")
print(object.size(rastGexp_hx), units = "Mb")

## extract sfc_POLYGON
sq <- colData(rastGexp_sq)$geometry
hx <- colData(rastGexp_hx)$geometry

## compute size of sfc_POLYGON
print(object.size(sq), units = "Kb")
print(object.size(hx), units = "Kb")

## compute size of one entry of sfc_POLYGON
print(object.size(sq[1]), units = "Kb")
print(object.size(hx[1]), units = "Kb")

## compute size of one sfg
print(object.size(sq[[1]]), units = "Kb")
print(object.size(hx[[1]]), units = "Kb")

## compute size of vertices matrix array
print(object.size(sq[[1]][[1]]), units = "Kb")
print(object.size(hx[[1]][[1]]), units = "Kb")

## check the distribution
df_sq <- do.call(rbind, lapply(seq(length(sq)), function(i) {
  return = data.frame(size = object.size(sq[[i]]))
}))
rownames(df_sq) <- colnames(rastGexp_sq)
unique(df_sq$size)

df_hx <- do.call(rbind, lapply(seq(length(hx)), function(i) {
  return = data.frame(size = object.size(hx[[i]]))
}))
rownames(df_hx) <- colnames(rastGexp_hx)
unique(df_hx$size)


# Test rotation -----------------------------------------------------------

data("merfish_mousePOA")
spe1 <- merfish_mousePOA
spe2 <- readRDS(file = here("outputs", "merfish_mouseBrain_preprocessed.RDS"))

ggplot(data.frame(spatialCoords(spe1)), aes(x = x, y = y)) +
  coord_fixed() +
  geom_point(size = 0.5, stroke = 0.01, alpha = 0.5) +
  theme_bw()

ggplot(data.frame(spatialCoords(spe2)), aes(x = x, y = y)) +
  coord_fixed() +
  geom_point(size = 0.5, stroke = 0.01, alpha = 0.5) +
  theme_bw()

# input = single SpatialExperiment
spe_list <- permutateByRotation(spe2, n_perm = 3)
names(spe_list)

df_plt <- do.call(rbind, lapply(seq_along(spe_list), function(i) {
  spe <- spe_list[[i]]
  return(data.frame(spatialCoords(spe), perm = as.character(i)))
}))

ggplot() +
  coord_fixed() +
  geom_point(data = df_plt, aes(x = x, y = y, col = perm, shape = perm), size = 0.5, stroke = 0.01, alpha = 0.5) +
  theme_bw()

# input = a list of SpatialExperiment objects
spe_list_orig <- list(spe1, spe2)
names(spe_list_orig) <- c("merfish_mousePOA", "merfish_mouseBrain")
# 
# input <- spe_list_orig
# pos_comb <- do.call(rbind, lapply(seq_along(input), function(i) {
#   pos <- SpatialExperiment::spatialCoords(input[[i]])
#   if (!is.null(names(input))) {
#     dataset <- names(input)[i]
#   } else {
#     dataset <- i
#   }
#   return(data.frame(dataset = dataset, x = pos[,1], y = pos[,2]))
# }))
# ## find the midrange point across combined x,y coordinates
# midrange_pt <- rearrr::midrange(pos_comb, cols = c("x", "y"))
# message(paste0("Rotating all datasets around (x, y) = (", midrange_pt$x, ", ", midrange_pt$y, ")"))
# 
# pos_orig <- as.data.frame(spatialCoords(spe1))
# angle <- 90
# pos_rotated <- rearrr::rotate_2d(data = pos_orig, degrees = angle, x_col = "x", y_col = "y", origin = midrange_pt, overwrite = TRUE)

spe_list <- permutateByRotation(spe_list_orig, n_perm = 3)
names(spe_list)

df_plt <- do.call(rbind, lapply(seq_along(spe_list), function(i) {
  spe <- spe_list[[i]]
  name <- strsplit(names(spe_list)[i], "_rotated_", fixed = TRUE)
  return(data.frame(spatialCoords(spe), dataset = name[[1]][1], angle = name[[1]][2], perm = as.character(i)))
}))

ggplot(df_plt, aes(x = x, y = y, col = angle, shape = dataset)) +
  coord_fixed() +
  geom_point(size = 0.5, stroke = 0.01, alpha = 0.5) +
  theme_bw()
