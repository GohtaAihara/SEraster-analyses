## This file integrates rasterization with nnSVG to analyze MERFISH mouse whole brain coronal section datasets. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(nnSVG)
library(here)
library(tidyr)
library(tibble)
library(dplyr)

par(mfrow=c(1,1))

dir <- here("outputs", "nnSVG_simulations")

dataset_name <- "nnSVG_simulations"


# Load dataset ------------------------------------------------------------

## filenames for saved simulation datasets
sim_names <- c(
  "sim_largeBandwidth_fullExpr",
  "sim_largeBandwidth_medExpr",
  "sim_largeBandwidth_lowExpr",
  "sim_medBandwidth_fullExpr",
  "sim_medBandwidth_medExpr",
  "sim_medBandwidth_lowExpr",
  "sim_smallBandwidth_fullExpr",
  "sim_smallBandwidth_medExpr",
  "sim_smallBandwidth_lowExpr"
)

# sim_names <- c(
#   "sim_largeBandwidth_fullExpr",
#   "sim_largeBandwidth_medExpr"
# )

## load each dataset when we analyze

# Run method -------------------------------------------------------------

## iterate over 1. dataset, 2. resolution, 3. rotation (save everything in one df for each dataset)
res_list <- c(list("singlecell"), as.list(seq(0.01, 0.09, by = 0.01)))
# res_list <- list("singlecell", 0.09)
n_rotation <- 10
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)

for (i in sim_names) {
  print(paste0("Dataset: ", i))
  
  ## load dataset,
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  
  ## iterate over resolutions
  nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
    print(paste0("Resolution: ", res))
    if (res == "singlecell") {
      ## iterate over rotations (single cell)
      df <- do.call(rbind, lapply(angle_deg_list, function(deg) {
        print(paste0("Rotation (degrees): ", deg))
        ## rotate xy coordinates
        spe_rotated <- SpatialExperiment::SpatialExperiment(
          assays = assays(spe),
          spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
        )
        
        num_points = dim(spe_rotated)[2]
        
        ## nnSVG
        spe_rotated <- nnSVG::nnSVG(
          spe_rotated,
          assay_name = "logcounts",
          BPPARAM = BiocParallel::MulticoreParam()
        )
        temp <- rownames_to_column(as.data.frame(rowData(spe_rotated)), var = "gene")
        temp <- cbind(rotation_deg = deg, num_points = num_points, temp)
        return(temp)
      }))
    } else {
      ## iterate over rotations (rasterization)
      df <- do.call(rbind, lapply(angle_deg_list, function(deg) {
        print(paste0("Rotation (degrees): ", deg))
        ## rotate xy coordinates
        spe_rotated <- SpatialExperiment::SpatialExperiment(
          assays = assays(spe),
          spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
        )
        
        ## rasterization
        spe_rast <- SEraster::rasterizeGeneExpression(spe_rotated, assay_name = "logcounts", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
        num_points = dim(spe_rast)[2]
        
        ## nnSVG
        spe_rast <- nnSVG::nnSVG(
          spe_rast,
          assay_name = "pixelval",
          BPPARAM = BiocParallel::MulticoreParam()
        )
        temp <- rownames_to_column(as.data.frame(rowData(spe_rast)), var = "gene")
        temp <- cbind(rotation_deg = deg, num_points = num_points, temp)
        return(temp)
      }))
    }
    return(data.frame(dataset = i, resolution = res, df))
  }))
  saveRDS(nnsvg_results, file = here("outputs", dataset_name, paste0(i, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
}


# Plot --------------------------------------------------------------------

## Figure x (single cell and rasterization visualizations)
res_list <- c(list("singlecell"), as.list(seq(0.01, 0.1, by = 0.01)))

for (i in sim_names) {
  ## load dataset
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  svg_example <- which(rowData(spe)$expressed)[1]
  
  for (res in res_list) {
    if (res == "singlecell") {
      ## plot
      df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], gene = assay(spe)[svg_example,])
      plt <- ggplot(df, aes(x = x, y = y, col = gene)) +
        coord_fixed() +
        geom_point(size = 1) +
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
      spe_rast <- SEraster::rasterizeGeneExpression(spe, resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
      
      ## plot
      df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], gene = assay(spe_rast)[svg_example,])
      plt <- ggplot(df, aes(x = x, y = y, fill = gene)) +
        coord_fixed() +
        geom_tile(width = res, height = res) +
        scale_fill_viridis_c() +
        labs(title = paste0("Resolution = ", res, " %")) +
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
    ggsave(plot = plt, filename = here("plots", dataset_name, paste0(i, "_gene", svg_example, "_",  res, ".pdf")))
  }
}

## Figure (performance comparison)
i <- sim_names[2]
n_rotation <- 2
df <- readRDS(file = here("outputs", dataset_name, paste0(i, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))

# Further exploration -----------------------------------------------------

## load dataset
i <- sim_names[1]
i
spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))

# ## single cell
# spe_nnsvg <- nnSVG::nnSVG(
#   spe,
#   BPPARAM = BiocParallel::MulticoreParam()
# )

# ## rasterize (debug)
# data <- assay(spe)
# pos <- spatialCoords(spe)
# bbox <- bbox <- sf::st_bbox(c(
#   xmin = floor(min(pos[,1])-resolution/2), 
#   xmax = ceiling(max(pos[,1])+resolution/2), 
#   ymin = floor(min(pos[,2])-resolution/2), 
#   ymax = ceiling(max(pos[,2])+resolution/2)
# ))
# fun <- "mean"
# 
# ## create grid for rasterization
# grid <- sf::st_make_grid(bbox, cellsize = resolution)
# 
# ## extract position
# pos_pixel <- sf::st_coordinates(sf::st_centroid(grid))
# colnames(pos_pixel) <- c("x", "y")
# rownames(pos_pixel) <- paste0("pixel",seq_along(pos_pixel[,1]))
# 
# ## convert pos to sf
# cells <- sf::st_as_sf(data.frame(pos), coords = c("x", "y"))
# ## get pixel ID for each cell
# ## since some cells are assigned to > 1 pixels if they are on the border, choose the 1st pixel id
# pixel_ids <- unlist(lapply(sf::st_intersects(cells, grid), function(sublist) sublist[1]))
# names(pixel_ids) <- rownames(pos)
# 
# ## store aggregated subsetted matrix data for each pixel, store cell IDs and number of cells for each pixel into a data frame
# out <- unlist(BiocParallel::bplapply(sort(unique(pixel_ids)), function(id){
#   ## get cell IDs for a particular pixel
#   cell_ids <- names(pixel_ids[pixel_ids == id])
#   
#   ## subset feature observation matrix
#   sub <- data[,cell_ids, drop = FALSE]
#   ## aggregate cell counts to create pixel value
#   if (fun == "mean") {
#     pixel_val <- rowMeans(sub)
#   } else if (fun == "sum") {
#     pixel_val <- rowSums(sub)
#   }
#   
#   ## store number of cells
#   meta_rast <- data.frame(num_cell = length(cell_ids))
#   ## store a list of cell IDs
#   meta_rast$cellID_list <- list(cell_ids)
#   
#   if (is.matrix(data)) {
#     return(list(pixel_val, meta_rast))
#   } else {
#     return(list(as(pixel_val, "CsparseMatrix"), meta_rast))
#   }
# }), recursive = FALSE)
# 
# ## extract rasterized sparse matrix
# data_rast <- do.call(cbind, out[seq(1,length(out),by=2)])
# 
# ## extract rasterized data frame
# meta_rast <- do.call(rbind, out[seq(1,length(out),by=2)+1])
# 
# ## set rownames/colnames for rasterized sparse matrix and rasterized data frame
# rownames(data_rast) <- rownames(data)
# colnames(data_rast) <- paste0("pixel", sort(unique(pixel_ids)))
# rownames(meta_rast) <- paste0("pixel", sort(unique(pixel_ids)))
# 
# ## subset rasterized pos
# pos_pixel <- pos_pixel[rownames(pos_pixel) %in% colnames(data_rast),]

## rasterize
resolution <- 0.1
spe_rast <- SEraster::rasterizeGeneExpression(spe, resolution = resolution, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

## plot rasterized data
svg_example <- which(rowData(spe)$expressed)[1]
df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], gene = assay(spe_rast)[svg_example,])

ggplot(df, aes(x = x, y = y, fill = gene)) +
  # geom_tile(width = resolution, height = resolution) +
  # geom_tile() +
  geom_point() +
  coord_fixed() +
  scale_fill_viridis_c(name = "svg_example") +
  labs(title = paste0("Resolution = ", resolution)) +
  theme_classic()

## see why geom_tile() doesn't work
a <- diff(spatialCoords(spe_rast))

spe_rast <- nnSVG::nnSVG(
  spe_rast,
  assay_name = "pixelval",
  BPPARAM = BiocParallel::MulticoreParam()
)
