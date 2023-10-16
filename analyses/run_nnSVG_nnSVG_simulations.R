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
# sim_names <- c(
#   "sim_largeBandwidth_fullExpr", 
#   "sim_largeBandwidth_medExpr", 
#   "sim_largeBandwidth_lowExpr", 
#   "sim_medBandwidth_fullExpr", 
#   "sim_medBandwidth_medExpr", 
#   "sim_medBandwidth_lowExpr", 
#   "sim_smallBandwidth_fullExpr", 
#   "sim_smallBandwidth_medExpr", 
#   "sim_smallBandwidth_lowExpr"
# )

sim_names <- c(
  "sim_largeBandwidth_fullExpr"
)

## load each dataset when we analyze

# Run method -------------------------------------------------------------

## iterate over 1. dataset, 2. resolution, 3. rotation (save everything in one df for each dataset)
res_list <- list("singlecell", 0.04, 0.08)
n_rotation <- 2
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)

for (i in sim_names) {
  ## load dataset,
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  
}


# Further exploration -----------------------------------------------------

## rasterize
resolution <- 0.07
spe_rast <- SEraster::rasterizeGeneExpression(spe, resolution = resolution, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

## plot rasterized data
svg_example <- which(rowData(spe)$expressed)[1]
df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], gene = assay(spe_rast)[svg_example,])

ggplot(df, aes(x = x, y = y, fill = gene)) +
  geom_tile(width = resolution, height = resolution) +
  coord_fixed() +
  scale_fill_viridis_c(name = "svg_example") +
  labs(title = paste0("Resolution = ", resolution)) +
  theme_classic()

spe_rast <- nnSVG::nnSVG(
  spe_rast,
  assay_name = "pixelval",
  n_threads = 14
)
