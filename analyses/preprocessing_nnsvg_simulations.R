## This file is for loading and modifyng spatial coordiantes of the simulated datasets from Weber, L.M., et al., 2023, Nature Communications. (Original files are taken from https://github.com/lmweber/nnSVG-analyses/blob/main/analyses/08_simulations/simulations_build_data.R)

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(rhdf5)
library(here)

par(mfrow=c(1,1))

dir_input <- here("inputs", "nnSVG_simulations")
dir_output <- here("outputs", "nnSVG_simulations")

# Load dataset ------------------------------------------------------------

## simulations for various bandwidth and expression
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

for (i in sim_names) {
  spe <- readRDS(file = here(dir_input, paste0("spe_", i, ".RDS")))
  plot(spatialCoords(spe), pch=".", main=i)
  
  ## modify
  ## add gene names
  rownames(spe) <- paste0("gene", seq(dim(spe)[1]))
  ## add cell IDs
  colnames(spe) <- paste0("cell", seq(dim(spe)[2]))
  rownames(spatialCoords(spe)) <- paste0("cell", seq(dim(spe)[2]))
  ## scale coordinates to 6,000 um by 6,000 um
  spatialCoords(spe) <- spatialCoords(spe) * 6000
  
  plot(spatialCoords(spe), pch=".", main=i)
  
  saveRDS(spe, file = here(dir_output, paste0("spe_", i, ".RDS")))
}

## simulations for shuffle based on medium bandwidth, medium expression
sim_names_shuffle <- c(
  "sim_shuffle00", 
  "sim_shuffle01", 
  "sim_shuffle02", 
  "sim_shuffle03", 
  "sim_shuffle04", 
  "sim_shuffle05", 
  "sim_shuffle06", 
  "sim_shuffle07", 
  "sim_shuffle08", 
  "sim_shuffle09", 
  "sim_shuffle10"
)

for (i in sim_names_shuffle) {
  spe <- readRDS(file = here(dir_input, paste0("spe_", i, ".RDS")))
  plot(spatialCoords(spe), pch=".", main=i)
  
  ## modify
  ## add gene names
  rownames(spe) <- paste0("gene", seq(dim(spe)[1]))
  ## add cell IDs
  colnames(spe) <- paste0("cell", seq(dim(spe)[2]))
  rownames(spatialCoords(spe)) <- paste0("cell", seq(dim(spe)[2]))
  ## scale coordinates to 6,000 um by 6,000 um
  spatialCoords(spe) <- spatialCoords(spe) * 6000
  
  plot(spatialCoords(spe), pch=".", main=i)
  
  saveRDS(spe, file = here(dir_output, paste0("spe_", i, ".RDS")))
}
