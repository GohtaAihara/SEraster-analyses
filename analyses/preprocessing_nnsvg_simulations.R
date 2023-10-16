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
  ## no modification
  
  saveRDS(spe, file = here(dir_output, paste0("spe_", i, ".RDS")))
}
