## This file is for loading and modifyng spatial coordiantes of the simulated datasets from Weber, L.M., et al., 2023, Nature Communications. (Original files are taken from https://github.com/lmweber/nnSVG-analyses/blob/main/analyses/08_simulations/simulations_build_data.R)

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(crawdad)
library(here)

par(mfrow=c(1,1))

dataset_name <- "crawdad_simulation"

# Load dataset ------------------------------------------------------------

data("sim")

## format into SpatialExperiment object
spe <- SpatialExperiment::SpatialExperiment(
  spatialCoords = as.matrix(sim[,c("x","y")]),
  colData = sim[,"celltypes", drop = FALSE]
)

saveRDS(spe, file = here("outputs", paste0(dataset_name, "_preprocessed.RDS")))