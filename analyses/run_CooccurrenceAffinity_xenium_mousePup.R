## This file integrates rasterization with CooccurrenceAffinity to analyze 10X Genomics Xenium mouse whole pup datasets. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(here)
library(CooccurrenceAffinity)

par(mfrow=c(1,1))

dataset_name = "xenium_mousePup"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = "outputs/xenium_mousePup_preprocessed.RDS")

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  geom_point(size = 0.1) +
  theme_classic()

ct_labels <- colData(spe)$cluster

# Run methods -------------------------------------------------------------

res <- 100

spe_rast <- SEraster::rasterizeCellType(spe, )
