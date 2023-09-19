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

par(mfrow=c(1,1))

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = "outputs/merfish_mouseBrain_preprocessed.RDS")

plot(spatialCoords(spe), pch=".")