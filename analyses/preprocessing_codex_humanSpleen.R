## This file is for loading and preprocessing MERFISH mouse whole brain datasets

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(stringr)

par(mfrow=c(1,1))

# Load dataset ------------------------------------------------------------

## use OneDrive directory for now
data <- readRDS('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/CRAWDAD/spleen_combined_volNorm_log10_mat.RDS')
dim(data)

calculateDensity(data)

## filter
cells <- rownames(data)
smpl_cells <- str_detect(cells, 'PKHL')
mtx <- data[smpl_cells, ]
dim(mtx)
calculateDensity(mtx)
tmp <- rownames(mtx)
rownames(mtx) <- str_replace(tmp, '.*?_', '')
mtx_id <- cbind(mtx, id = rownames(mtx))
