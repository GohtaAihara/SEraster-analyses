## This file is for loading and preprocessing seqFISH mouse kidney dataset

## https://spatialgenomics.com/data/

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

dataset_name <- "seqfish_mouseKidney"

# # Load dataset ------------------------------------------------------------
# 
# ## use OneDrive directory for now
# data <- readRDS('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/seqFISH_mouseKidney/seqfish_mouseKidney_norm_data.rds')
# 
# ## genes-by-cells matrix
# gexp <- data$all.norms$nonorm
# dim(gexp)
# 
# ## x,y coordinates
# pos <- data.frame(x = data$meta$center_x, y = data$meta$center_y)
# 
# ## meta
# meta <- data.frame(region = data$meta$region)
# 

# Load dataset ------------------------------------------------------------

section <- 1

## genes-by-cells matrix
gexp <- read.csv(paste0('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/seqFISH_mouseKidney/SG_MouseKidneyDataRelease_CxG_section', section, '.csv'))
class(gexp)
dim(gexp)

## meta data
meta <- read.csv(paste0('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/seqFISH_mouseKidney/SG_MouseKidneyDataRelease_CellCoordinates_section', section, '.csv'))
class(meta)
dim(meta)
