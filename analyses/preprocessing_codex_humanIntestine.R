## This file is for loading and preprocessing CODEX human intestine datasets

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)

par(mfrow=c(1,1))

# Load dataset ------------------------------------------------------------

## use OneDrive directory for now
data <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/CODEX_humanIntestine/23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv', row.names = 1)

tissue_locations <- unique(data$Tissue_location)
tissues <- unique(data$tissue)
donors <- unique(data$donor)
unique_regions <- unique(data$unique_region)

genes <- colnames(data[1:which(colnames(data) == "x")-1])

id <- 20
unique_regions[id]

data_sub <- data[data$unique_region == unique_regions[id],]
gexp <- data_sub[,genes]
meta <- data_sub[,which(colnames(data_sub) == "Tissue_location"):length(colnames(data_sub))]
pos <- data_sub[,c("Xcorr", "Ycorr")]
colnames(pos) <- c("x", "y")

## make all coordinates positive
pos[,1] <- pos[,1] - min(pos[,1])
pos[,2] <- pos[,2] - min(pos[,2])

## plot
plot(pos, pch=".")

## density of feature x observation matrix
calculateDensity(gexp)
