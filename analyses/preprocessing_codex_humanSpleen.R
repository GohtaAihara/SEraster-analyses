## This file is for loading and preprocessing MERFISH mouse whole brain datasets

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(crawdad)
library(Matrix)
library(ggplot2)
library(stringr)
library(here)

par(mfrow=c(1,1))

dataset_name <- "codex_humanSpleen"

# Load dataset ------------------------------------------------------------

## use OneDrive directory for now
data <- readRDS('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/CRAWDAD/spleen_combined_volNorm_log10_mat.RDS')
dim(data)

calculateDensity(data)

## filter
cells <- rownames(data)
smpl_cells <- str_detect(cells, 'PKHL')
gexp <- data[smpl_cells, ]
dim(gexp)
temp <- rownames(gexp)
rownames(gexp) <- str_replace(temp, '.*?_', '')
gexp <- t(gexp)

## get position and celltype labels from CRAWDAD dataset
data("pkhl")
meta <- pkhl
pos <- cbind(meta$x, meta$y)
colnames(pos) <- c("x", "y")

colnames(gexp) <- rownames(meta) <- rownames(pos) <- paste0("cell-", colnames(gexp))

## make all coordinates positive
pos[,1] <- pos[,1] - min(pos[,1])
pos[,2] <- pos[,2] - min(pos[,2])

## plot
plot(pos, pch=".")

## density of feature x observation matrix (store as dense matrix)
calculateDensity(gexp)

# format into SpatialExperiment class -------------------------------------

coldata <- meta[,"celltypes", drop = FALSE]

spe <- SpatialExperiment::SpatialExperiment(
  assays = list(lognorm = gexp),
  spatialCoords = as.matrix(pos),
  colData = coldata
)

saveRDS(spe, file = here("outputs", paste0(dataset_name, "_preprocessed.RDS")))
