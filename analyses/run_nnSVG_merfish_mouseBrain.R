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

# Load dataset ------------------------------------------------------------

## use OneDrive directory for now
gexp <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_cortex_For_Daniel/raw/datasets_mouse_brain_map_BrainReceptorShowcase_Slice2_Replicate1_cell_by_gene_S2R1.csv.gz', row.names = 1)
meta <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_cortex_For_Daniel/raw/datasets_mouse_brain_map_BrainReceptorShowcase_Slice2_Replicate1_cell_metadata_S2R1.csv.gz', row.names = 1)
pos <- meta[, c('center_x', 'center_y')]
dim(gexp)

## flip Y coordinates
pos[,2] <- -pos[,2]
## make all coordinates positive
pos[,1] <- pos[,1] - min(pos[,1])
pos[,2] <- pos[,2] - min(pos[,2])

## plot
plot(pos, pch=".")

## density of feature x observation matrix
calculateDensity(gexp)

# Preprocessing -----------------------------------------------------------

## filter genes
bad_genes <- colnames(gexp)[grepl('Blank', colnames(gexp))]
gexp <- gexp[, !(colnames(gexp) %in% bad_genes)]

## filter cells
par(mfrow=c(2,3))
hist(log10(rowSums(gexp)+1))
hist(log10(rowSums(gexp>0)+1))
hist(log10(rowSums(gexp>1)+1))

good.cells <- rownames(gexp)[rowSums(gexp>2) > 0]
gexp <- gexp[good.cells,]
hist(log10(rowSums(gexp)+1))
hist(log10(rowSums(gexp>0)+1))
hist(log10(rowSums(gexp>4)+1))

pos1 <- pos1[intersect(rownames(pos1), good.cells),]
pos2 <- pos2[intersect(rownames(pos2), good.cells),]
pos3 <- pos3[intersect(rownames(pos3), good.cells),]
par(mfrow=c(1,3))
plot(pos1, pch=".")
plot(pos2, pch=".")
plot(pos3, pch=".")

calculateDensity(gexp)