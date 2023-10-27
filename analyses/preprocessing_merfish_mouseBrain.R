## This file is for loading and preprocessing MERFISH mouse whole brain datasets

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(rhdf5)

par(mfrow=c(1,1))

# Load dataset ------------------------------------------------------------

## use OneDrive directory for now
## load feature x gene counts matrix
gexp <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_cortex_For_Daniel/raw/datasets_mouse_brain_map_BrainReceptorShowcase_Slice2_Replicate1_cell_by_gene_S2R1.csv.gz', row.names = 1)

## convert to sparse matrix
gexp <- as(t(gexp), "CsparseMatrix")

## load meta data
meta <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_cortex_For_Daniel/raw/datasets_mouse_brain_map_BrainReceptorShowcase_Slice2_Replicate1_cell_metadata_S2R1.csv.gz', row.names = 1)
pos <- meta[, c('center_x', 'center_y')]

colnames(pos) <- c("x", "y")
colnames(gexp) <- rownames(meta) <- rownames(pos) <- paste0("cell-", rownames(meta))

## flip Y coordinates
pos[,2] <- -pos[,2]
## make all coordinates positive
pos[,1] <- pos[,1] - min(pos[,1])
pos[,2] <- pos[,2] - min(pos[,2])

## plot
plot(pos, pch=".")

## density of feature x observation matrix
calculateDensity(gexp)

## load another metadata
meta2 <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_cortex_For_Daniel/s1r1_metadata.csv.gz', row.names = 1)

# Preprocessing -----------------------------------------------------------

## filter genes
bad_genes <- rownames(gexp)[grepl('Blank', rownames(gexp))]
gexp <- gexp[!(rownames(gexp) %in% bad_genes),]

## filter cells
par(mfrow=c(2,3))
hist(log10(colSums(gexp)+1))
hist(log10(colSums(gexp>0)+1))
hist(log10(colSums(gexp>1)+1))

good_cells <- colnames(gexp)[colSums(gexp>2) > 0]
gexp <- gexp[,colnames(gexp) %in% good_cells]
hist(log10(colSums(gexp)+1))
hist(log10(colSums(gexp>0)+1))
hist(log10(colSums(gexp>1)+1))

pos <- pos[rownames(pos) %in% good_cells,]

par(mfrow=c(1,1))
plot(pos, pch=".")

calculateDensity(gexp)

## normalization

vol <- meta$volume
names(vol) <- rownames(meta)
vol <- vol[good_cells]

hist(log10(vol))
plot(vol, colSums(gexp), pch = ".")

par(mfrow=c(2,1))
hist(log10(colSums(gexp)+1))
hist(log10(colSums(gexp/vol*mean(vol))+1))

gexp_lognorm <- as(log10(gexp/vol*mean(vol) + 1), "CsparseMatrix")

calculateDensity(gexp_lognorm)

# format into SpatialExperiment class -------------------------------------

spe <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = gexp, lognorm = gexp_lognorm),
  spatialCoords = as.matrix(pos)
)

saveRDS(spe, file = "outputs/merfish_mouseBrain_preprocessed.RDS")
