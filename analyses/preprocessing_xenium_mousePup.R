## This file is for loading and preprocessing 10X Genomics Xenium human breast cancer dataset

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

## load feature x observation count matrix
file <- "~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_mousePup/cell_feature_matrix.h5"
rhdf5::h5ls(file)
data <- rhdf5::h5read("~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_mousePup/cell_feature_matrix.h5", "matrix")
rhdf5::h5closeAll()
names(data)

counts <- data$data
indices <- data$indices
indptr <- data$indptr
shp <- data$shape
features <- data$features
barcodes <- data$barcodes

## convert to a sparse matrix
gexp <- sparseMatrix(i = indices[] + 1, p = indptr[], 
                   x = as.numeric(x = counts[]), dims = shp[], repr = "T")
rownames(gexp) <- features$name
colnames(gexp) <- barcodes
class(gexp)
gexp <- as(gexp, "CsparseMatrix")
dim(gexp)

## load metadata
meta <- read.csv("~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_mousePup/cells.csv.gz")
head(meta)
pos <- cbind(meta$x_centroid, meta$y_centroid)
colnames(pos) <- c("x", "y")
rownames(pos) <- meta$cell_id
plot(meta$cell_area, meta$total_counts, pch=".")
plot(meta$nucleus_area, meta$total_counts, pch=".")
plot(meta$cell_area, meta$nucleus_area, pch=".")

## load analysis clusters
clusters <- read.csv("~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_mousePup/analysis/clustering/gene_expression_graphclust/clusters.csv")
rownames(clusters) <- clusters$Barcode

rownames(meta) <- rownames(pos) <- colnames(gexp)

calculateDensity(gexp)

# Preprocessing -----------------------------------------------------------

## remove genes that are not in the gene panel
gene_panel <- read.csv("~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_mousePup/Xenium_mMulti_v1_metadata_v2.csv")
gene_excluded <- setdiff(features$name, gene_panel$Gene)

gexp <- gexp[gene_panel$Gene,]
dim(gexp)
calculateDensity(gexp)

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
cell_area <- meta$cell_area
nuc_area <- meta$nucleus_area
names(cell_area) <- names(nuc_area) <- rownames(meta)
cell_area <- cell_area[good_cells]
nuc_area <- nuc_area[good_cells]

hist(cell_area)
hist(nuc_area)

par(mfrow=c(3,1))
hist(log10(colSums(gexp)+1))
hist(log10(colSums(gexp/cell_area*mean(cell_area))+1))
hist(log10(colSums(gexp/nuc_area*mean(nuc_area))+1))

## normalize by nucleus area
gexp_lognorm <- as(log10(gexp/nuc_area*mean(nuc_area) + 1), "CsparseMatrix")

calculateDensity(gexp_lognorm)

# Format into SpatialExperiment class -------------------------------------

coldata <- clusters[rownames(clusters) %in% good_cells,"Cluster", drop = FALSE]
colnames(coldata) <- c("cluster")
coldata$cluster <- as.factor(coldata$cluster)

spe <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = gexp, lognorm = gexp_lognorm),
  spatialCoords = as.matrix(pos),
  colData = coldata
)

saveRDS(spe, file = "outputs/xenium_mousePup_preprocessed.RDS")
