# test drosophila 3D data

# single-cell 3D spatiotemporal Drosophila embryogenesis to metamorphosis data from https://doi.org/10.1101/2024.02.06.577903


# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(anndata)
library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(rhdf5)
library(here)
library(tibble)
library(dplyr)
library(reshape2)
library(pheatmap)

par(mfrow=c(1,1))

dataset_name = "stereoseq_drosophilaEmbryogenesis"


# Load dataset ------------------------------------------------------------

## use anndata
ad <- anndata::read_h5ad('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Stereo-seq_drosophilaEmbryogenesis/E14-16h_a_count_normal_stereoseq.h5ad')

## directly load h5ad file
# E14-16h
# file <- '~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Stereo-seq_drosophilaEmbryogenesis/E14-16h_a_count_normal_stereoseq.h5ad'
# L2
file <- '~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Stereo-seq_drosophilaEmbryogenesis/L2_a_count_normal_stereoseq.h5ad'
rhdf5::h5ls(file)
# genes x cells (normalized?)
x <- rhdf5::h5read(file, "X")
# layers (raw counts)
layers <- rhdf5::h5read(file, "layers")
# obs
obs <- rhdf5::h5read(file, "obs")
rhdf5::h5closeAll()

## extract info
pos <- data.frame(obs$new_x, obs$new_y, obs$new_z)
rownames(pos) <- obs$`_index`
colnames(pos) <- c("x","y","z")
dim(pos)

## plot
plotly::plot_ly(data = pos, x = ~x, y = ~y, z = ~z, size = 2)
