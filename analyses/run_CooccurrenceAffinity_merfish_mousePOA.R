## This file integrates rasterization with CooccurrenceAffinity to analyze MERFISH mouse preoptic area (mousePOA) datasets. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::document()
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(ggrastr)
library(gridExtra)
library(here)
library(CooccurrenceAffinity)
library(tidyr)
library(tibble)
library(dplyr)
library(reshape)
library(DescTools)

par(mfrow=c(1,1))

dataset_name <- "merfish_mousePOA"
method <- "CooccurrenceAffinity"

# Load dataset ------------------------------------------------------------

# specify
animal <- 1
sex <- "Female"
behavior <- "Naive"
bregma <- "-0.29"

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_preprocessed.RDS")))

plot(spatialCoords(spe), pch=".", asp=1)

ct_labels <- as.factor(colData(spe)$celltype)

# Run method --------------------------------------------------------------

## modify code from "run_CooccurrenceAffinity_codex_humanSpleen.R"


# Plot --------------------------------------------------------------

## modify code from "run_CooccurrenceAffinity_codex_humanSpleen.R" (hiearchical clustering of heatmap)


# Questions ---------------------------------------------------------------

## can we identify cell-type cooccurrence patterns in mPOA? if so, can they be validated by visualizing them at single-cell resolution?
## are cell-type cooccurrence patterns consistent across mice for each bregma?
