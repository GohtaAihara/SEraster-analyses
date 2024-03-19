## This file integrates rasterization with CooccurrenceAffinity to analyze CODEX human intestine datasets. Produce figures for 

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

dataset_name = "codex_humanIntestine"
method = "CooccurrenceAffinity"

# Load dataset ------------------------------------------------------------

donors <- c("B009", "B010", "B011", "B012")
tissue_locations <- c("Transverse", "Proximal Jejunum", "Duodenum", "Ascending", "Ileum", "Mid-jejunum", "Descending", "Descending - Sigmoid")

donor <- donors[[1]]
tissue_location <- tissue_locations[[1]]

