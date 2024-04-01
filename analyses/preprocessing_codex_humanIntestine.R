## This file is for loading and preprocessing CODEX human intestine datasets

## Dryad information: https://datadryad.org/stash/dataset/doi:10.5061/dryad.pk0p2ngrf

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(here)

par(mfrow=c(1,1))

dataset_name <- "codex_humanIntestine"

# Load dataset ------------------------------------------------------------

## use OneDrive directory for now
data <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Hickey_Colon_CODEX/single_cells.csv.gz', row.names = 1)

donors <- c("B009", "B010", "B011", "B012")
tissue_locations <- c("Transverse", "Proximal Jejunum", "Duodenum", "Ascending", "Ileum", "Mid-jejunum", "Descending", "Descending - Sigmoid")

for (donor in donors) {
  for (tissue_location in tissue_locations) {
    print(donor)
    print(tissue_location)
    
    ## subset to tissue array
    data_sub <- data[data$donor == donor & data$Tissue_location == tissue_location,]
    
    if (nrow(data_sub) > 0) {
      ## SpatialCoordinates
      # extract
      pos <- data_sub[,c("Xcorr", "Ycorr")]
      colnames(pos) <- c("x", "y")
      # plot
      plot(pos)
      
      ## colData
      coldata <- data_sub[,c("size", "array", "donor", "cell_type", "Neighborhood", "Community", "Tissue.Unit")]
      
      rownames(pos) <- rownames(coldata) <- paste0("cell", seq(dim(pos)[1]))
      
      ## format into SpatialExperiment class
      spe <- SpatialExperiment::SpatialExperiment(
        spatialCoords = as.matrix(pos),
        colData = coldata
      )
      
      ## save
      saveRDS(spe, file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_loc_", tissue_location, "_preprocessed.RDS")))
    }
  }
}


# temporary ---------------------------------------------------------------

data <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/CODEX_humanIntestine/23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv', row.names = 1)

donor <- "B006"
# tissue_location <- "Proximal Jejunum"
tissue_location <- "Ascending"

celltypes <- paste0(unique(data$Cell.Type), collapse = ", ")

## subset to tissue array
data_sub <- data[data$donor == donor & data$Tissue_location == tissue_location,]

if (nrow(data_sub) > 0) {
  ## SpatialCoordinates
  # extract
  pos <- data_sub[,c("Xcorr", "Ycorr")]
  colnames(pos) <- c("x", "y")
  # plot
  plot(pos)
  
  ## colData
  coldata <- data_sub[,c("array", "region", "donor", "Cell.Type", "Neighborhood", "Community", "Tissue.Unit")]
  
  rownames(pos) <- rownames(coldata) <- paste0("cell", seq(dim(pos)[1]))
  
  ## format into SpatialExperiment class
  spe <- SpatialExperiment::SpatialExperiment(
    spatialCoords = as.matrix(pos),
    colData = coldata
  )
  
  ## save
  saveRDS(spe, file = here("outputs", paste0(dataset_name, "_donor_", donor, "_tissue_loc_", tissue_location, "_preprocessed.RDS")))
}




# Format into SpatialExperiment class ----------------------------------

spe <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = gexp, lognorm = gexp_lognorm),
  spatialCoords = as.matrix(pos),
  colData = coldata
)

saveRDS(spe, file = here("outputs", paste0(dataset_name, "_preprocessed.RDS")))

# ## use OneDrive directory for now
# data <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/CODEX_humanIntestine/23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv', row.names = 1)
# 
# tissue_locations <- unique(data$Tissue_location)
# tissues <- unique(data$tissue)
# donors <- unique(data$donor)
# unique_regions <- unique(data$unique_region)
# 
# genes <- colnames(data[1:which(colnames(data) == "x")-1])
# 
# id <- 20
# unique_regions[id]
# 
# data_sub <- data[data$unique_region == unique_regions[id],]
# gexp <- data_sub[,genes]
# meta <- data_sub[,which(colnames(data_sub) == "Tissue_location"):length(colnames(data_sub))]
# pos <- data_sub[,c("Xcorr", "Ycorr")]
# colnames(pos) <- c("x", "y")
# 
# ## make all coordinates positive
# pos[,1] <- pos[,1] - min(pos[,1])
# pos[,2] <- pos[,2] - min(pos[,2])
# 
# ## plot
# plot(pos, pch=".")
# 
# ## density of feature x observation matrix
# calculateDensity(gexp)
