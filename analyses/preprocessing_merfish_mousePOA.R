## This file is for loading and preprocessing MERFISH mouse preoptic area dataset

## https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248
## "Expression values for the 135 genes measured in the combinatorial smFISH run 
## were determined as the total counts per cell divided by the cell volume and scaled by 1000"

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

dataset_name <- "merfish_mousePOA"

# Load dataset ------------------------------------------------------------

## use OneDrive directory for now
# data <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_mousePOA/merfish_mousePOA_all_cells.csv')

## mayling's local directory
data <- read.csv('~/Downloads/lab/data/merfish_mousePOA_all_cells.csv')

## look at all the conditions
conditions <- paste0(data$Animal_ID, "_", data$Animal_sex, "_", data$Behavior, "_", data$Bregma)
unique(conditions)

## subset
animal <- 1
sex <- "Female"
behavior <- "Naive"
bregma <- "-0.29"
data_sub <- data[(data$Animal_ID == animal & data$Animal_sex == sex & data$Behavior == behavior & data$Bregma == bregma),]
dim(data_sub)

# animals <- seq(1, 11, by = 1)
# sexes <- c("Female", "Male")
# bregmas <- seq(-0.29, 0.26, by = 0.05)

animals <- unique(data$Animal_ID)
sexes <- unique(data$Animal_sex)
bregmas <- unique(data$Bregma)

count = 0


for (animal in animals) {
  for (sex in sexes) {
    for (bregma in bregmas) {
      data_sub <- data[(data$Animal_ID == animal & data$Animal_sex == sex & data$Behavior == behavior & data$Bregma == bregma),]
      current <- paste(animal, sex, behavior, bregma, sep = "_")
      # checks for unique conditions and valid number of cells
      
      if(current %in% unique(conditions) & nrow(data_sub) > 0) {
        print(paste0(animal, "_", sex, "_", behavior, "_",  bregma))
        data_sub <- data[(data$Animal_ID == animal & data$Animal_sex == sex & data$Behavior == behavior & data$Bregma == bregma),]
        ## extract features x observations matrix, spatial coordinates, meta data
        ## genes x cells matrix ("total counts per cell divided by the cell volume and scaled by 1000")
        mat <- as(t(data_sub[,10:ncol(data_sub)]), "CsparseMatrix")
        blanks <- rownames(mat)[grepl("Blank", rownames(mat))]
        mat <- mat[setdiff(rownames(mat),blanks),]
        
        ## spatial coordinates
        pos <- data_sub[,c("Centroid_X", "Centroid_Y")]
        colnames(pos) <- c("x","y")
        ## make x,y coordinates positive
        pos[,1] <- pos[,1] - min(pos[,1])
        pos[,2] <- pos[,2] - min(pos[,2])
        
        ## meta data
        meta <- data_sub[,c("Bregma", "Cell_class", "Neuron_cluster_ID")]
        colnames(meta) <- c("bregma", "celltype", "neurontype")
        
        colnames(mat) <- rownames(pos) <- rownames(meta) <- data_sub$Cell_ID
        
        ## filter genes with NaN values
        bad_genes <- names(which(rowSums(is.nan(mat)) > 0))
        mat <- mat[setdiff(rownames(mat),bad_genes),]
        dim(mat)
        
        ## filter cells with NaN values
        bad_cells <- names(which(colSums(is.nan(mat)) > 0))
        mat <- mat[,setdiff(colnames(mat),bad_cells)]
        pos <- pos[setdiff(rownames(pos),bad_cells),]
        meta <- meta[setdiff(rownames(pos),bad_cells),]
        
        ## log transformation
        par(mfrow=c(2,1))
        hist(colSums(mat))
        hist(log10(colSums(mat) + 1))
        
        mat_lognorm <- as(log10(mat + 1), "CsparseMatrix")
        
        calculateDensity(mat_lognorm)
        
        # format into SpatialExperiment class -------------------------------------
        
        spe <- SpatialExperiment::SpatialExperiment(
          assays = list(volnorm = mat, lognorm = mat_lognorm),
          spatialCoords = as.matrix(pos),
          colData = meta
        )
        count = count + 1
        #print(animals, "_", sexes, "_", behavior, "_",  bregmas)
        #saveRDS(spe, file = here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_preprocessed.RDS")))
        
      }
    }
  }
}
print(count)

  