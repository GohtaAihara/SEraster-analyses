## This file integrates rasterization with nnSVG to analyze 10X Genomics Xenium mouse whole pup datasets. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(here)

par(mfrow=c(1,1))

dataset_name = "xenium_mousePup"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = "outputs/xenium_mousePup_preprocessed.RDS")

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  geom_point(size = 0.1) +
  theme_classic()

ct_labels <- colData(spe)$cluster

# Run methods -------------------------------------------------------------

res_list <- list(100)

## Global nnSVG
nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
  print(paste0("Resolution: ", res))
  if (res == "singlecell") {
    num_points = dim(spe)[2]
    
    ## nnSVG
    spe <- nnSVG::nnSVG(
      spe,
      assay_name = "lognorm",
      BPPARAM = BiocParallel::MulticoreParam()
    )
    df <- tibble::rownames_to_column(as.data.frame(rowData(spe)), var = "gene")
  } else {
    ## rasterization
    spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
    num_points = dim(spe_rast)[2]
    
    ## nnSVG
    spe_rast <- nnSVG::nnSVG(
      spe_rast,
      assay_name = "pixelval",
      BPPARAM = BiocParallel::MulticoreParam()
    )
    df <- tibble::rownames_to_column(as.data.frame(rowData(spe_rast)), var = "gene")
  }
  return(data.frame(dataset = dataset_name, resolution = res, num_cells = dim(spe)[2], num_pixels = dim(spe_rast)[2], df))
}))
## save results
saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global.RDS")))

## cluster-specific nnSVG
nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
  out <- do.call(rbind, lapply(levels(ct_labels), function(ct_label) {
    ## subset cluster of interest
    spe_sub <- spe[,spe$cluster == ct_label]
    print(paste0("# of cells in cluster ", ct_label, " : ", dim(spe_sub)[2]))
    
    ## rasterization
    spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
    
    ## nnSVG
    ## using try() to handle error
    spe_sub_rast_nnsvg <- try({nnSVG::nnSVG(
      spe_sub_rast,
      assay_name = "pixelval",
      BPPARAM = BiocParallel::MulticoreParam()
    )})
    
    if (class(spe_sub_rast_nnsvg) == "try-error") {
      ## do not save anything for clusters that caused error in nnSVG
      return(NULL)
    } else {
      df <- tibble::rownames_to_column(as.data.frame(rowData(spe_sub_rast_nnsvg)), var = "gene")
      return(data.frame(
        cluster = ct_label, 
        num_cells = dim(spe_sub)[2], 
        num_pixels = dim(spe_sub_rast)[2], 
        df))
    }
  }))
  return(cbind(dataset = dataset_name, resolution = res, out))
}))
## save results
saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_ct_specific.RDS")))


# Plot --------------------------------------------------------------------



# Further exploration -----------------------------------------------------

ct_label <- 2
spe_sub <- spe[,spe$cluster == ct_label]

df <- data.frame(spatialCoords(spe_sub), colData(spe_sub))
ggplot(df, aes(x = x, y = y, col = cluster)) +
  geom_point(size = 0.1) +
  theme_classic()

## rasterization
res <- 100
spe_sub_rast <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

df <- data.frame(spatialCoords(spe_sub_rast), colData(spe_sub_rast), transcripts = colMeans(assay(spe_sub_rast)))
ggplot(df, aes(x = x, y = y, fill = transcripts)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_classic()

## nnSVG
## using try() to handle error
start <- Sys.time()
spe_sub_rast_nnsvg <- try({nnSVG::nnSVG(
  spe_sub_rast,
  assay_name = "pixelval",
  BPPARAM = BiocParallel::MulticoreParam()
)})
difftime(Sys.time(), start)

rowData(spe_sub_rast_nnsvg)
