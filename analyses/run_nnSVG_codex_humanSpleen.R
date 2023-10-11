## This file integrates rasterization with nnSVG to analyze MERFISH mouse whole brain coronal section datasets. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::document()
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(nnSVG)
library(here)
library(tidyr)
library(tibble)
library(dplyr)

par(mfrow=c(1,1))

dataset_name <- "codex_humanSpleen"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_preprocessed.RDS")))

plot(spatialCoords(spe), pch=".", asp=1)

# Run method --------------------------------------------------------------

res_list <- list("singlecell", 50, 100, 200, 400)

## Rotate dataset, rasterize, run nnSVG for each resolution
n_rotation <- 10
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)

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
    df <- rownames_to_column(as.data.frame(rowData(spe)), var = "gene")
    df <- cbind(rotation_deg = NA, num_points = num_points, df)
    print(df)
    return(df)
  } else {
    df <- do.call(rbind, lapply(angle_deg_list, function(deg) {
      print(paste0("Rotation (degrees): ", deg))
      ## rotate xy coordinates
      spe_rotated <- SpatialExperiment::SpatialExperiment(
        assays = assays(spe),
        spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
      )
      
      ## rasterization
      spe_rast <- SEraster::rasterizeGeneExpression(spe_rotated, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
      num_points = dim(spe_rast)[2]
      
      ## nnSVG
      spe_rast <- nnSVG::nnSVG(
        spe_rast,
        assay_name = "pixelval",
        BPPARAM = BiocParallel::MulticoreParam()
      )
      temp <- rownames_to_column(as.data.frame(rowData(spe_rast)), var = "gene")
      temp <- cbind(rotation_deg = deg, num_points = num_points, temp)
      print(temp)
      return(temp)
    }))
  }
  print(df)
  return(data.frame(dataset = dataset_name, resolution = res, df))
}))
saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))

## Measure runtime for each resolution
n_itr <- 5

runtime_results <- do.call(rbind, lapply(res_list, function(res) {
  out <- do.call(rbind, lapply(seq(n_itr), function(i) {
    print(paste0("Resolution: ", res, ", trial: ", i))
    if (res == "singlecell") {
      num_points = dim(spe)[2]
      
      start <- Sys.time()
      nnSVG::nnSVG(
        spe,
        assay_name = "lognorm",
        BPPARAM = BiocParallel::MulticoreParam()
      )
      runtime <- difftime(Sys.time(), start, units = "secs")
      return(data.frame(trial = i, num_points = num_points, runtime_rast = NA, runtime_nnsvg = runtime, runtime_total = runtime))
    } else {
      start1 <- Sys.time()
      spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
      runtime_rast <- difftime(Sys.time(), start1, units = "secs")
      
      num_points = dim(spe_rast)[2]
      
      start2 <- Sys.time()
      nnSVG::nnSVG(
        spe_rast,
        assay_name = "pixelval",
        BPPARAM = BiocParallel::MulticoreParam()
      )
      end <- Sys.time()
      runtime_nnsvg <- difftime(end, start2, units = "secs")
      runtime_total <- difftime(end, start1, units = "secs")
      return(data.frame(trial = i, num_points = num_points, runtime_rast = runtime_rast, runtime_nnsvg = runtime_nnsvg, runtime_total = runtime_total))
    }
  }))
  return(data.frame(dataset = dataset_name, resolution = res, out))
}))
saveRDS(runtime_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global_runtime.RDS")))


# Plot --------------------------------------------------------------------

## visualize cell type
df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], celltype = colData(spe)[,"celltypes"])
ggplot(df, aes(x = x, y = y, col = celltype)) +
  coord_fixed() +
  geom_point(size = 0.1) +
  theme_classic()

## Figure 1a (spatial plots)
# res <- list("singlecell", 50, 100, 200, 400)
res_list <- list(200, 400)

## total transcripts
for (res in res_list) {
  if (res == "singlecell") {
    ## plot
    df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], transcripts = colSums(assay(spe, "counts")))
    plt <- ggplot(df, aes(x = x, y = y, col = transcripts)) +
      coord_fixed() +
      geom_point(size = 0.1) +
      scale_color_viridis_c() +
      labs(title = "Single cell") +
      theme_bw() +
      theme(
        legend.position="none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
      )
  } else {
    ## rasterize
    spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "counts", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
    
    ## plot
    df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], celltype = colData(spe)[,"celltypes"])
    plt <- ggplot(df, aes(x = x, y = y, col = transcripts)) +
      coord_fixed() +
      geom_tile() +
      scale_fill_viridis_c() +
      labs(title = paste0("Resolution = ", i, " um")) +
      theme_bw() +
      theme(
        legend.position="none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
      )
  }
  ## save plot
  ggsave(plot = plt, filename = here("plots", dataset_name, paste0(dataset_name, "_tot_counts_", i, ".pdf")))
}

## celltypes
# pie charts?


# Further exploration -----------------------------------------------------

## nnSVG at single cell
spe <- nnSVG::nnSVG(
  spe,
  assay_name = "lognorm",
  BPPARAM = BiocParallel::MulticoreParam()
)
View(rowData(spe))

## rasterization
res <- 100
spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
dim(spe_rast)[2]
protein <- "CD3e"
df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], protein = assay(spe_rast)[protein,])
ggplot(df, aes(x = x, y = y, fill = protein)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_classic()

## nnSVG
spe_rast <- nnSVG::nnSVG(
  spe_rast,
  assay_name = "pixelval"
)
