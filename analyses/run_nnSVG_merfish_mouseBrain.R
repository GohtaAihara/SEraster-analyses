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

par(mfrow=c(1,1))

dataset_name <- "merfish_mouseBrain"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_preprocessed.RDS")))

plot(spatialCoords(spe), pch=".")


# Run method --------------------------------------------------------------

res_list <- list("singlecell", 50, 100, 200, 400)

## Run nnSVG once for each resolution
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
  return(data.frame(dataset = dataset_name, resolution = res, num_points = num_points, df))
}))
saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global.RDS")))


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


## Measure memory for each resolution (serial computing due to bench package)
res_list <- list(100, 400)
n_itr <- 2

memory_results <- do.call(rbind, lapply(res_list, function(res) {
  out <- do.call(rbind, lapply(seq(n_itr), function(i) {
    print(paste0("Resolution: ", res, ", trial: ", i))
    if (res == "singlecell") {
      num_points = dim(spe)[2]
      
      mem <- bench::bench_memory(
        nnSVG::nnSVG(
          spe,
          assay_name = "lognorm",
          n_threads = 1
        )
      )
      return(data.frame(trial = i, num_points = num_points, memory_rast = NA, memory_nnsvg = mem$mem_alloc, memory_total = mem$mem_alloc))
    } else {
      mem_rast <- bench::bench_memory(
        spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", n_threads = 1)
      )
      
      num_points = dim(spe_rast)[2]
      
      mem_nnsvg <- bench::bench_memory(
        nnSVG::nnSVG(
          spe_rast,
          assay_name = "pixelval",
          n_threads = 1
        )
      )
      return(data.frame(trial = i, num_points = num_points, memory_rast = mem_rast$mem_alloc, memory_nnsvg = mem_nnsvg$mem_alloc, memory_total = mem_rast$mem_alloc + mem_nnsvg$mem_alloc))
    }
  }))
  return(data.frame(dataset = dataset_name, resolution = res, out))
}))
saveRDS(memory_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global_memory.RDS")))




res <- c(100)
spe_rast1 <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = 100, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
spe_rast2 <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = 200, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

memory_results <- bench::mark(
  nnSVG::nnSVG(
    spe_rast1,
    assay_name = "pixelval",
    n_threads = 1
  ),
  nnSVG::nnSVG(
    spe_rast2,
    assay_name = "pixelval",
    n_threads = 1
  )
)

x <- 1:1e6
result <- bench::mark(
  is.numeric(sqrt(x)),
  is.numeric(x^2)
)

test <- bench::bench_memory(
  spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = 400, fun = "mean", BPPARAM = BiocParallel::MulticoreParam(workers = 1))
  # sqrt(x)
)

test <- bench::bench_memory(
  spe_rast <- SEraster::rasterizeSparseMatrix2(assay(spe, "lognorm"), spatialCoords(spe), resolution = 400, fun = "mean", n_threads = 22)
)

spe_rast <- SEraster::rasterizeSparseMatrix2(assay(spe, "lognorm"), spatialCoords(spe), resolution = 400, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
spe_rast2 <- SEraster::rasterizeSparseMatrix2(assay(spe, "lognorm"), spatialCoords(spe), resolution = 400, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

res <- 50
test <- bench::mark(
  class(SEraster::rasterizeSparseMatrix(assay(spe, "lognorm"), spatialCoords(spe), resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())),
  class(SEraster::rasterizeSparseMatrix2(assay(spe, "lognorm"), spatialCoords(spe), resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())),
  memory = FALSE
)

# Plot --------------------------------------------------------------------

## Figure 1a (spatial plots)
# res <- list("singlecell", 50, 100, 200, 400)
res_list <- list(200, 400)

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
    df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], transcripts = colSums(assay(spe_rast)))
    plt <- ggplot(df, aes(x = x, y = y, fill = transcripts)) +
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

## Figure 1b (runtime and memory comparison)

## Figure 1c (performance comparison)

## Figure 1d (nnSVG results correlation)
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global.RDS")))

# Further exploration -----------------------------------------------------

res <- 400
spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
dim(spe_rast)

df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], transcripts = colMeans(assay(spe_rast, "pixelval")))
ggplot(df, aes(x = x, y = y, fill = transcripts)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean transcripts/bin") +
  theme_classic()

## nnSVG
spe_rast <- nnSVG::nnSVG(
  spe_rast,
  assay_name = "pixelval",
  BPPARAM = BiocParallel::MulticoreParam()
)
View(as.data.frame(rowData(spe_rast)))
