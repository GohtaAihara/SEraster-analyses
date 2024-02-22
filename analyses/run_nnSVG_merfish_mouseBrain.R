## This file integrates rasterization with nnSVG to analyze MERFISH mouse whole brain coronal section datasets. Produce figures for 

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::document()
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
# library(sf)
library(ggplot2)
library(ggrastr)
library(nnSVG)
library(here)
library(tidyr)
library(tibble)
library(dplyr)
library(pryr)
library(R.utils)
library(Seurat)

par(mfrow=c(1,1))

dataset_name <- "merfish_mouseBrain"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_preprocessed.RDS")))

plot(spatialCoords(spe), pch=".", asp=1)


# Run method --------------------------------------------------------------

res_list <- list("singlecell", 50, 100, 200, 400)
res_list
# res_list <- list(50, 100, 200, 400)

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

## Rotate dataset, rasterize, run nnSVG for each resolution
start_res <- 50
end_res <- 400
interval_res <- 10
res_list <- c("singlecell", as.list(seq(start_res, end_res, by = interval_res)))
# res_list <- c(as.list(seq(start_res, end_res, by = interval_res)))
res_list

n_rotation <- 10
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)
# angle_deg_list <- c(216, 252, 288, 324)

nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
  print(paste0("Resolution: ", res))
  if (res == "singlecell") {
    
    ## nnSVG
    spe <- nnSVG::nnSVG(
      spe,
      assay_name = "lognorm",
      BPPARAM = BiocParallel::MulticoreParam()
    )
    df <- rownames_to_column(as.data.frame(rowData(spe)), var = "gene")
    df <- cbind(rotation_deg = NA, num_points = dim(spe)[2], df)
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
      
      ## nnSVG
      ## using try() to handle error and withTimeout() to handle > 30 mins runtime
      # spe_rast <- try({
      #   withTimeout(
      #     {
      #       nnSVG::nnSVG(
      #         spe_rast,
      #         assay_name = "pixelval",
      #         BPPARAM = BiocParallel::MulticoreParam()
      #       )
      #     },
      #     timeout = 1800,
      #     onTimeout = c("error")
      #   )
      # })
      
      ## using try() to handle error
      spe_rast <- try({
        nnSVG::nnSVG(
          spe_rast,
          assay_name = "pixelval",
          BPPARAM = BiocParallel::MulticoreParam()
        )
      })
      
      # spe_rast <- nnSVG::nnSVG(
      #   spe_rast,
      #   assay_name = "pixelval",
      #   BPPARAM = BiocParallel::MulticoreParam()
      # )
      
      if (class(spe_rast) == "try-error") {
        ## do not save anything for clusters that caused error in nnSVG
        return(NULL)
        
      } else {
        temp <- rownames_to_column(as.data.frame(rowData(spe_rast)), var = "gene")
        temp <- cbind(rotation_deg = deg, num_points = dim(spe_rast)[2], temp)
        return(temp)
      }
    }))
  }
  return(data.frame(dataset = dataset_name, resolution = res, df))
}))
# saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global_", "n_rotation_", n_rotation, "_", start_res, "-", end_res, "-by-", interval_res, ".RDS")))

## Measure runtime for each resolution
device <- "macbookpro"
n_itr <- 5

if (device == "macbookpro") {
  bpparam <- BiocParallel::MulticoreParam(workers = parallel::detectCores(logical = FALSE))
  bpparam
}

runtime_results <- do.call(rbind, lapply(res_list, function(res) {
  out <- do.call(rbind, lapply(seq(n_itr), function(i) {
    print(paste0("Resolution: ", res, ", trial: ", i))
    if (res == "singlecell") {
      num_points = dim(spe)[2]
      
      start <- Sys.time()
      nnSVG::nnSVG(
        spe,
        assay_name = "lognorm",
        BPPARAM = bpparam
      )
      runtime <- difftime(Sys.time(), start, units = "secs")
      return(data.frame(trial = i, num_points = num_points, runtime_rast = NA, runtime_nnsvg = runtime, runtime_total = runtime))
    } else {
      start1 <- Sys.time()
      spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = bpparam)
      runtime_rast <- difftime(Sys.time(), start1, units = "secs")
      
      num_points = dim(spe_rast)[2]
      
      start2 <- Sys.time()
      nnSVG::nnSVG(
        spe_rast,
        assay_name = "pixelval",
        BPPARAM = bpparam
      )
      end <- Sys.time()
      runtime_nnsvg <- difftime(end, start2, units = "secs")
      runtime_total <- difftime(end, start1, units = "secs")
      return(data.frame(trial = i, num_points = num_points, runtime_rast = runtime_rast, runtime_nnsvg = runtime_nnsvg, runtime_total = runtime_total))
    }
  }))
  return(data.frame(dataset = dataset_name, resolution = res, out))
}))
saveRDS(runtime_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global_runtime_", device, ".RDS")))


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

## test
library(pryr)

test <- mem_change(
  nnSVG::nnSVG(
    spe,
    assay_name = "lognorm",
    n_threads = 1
  )
)
test

data("merfish_mousePOA")

start <- Sys.time()
mem <- mem_change(
  rast <- SEraster::rasterizeGeneExpression(merfish_mousePOA)
)
runtime <- difftime(Sys.time(), start, units = "secs")
mem
runtime
rast

mem_change(
  nnSVG::nnSVG(
    merfish_mousePOA,
    assay_name = "volnorm",
    n_threads = 1
  )
)

library(profvis)
profvis(
  nnSVG::nnSVG(
    merfish_mousePOA,
    assay_name = "volnorm",
    n_threads = 1
  )
)

setwd("~/Desktop/SEraster")
source("R/main.R")

profvis(rasterizeGeneExpression(merfish_mousePOA))

resolution = 100
data <- assay(merfish_mousePOA)
pos <- spatialCoords(merfish_mousePOA)
bbox <- sf::st_bbox(c(
  xmin = floor(min(pos[,1])-resolution/2), 
  xmax = ceiling(max(pos[,1])+resolution/2), 
  ymin = floor(min(pos[,2])-resolution/2), 
  ymax = ceiling(max(pos[,2])+resolution/2)
))

profvis(rasterizeMatrix(data, pos, bbox))

setwd("~/Desktop/SEraster-analyses/")

## Measure runtime (Sys.time) and memory (pryr::mem_change)

res_list <- list("singlecell", 50, 100, 200, 400)

device <- "MacStudio"
n_cores <- 1

if (device == "MacStudio") {
  bpparam <- BiocParallel::MulticoreParam(workers = n_cores)
  bpparam
}

n_itr <- 5

resource_results <- do.call(rbind, lapply(res_list, function(res) {
  out <- do.call(rbind, lapply(seq(n_itr), function(i) {
    print(paste0("Resolution: ", res, ", trial: ", i))
    if (res == "singlecell") {
      num_points = dim(spe)[2]
      
      start <- Sys.time()
      mem <- mem_change(
        spe <- nnSVG::nnSVG(
          spe,
          assay_name = "lognorm",
          BPPARAM = bpparam
        )
      )
      runtime <- difftime(Sys.time(), start, units = "secs")
      return(data.frame(trial = i, num_points = num_points, runtime_rast = NA, runtime_nnsvg = runtime, runtime_total = runtime, mem_rast = NA, mem_nnsvg = mem, mem_total = mem))
    } else {
      start1 <- Sys.time()
      mem_rast <- mem_change(
        spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = bpparam)
      )
      runtime_rast <- difftime(Sys.time(), start1, units = "secs")
      
      num_points = dim(spe_rast)[2]
      
      start2 <- Sys.time()
      mem_nnsvg <- mem_change(
        spe_rast <- nnSVG::nnSVG(
          spe_rast,
          assay_name = "pixelval",
          BPPARAM = bpparam
        )
      )
      end <- Sys.time()
      runtime_nnsvg <- difftime(end, start2, units = "secs")
      runtime_total <- difftime(end, start1, units = "secs")
      return(data.frame(trial = i, num_points = num_points, runtime_rast = runtime_rast, runtime_nnsvg = runtime_nnsvg, runtime_total = runtime_total, mem_rast = mem_rast, mem_nnsvg = mem_nnsvg, mem_total = mem_rast + mem_nnsvg))
    }
  }))
  return(data.frame(dataset = dataset_name, resolution = res, out))
}))
saveRDS(resource_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global_runtime_memory_", device, "_n=", n_itr, ".RDS")))

## compare the performance with other preprocessing methods (uniform sampling, geometric sketch)
## for each method, subset to the number of pixels for each resolution

# get number of pixels for each resolution
# v1
# res_list <- list(50, 100, 200, 400)
# v2
start_res <- 50
end_res <- 400
interval_res <- 10
res_list <- as.list(seq(start_res, end_res, by = interval_res))

n_rotation <- 10
# v1
# df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
# v2
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_", "n_rotation_", n_rotation, "_", start_res, "-", end_res, "-by-", interval_res, ".RDS")))

deg <- 0

num_pixels <- do.call(rbind, lapply(unique(df$resolution), function(res) {
  if (res != "singlecell") {
    df_sub <- df[df$resolution == res,]
    return(data.frame(resolution = res, num_pixels = round(mean(unique(df_sub$num_points)))))
  }
}))
num_pixels

seeds <- seq(1,10)

# geometric sketching
nnsvg_results_sketch <- do.call(rbind, lapply(seq_along(res_list), function(i) {
  out <- do.call(rbind, lapply(seeds, function(seed) {
    print(paste0("Resolution: ", res_list[[i]], ", # of pixels: ", num_pixels[i,2], ", Seed: ", seed))
    
    # create a Seurat object from raw counts --> "counts" layer in RNA assay
    obj <- Seurat::CreateSeuratObject(assay(spe, "counts"))
    # add log-normalized values as "data" layer in RNA assay
    obj[["RNA"]]$data <- assay(spe, "lognorm")
    # find highly variable genes in the dataset (all genes in this case)
    obj <- Seurat::FindVariableFeatures(obj)

    # perform geometry sketching
    obj <- Seurat::SketchData(
      object = obj,
      assay = "RNA",
      ncells = num_pixels[i,2],
      method = "LeverageScore",
      sketched.assay = "sketch",
      seed = seed
    )

    # get cell IDs from sampled dataset
    selected_cells <- colnames(obj[["sketch"]]$counts)

    # subset the original SpatialExperiment data
    spe_sketch <- spe[,colnames(spe) %in% selected_cells]
    stopifnot(dim(spe_sketch)[2] == num_pixels[i,2])

    # sanity check
    pos_sketch <- spatialCoords(spe_sketch)
    meta_sketch <- colData(spe_sketch)
    gexp_sketch <- assay(spe_sketch, "lognorm")
    stopifnot(identical(rownames(pos_sketch), colnames(gexp_sketch)))
    stopifnot(identical(rownames(meta_sketch), colnames(gexp_sketch)))
    
    # perform nnSVG
    # using try() to handle error
    spe_sketch <- try({
      nnSVG::nnSVG(
        spe_sketch,
        assay_name = "lognorm",
        BPPARAM = BiocParallel::MulticoreParam()
      )
    })
    
    if (class(spe_sketch) == "try-error") {
      ## do not save anything for clusters that caused error in nnSVG
      return(NULL)
    } else {
      temp <- rownames_to_column(as.data.frame(rowData(spe_sketch)), var = "gene")
      temp <- cbind(seed = seed, num_points = dim(spe_sketch)[2], temp)
      return(temp)
    }
  }))
  
  if (is.null(out)) {
    print(paste0("All conditions for resolution ", res_list[[i]], " failed"))
    return(NULL)
  } else {
    return(data.frame(dataset = dataset_name, resolution = res_list[[i]], method = "geometric_sketching", out))
  }
  
}))
saveRDS(nnsvg_results_sketch, file = here("outputs", paste0(dataset_name, "_nnsvg_global_geometric_sketching_v2.RDS")))

# uniform sampling
nnsvg_results_uniform <- do.call(rbind, lapply(seq_along(res_list), function(i) {
  out <- do.call(rbind, lapply(seeds, function(seed) {
    print(paste0("Resolution: ", res_list[[i]], ", # of pixels: ", num_pixels[i,2], ", Seed: ", seed))
    
    # randomly sample cell IDs from the original SpatialExperiment data
    set.seed(seed)
    selected_cells <- sample(colnames(spe), num_pixels[i,2])
    
    # subset the original SpatialExperiment data
    spe_uniform <- spe[,colnames(spe) %in% selected_cells]
    stopifnot(dim(spe_uniform)[2] == num_pixels[i,2])
    
    # sanity check
    pos_uniform <- spatialCoords(spe_uniform)
    meta_uniform <- colData(spe_uniform)
    gexp_uniform <- assay(spe_uniform, "lognorm")
    stopifnot(identical(rownames(pos_uniform), colnames(gexp_uniform)))
    stopifnot(identical(rownames(meta_uniform), colnames(gexp_uniform)))
    
    # perform nnSVG
    # using try() to handle error
    spe_uniform <- try({
      nnSVG::nnSVG(
        spe_uniform,
        assay_name = "lognorm",
        BPPARAM = BiocParallel::MulticoreParam()
      )
    })
    
    if (class(spe_uniform) == "try-error") {
      ## do not save anything for clusters that caused error in nnSVG
      return(NULL)
    } else {
      temp <- rownames_to_column(as.data.frame(rowData(spe_uniform)), var = "gene")
      temp <- cbind(seed = seed, num_points = dim(spe_uniform)[2], temp)
      return(temp)
    }
  }))
  
  if (is.null(out)) {
    print(paste0("All conditions for resolution ", res_list[[i]], " failed"))
    return(NULL)
  } else {
    return(data.frame(dataset = dataset_name, resolution = res_list[[i]], method = "uniform", out))
  }
  
}))
saveRDS(nnsvg_results_uniform, file = here("outputs", paste0(dataset_name, "_nnsvg_global_uniform_v2.RDS")))

# Plot --------------------------------------------------------------------

col_celltype <- gg_color_hue(length(levels(colData(spe)$celltype)))
names(col_celltype) <- levels(colData(spe)$celltype)

## Figure 1 (schematics)
selected_genes <- c("Baiap2", "Slc17a6", "Gpr151")
selected_celltypes <- c("Excitatory Neurons", "Inhibitory interneurons", "Oligodendrocytes")
resolution <- 100
spe_rast_gexp <- SEraster::rasterizeGeneExpression(spe, assay_name = "counts", resolution = resolution, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
spe_rast_ct <- SEraster::rasterizeCellType(spe, col_name = "celltype", resolution = resolution, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())

## create bbox
pos <- spatialCoords(spe)
bbox <- sf::st_bbox(c(
  xmin = floor(min(pos[,1])-resolution/2), 
  xmax = ceiling(max(pos[,1])+resolution/2), 
  ymin = floor(min(pos[,2])-resolution/2), 
  ymax = ceiling(max(pos[,2])+resolution/2)
))

## create grid for rasterization
grid <- sf::st_make_grid(bbox, cellsize = resolution)
grid_coord <- st_coordinates(grid)

## single cell total transcripts
df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], transcripts = colSums(assay(spe, "counts")), celltype = colData(spe)$celltype)
ggplot(df, aes(x = x, y = y, col = transcripts)) +
  coord_fixed() +
  rasterise(geom_point(size = 1, stroke = 0), dpi = 300) +
  scale_color_viridis_c() +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_sc_global_tot_counts.pdf")))

## single cell celltypes
ggplot(df, aes(x = x, y = y, col = celltype)) +
  coord_fixed() +
  rasterise(geom_point(size = 1, stroke = 0), dpi = 300) +
  scale_color_manual(name = "Cell type", values = col_celltype) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(
    legend.position="right",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_sc_global_celltype_with_legends.pdf")))

ggplot(df, aes(x = x, y = y, col = celltype)) +
  coord_fixed() +
  rasterise(geom_point(size = 1, stroke = 0), dpi = 300) +
  scale_color_manual(name = "Cell type", values = col_celltype) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_sc_global_celltype_without_legends.pdf")))

## global gene expression
for (gene in selected_genes) {
  ## plot single cell
  df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], gene = assay(spe, "counts")[gene,])
  ggplot(df, aes(x = x, y = y, col = gene)) +
    coord_fixed() +
    rasterise(geom_point(size = 1, stroke = 0), dpi = 300) +
    scale_color_viridis_c() +
    geom_hline(yintercept = grid_coord[,2], linetype = "solid", color = "gray") +
    geom_vline(xintercept = grid_coord[,1], linetype = "solid",color = "gray") +
    theme_bw() +
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_sc_global_gene_", gene, ".pdf")))
  
  ## plot rasterized
  df <- data.frame(x = spatialCoords(spe_rast_gexp)[,1], y = spatialCoords(spe_rast_gexp)[,2], gene = assay(spe_rast_gexp, "pixelval")[gene,])
  ggplot(df, aes(x = x, y = y, fill = gene)) +
    coord_fixed() +
    geom_tile(color = "gray") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_rast_global_gene_", gene, ".pdf")))
}

## cell type-specific gene expression
## subset by cell type
ct_label <- selected_celltypes[[1]]
spe_sub <- spe[,spe$celltype %in% ct_label]
spe_sub_rast_gexp <- SEraster::rasterizeGeneExpression(spe_sub, assay_name = "counts", resolution = resolution, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
for (gene in selected_genes) {
  ## plot single cell
  df <- data.frame(x = spatialCoords(spe_sub)[,1], y = spatialCoords(spe_sub)[,2], gene = assay(spe_sub, "counts")[gene,])
  ggplot(df, aes(x = x, y = y, col = gene)) +
    coord_fixed() +
    rasterise(geom_point(size = 1, stroke = 0), dpi = 300) +
    scale_color_viridis_c() +
    geom_hline(yintercept = grid_coord[,2], linetype = "solid", color = "gray") +
    geom_vline(xintercept = grid_coord[,1], linetype = "solid",color = "gray") +
    theme_bw() +
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_sc_celltype_", ct_label, "_gene_", gene, ".pdf")))
  
  ## plot rasterized
  df <- data.frame(x = spatialCoords(spe_sub_rast_gexp)[,1], y = spatialCoords(spe_sub_rast_gexp)[,2], gene = assay(spe_sub_rast_gexp, "pixelval")[gene,])
  ggplot(df, aes(x = x, y = y, fill = gene)) +
    coord_fixed() +
    geom_tile(color = "gray") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_rast_celltype_", ct_label, "_gene_", gene, ".pdf")))
}

## cell type
for (ct_label in selected_celltypes) {
  ## plot single cell
  ## subset by cell type
  spe_sub <- spe[,spe$celltype %in% ct_label]
  ggplot(data.frame(spatialCoords(spe)), aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 1, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = data.frame(spatialCoords(spe_sub)), aes(x = x, y = y), color = col_celltype[[ct_label]], size = 1, stroke = 0), dpi = 300) +
    geom_hline(yintercept = grid_coord[,2], linetype = "solid", color = "gray") +
    geom_vline(xintercept = grid_coord[,1], linetype = "solid",color = "gray") +
    theme_bw() +
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_sc_celltype_", ct_label, ".pdf")))
  
  ## plot rasterized
  df <- data.frame(x = spatialCoords(spe_rast_ct)[,1], y = spatialCoords(spe_rast_ct)[,2], celltype = assay(spe_rast_ct, "pixelval")[ct_label,])
  ggplot(df, aes(x = x, y = y, fill = celltype)) +
    coord_fixed() +
    geom_tile(color = "gray") +
    scale_fill_viridis_c(option = "inferno") +
    theme_bw() +
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_rast_celltype_", ct_label, ".pdf")))
}

## plot inferno legend
ggplot(df, aes(x = x, y = y, fill = celltype)) +
  coord_fixed() +
  geom_tile(color = "gray") +
  scale_fill_viridis_c(option = "inferno") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
ggsave(file = here("plots", dataset_name, paste0(dataset_name, "_schematic_rast_global_celltype_", ct_label, "with_legend.pdf")))

## plot various resolutions
res_list <- list(50, 100, 200, 400)

for (res in res_list) {
  ## rasterize
  spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "counts", resolution = res, fun = "sum", BPPARAM = BiocParallel::MulticoreParam())
  
  ## plot
  df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], transcripts = colSums(assay(spe_rast)))
  plt <- ggplot(df, aes(x = x, y = y, fill = transcripts)) +
    coord_fixed() +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  ## save plot
  ggsave(plot = plt, filename = here("plots", dataset_name, paste0(dataset_name, "_schematic_rast_resolution_", res, ".pdf")))
}

## Figure 2a (spatial plots)
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
      labs(title = paste0("Resolution = ", res, " um")) +
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
  ggsave(plot = plt, filename = here("plots", dataset_name, paste0(dataset_name, "_tot_counts_", res, ".pdf")))
}


## Figure 2c (runtime and memory comparison)
# df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_runtime.RDS")))
# df$resolution <- factor(df$resolution, levels = c("singlecell", "50", "100", "200", "400"))
# 
# col_res <- c("#666666", gg_color_hue(4))
# 
# # compute fold change for total time across resolution
# df_summary <- df %>%
#   group_by(resolution) %>%
#   summarize(
#     avg_runtime_rast = mean(runtime_rast, na.rm = TRUE),
#     avg_runtime_nnsvg = mean(runtime_nnsvg, na.rm = TRUE),
#     avg_runtime_total = mean(runtime_total, na.rm = TRUE)
#   )
# singlecell_time <- df_summary %>%
#   filter(resolution == "singlecell") %>%
#   select(avg_runtime_total) %>%
#   unlist()
# df_summary <- df_summary %>%
#   mutate(
#     avg_runtime_total = as.numeric(avg_runtime_total),
#     fold_change = singlecell_time / avg_runtime_total
#   )
# 
# # total runtime
# ggplot(df, aes(x = num_points, y = as.numeric(runtime_total), col = resolution)) +
#   geom_boxplot(width = 6000, lwd = 0.75, outlier.shape = NA) +
#   geom_jitter(width = 1000, size = 1, alpha = 0.75) +
#   scale_x_continuous(breaks = unique(df$num_points)) + 
#   scale_color_manual(values = col_res) +
#   labs(title = "Total runtime",
#        x = "Number of spatial points",
#        y = "Runtime (secs)",
#        col = "Resolution") +
#   theme_bw() +
#   theme(panel.grid.minor.x = element_blank(), 
#         panel.grid.minor.y = element_blank(), 
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_total.pdf")), width = 6, heigh = 5, dpi = 300)
# 
# ggplot(df, aes(x = resolution, y = as.numeric(runtime_total), col = resolution)) +
#   geom_boxplot(lwd = 0.5, outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.75) +
#   scale_color_manual(values = col_res) +
#   labs(title = "Total runtime",
#        x = "Resolution",
#        y = "Runtime (secs)",
#        col = "Resolution") +
#   theme_bw() +
#   theme(panel.grid.minor.x = element_blank(), 
#         panel.grid.minor.y = element_blank(),
#         legend.position = "none")
# ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_total_v2.pdf")), width = 4, height = 8, dpi = 300)
# 
# # nnSVG runtime
# ggplot(df, aes(x = num_points, y = as.numeric(runtime_nnsvg), col = resolution)) +
#   geom_boxplot(width = 6000, lwd = 0.75, outlier.shape = NA) +
#   geom_jitter(width = 1000, size = 1, alpha = 0.75) +
#   scale_x_continuous(breaks = unique(df$num_points)) + 
#   scale_color_manual(values = col_res) +
#   labs(title = "nnSVG runtime",
#        x = "Number of spatial points",
#        y = "Runtime (secs)",
#        col = "Resolution") +
#   theme_bw() +
#   theme(panel.grid.minor.x = element_blank(), 
#         panel.grid.minor.y = element_blank(), 
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_nnsvg.pdf")), width = 6, heigh = 5, dpi = 300)
# 
# ggplot(df, aes(x = resolution, y = as.numeric(runtime_nnsvg), col = resolution)) +
#   geom_boxplot(lwd = 0.75, outlier.shape = NA) +
#   geom_jitter(width = 0.1, alpha = 0.75) +
#   scale_color_manual(values = col_res) +
#   labs(title = "nnSVG runtime",
#        x = "Resolution",
#        y = "Runtime (secs)",
#        col = "Resolution") +
#   theme_bw() +
#   theme(panel.grid.minor.x = element_blank(), 
#         panel.grid.minor.y = element_blank(),
#         legend.position = "none")
# ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_nnsvg_v2.pdf")), width = 5, heigh = 5, dpi = 300)
# 
# # SEraster runtime
# ggplot(df[df$resolution != "singlecell",], aes(x = num_points, y = as.numeric(runtime_rast), col = resolution)) +
#   geom_boxplot(width = 2000, lwd = 0.75, outlier.shape = NA) +
#   geom_jitter(width = 500, size = 1, alpha = 0.75) +
#   scale_x_continuous(breaks = unique(df$num_points)) + 
#   scale_color_manual(values = col_res) +
#   labs(title = "SEraster runtime",
#        x = "Number of spatial points",
#        y = "Runtime (secs)",
#        col = "Resolution") +
#   theme_bw() +
#   theme(panel.grid.minor.x = element_blank(), 
#         panel.grid.minor.y = element_blank(), 
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_rast.pdf")), width = 6, heigh = 5, dpi = 300)
# 
# ggplot(df[df$resolution != "singlecell",], aes(x = resolution, y = as.numeric(runtime_rast), col = resolution)) +
#   geom_boxplot(lwd = 0.75, outlier.shape = NA) +
#   geom_jitter(width = 0.1, alpha = 0.75) +
#   scale_color_manual(values = col_res) +
#   labs(title = "SEraster runtime",
#        x = "Resolution",
#        y = "Runtime (secs)",
#        col = "Resolution") +
#   theme_bw() +
#   theme(panel.grid.minor.x = element_blank(), 
#         panel.grid.minor.y = element_blank(),
#         legend.position = "none")
# ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_rast_v2.pdf")), width = 5, heigh = 5, dpi = 300)

## Figure 2c new (runtime and memory comparison)
device <- "MacStudio"
n_itr <- 5
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_runtime_memory_", device, "_n=", n_itr, ".RDS")))

col_res <- c("#666666", gg_color_hue(4))

# compute fold change for total time across resolution
df_summary <- df %>%
  group_by(resolution) %>%
  summarize(
    avg_runtime_rast = mean(runtime_rast, na.rm = TRUE),
    avg_runtime_nnsvg = mean(runtime_nnsvg, na.rm = TRUE),
    avg_runtime_total = mean(runtime_total, na.rm = TRUE),
    avg_mem_rast = mean(mem_rast, na.rm = TRUE),
    avg_mem_nnsvg = mean(mem_nnsvg, na.rm = TRUE),
    avg_mem_total = mean(mem_total, na.rm = TRUE)
  )
singlecell_time <- df_summary %>%
  filter(resolution == "singlecell") %>%
  select(avg_runtime_total) %>%
  unlist()
df_summary <- df_summary %>%
  mutate(
    avg_runtime_total = as.numeric(avg_runtime_total),
    percentage = avg_runtime_total / singlecell_time * 100
  )

df$resolution <- factor(df$resolution, levels = c("singlecell", 50, 100, 200, 400))

# total runtime
set.seed(0)
ggplot(df, aes(x = resolution, y = as.numeric(runtime_total)/60, col = resolution)) +
  geom_boxplot(lwd = 0.5, outlier.shape = NA, linewidth = 1) +
  geom_jitter(width = 0.3, alpha = 0.75, size = 2.5, stroke = 0) +
  geom_hline(yintercept = df_summary[df_summary$resolution == "singlecell",]$avg_runtime_total/60, color = "black", linetype = "dashed") +
  geom_text(data = df_summary, aes(x = resolution, y = as.numeric(avg_runtime_total)/60, label = paste0(round(percentage, digits = 1), "%")), vjust = -1.5, size = 5) +
  scale_color_manual(values = col_res) +
  # scale_y_continuous(expand = expansion(mult = 0.1)) +
  ylim(0,150) +
  labs(title = "Total runtime",
       x = "Rasterization Resolution (µm)",
       y = "Runtime (minutes)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        legend.position = "none")
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_total_", device, "_n=", n_itr, ".pdf")), width = 4, heigh = 8, dpi = 300)

# total memory
ggplot(df, aes(x = resolution, y = as.numeric(mem_total)*1e-6, col = resolution)) +
  geom_boxplot(lwd = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.75) +
  scale_color_manual(values = col_res) +
  labs(title = "Total memory",
       x = "Resolution",
       y = "Memory (MB)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        legend.position = "none")

# SEraster memory
ggplot(df, aes(x = resolution, y = as.numeric(mem_rast)*1e-6, col = resolution)) +
  geom_boxplot(lwd = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.75) +
  scale_color_manual(values = col_res) +
  labs(title = "SEraster memory",
       x = "Resolution",
       y = "Memory (MB)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        legend.position = "none")

# nnSVG memory
ggplot(df, aes(x = resolution, y = as.numeric(mem_nnsvg)*1e-6, col = resolution)) +
  geom_boxplot(lwd = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.75) +
  scale_color_manual(values = col_res) +
  labs(title = "nnSVG memory",
       x = "Resolution",
       y = "Memory (MB)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        legend.position = "none")

## Figure 2d (performance comparison)
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global.RDS")))
# set a threshold p value
alpha <- 0.05
df_perf <- do.call(rbind, lapply(unique(df$resolution), function(res) {
  if (res != "singlecell") {
    sc <- df[df$resolution == "singlecell",]
    rast <- df[df$resolution == res,]
    results_sig <- do.call(rbind, lapply(rast$gene, function(gene) {
      return(data.frame(gene = gene, pred = rast[rast$gene == gene, "padj"] <= alpha, obs = sc[sc$gene == gene, "padj"] <= alpha))
    }))
    out <- calculatePerformanceMetrics(results_sig)
    return(data.frame(resolution = res, out))
  }
}))
df_perf1 <- df_perf %>%
  mutate(resolution = factor(resolution, levels = c("singlecell", "50", "100", "200", "400"))) %>%
  select(resolution, TPR, FPR, PPV, F1, ACC) %>%
  pivot_longer(!resolution, names_to = "metrics", values_to = "values")

ggplot(df_perf1, aes(x = resolution, y = values, fill = resolution)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Performance metrics comparison",
       x = "Resolution",
       y = "Value",
       fill = "Resolution") +
  facet_wrap(~metrics) +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_perf_metric.pdf")), width = 6, heigh = 5, dpi = 300)

df_perf2 <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, TPR, TNR, PPV, F1, ACC) %>%
  pivot_longer(!resolution, names_to = "metrics", values_to = "values")

ggplot(df_perf2, aes(x = resolution, y = values, col = metrics)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = unique(df_perf2$resolution)) + 
  # ylim(0,1) +
  labs(title = "Performance metrics comparison",
       x = "Resolution",
       y = "Performance",
       col = "Metric") +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_perf_metric_v2.pdf")), width = 6, heigh = 5, dpi = 300)

## performance comparison with rotations
n_rotation <- 10
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
# set a threshold p value
alpha <- 0.05
df_perf <- do.call(rbind, lapply(unique(df$resolution), function(res) {
  if (res != "singlecell") {
    sc <- df[df$resolution == "singlecell",]
    out <- do.call(rbind, lapply(angle_deg_list, function(deg) {
      rast <- df[df$resolution == res & df$rotation_deg == deg,]
      results_sig <- do.call(rbind, lapply(rast$gene, function(gene) {
        return(data.frame(gene = gene, pred = rast[rast$gene == gene, "padj"] <= alpha, obs = sc[sc$gene == gene, "padj"] <= alpha))
      }))
      out <- calculatePerformanceMetrics(results_sig)
      return(data.frame(rotation_deg = deg, out))
    }))
    return(data.frame(resolution = res, out))
  }
}))

df_perf_raw <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, rotation_deg, TPR, TNR, PPV) %>%
  pivot_longer(!c(resolution, rotation_deg), names_to = "metrics", values_to = "values")

df_perf_summary <- do.call(rbind, lapply(unique(df_perf$resolution), function(res) {
  out <- do.call(rbind, lapply(c("TPR", "TNR", "PPV"), function(metric) {
    temp <- df_perf[df_perf$resolution == res, metric]
    return(data.frame(metrics = metric, mean = mean(temp), sd = sd(temp)))
  }))
  return(data.frame(resolution = as.numeric(res), out))
}))

df_perf2 <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, TPR, TNR, PPV) %>%
  pivot_longer(!resolution, names_to = "metrics", values_to = "values")

set.seed(0)
ggplot(df_perf2, aes(x = resolution, y = values, col = metrics)) +
  geom_jitter(width = 10, alpha = 0.3, size = 2, stroke = 0) +
  geom_line(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics)) +
  geom_point(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics), size = 1) +
  geom_errorbar(data = df_perf_summary, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics), width = 10) +
  scale_x_continuous(breaks = unique(df_perf2$resolution)) + 
  ylim(0,1) +
  labs(title = "Performance",
       x = "Rasterization Resolution",
       y = "Performance",
       col = "Metric") +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_perf_metric_summary.pdf")), width = 6, heigh = 5, dpi = 300)

# plot TP, FP, TN, FN
df_perf3 <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, TP, FP, TN, FN)
df_perf3$resolution_label <- factor(df_perf3$resolution, levels = c(50, 100, 200, 400))

col_res <- c(gg_color_hue(4))

p1 <- ggplot(df_perf3, aes(x = resolution, y = TP, col = resolution_label)) +
  geom_boxplot(aes(group = resolution)) +
  geom_jitter() +
  scale_x_continuous(breaks = unique(df_perf3$resolution)) + 
  scale_color_manual(name = "Resolution", values = col_res) +
  labs(x = "Rasterization Resolution (μm)",
       y = "TP") +
  theme_bw() +
  theme(legend.position = "none")
p2 <- ggplot(df_perf3, aes(x = resolution, y = FP, col = resolution_label)) +
  geom_boxplot(aes(group = resolution)) +
  geom_jitter() +
  scale_x_continuous(breaks = unique(df_perf3$resolution)) + 
  scale_color_manual(name = "Resolution", values = col_res) +
  labs(x = "Rasterization Resolution (μm)",
       y = "FP") +
  theme_bw() +
  theme(legend.position = "none")
p3 <- ggplot(df_perf3, aes(x = resolution, y = TN, col = resolution_label)) +
  geom_boxplot(aes(group = resolution)) +
  geom_jitter() +
  scale_x_continuous(breaks = unique(df_perf3$resolution)) + 
  scale_color_manual(name = "Resolution", values = col_res) +
  labs(x = "Rasterization Resolution (μm)",
       y = "TN") +
  theme_bw() +
  theme(legend.position = "none")
p4 <- ggplot(df_perf3, aes(x = resolution, y = FN, col = resolution_label)) +
  geom_boxplot(aes(group = resolution)) +
  geom_jitter() +
  scale_x_continuous(breaks = unique(df_perf3$resolution)) + 
  scale_color_manual(name = "Resolution", values = col_res) +
  labs(x = "Rasterization Resolution (μm)",
       y = "FN") +
  theme_bw() +
  theme(legend.position = "none")
grid.arrange(p1, p2, p3, p4, ncol = 2)

# genes that contribute to the difference
alpha <- 0.05
# get svgs, non-svgs
svgs <- df[(df$resolution == "singlecell" & df$padj <= alpha),]$gene # 401
non_svgs <- df[(df$resolution == "singlecell" & df$padj > alpha),]$gene # 82
# label each gene as TP, TN, FP, FN
df_perf4 <- df_sub %>%
  mutate(confusion_matrix = case_when(
    gene %in% svgs & padj <= alpha ~ "TP",
    gene %in% non_svgs & padj > alpha ~ "TN",
    gene %in% non_svgs & padj <= alpha ~ "FP",
    gene %in% svgs & padj > alpha ~ "FN"
  )) %>%
  select(resolution, rotation_deg, gene, padj, confusion_matrix)
# sanity check
test <- df_perf4[df_perf4$resolution == "singlecell",]
test_svgs <- test[test$confusion_matrix %in% c("TP", "FN"),]$gene
length(test_svgs)
setdiff(test_svgs, svgs)
test_non_svgs <- test[test$confusion_matrix %in% c("FP", "TN"),]$gene
length(test_non_svgs)
setdiff(test_non_svgs, non_svgs)
# subset to relevant resolutions and permutations

##WIP

# what if we make the number of SVG and non-SVG the same? (balanced dataset)
svgs <- df[(df$resolution == "singlecell" & df$padj <= 0.05),]$gene # 401
non_svgs <- df[(df$resolution == "singlecell" & df$padj > 0.05),]$gene # 82
# randomly select 82 svgs
set.seed(0)
idx <- sample(1:length(svgs), length(non_svgs), replace = FALSE)
genes_selected <- c(svgs[idx], non_svgs)
# subset data
df_selected <- df[df$gene %in% genes_selected,]
# sanity check
table(df_selected[df_selected$resolution == "singlecell",]$padj <= 0.05)
# re-run perf analysis and plot
df_perf <- do.call(rbind, lapply(unique(df_selected$resolution), function(res) {
  if (res != "singlecell") {
    sc <- df_selected[df_selected$resolution == "singlecell",]
    out <- do.call(rbind, lapply(angle_deg_list, function(deg) {
      rast <- df_selected[df_selected$resolution == res & df_selected$rotation_deg == deg,]
      results_sig <- do.call(rbind, lapply(rast$gene, function(gene) {
        return(data.frame(gene = gene, pred = rast[rast$gene == gene, "padj"] <= alpha, obs = sc[sc$gene == gene, "padj"] <= alpha))
      }))
      out <- calculatePerformanceMetrics(results_sig)
      return(data.frame(rotation_deg = deg, out))
    }))
    return(data.frame(resolution = res, out))
  }
}))

df_perf_raw <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, rotation_deg, TPR, TNR, PPV) %>%
  pivot_longer(!c(resolution, rotation_deg), names_to = "metrics", values_to = "values")

df_perf_summary <- do.call(rbind, lapply(unique(df_perf$resolution), function(res) {
  out <- do.call(rbind, lapply(c("TPR", "TNR", "PPV"), function(metric) {
    temp <- df_perf[df_perf$resolution == res, metric]
    return(data.frame(metrics = metric, mean = mean(temp), sd = sd(temp)))
  }))
  return(data.frame(resolution = as.numeric(res), out))
}))

df_perf2 <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, TPR, TNR, PPV) %>%
  pivot_longer(!resolution, names_to = "metrics", values_to = "values")

ggplot(df_perf2, aes(x = resolution, y = values, col = metrics)) +
  geom_jitter(width = 10, alpha = 0.3) +
  geom_line(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics)) +
  geom_point(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics), size = 1) +
  geom_errorbar(data = df_perf_summary, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics), width = 10) +
  scale_x_continuous(breaks = unique(df_perf2$resolution)) + 
  ylim(0,1) +
  labs(title = "Performance",
       x = "Rasterization Resolution",
       y = "Performance",
       col = "Metric") +
  theme_bw()

# plot TP, FP, TN, FN
df_perf3 <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, TP, FP, TN, FN)
df_perf3$resolution_label <- factor(df_perf3$resolution, levels = c(50, 100, 200, 400))

col_res <- c(gg_color_hue(4))

p1 <- ggplot(df_perf3, aes(x = resolution, y = TP, col = resolution_label)) +
  geom_boxplot(aes(group = resolution)) +
  geom_jitter() +
  scale_x_continuous(breaks = unique(df_perf3$resolution)) + 
  scale_color_manual(name = "Resolution", values = col_res) +
  labs(x = "Rasterization Resolution (μm)",
       y = "TP") +
  theme_bw() +
  theme(legend.position = "none")
p2 <- ggplot(df_perf3, aes(x = resolution, y = FP, col = resolution_label)) +
  geom_boxplot(aes(group = resolution)) +
  geom_jitter() +
  scale_x_continuous(breaks = unique(df_perf3$resolution)) + 
  scale_color_manual(name = "Resolution", values = col_res) +
  labs(x = "Rasterization Resolution (μm)",
       y = "FP") +
  theme_bw() +
  theme(legend.position = "none")
p3 <- ggplot(df_perf3, aes(x = resolution, y = TN, col = resolution_label)) +
  geom_boxplot(aes(group = resolution)) +
  geom_jitter() +
  scale_x_continuous(breaks = unique(df_perf3$resolution)) + 
  scale_color_manual(name = "Resolution", values = col_res) +
  labs(x = "Rasterization Resolution (μm)",
       y = "TN") +
  theme_bw() +
  theme(legend.position = "none")
p4 <- ggplot(df_perf3, aes(x = resolution, y = FN, col = resolution_label)) +
  geom_boxplot(aes(group = resolution)) +
  geom_jitter() +
  scale_x_continuous(breaks = unique(df_perf3$resolution)) + 
  scale_color_manual(name = "Resolution", values = col_res) +
  labs(x = "Rasterization Resolution (μm)",
       y = "FN") +
  theme_bw() +
  theme(legend.position = "none")
grid.arrange(p1, p2, p3, p4, ncol = 2)

## more resolutions
start_res <- 50
end_res <- 400
interval_res <- 10
n_rotation <- 10
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_", "n_rotation_", n_rotation, "_", start_res, "-", end_res, "-by-", interval_res, ".RDS")))

# set a threshold p value
alpha <- 0.05
df_perf <- do.call(rbind, lapply(unique(df$resolution), function(res) {
  if (res != "singlecell") {
    sc <- df[df$resolution == "singlecell",]
    df_res <- df[df$resolution == res,]
    out <- do.call(rbind, lapply(unique(df_res$rotation_deg), function(deg) {
      print(paste0("Resolution = ", res, ", Angle = ", deg))
      rast <- df[df$resolution == res & df$rotation_deg == deg,]
      results_sig <- do.call(rbind, lapply(rast$gene, function(gene) {
        return(data.frame(gene = gene, pred = rast[rast$gene == gene, "padj"] <= alpha, obs = sc[sc$gene == gene, "padj"] <= alpha))
      }))
      out <- calculatePerformanceMetrics(results_sig)
      return(data.frame(rotation_deg = deg, out))
    }))
    return(data.frame(resolution = res, out))
  }
}))

df_perf_raw <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, rotation_deg, TPR, TNR, PPV) %>%
  pivot_longer(!c(resolution, rotation_deg), names_to = "metrics", values_to = "values")

df_perf_summary <- do.call(rbind, lapply(unique(df_perf$resolution), function(res) {
  out <- do.call(rbind, lapply(c("TPR", "TNR", "PPV"), function(metric) {
    temp <- df_perf[df_perf$resolution == res, metric]
    return(data.frame(metrics = metric, mean = mean(temp), sd = sd(temp)))
  }))
  return(data.frame(resolution = as.numeric(res), out))
}))

df_perf2 <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, TPR, TNR, PPV) %>%
  pivot_longer(!resolution, names_to = "metrics", values_to = "values")

set.seed(0)
ggplot(df_perf2, aes(x = resolution, y = values, col = metrics)) +
  geom_jitter(width = 1, alpha = 0.3, size = 2, stroke = 0) +
  geom_line(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics)) +
  geom_point(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics), size = 1) +
  geom_errorbar(data = df_perf_summary, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics), width = 10) +
  scale_x_continuous(breaks = unique(df_perf2$resolution)) + 
  ylim(0,1) +
  labs(title = "Performance",
       x = "Rasterization Resolution",
       y = "Performance",
       col = "Metric") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_perf_metric_summary_n_rotation_", n_rotation, "_", start_res, "-", end_res, "-by-", interval_res, ".pdf")), width = 6, heigh = 5, dpi = 300)

## Figure 1E (nnSVG results comparison)
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global.RDS")))
df <- df %>%
  mutate(resolution = factor(resolution, levels = c("singlecell", "50", "100", "200", "400"))) %>%
  select(gene, resolution, phi, prop_sv, LR_stat, rank, padj)
df <- data.frame(df[df$resolution == "singlecell",], df[df$resolution != "singlecell",])

# compute Spearman correlations
corr_results <- do.call(rbind, lapply(unique(df$resolution.1), function(res) {
  temp <- df[df$resolution.1 == res,]
  return(data.frame(
    resolution = res,
    corr_LR_stat = cor(temp$LR_stat, temp$LR_stat.1, method = "spearman"),
    corr_rank = cor(temp$rank, temp$rank.1, method = "spearman"),
    corr_prop_sv = cor(temp$prop_sv, temp$prop_sv.1, method = "spearman"),
    corr_ls = cor(1/temp$phi, 1/temp$phi.1, method = "spearman"),
    corr_padj = cor(temp$padj, temp$padj.1, method = "spearman")
  ))
}))

# create Spearman correlation labels
text_LR_stat <- data.frame(
  x = 7500,
  y = 50000,
  resolution.1 = corr_results$resolution,
  label = paste0("cor = ", round(corr_results$corr_LR_stat, 3))
)

text_rank <- data.frame(
  x = 150,
  y = 400,
  resolution.1 = corr_results$resolution,
  label = paste0("cor = ", round(corr_results$corr_rank, 3))
)

text_prop_sv <- data.frame(
  x = 0.25,
  y = 0.75,
  resolution.1 = corr_results$resolution,
  label = paste0("cor = ", round(corr_results$corr_prop_sv, 3))
)

text_ls <- data.frame(
  x = 0,
  y = 6,
  resolution.1 = corr_results$resolution,
  label = paste0("cor = ", round(corr_results$corr_ls, 3))
)

text_padj <- data.frame(
  x = 5,
  y = 12.5,
  resolution.1 = corr_results$resolution,
  label = paste0("cor = ", round(corr_results$corr_padj, 3))
)

# LR statistic
ggplot(df, aes(x = LR_stat.1, y = LR_stat, col = resolution.1)) +
  geom_point(size = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(data = text_LR_stat, aes(x = x, y = y, label = label), size = 5, color = "black") +
  labs(title = "LR statistic single cell vs. rasterization",
       x = "LR statistic rasterization",
       y = "LR statistic single cell",
       col = "Resolution") +
  facet_wrap(~ resolution.1) +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_LR_stat.pdf")), width = 6, heigh = 5, dpi = 300)

# rank
ggplot(df, aes(x = rank.1, y = rank, col = resolution.1)) +
  geom_point(size = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(data = text_rank, aes(x = x, y = y, label = label), size = 5, color = "black") +
  labs(title = "Rank single cell vs. rasterization",
       x = "Rank rasterization",
       y = "Rank single cell",
       col = "Resolution") +
  facet_wrap(~ resolution.1) +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_rank.pdf")), width = 6, heigh = 5, dpi = 300)

# proportion of spatial variance (effect size)
ggplot(df, aes(x = prop_sv.1, y = prop_sv, col = resolution.1)) +
  geom_point(size = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(data = text_prop_sv, aes(x = x, y = y, label = label), size = 5, color = "black") +
  labs(title = "Effect size single cell vs. rasterization",
       x = "Effect size rasterization",
       y = "Effect size single cell",
       col = "Resolution") +
  facet_wrap(~ resolution.1) +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_effect_size.pdf")), width = 6, heigh = 5, dpi = 300)

# length scale
ggplot(df, aes(x = log10(1/phi.1), y = log10(1/phi), col = resolution.1)) +
  geom_point(size = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(data = text_ls, aes(x = x, y = y, label = label), size = 5, color = "black") +
  labs(title = "Length scale single cell vs. rasterization",
       x = "log10(length scale rasterization)",
       y = "log10(length scale single cell)",
       col = "Resolution") +
  facet_wrap(~ resolution.1) +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_log10_length_scale.pdf")), width = 6, heigh = 5, dpi = 300)

# adjusted p value
ggplot(df, aes(x = -log10(padj.1), y = -log10(padj), col = resolution.1)) +
  geom_point(size = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(data = text_padj, aes(x = x, y = y, label = label), size = 5, color = "black") +
  labs(title = "Adjusted p value single cell vs. rasterization",
       x = "-log10(padj) rasterization",
       y = "-log10(padj) single cell",
       col = "Resolution") +
  facet_wrap(~ resolution.1) +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_adjusted_pval.pdf")), width = 6, heigh = 5, dpi = 300)

## Supplementary Figure 1 (number of spatial points for various resolutions)
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_runtime.RDS"))) %>%
  mutate(resolution = factor(resolution, levels = c("singlecell", "50", "100", "200", "400"))) %>%
  filter(trial == 1)

col_res <- c("#666666", gg_color_hue(4))

ggplot(df, aes(x = resolution, y = num_points, col = resolution, label = num_points)) +
  geom_point() +
  geom_text(vjust = -1, size = 5) +
  scale_color_manual(values = col_res) +
  ylim(c(0,9e+4)) +
  labs(x = "Rasterization Resolution",
       y = "Number of spatial points",
       col = "Resolution") +
  theme_bw() +
  theme(
    legend.position="none"
  )
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_num_points.pdf")), width = 6, heigh = 5, dpi = 300)

## Supplementary Figure x (performance comparison between SEraster, geometric sketching, and uniform sampling)
# load SEraster results
n_rotation <- 10
df_rast <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
# load geometric sketching
df_sketch <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_geometric_sketching.RDS")))
# load uniform sampling
df_uniform <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_uniform.RDS")))

# set a threshold p value
alpha <- 0.05

# extract single-cell data from SEraster results
sc <- df_rast[df_rast$resolution == "singlecell",]
df_rast <- df_rast[df_rast$resolution != "singlecell",]

# set resolutions
res_list <- list(50, 100, 200, 400)

# set metrics
metrics <- c("TPR", "TNR", "PPV")

# compute performance metrics
# SEraster (permutation based on rotation)
perf_rast <- do.call(rbind, lapply(res_list, function(res) {
  # skip if the dataset does not contain this resolution
  if (!res %in% unique(df_rast$resolution)) {
    return(NULL)
  } else {
    out <- do.call(rbind, lapply(seq_along(unique(df_rast$rotation_deg)), function(i) {
      deg <- unique(df_rast$rotation_deg)[i]
      print(paste0("SEraster, Resolution = ", res, ", Angle = ", deg))
      sampled <- df_rast[df_rast$resolution == res & df_rast$rotation_deg == deg,]
      results_sig <- do.call(rbind, lapply(sampled$gene, function(gene) {
        return(data.frame(gene = gene, pred = sampled[sampled$gene == gene, "padj"] <= alpha, obs = sc[sc$gene == gene, "padj"] <= alpha))
      }))
      out <- calculatePerformanceMetrics(results_sig)
      return(data.frame(permutation = i, out))
    }))
    return(data.frame(method = "SEraster", resolution = res, out))
  }
}))

# geometric sketching (permutation based on seed)
perf_sketch <- do.call(rbind, lapply(res_list, function(res) {
  # skip if the dataset does not contain this resolution
  if (!res %in% unique(df_sketch$resolution)) {
    return(NULL)
  } else {
    out <- do.call(rbind, lapply(seq_along(unique(df_sketch$seed)), function(i) {
      seed <- unique(df_sketch$seed)[i]
      print(paste0("Geometric sketching, Resolution = ", res, ", Seed = ", seed))
      sampled <- df_sketch[df_sketch$resolution == res & df_sketch$seed == seed,]
      results_sig <- do.call(rbind, lapply(sampled$gene, function(gene) {
        return(data.frame(gene = gene, pred = sampled[sampled$gene == gene, "padj"] <= alpha, obs = sc[sc$gene == gene, "padj"] <= alpha))
      }))
      out <- calculatePerformanceMetrics(results_sig)
      return(data.frame(permutation = i, out))
    }))
    return(data.frame(method = "geometric_sketching", resolution = res, out))
  }
}))

# uniform sampling (permutation based on seed)
perf_uniform <- do.call(rbind, lapply(res_list, function(res) {
  # skip if the dataset does not contain this resolution
  if (!res %in% unique(df_uniform$resolution)) {
    return(NULL)
  } else {
    out <- do.call(rbind, lapply(seq_along(unique(df_uniform$seed)), function(i) {
      seed <- unique(df_uniform$seed)[i]
      print(paste0("Uniform sampling, Resolution = ", res, ", Seed = ", seed))
      sampled <- df_uniform[df_uniform$resolution == res & df_uniform$seed == seed,]
      results_sig <- do.call(rbind, lapply(sampled$gene, function(gene) {
        return(data.frame(gene = gene, pred = sampled[sampled$gene == gene, "padj"] <= alpha, obs = sc[sc$gene == gene, "padj"] <= alpha))
      }))
      out <- calculatePerformanceMetrics(results_sig)
      return(data.frame(permutation = i, out))
    }))
    return(data.frame(method = "uniform", resolution = res, out))
  }
}))

# combine
perf_comb <- rbind(perf_rast, perf_sketch, perf_uniform)

for (metric in metrics) {
  print(paste0("Plotting ", metric))
  
  df <- data.frame(perf_comb[,c("method", "resolution", "permutation")], values = perf_comb[,metric])
  
  df_summary <- df %>%
    group_by(method, resolution) %>%
    summarize(mean = mean(values), sd = sd(values))
  
  # plot
  set.seed(0)
  ggplot(df, aes(x = resolution, y = values, col = method)) +
    geom_jitter(width = 10, alpha = 0.3, size = 2, stroke = 0) +
    geom_line(data = df_summary, aes(x = resolution, y = mean, col = method)) +
    geom_point(data = df_summary, aes(x = resolution, y = mean, col = method), size = 1) +
    geom_errorbar(data = df_summary, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = method), width = 10) +
    scale_x_continuous(breaks = unique(df$resolution)) +
    ylim(0,1) +
    labs(title = metric,
         x = "Rasterization Resolution",
         y = "Performance",
         col = "Sampling Method") +
    theme_bw()
  
  # save plot
  ggsave(filename = here("plots", paste0(dataset_name, "_perf_comparison_across_sampling_methods_", metric, ".pdf")))
}

# Further exploration -----------------------------------------------------

## mean-variance relationship
resolution <- 200
spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = resolution, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
df_sc_counts <- data.frame(mean = rowMeans(assay(spe, "counts")), var = apply(assay(spe, "counts"), 1, var), dataset = "single cell (counts)")
df_sc_counts <- rownames_to_column(df_sc_counts,  var = "gene")
df_sc_lognorm <- data.frame(mean = rowMeans(assay(spe, "lognorm")), var = apply(assay(spe, "lognorm"), 1, var), dataset = "single cell (lognorm)")
df_sc_lognorm <- rownames_to_column(df_sc_lognorm,  var = "gene")
res_list <- c(50, 100, 200, 400)
df_rast <- do.call(rbind, lapply(res_list, function(res) {
  spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
  out <- data.frame(mean = rowMeans(assay(spe_rast, "pixelval")), var = apply(assay(spe_rast, "pixelval"), 1, var), dataset = paste0("resolution = ", res))
  out <- rownames_to_column(out, var = "gene")
  return(out)
}))
df_comb <- rbind(df_sc_counts, df_sc_lognorm, df_rast) %>%
  mutate(dataset = factor(dataset, levels = c("single cell (counts)", "single cell (lognorm)", "resolution = 50", "resolution = 100", "resolution = 200", "resolution = 400")))
ggplot(df_comb, aes(x = log10(mean), y = log10(var), col = dataset)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  labs(title = "Comparison of mean-variance relationship") +
  theme_classic()

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

## bugfix
n_rotation <- 20
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)
res_list <- list(50, 100, 200, 400)

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
    return(cbind(rotation_deg = NA, num_points = num_points, df))
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
      return(cbind(rotation_deg = deg, num_points = num_points, temp))
    }))
  }
  return(data.frame(dataset = dataset_name, resolution = res, df))
}))

## nnSVG failed at 270 degrees
## rotate xy coordinates
deg <- angle_deg_list[18]
deg
spe_rotated <- SpatialExperiment::SpatialExperiment(
  assays = assays(spe),
  spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
)
plot(spatialCoords(spe_rotated), pch=".", asp=1)

## rasterization
res <- 50
spe_rast <- SEraster::rasterizeGeneExpression(spe_rotated, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
num_points = dim(spe_rast)[2]
gene <- "Baiap2"
df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], gene = assay(spe_rast)[gene,])
ggplot(df, aes(x = x, y = y, fill = gene)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_classic()

## nnSVG
spe_rast <- nnSVG::nnSVG(
  spe_rast,
  assay_name = "pixelval",
  BPPARAM = BiocParallel::MulticoreParam()
)
