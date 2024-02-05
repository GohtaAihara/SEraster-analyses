## This file integrates rasterization with nnSVG to analyze MERFISH mouse whole brain coronal section datasets.

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(sf)
library(ggplot2)
library(ggrastr)
library(nnSVG)
library(here)
library(tidyr)
library(tibble)
library(dplyr)

par(mfrow=c(1,1))

dataset_name <- "merfish_mousePOA"

# Load dataset ------------------------------------------------------------

# specify
animal <- 1
sex <- "Female"
behavior <- "Naive"
bregma <- "-0.29"

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_preprocessed.RDS")))

plot(spatialCoords(spe), pch=".", asp=1)


# Run method --------------------------------------------------------------

## test if each resolution can be run (it seems like 200 um works but not 400 um)
res <- 200
## SEraster
spe_rast <- SEraster::rasterizeGeneExpression(spe, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
# assays(0)

dim(spe_rast)[2] # number of spatial points

## nnSVG
spe_rast <- nnSVG::nnSVG(
  spe_rast,
  assay_name = "pixelval",
  BPPARAM = BiocParallel::MulticoreParam()
)

## Rotate dataset, rasterize, run nnSVG for each resolution
res_list <- list("singlecell", 50, 100, 200) # feel free to add more resolutions
n_rotation <- 10 # use 10 for the actual figure!
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)

## set BiocParallel parameters
# the number of workers is set to be the number of physical cores
bpparam <- BiocParallel::MulticoreParam()
# bpparam <- BiocParallel::MulticoreParam(workers = parallel::detectCores(logical = FALSE))

nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
  print(paste0("Resolution: ", res))
  if (res == "singlecell") {
    num_points = dim(spe)[2]
    
    ## nnSVG
    spe <- nnSVG::nnSVG(
      spe,
      assay_name = "lognorm",
      BPPARAM = bpparam
    )
    df <- rownames_to_column(as.data.frame(rowData(spe)), var = "gene")
    df <- cbind(rotation_deg = NA, num_points = num_points, df)
  } else {
    df <- do.call(rbind, lapply(angle_deg_list, function(deg) {
      print(paste0("Rotation (degrees): ", deg))
      ## rotate xy coordinates
      spe_rotated <- SpatialExperiment::SpatialExperiment(
        assays = assays(spe),
        spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
      )
      
      ## rasterization
      spe_rast <- SEraster::rasterizeGeneExpression(spe_rotated, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = bpparam)
      num_points = dim(spe_rast)[2]
      
      ## nnSVG (might want to use try(), see "run_nnSVG_nnSVG_simulations.R")
      spe_rast <- nnSVG::nnSVG(
        spe_rast,
        assay_name = "pixelval",
        BPPARAM = bpparam
      )
      temp <- rownames_to_column(as.data.frame(rowData(spe_rast)), var = "gene")
      temp <- cbind(rotation_deg = deg, num_points = num_points, temp)
      return(temp)
    }))
  }
  return(data.frame(dataset = dataset_name, resolution = res, df))
}))
saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))

# Plot --------------------------------------------------------------------

## Figure x (performance comparison)
n_rotation <- 1
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
  select(resolution, rotation_deg, TPR, specificity, PPV) %>%
  pivot_longer(!c(resolution, rotation_deg), names_to = "metrics", values_to = "values")

df_perf_summary <- do.call(rbind, lapply(unique(df_perf$resolution), function(res) {
  out <- do.call(rbind, lapply(c("TPR", "specificity", "PPV"), function(metric) {
    temp <- df_perf[df_perf$resolution == res, metric]
    return(data.frame(metrics = metric, mean = mean(temp), sd = sd(temp)))
  }))
  return(data.frame(resolution = as.numeric(res), out))
}))

df_perf2 <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, TPR, specificity, PPV) %>%
  pivot_longer(!resolution, names_to = "metrics", values_to = "values")

ggplot(df_perf2, aes(x = resolution, y = values, col = metrics)) +
  # geom_jitter(width = 10, alpha = 0.3) +
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


# Questions ---------------------------------------------------------------

# Is the performance consistent across bregma sections for each mouse OR across mice for each bregma section?
# How does the correlation of nnSVG outputs look like single-cell vs. rasterized resolutions (e.g. gene ranking)?
