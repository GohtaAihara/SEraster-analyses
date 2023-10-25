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
library(tidyr)
library(tibble)
library(dplyr)

par(mfrow=c(1,1))

dir <- here("outputs", "nnSVG_simulations")

dataset_name <- "nnSVG_simulations"


# Load dataset ------------------------------------------------------------

## simulations for various bandwidth and expression
## filenames for saved simulation datasets
sim_names <- c(
  "sim_largeBandwidth_fullExpr",
  "sim_largeBandwidth_medExpr",
  "sim_largeBandwidth_lowExpr",
  "sim_medBandwidth_fullExpr",
  "sim_medBandwidth_medExpr",
  "sim_medBandwidth_lowExpr",
  "sim_smallBandwidth_fullExpr",
  "sim_smallBandwidth_medExpr",
  "sim_smallBandwidth_lowExpr"
)

# sim_names <- c(
#   "sim_largeBandwidth_fullExpr",
#   "sim_largeBandwidth_medExpr"
# )

## simulations for shuffle based on medium bandwidth, medium expression
## filenames for saved simulation datasets
sim_names_shuffle <- c(
  "sim_shuffle00", 
  "sim_shuffle01", 
  "sim_shuffle02", 
  "sim_shuffle03", 
  "sim_shuffle04", 
  "sim_shuffle05", 
  "sim_shuffle06", 
  "sim_shuffle07", 
  "sim_shuffle08", 
  "sim_shuffle09", 
  "sim_shuffle10"
)

## load each dataset when we analyze

# Run method -------------------------------------------------------------

## simulations for various bandwidth and expression
## iterate over 1. dataset, 2. resolution, 3. rotation (save everything in one df for each dataset)
res_list <- c(list("singlecell"), as.list(seq(0.01, 0.1, by = 0.01)*6000))
n_rotation <- 10
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)

for (i in sim_names) {
  print(paste0("Dataset: ", i))
  
  ## load dataset,
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  
  ## iterate over resolutions
  nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
    print(paste0("Resolution: ", res))
    if (res == "singlecell") {
      ## iterate over rotations (single cell)
      df <- do.call(rbind, lapply(angle_deg_list, function(deg) {
        print(paste0("Rotation (degrees): ", deg))
        ## rotate xy coordinates
        spe_rotated <- SpatialExperiment::SpatialExperiment(
          assays = assays(spe),
          spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
        )
        
        ## nnSVG
        spe_rotated_nnsvg <- try({nnSVG::nnSVG(
          spe_rotated,
          assay_name = "logcounts",
          BPPARAM = BiocParallel::MulticoreParam()
        )})
        
        if (class(spe_rotated_nnsvg) == "try-error") {
          ## do not save anything for dataset/resolution/rotation that caused error in nnSVG
          return(NULL)
        } else {
          temp <- rownames_to_column(as.data.frame(rowData(spe_rotated_nnsvg)), var = "gene")
          temp <- cbind(rotation_deg = deg, num_pixels = dim(spe_rotated)[2], temp)
          return(temp)
        }
      }))
    } else {
      ## iterate over rotations (rasterization)
      df <- do.call(rbind, lapply(angle_deg_list, function(deg) {
        print(paste0("Rotation (degrees): ", deg))
        ## rotate xy coordinates
        spe_rotated <- SpatialExperiment::SpatialExperiment(
          assays = assays(spe),
          spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
        )
        
        ## rasterization
        spe_rast <- SEraster::rasterizeGeneExpression(spe_rotated, assay_name = "logcounts", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
        
        ## nnSVG
        spe_rast_nnsvg <- try({nnSVG::nnSVG(
          spe_rast,
          assay_name = "pixelval",
          BPPARAM = BiocParallel::MulticoreParam()
        )})
        
        if (class(spe_rast_nnsvg) == "try-error") {
          ## do not save anything for dataset/resolution/rotation that caused error in nnSVG
          return(NULL)
        } else {
          temp <- rownames_to_column(as.data.frame(rowData(spe_rast_nnsvg)), var = "gene")
          temp <- cbind(rotation_deg = deg, num_pixels = dim(spe_rast)[2], temp)
          return(temp)
        }
      }))
    }
    return(data.frame(dataset = i, resolution = res, df))
  }))
  saveRDS(nnsvg_results, file = here("outputs", dataset_name, paste0(i, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
}


## simulations for shuffle based on medium bandwidth, medium expression
## iterate over 1. dataset, 2. resolution, 3. rotation (save everything in one df for each dataset)
res_list <- c(list("singlecell"), as.list(seq(0.01, 0.1, by = 0.01)*6000))
n_rotation <- 10
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)

for (i in sim_names_shuffle) {
  print(paste0("Dataset: ", i))
  
  ## load dataset,
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  
  ## iterate over resolutions
  nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
    print(paste0("Resolution: ", res))
    if (res == "singlecell") {
      ## iterate over rotations (single cell)
      df <- do.call(rbind, lapply(angle_deg_list, function(deg) {
        print(paste0("Rotation (degrees): ", deg))
        ## rotate xy coordinates
        spe_rotated <- SpatialExperiment::SpatialExperiment(
          assays = assays(spe),
          spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
        )
        
        ## nnSVG
        spe_rotated_nnsvg <- try({nnSVG::nnSVG(
          spe_rotated,
          assay_name = "logcounts",
          BPPARAM = BiocParallel::MulticoreParam()
        )})
        
        if (class(spe_rotated_nnsvg) == "try-error") {
          ## do not save anything for dataset/resolution/rotation that caused error in nnSVG
          return(NULL)
        } else {
          temp <- rownames_to_column(as.data.frame(rowData(spe_rotated_nnsvg)), var = "gene")
          temp <- cbind(rotation_deg = deg, num_pixels = dim(spe_rotated)[2], temp)
          return(temp)
        }
      }))
    } else {
      ## iterate over rotations (rasterization)
      df <- do.call(rbind, lapply(angle_deg_list, function(deg) {
        print(paste0("Rotation (degrees): ", deg))
        ## rotate xy coordinates
        spe_rotated <- SpatialExperiment::SpatialExperiment(
          assays = assays(spe),
          spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
        )
        
        ## rasterization
        spe_rast <- SEraster::rasterizeGeneExpression(spe_rotated, assay_name = "logcounts", resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
        
        ## nnSVG
        spe_rast_nnsvg <- try({nnSVG::nnSVG(
          spe_rast,
          assay_name = "pixelval",
          BPPARAM = BiocParallel::MulticoreParam()
        )})
        
        if (class(spe_rast_nnsvg) == "try-error") {
          ## do not save anything for dataset/resolution/rotation that caused error in nnSVG
          return(NULL)
        } else {
          temp <- rownames_to_column(as.data.frame(rowData(spe_rast_nnsvg)), var = "gene")
          temp <- cbind(rotation_deg = deg, num_pixels = dim(spe_rast)[2], temp)
          return(temp)
        }
      }))
    }
    return(data.frame(dataset = i, resolution = res, df))
  }))
  saveRDS(nnsvg_results, file = here("outputs", dataset_name, paste0(i, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
}

# Plot --------------------------------------------------------------------

## Figure x (single cell)
## simulations for various bandwidth and expression
sim_names <- c(
  "sim_largeBandwidth_fullExpr",
  "sim_largeBandwidth_medExpr",
  "sim_largeBandwidth_lowExpr",
  "sim_medBandwidth_fullExpr",
  "sim_medBandwidth_medExpr",
  "sim_medBandwidth_lowExpr",
  "sim_smallBandwidth_fullExpr",
  "sim_smallBandwidth_medExpr",
  "sim_smallBandwidth_lowExpr"
)

df <- do.call(rbind, lapply(sim_names, function(i) {
  ## load dataset
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  svg_example <- which(rowData(spe)$expressed)[1]
  
  meta <- unlist(strsplit(gsub("^sim_(.*?)Bandwidth_(.*?)Expr$", "\\1 \\2", i), " "))
  
  return(data.frame(dataset = i, bandwidth = meta[1], expression = meta[2], x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], gene = assay(spe)[svg_example,]))
}))

df <- df %>%
  mutate(bandwidth = factor(bandwidth, levels = c("large", "med", "small")),
         expression = factor(expression, levels = c("full", "med", "low")))

ggplot(df, aes(x = x, y = y, col = gene)) +
  facet_grid(bandwidth ~ expression) +
  coord_fixed() +
  geom_point(size = 0.7) +
  scale_color_viridis_c(name = "Log-normalized\nexpression") +
  labs(title = "Simulated dataset",
       x = "x (um)",
       y = "y (um)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )
ggsave(filename = here("plots", dataset_name, paste0("nnsvg_sim_singlecell.pdf")), width = 12, height = 10, dpi = 300)

## simulations for shuffle based on medium bandwidth, medium expression
sim_names_shuffle <- c(
  "sim_shuffle00", 
  "sim_shuffle01", 
  "sim_shuffle02", 
  "sim_shuffle03", 
  "sim_shuffle04", 
  "sim_shuffle05", 
  "sim_shuffle06", 
  "sim_shuffle07", 
  "sim_shuffle08", 
  "sim_shuffle09", 
  "sim_shuffle10"
)

df <- do.call(rbind, lapply(sim_names_shuffle, function(i) {
  ## load dataset
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  svg_example <- which(rowData(spe)$expressed)[1]
  
  meta <- as.numeric(gsub("[^0-9]", "", i))
  
  return(data.frame(dataset = i, shuffle = paste0(meta*10, "% shuffled"), x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], gene = assay(spe)[svg_example,]))
}))

df <- df %>%
  mutate(shuffle = factor(shuffle, levels = paste0(seq(0, 100, by = 10), "% shuffled")))

ggplot(df, aes(x = x, y = y, col = gene)) +
  facet_wrap(~shuffle) +
  coord_fixed() +
  geom_point(size = 0.7) +
  scale_color_viridis_c(name = "Log-normalized\nexpression") +
  labs(title = "Simulated dataset (shuffle)",
       x = "x (um)",
       y = "y (um)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )
ggsave(filename = here("plots", dataset_name, paste0("nnsvg_sim_shuffle_singlecell.pdf")), width = 15, height = 10, dpi = 300)

## Figure x (single cell and rasterization visualizations)
res_list <- c(list("singlecell"), as.list(seq(0.01, 0.1, by = 0.01)*6000))

for (i in sim_names) {
  ## load dataset
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  svg_example <- which(rowData(spe)$expressed)[1]
  
  for (res in res_list) {
    if (res == "singlecell") {
      ## plot
      df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], gene = assay(spe)[svg_example,])
      plt <- ggplot(df, aes(x = x, y = y, col = gene)) +
        coord_fixed() +
        geom_point(size = 1) +
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
      spe_rast <- SEraster::rasterizeGeneExpression(spe, resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
      
      ## plot
      df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], gene = assay(spe_rast)[svg_example,])
      plt <- ggplot(df, aes(x = x, y = y, fill = gene)) +
        coord_fixed() +
        geom_tile(width = res, height = res) +
        scale_fill_viridis_c() +
        labs(title = paste0("Resolution = ", res, " %")) +
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
    ggsave(plot = plt, filename = here("plots", dataset_name, paste0(i, "_gene", svg_example, "_",  res, ".pdf")))
  }
}

## Figure (performance comparison)
res_list <- c(list("singlecell"), as.list(seq(0.01, 0.1, by = 0.01)*6000))

## simulations for various bandwidth and expression
sim_names <- c(
  "sim_largeBandwidth_fullExpr",
  "sim_largeBandwidth_medExpr",
  "sim_largeBandwidth_lowExpr",
  "sim_medBandwidth_fullExpr",
  "sim_medBandwidth_medExpr",
  "sim_medBandwidth_lowExpr",
  "sim_smallBandwidth_fullExpr",
  "sim_smallBandwidth_medExpr",
  "sim_smallBandwidth_lowExpr"
)

# sim_names <- c(
#   "sim_largeBandwidth_fullExpr",
#   "sim_largeBandwidth_medExpr",
#   "sim_largeBandwidth_lowExpr"
# )
i <- sim_names[1]
n_rotation <- 10
spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
df <- readRDS(file = here("outputs", dataset_name, paste0(i, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
temp <- df[df$resolution == "0.01" & df$rotation_deg == 36,]
results_sig <- data.frame(gene = temp$gene, pred = temp$padj <= alpha, obs = rowData(spe)$expressed)
calculatePerformanceMetrics(results_sig)

# set a threshold p value
alpha <- 0.05
n_rotation <- 10
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)
df_perf <- do.call(rbind, lapply(sim_names, function(i) {
  ## load original spe to get T/F labels
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  ## load nnSVG results
  df <- readRDS(file = here("outputs", dataset_name, paste0(i, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
  
  out <- do.call(rbind, lapply(unique(df$resolution), function(res) {
    temp <- df[df$resolution == res,]
    out2 <- do.call(rbind, lapply(unique(temp$rotation_deg), function(deg) {
      temp2 <- temp[temp$rotation_deg == deg,]
      results_sig <- data.frame(gene = temp2$gene, pred = temp2$padj <= alpha, obs = rowData(spe)$expressed)
      perf <- calculatePerformanceMetrics(results_sig)
      return(data.frame(rotation_deg = deg, perf))
    }))
    return(data.frame(resolution = res, out2))
  }))
  
  meta <- unlist(strsplit(gsub("^sim_(.*?)Bandwidth_(.*?)Expr$", "\\1 \\2", i), " "))
  
  return(data.frame(dataset = i, bandwidth = meta[1], expression = meta[2], out))
}))

df_perf_raw <- df_perf %>%
  mutate(resolution = factor(resolution, levels = res_list),
         bandwidth = factor(bandwidth, levels = c("large", "med", "small")),
         expression = factor(expression, levels = c("full", "med", "low"))) %>%
  select(dataset, bandwidth, expression, resolution, rotation_deg, TPR, specificity, PPV, F1, ACC) %>%
  pivot_longer(!c(dataset, bandwidth, expression, resolution, rotation_deg), names_to = "metrics", values_to = "values")

df_perf_summary <- do.call(rbind, lapply(unique(df_perf$resolution), function(res) {
  out <- do.call(rbind, lapply(unique(df_perf$dataset), function(dataset) {
    out2 <- do.call(rbind, lapply(c("TPR", "specificity", "PPV", "F1", "ACC"), function(metric) {
      temp <- df_perf[df_perf$dataset == dataset & df_perf$resolution == res, metric]
      return(data.frame(metrics = metric, mean = mean(temp), sd = sd(temp)))
    }))
    
    meta <- unlist(strsplit(gsub("^sim_(.*?)Bandwidth_(.*?)Expr$", "\\1 \\2", dataset), " "))
    
    return(data.frame(dataset = dataset, bandwidth = meta[1], expression = meta[2], resolution = factor(res, levels = res_list), out2))
  }))
}))
df_perf_summary <- df_perf_summary %>%
  mutate(resolution = factor(resolution, levels = res_list),
         bandwidth = factor(bandwidth, levels = c("large", "med", "small")),
         expression = factor(expression, levels = c("full", "med", "low")))

df_perf_summary2 <- df_perf_summary
df_perf_summary2$resolution <- sub("singlecell", 0, df_perf_summary2$resolution)
df_perf_summary2$resolution <- as.numeric(df_perf_summary2$resolution)

ggplot(df_perf_raw, aes(x = resolution, y = values, col = metrics)) +
  facet_grid(bandwidth ~ expression) +
  # geom_jitter(width = 10, alpha = 0.3) +
  # geom_line(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics)) +
  geom_point(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics), size = 1) +
  geom_errorbar(data = df_perf_summary, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics)) +
  # scale_x_discrete(breaks = unique(df_perf_raw$resolution)) + 
  ylim(0,1) +
  labs(title = "Performance metrics comparison",
       x = "Resolution",
       y = "Performance",
       col = "Metric") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_perf_metric_summary.pdf")), width = 8, heigh = 8, dpi = 300)

ggplot(df_perf_raw, aes(x = resolution, y = values, col = metrics)) +
  facet_grid(bandwidth ~ expression) +
  # geom_jitter(width = 10, alpha = 0.3) +
  geom_line(data = df_perf_summary2, aes(x = resolution, y = mean, col = metrics)) +
  geom_point(data = df_perf_summary2, aes(x = resolution, y = mean, col = metrics), size = 1) +
  geom_errorbar(data = df_perf_summary2, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics)) +
  scale_x_continuous(breaks = unique(df_perf_summary2$resolution)) +
  ylim(0,1) +
  labs(title = "Performance metrics comparison",
       x = "Resolution",
       y = "Performance",
       col = "Metric") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_perf_metric_summary_v2.pdf")), width = 12, heigh = 10, dpi = 300)


## simulations for shuffle based on medium bandwidth, medium expression
sim_names_shuffle <- c(
  "sim_shuffle00", 
  "sim_shuffle01", 
  "sim_shuffle02", 
  "sim_shuffle03", 
  "sim_shuffle04", 
  "sim_shuffle05", 
  "sim_shuffle06", 
  "sim_shuffle07", 
  "sim_shuffle08", 
  "sim_shuffle09", 
  "sim_shuffle10"
)

# set a threshold p value
alpha <- 0.05
n_rotation <- 10
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)
df_perf <- do.call(rbind, lapply(sim_names_shuffle, function(i) {
  ## load original spe to get T/F labels
  spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))
  ## load nnSVG results
  df <- readRDS(file = here("outputs", dataset_name, paste0(i, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
  
  out <- do.call(rbind, lapply(unique(df$resolution), function(res) {
    temp <- df[df$resolution == res,]
    out2 <- do.call(rbind, lapply(unique(temp$rotation_deg), function(deg) {
      temp2 <- temp[temp$rotation_deg == deg,]
      results_sig <- data.frame(gene = temp2$gene, pred = temp2$padj <= alpha, obs = rowData(spe)$expressed)
      perf <- calculatePerformanceMetrics(results_sig)
      return(data.frame(rotation_deg = deg, perf))
    }))
    return(data.frame(resolution = res, out2))
  }))
  
  meta <- as.numeric(gsub("[^0-9]", "", i))
  
  return(data.frame(dataset = i, shuffle_label = paste0(meta*10, "% shuffled"), shuffle_num = meta*10, out))
}))

df_perf_raw <- df_perf %>%
  mutate(resolution = factor(resolution, levels = res_list),
         shuffle_label = factor(shuffle_label, levels = paste0(seq(0, 100, by = 10), "% shuffled")),
         shuffle_num = factor(shuffle_num, levels = seq(0, 100, by = 10))) %>%
  select(dataset, shuffle_label, shuffle_num, resolution, rotation_deg, TPR, specificity, PPV, F1, ACC) %>%
  pivot_longer(!c(dataset, shuffle_label, shuffle_num, resolution, rotation_deg), names_to = "metrics", values_to = "values")

df_perf_summary <- do.call(rbind, lapply(unique(df_perf$resolution), function(res) {
  out <- do.call(rbind, lapply(unique(df_perf$dataset), function(dataset) {
    out2 <- do.call(rbind, lapply(c("TPR", "specificity", "PPV", "F1", "ACC"), function(metric) {
      temp <- df_perf[df_perf$dataset == dataset & df_perf$resolution == res, metric]
      return(data.frame(metrics = metric, mean = mean(temp), sd = sd(temp)))
    }))
    
    meta <- as.numeric(gsub("[^0-9]", "", dataset))
    
    return(data.frame(dataset = dataset, shuffle_label = paste0(meta*10, "% shuffled"), shuffle_num = meta*10, resolution = factor(res, levels = res_list), out2))
  }))
}))
df_perf_summary <- df_perf_summary %>%
  mutate(resolution = factor(resolution, levels = res_list),
         shuffle_label = factor(shuffle_label, levels = paste0(seq(0, 100, by = 10), "% shuffled")))

ggplot(df_perf_raw, aes(x = shuffle_num, y = values, col = metrics)) +
  facet_wrap(~resolution) +
  # geom_jitter(width = 10, alpha = 0.3) +
  geom_line(data = df_perf_summary, aes(x = shuffle_num, y = mean, col = metrics)) +
  geom_point(data = df_perf_summary, aes(x = shuffle_num, y = mean, col = metrics), size = 1) +
  geom_errorbar(data = df_perf_summary, aes(x = shuffle_num, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics)) +
  # scale_x_continuous(breaks = unique(df_perf_summary$shuffle_num)) +
  ylim(0,1) +
  labs(title = "Performance metrics comparison",
       x = "% shuffled",
       y = "Performance",
       col = "Metric") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_shuffle_perf_metric_summary.pdf")), width = 12, heigh = 8, dpi = 300)

# Further exploration -----------------------------------------------------

## load dataset
i <- sim_names[1]
i
spe <- readRDS(file = here(dir, paste0("spe_", i, ".RDS")))

# ## single cell
# spe_nnsvg <- nnSVG::nnSVG(
#   spe,
#   BPPARAM = BiocParallel::MulticoreParam()
# )

# ## rasterize (debug)
# data <- assay(spe)
# pos <- spatialCoords(spe)
# bbox <- bbox <- sf::st_bbox(c(
#   xmin = floor(min(pos[,1])-resolution/2), 
#   xmax = ceiling(max(pos[,1])+resolution/2), 
#   ymin = floor(min(pos[,2])-resolution/2), 
#   ymax = ceiling(max(pos[,2])+resolution/2)
# ))
# fun <- "mean"
# 
# ## create grid for rasterization
# grid <- sf::st_make_grid(bbox, cellsize = resolution)
# 
# ## extract position
# pos_pixel <- sf::st_coordinates(sf::st_centroid(grid))
# colnames(pos_pixel) <- c("x", "y")
# rownames(pos_pixel) <- paste0("pixel",seq_along(pos_pixel[,1]))
# 
# ## convert pos to sf
# cells <- sf::st_as_sf(data.frame(pos), coords = c("x", "y"))
# ## get pixel ID for each cell
# ## since some cells are assigned to > 1 pixels if they are on the border, choose the 1st pixel id
# pixel_ids <- unlist(lapply(sf::st_intersects(cells, grid), function(sublist) sublist[1]))
# names(pixel_ids) <- rownames(pos)
# 
# ## store aggregated subsetted matrix data for each pixel, store cell IDs and number of cells for each pixel into a data frame
# out <- unlist(BiocParallel::bplapply(sort(unique(pixel_ids)), function(id){
#   ## get cell IDs for a particular pixel
#   cell_ids <- names(pixel_ids[pixel_ids == id])
#   
#   ## subset feature observation matrix
#   sub <- data[,cell_ids, drop = FALSE]
#   ## aggregate cell counts to create pixel value
#   if (fun == "mean") {
#     pixel_val <- rowMeans(sub)
#   } else if (fun == "sum") {
#     pixel_val <- rowSums(sub)
#   }
#   
#   ## store number of cells
#   meta_rast <- data.frame(num_cell = length(cell_ids))
#   ## store a list of cell IDs
#   meta_rast$cellID_list <- list(cell_ids)
#   
#   if (is.matrix(data)) {
#     return(list(pixel_val, meta_rast))
#   } else {
#     return(list(as(pixel_val, "CsparseMatrix"), meta_rast))
#   }
# }), recursive = FALSE)
# 
# ## extract rasterized sparse matrix
# data_rast <- do.call(cbind, out[seq(1,length(out),by=2)])
# 
# ## extract rasterized data frame
# meta_rast <- do.call(rbind, out[seq(1,length(out),by=2)+1])
# 
# ## set rownames/colnames for rasterized sparse matrix and rasterized data frame
# rownames(data_rast) <- rownames(data)
# colnames(data_rast) <- paste0("pixel", sort(unique(pixel_ids)))
# rownames(meta_rast) <- paste0("pixel", sort(unique(pixel_ids)))
# 
# ## subset rasterized pos
# pos_pixel <- pos_pixel[rownames(pos_pixel) %in% colnames(data_rast),]

## rasterize
resolution <- 0.1
spe_rast <- SEraster::rasterizeGeneExpression(spe, resolution = resolution, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())

## plot rasterized data
svg_example <- which(rowData(spe)$expressed)[1]
df <- data.frame(x = spatialCoords(spe_rast)[,1], y = spatialCoords(spe_rast)[,2], gene = assay(spe_rast)[svg_example,])

ggplot(df, aes(x = x, y = y, fill = gene)) +
  # geom_tile(width = resolution, height = resolution) +
  # geom_tile() +
  geom_point() +
  coord_fixed() +
  scale_fill_viridis_c(name = "svg_example") +
  labs(title = paste0("Resolution = ", resolution)) +
  theme_classic()

## see why geom_tile() doesn't work
a <- diff(spatialCoords(spe_rast))

spe_rast <- nnSVG::nnSVG(
  spe_rast,
  assay_name = "pixelval",
  BPPARAM = BiocParallel::MulticoreParam()
)
