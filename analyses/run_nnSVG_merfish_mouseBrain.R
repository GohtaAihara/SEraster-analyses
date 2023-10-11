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

dataset_name <- "merfish_mouseBrain"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_preprocessed.RDS")))

plot(spatialCoords(spe), pch=".", asp=1)


# Run method --------------------------------------------------------------

res_list <- list("singlecell", 50, 100, 200, 400)
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
      return(temp)
    }))
  }
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
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global_runtime.RDS")))
df$resolution <- factor(df$resolution, levels = c("singlecell", "50", "100", "200", "400"))

# total runtime
ggplot(df, aes(x = num_points, y = as.numeric(runtime_total), col = resolution)) +
  geom_boxplot(width = 6000, lwd = 0.75, outlier.shape = NA) +
  geom_jitter(width = 1000, size = 1, alpha = 0.75) +
  scale_x_continuous(breaks = unique(df$num_points)) + 
  labs(title = "Total runtime",
       x = "Number of spatial points",
       y = "Runtime (secs)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_total.pdf")), width = 6, heigh = 5, dpi = 300)

ggplot(df, aes(x = resolution, y = as.numeric(runtime_total), col = resolution)) +
  geom_boxplot(lwd = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  labs(title = "Total runtime",
       x = "Resolution",
       y = "Runtime (secs)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        legend.position = "none")
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_total_v2.pdf")), width = 5, heigh = 5, dpi = 300)

# nnSVG runtime
ggplot(df, aes(x = num_points, y = as.numeric(runtime_nnsvg), col = resolution)) +
  geom_boxplot(width = 6000, lwd = 0.75, outlier.shape = NA) +
  geom_jitter(width = 1000, size = 1, alpha = 0.75) +
  scale_x_continuous(breaks = unique(df$num_points)) + 
  labs(title = "nnSVG runtime",
       x = "Number of spatial points",
       y = "Runtime (secs)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_nnsvg.pdf")), width = 6, heigh = 5, dpi = 300)

ggplot(df, aes(x = resolution, y = as.numeric(runtime_nnsvg), col = resolution)) +
  geom_boxplot(lwd = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  labs(title = "nnSVG runtime",
       x = "Resolution",
       y = "Runtime (secs)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        legend.position = "none")
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_nnsvg_v2.pdf")), width = 5, heigh = 5, dpi = 300)

# SEraster runtime
ggplot(df[df$resolution != "singlecell",], aes(x = num_points, y = as.numeric(runtime_rast), col = resolution)) +
  geom_boxplot(width = 2000, lwd = 0.75, outlier.shape = NA) +
  geom_jitter(width = 500, size = 1, alpha = 0.75) +
  scale_x_continuous(breaks = unique(df$num_points)) + 
  labs(title = "SEraster runtime",
       x = "Number of spatial points",
       y = "Runtime (secs)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_rast.pdf")), width = 6, heigh = 5, dpi = 300)

ggplot(df[df$resolution != "singlecell",], aes(x = resolution, y = as.numeric(runtime_rast), col = resolution)) +
  geom_boxplot(lwd = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  labs(title = "SEraster runtime",
       x = "Resolution",
       y = "Runtime (secs)",
       col = "Resolution") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        legend.position = "none")
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_runtime_rast_v2.pdf")), width = 5, heigh = 5, dpi = 300)

## Figure 1c (performance comparison)
df <- readRDS(file = here("outputs", paste0(dataset_name, "_nnsvg_global.RDS")))
# set a threshold p value
alpha <- 0.05
df_perf <- do.call(rbind, lapply(unique(df$resolution), function(res) {
  if (res != "singlecell") {
    sc <- df[df$resolution == "singlecell",]
    rast <- df[df$resolution == res,]
    results_sig <- do.call(rbind, lapply(rast$gene, function(gene) {
      return(data.frame(gene = gene, pixel = rast[rast$gene == gene, "padj"] <= alpha, singlecell = sc[sc$gene == gene, "padj"] <= alpha))
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
  select(resolution, TPR, specificity, PPV, F1, ACC) %>%
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
        return(data.frame(gene = gene, pixel = rast[rast$gene == gene, "padj"] <= alpha, singlecell = sc[sc$gene == gene, "padj"] <= alpha))
      }))
      out <- calculatePerformanceMetrics(results_sig)
      return(data.frame(rotation_deg = deg, out))
    }))
    return(data.frame(resolution = res, out))
  }
}))

df_perf_raw <- df_perf %>%
  mutate(resolution = as.numeric(resolution)) %>%
  select(resolution, rotation_deg, TPR, specificity, PPV, F1, ACC) %>%
  pivot_longer(!c(resolution, rotation_deg), names_to = "metrics", values_to = "values")

df_perf_summary <- do.call(rbind, lapply(unique(df_perf$resolution), function(res) {
  out <- do.call(rbind, lapply(c("TPR", "specificity", "PPV", "F1", "ACC"), function(metric) {
    temp <- df_perf[df_perf$resolution == res, metric]
    return(data.frame(metrics = metric, mean = mean(temp), sd = sd(temp)))
  }))
  return(data.frame(resolution = as.numeric(res), out))
}))

ggplot(df_perf2, aes(x = resolution, y = values, col = metrics)) +
  # geom_jitter(width = 10, alpha = 0.3) +
  geom_line(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics)) +
  geom_point(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics), size = 1) +
  geom_errorbar(data = df_perf_summary, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics), width = 10) +
  scale_x_continuous(breaks = unique(df_perf2$resolution)) + 
  ylim(0,1) +
  labs(title = "Performance metrics comparison",
       x = "Resolution",
       y = "Performance",
       col = "Metric") +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_perf_metric_summary.pdf")), width = 6, heigh = 5, dpi = 300)

## Figure 1d (nnSVG results comparison)
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
