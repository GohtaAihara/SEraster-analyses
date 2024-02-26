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
library(rearrr)
library(UpSetR)

par(mfrow=c(1,1))

dataset_name <- "merfish_mousePOA"

# Load dataset ------------------------------------------------------------

# data <- read.csv('~/Downloads/lab/data/merfish_mousePOA_all_cells.csv')

data <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_mousePOA/merfish_mousePOA_all_cells.csv')

animal <- 1
sex <- "Female"
behavior <- "Naive"
bregma <- "-0.19"

animals <- unique(data$Animal_ID)
sexes <- unique(data$Animal_sex)
bregmas <- unique(data$Bregma)

count = 0

## set resolution parameters
res_list <- list("singlecell", 50, 100, 200) # feel free to add more resolutions

## set rotation/permutation parameters
n_rotation <- 10 # use 10 for the actual figure!
angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)

## set BiocParallel parameters
# the number of workers is set to be the number of physical cores
# bpparam <- BiocParallel::MulticoreParam()
# bpparam <- BiocParallel::MulticoreParam(workers = parallel::detectCores(logical = FALSE))
bpparam <- BiocParallel::MulticoreParam(workers = 10)
bpparam

#run on mac studio
behavior <- "Naive"
for (animal in animals) {
  for (sex in sexes) {
    for (bregma in bregmas) {
      if (file.exists(here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_preprocessed.RDS")))) {
        count = count + 1
        print(paste0("Animal: ", animal, ", Sex: ", sex, ", Bregma: ", bregma))
        spe <- readRDS(file = here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_preprocessed.RDS")))
        plot(spatialCoords(spe), pch=".", asp=1)
        
        ## Rotate dataset, rasterize, run nnSVG for each resolution
        nnsvg_results <- do.call(rbind, lapply(res_list, function(res) {
          print(paste0("Resolution: ", res))
          if (res == "singlecell") {
            num_points = dim(spe)[2]
            
            ## nnSVG
            spe <- try({
              nnSVG::nnSVG(
                spe,
                assay_name = "lognorm",
                BPPARAM = bpparam
              )
            })
            if (class(spe) == "try-error") {
              ## do not save anything for clusters that caused error in nnSVG
              return(NULL)
              
            } else {
              df <- rownames_to_column(as.data.frame(rowData(spe)), var = "gene")
              df <- cbind(rotation_deg = NA, num_points = num_points, df)
            }
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
              spe_rast <- try({
                nnSVG::nnSVG(
                  spe_rast,
                  assay_name = "pixelval",
                  BPPARAM = bpparam
                )
              })
              if (class(spe_rast) == "try-error") {
                ## do not save anything for clusters that caused error in nnSVG
                return(NULL)
                
              } else {
                temp <- rownames_to_column(as.data.frame(rowData(spe_rast)), var = "gene")
                temp <- cbind(rotation_deg = deg, num_points = num_points, temp)
                return(temp)
              }
            }))
          }
          if (is.null(df)) {
            print(paste0("All permutations (rotations) for resolution ", res, " failed"))
            return(NULL)
          } else {
            return(data.frame(dataset = dataset_name, resolution = res, df))
          }
        }))
        saveRDS(nnsvg_results, file = here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
      }
    }
  }
}


# Bug fix -----------------------------------------------------------------

behavior <- "Naive"

animal <- 2
sex <- "Female"
bregma <- -0.29

bpparam <- BiocParallel::MulticoreParam(workers = 10)

spe <- readRDS(file = here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_preprocessed.RDS")))

plot(spatialCoords(spe), pch=".", asp=1)

## single-cell
spe <- nnSVG::nnSVG(
  spe,
  assay_name = "lognorm",
  BPPARAM = bpparam
)

## rasterization
res <- 200
deg <- 36
## rotate xy coordinates
spe_rotated <- SpatialExperiment::SpatialExperiment(
  assays = assays(spe),
  spatialCoords = rotateAroundCenter(spatialCoords(spe), deg)
)

## rasterization
spe_rast <- SEraster::rasterizeGeneExpression(spe_rotated, assay_name = "lognorm", resolution = res, fun = "mean", BPPARAM = bpparam)
num_points = dim(spe_rast)[2]

## nnSVG (might want to use try(), see "run_nnSVG_nnSVG_simulations.R")
spe_rast <- try({
  nnSVG::nnSVG(
    spe_rast,
    assay_name = "pixelval",
    BPPARAM = bpparam
  )
})

# Plot --------------------------------------------------------------------

for (animal in animals) {
  for (sex in sexes) {
    for (bregma in bregmas) {
      if (file.exists(here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_preprocessed.RDS")))) {
        ## Figure x (performance comparison)
        n_rotation <- 1
        angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)
        df <- readRDS(file = here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))
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
        
        ggplot(df_perf2, aes(x = resolution, y = values, col = metrics)) +
          geom_jitter(width = 10, alpha = 0.3) +
          geom_line(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics)) +
          geom_point(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics), size = 1) +
          geom_errorbar(data = df_perf_summary, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics), width = 10) +
          scale_x_continuous(breaks = unique(df_perf2$resolution)) + 
          ylim(0,1) +
          labs(title = paste("mPOA Performance for", animal, sex, behavior, bregma),
               x = "Rasterization Resolution",
               y = "Performance",
               col = "Metric") +
          theme_bw()
        ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_perf_metric_summary.pdf")), width = 6, heigh = 5, dpi = 300)
      }
    }
  }
}
        
## visualize intersection of SVGs across different resolutions
u_df <- data.frame(gene = df$gene,
                       resolution = df$resolution, 
                       identified = ifelse(df$pval < 0.05, 1, 0))

sc_df <- u_df[u_df$resolution == "singlecell", ]
res_50 <- u_df[u_df$resolution == 50, ]
res_100 <- u_df[u_df$resolution == 100, ]
res_200 <- u_df[u_df$resolution == 200, ]

upset_df <- data.frame(gene = sc_df$gene,
                   sc = sc_df$identified,
                   r50 = res_50$identified,
                   r100 = res_100$identified,
                   r200 = res_200$identified)

upset(upset_df)

## identify which genes + visualize their patterns

## plot single cell
df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], gene = assay(spe, "counts")[gene,])
ggplot(df, aes(x = x, y = y, col = gene)) +
  coord_fixed() +
  geom_point(size = 1, stroke = 0) +
  scale_color_viridis_c() +
  theme_bw() +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
identified_sc <- u_df[u_df$resolution == "singlecell" & u_df$pval < 0.05, ]

# identified_200um <- upset_df[upset_df$resolution == 200 & upset_df$pval < 0.05, ]
# not_identified_200um <- upset_df[upset_df$resolution == 200 & upset_df$pval > 0.05, ]
# res_200um <- upset_df[upset_df$resolution == 200, ]
# res_200um$pval <- ifelse(res_200um$pval < 0.05, 1, 0)
# res_200um <- data.frame(as.list(res_200um))



# Questions ---------------------------------------------------------------

# Is the performance consistent across bregma sections for each mouse OR across mice for each bregma section?
# only tested for 1_Female_Naive_0.29
# PPV slightly higher, TPR looks the correct

# How does the correlation of nnSVG outputs look like single-cell vs. rasterized resolutions (e.g. gene ranking)?


# UpSetR
# Jake R Conway, Alexander Lex, Nils Gehlenborg UpSetR: An R Package for the Visualization of Intersecting Sets and their Properties doi: 
# https://doi.org/10.1093/bioinformatics/btx364
        
        
        
        
        
        
        