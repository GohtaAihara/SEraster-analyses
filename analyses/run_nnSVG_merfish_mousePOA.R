## This file integrates rasterization with nnSVG to analyze MERFISH mPOA datasets.

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

data <- read.csv('~/Downloads/lab/data/merfish_mousePOA_all_cells.csv')

#data <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_mousePOA/merfish_mousePOA_all_cells.csv')

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
#df <- readRDS(file = here("outputs", paste0(dataset_name, "_animal1", "_sex", sex, "_behavior", behavior, "_bregma_-0.19", "_nnsvg_global_", "n_rotation_", n_rotation, ".RDS")))

animal <- 1
sex <- "Female"
behavior <- "Naive"
bregma <- "0.26"

animals <- unique(data$Animal_ID)
sexes <- unique(data$Animal_sex)
bregmas <- unique(data$Bregma)

df_perf_all <- data.frame()
df_perf_all2 <- data.frame()
 
# there are 83 unique conditions       
for (animal in animals) {
  for (sex in sexes) {
    for (bregma in bregmas) {
      if (file.exists(here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_nnsvg_global_n_rotation_10", ".RDS")))) {
          ## Figure x (performance comparison)
          n_rotation <- 10
          angle_deg_list <- seq(0, 360-0.1, by = 360/n_rotation)
          df <- readRDS(file = here("outputs", paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_nnsvg_global_n_rotation_10", ".RDS")))
          # set a threshold p value
          alpha <- 0.05
          if ("singlecell" %in% unique(df$resolution)) {
            df_perf <- do.call(rbind, lapply(unique(df$resolution), function(res) {
              df_sub <- df[df$resolution == res,]
              if (res != "singlecell") {
                sc <- df[df$resolution == "singlecell",]
                out <- do.call(rbind, lapply(unique(df_sub$rotation_deg), function(deg) {
                  rast <- df_sub[df_sub$rotation_deg == deg,]
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
            
            #df_perf_all <- do.call(rbind, df_perf_summary)
            df_perf_all <- rbind(df_perf_all, df_perf_summary)
            
            df_perf2 <- df_perf %>%
              mutate(resolution = as.numeric(resolution)) %>%
              select(resolution, TPR, TNR, PPV) %>%
              pivot_longer(!resolution, names_to = "metrics", values_to = "values")
            
            #df_perf_all2 <- do.call(rbind, df_perf2)
            df_perf_all2 <- rbind(df_perf_all2, df_perf2)
            #df_perf_all2 <-df_perf %>%
              # mutate(resolution = as.numeric(resolution)) %>%
              # select(resolution, TPR, TNR, PPV) %>%
              # pivot_longer(!resolution, names_to = "metrics", values_to = "values")
            
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
            
           
            #ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_perf_metric_summary.pdf")), width = 6, heigh = 5, dpi = 300)
            
            
            # ggplot(df_perf2, aes(x = resolution, y = values, col = metrics)) +
            #   geom_jitter(width = 10, alpha = 0.3) +
            #   geom_line(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics)) +
            #   geom_point(data = df_perf_summary, aes(x = resolution, y = mean, col = metrics), size = 1) +
            #   geom_errorbar(data = df_perf_summary, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics), width = 10) +
            #   scale_x_continuous(breaks = unique(df_perf2$resolution)) + 
            #   ylim(0,1) +
            #   labs(title = paste("mPOA Performance for", animal, sex, behavior, bregma),
            #        x = "Rasterization Resolution",
            #        y = "Performance",
            #        col = "Metric") +
            #   theme_bw()
            # ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma, "_perf_metric_summary.pdf")), width = 6, heigh = 5, dpi = 300)
           
          } else {
              paste0("Datasets that didn't have sc resolution: ", dataset_name, "_animal", animal, "_sex", sex, "_behavior", behavior, "_bregma", bregma)
          }
          
      }
    }
  }
}
  
df_perf_all <- na.omit(df_perf_all)
#df_perf_all2 <- na.omit(df_perf_all2)

df_perf_all2
df_50TPR2 <- df_perf_all2[df_perf_all2$resolution == 50 & df_perf_all2$metrics == "TPR", ]

# 50TPR ######
df_50TPR <- df_perf_all[df_perf_all$resolution == 50 & df_perf_all$metrics == "TPR", ]
mean_50TPR <- mean(df_50TPR$mean)
sd_50TPR <- sd(df_50TPR$mean)

# 50TNR
df_50TNR <- df_perf_all[df_perf_all$resolution == 50 & df_perf_all$metrics == "TNR", ]
mean_50TNR <- mean(df_50TNR$mean)
sd_50TNR <- sd(df_50TNR$mean)

# 50PPV 
df_50PPV <- df_perf_all[df_perf_all$resolution == 50 & df_perf_all$metrics == "PPV", ]
mean_50PPV <- mean(df_50PPV$mean)
sd_50PPV <- sd(df_50PPV$mean)

# 100TPR #####
df_100TPR <- df_perf_all[df_perf_all$resolution == 100 & df_perf_all$metrics == "TPR", ]
mean_100TPR <- mean(df_100TPR$mean)
sd_100TPR <- sd(df_100TPR$mean)

# 100TPR
df_100TNR <- df_perf_all[df_perf_all$resolution == 100 & df_perf_all$metrics == "TNR", ]
mean_100TNR <- mean(df_100TNR$mean)
sd_100TNR <- sd(df_100TNR$mean)

# 100PPV 
df_100PPV <- df_perf_all[df_perf_all$resolution == 100 & df_perf_all$metrics == "PPV", ]
mean_100PPV <- mean(df_100PPV$mean)
sd_100PPV <- sd(df_100PPV$mean)

# 200TPR #####
df_200TPR <- df_perf_all[df_perf_all$resolution == 200 & df_perf_all$metrics == "TPR", ]
mean_200TPR <- mean(df_200TPR$mean)
sd_200TPR <- sd(df_200TPR$mean)

# 200TPR
df_200TNR <- df_perf_all[df_perf_all$resolution == 200 & df_perf_all$metrics == "TNR", ]
mean_200TNR <- mean(df_200TNR$mean)
sd_200TNR <- sd(df_200TNR$mean)

# 200PPV 
df_200PPV <- df_perf_all[df_perf_all$resolution == 200 & df_perf_all$metrics == "PPV", ]
mean_200PPV <- mean(df_200PPV$mean)
sd_200PPV <- sd(df_200PPV$mean)

# combine performance for biological replicates
df_perf_all_summary <- data.frame(
                resolution = c(50, 50, 50, 100, 100, 100, 200, 200, 200),
                metrics = c("TPR", "TNR", "PPV", "TPR", "TNR", "PPV", "TPR", "TNR", "PPV"),
                mean = c(mean_50TPR, mean_50TNR, mean_50PPV, 
                         mean_100TPR, mean_100TNR, mean_100PPV,
                         mean_200TPR, mean_200TNR, mean_200PPV),
                sd = c(sd_50TPR, sd_50TNR, sd_50PPV, 
                       sd_100TPR, sd_100TNR, sd_100PPV,
                       sd_200TPR, sd_200TNR, sd_200PPV))

# df_perf_all_summary <- df_perf_all %>%
#   group_by(resolution, metrics) %>%
#   summarise(mean = mean(mean), sd = sd(mean))

## plot all biological replicate results
# df_perf2 is not the correct df
ggplot(df_perf_all, aes(x = resolution, y = mean, col = metrics)) +
  geom_jitter(width = 10, alpha = 0.3, size = 2, stroke = 0) +
  geom_line(data = df_perf_all_summary, aes(x = resolution, y = mean, col = metrics)) +
  geom_point(data = df_perf_all_summary, aes(x = resolution, y = mean, col = metrics), size = 1) +
  geom_errorbar(data = df_perf_all_summary, aes(x = resolution, y = mean, ymin = mean-sd, ymax = mean+sd, col = metrics), width = 10) +
  scale_x_continuous(breaks = unique(df_perf_all$resolution)) + 
  ylim(0,1) +
  labs(title = "Performance",
       x = "Rasterization Resolution",
       y = "Performance",
       col = "Metric") +
  theme_bw()
        
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "biological_replicates_perf_metric_summary.pdf")), width = 6, heigh = 5, dpi = 300)


# upset plot --------------------------------------------------------------

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
intersect(sc_df, res_50)

## identify which genes + visualize their patterns
# at 50um
# Rgs5 - only FN
# Avpr2 - only FP

# plot different resolutions ----------------------------------------------

gene <- rownames(spe)

df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], gene = colSums(assay(spe, "lognorm")[gene,]))

df <- data.frame(x = spatialCoords(spe)[,1], y = spatialCoords(spe)[,2], gene = (assay(spe, "lognorm")["Rgs5",]))

plot_sc <- ggplot(df, aes(x = x, y = y, col = gene)) +
      coord_fixed() +
      geom_point(size = 1, stroke = 0) +
      scale_color_viridis_c() +
      theme_bw() +
        theme(
           #legend.position="none",
           panel.grid = element_blank(),
           axis.title = element_blank(),
           axis.text = element_blank(),
           axis.ticks = element_blank(),
        )

# 50 
rastGexp50 <- SEraster::rasterizeGeneExpression(merfish_mousePOA, assay_name="volnorm", resolution = 50)
rastGexp100 <- SEraster::rasterizeGeneExpression(merfish_mousePOA, assay_name="volnorm", resolution = 100)
rastGexp200 <- SEraster::rasterizeGeneExpression(merfish_mousePOA, assay_name="volnorm", resolution = 200)

only_FN_at_50 <- "Rgs5"

# only FN at 50
plot_50 <- SEraster::plotRaster(rastGexp50, feature_name = only_FN_at_50, name = only_FN_at_50)
plot_100 <- SEraster::plotRaster(rastGexp100, feature_name = only_FN_at_50, name = only_FN_at_50)
plot_200 <- SEraster::plotRaster(rastGexp200, feature_name = only_FN_at_50, name = only_FN_at_50)

ggarrange(plot_sc, 
          plot_50, 
          plot_100, 
          plot_200, 
          labels = c("Single Cell", "50 um", "100 um", "200 um"), 
          ncol = 2, 
          nrow = 2) #,
          #top = "Gene expression pattern of Rgs5 across resolutions")

# only FP at 50
SEraster::plotRaster(rastGexp50, feature_name = "Avpr2", name = "Avpr2")
SEraster::plotRaster(rastGexp100, feature_name = "Avpr2", name = "Avpr2")
SEraster::plotRaster(rastGexp200, feature_name = "Avpr2", name = "Avpr2")

# only TP at 200 um
SEraster::plotRaster(rastGexp50, feature_name = "Klf4", name = "Klf4")
SEraster::plotRaster(rastGexp100, feature_name = "Klf4", name = "Klf4")
SEraster::plotRaster(rastGexp200, feature_name = "Klf4", name = "Klf4")


# plot performance for biological replicates ------------------------------


# Questions ---------------------------------------------------------------

# Is the performance consistent across bregma sections for each mouse OR across mice for each bregma section?
# only tested for 1_Female_Naive_0.29
# PPV slightly higher, TPR looks the correct

# How does the correlation of nnSVG outputs look like single-cell vs. rasterized resolutions (e.g. gene ranking)?


# UpSetR
# Jake R Conway, Alexander Lex, Nils Gehlenborg UpSetR: An R Package for the Visualization of Intersecting Sets and their Properties doi: 
# https://doi.org/10.1093/bioinformatics/btx364
        
        
        
        
        
        
        