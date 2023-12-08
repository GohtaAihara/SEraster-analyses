## simulate true positives and negatives

## Trevor's recommendation, simulate two spatially autocorrelated but independent cell-types
## https://hastie.su.domains/Papers/biodiversity/Biodiversity.pdf

# Set up ------------------------------------------------------------------

setwd("~/Desktop/SEraster")
devtools::document()
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

source("analyses/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(here)
library(dplyr)

par(mfrow=c(1,1))

dataset_name <- "autocorrelation_simulations"

# Generate simulated data -------------------------------------------------

## RandomFields not available anymore
## write own code / copy from their source

# Set parameters
set.seed(10)
N <- 2000  # Number of locations
nugget_variance <- 0.1  # Nugget variance (τ^2)
range_parameter <- 0.5  # Range parameter (κ)
smoothness_parameter <- 0.3  # Smoothness parameter (φ)

# Define locations in [0, 1] × [0, 1]
locations <- cbind(runif(N), runif(N)) 
colnames(locations) <- c('x', 'y')

# Function to calculate Matern covariance
matern_covariance <- function(h, range, smoothness) {
  nu <- smoothness
  term1 <- 2^(1 - nu)/gamma(nu)
  term2 <- (sqrt(2 * nu) * h / range)
  term3 <- besselK(term2, nu)
  return(term1 * (term2^nu) * term3)
}

# Function to calculate the covariance matrix
calculate_covariance_matrix <- function(locations, range, smoothness) {
  N <- nrow(locations)
  cov_matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      h <- sqrt(sum((locations[i,] - locations[j,])^2))
      cov_matrix[i, j] <- matern_covariance(h, range, smoothness)
    }
  }
  return(cov_matrix)
}

# Generate covariance matrix for the Gaussian random field
cov_matrix <- calculate_covariance_matrix(locations, range_parameter, smoothness_parameter)
diag(cov_matrix) <- 1

# Generate realizations of the Gaussian random field
W1 <- MASS::mvrnorm(1, rep(0, N), cov_matrix)
W2 <- MASS::mvrnorm(1, rep(0, N), cov_matrix)

# Generate independent and identically distributed errors
Z1 <- rnorm(N, mean = 0, sd = sqrt(nugget_variance))
Z2 <- rnorm(N, mean = 0, sd = sqrt(nugget_variance))

# Generate realizations of the process Xsi
Xs1 <- W1 + Z1 
Xs2 <- W2 + Z2
Xs3 <- W1 + Z2 # should correlate with Xs1

# Print or plot the results as needed
print(Xs1)
print(Xs2) 

# Set names
rownames(locations) <- names(Xs1) <- names(Xs2) <- paste0('cell', 1:N)

# Confirm Rob and Trevor's original concerns with pixel data
plot(Xs1, Xs2)
cor(Xs1, Xs2) ## indeed not 0
cor.test(Xs1, Xs2) ## indeed not 0
plot(Xs1, Xs3)
cor(Xs1, Xs3) ## high
cor.test(Xs1, Xs3) ## high

par(mfrow=c(1,2), mar=rep(2,4))
MERINGUE::plotEmbedding(locations, col=Xs1, cex=1) ## expression of gene A
MERINGUE::plotEmbedding(locations, col=Xs2, cex=1)  ## expression of gene B

## format into SpatialExperiment
gexp <- rbind(Xs1,Xs2,Xs3)
rownames(gexp) <- paste0("gene",1:dim(gexp)[1])
spe_gexp <- SpatialExperiment::SpatialExperiment(
  assays = list(lognorm = gexp),
  spatialCoords = locations
)

## cell type (each point stores the number of cells)
Ys1 <- round((Xs1 - min(Xs1))*10)
Ys2 <- round((Xs2 - min(Xs2))*10)
Ys3 <- round((Xs3 - min(Xs3))*10)
MERINGUE::plotEmbedding(locations, col=Ys1, cex=1) 

## format into SpatialExperiment
num_ct <- rbind(Ys1,Ys2,Ys3)
rownames(num_ct) <- paste0("celltype",1:dim(num_ct)[1])
spe_num_ct <- SpatialExperiment::SpatialExperiment(
  assays = list(num_ct = num_ct),
  spatialCoords = locations
)
ct_prop <- rowSums(num_ct)/sum(num_ct)

# Run methods -------------------------------------------------------------

## pixel-wise autocorrelation for gene expression
res_list <- seq(0,1-1e-10,by = 0.01)
autocorrelation_results <- do.call(rbind, lapply(res_list, function(res) {
  print(paste0("Resolution = ", res))
  if (res == 0) {
    out_12 <- cor.test(assay(spe_gexp)[1,], assay(spe_gexp)[2,])
    out_13 <- cor.test(assay(spe_gexp)[1,], assay(spe_gexp)[3,])
  } else {
    spe_rast <- SEraster::rasterizeGeneExpression(spe_gexp, resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
    out_12 <- cor.test(assay(spe_rast)[1,], assay(spe_rast)[2,])
    out_13 <- cor.test(assay(spe_rast)[1,], assay(spe_rast)[3,])
  }
  df_12 <- data.frame(
    resolution = res,
    condition = "independent",
    cor_estimate = out_12$estimate,
    cor_pval = out_12$p.value,
    cor_ci_min = out_12$conf.int[1],
    cor_ci_max = out_12$conf.int[2]
  )
  df_13 <- data.frame(
    resolution = res,
    condition = "correlated",
    cor_estimate = out_13$estimate,
    cor_pval = out_13$p.value,
    cor_ci_min = out_13$conf.int[1],
    cor_ci_max = out_13$conf.int[2]
  )
  return(rbind(df_12,df_13))
}))
saveRDS(autocorrelation_results, file = here("outputs", paste0(dataset_name, "_gexp.RDS")))

## cell type colocalization with CooccurrenceAffinity
res_list <- seq(0,1-1e-10,by = 0.01)

## set confidence interval method ("CP", "Blaker", "midQ", or "midP")
CI_method <- "Blaker"
## create a dictionary and index of CI method
CI_method_dict <- c("CP" = 6, "Blaker" = 7, "midQ" = 8, "midP" = 9)
## set confidence interval level (default = 0.95)
CI_lev <- 0.95
## set pval method ("Blaker" or "midP", default = "Blaker")
pval_method <- "Blaker"

cooccurrence_results <- do.call(rbind, lapply(res_list, function(res) {
  print(paste0("Resolution = ", res))
  if (res == 0) {
    mat <- assay(spe_num_ct)
    mat_re <- do.call(rbind, lapply(rownames(spe_num_ct), function(ct_label) {
      ## relative enrichment = celltype observed / celltype expected = celltype observed / (celltype frequency * total # of cells in the pixel)
      mat[ct_label,] / (sum(mat[ct_label,]) / sum(mat) * colSums(mat))
    }))
    rownames(mat_re) <- rownames(mat)
    
    ## binarize (1 if RE >= 1, 0 if RE < 1)
    mat_bin <- ifelse(mat_re >= 1, 1, 0)
    
  } else {
    spe_rast <- SEraster::rasterizeGeneExpression(spe_num_ct, resolution = res, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())
    mat <- assay(spe_rast)
    mat_re <- do.call(rbind, lapply(rownames(spe_rast), function(ct_label) {
      ## relative enrichment = celltype observed / celltype expected = celltype observed / (celltype frequency * total # of cells in the pixel)
      mat[ct_label,] / (sum(mat[ct_label,]) / sum(mat) * colSums(mat))
    }))
    rownames(mat_re) <- rownames(mat)
    
    ## binarize (1 if RE >= 1, 0 if RE < 1)
    mat_bin <- ifelse(mat_re >= 1, 1, 0)
  }
  
  ## get pair combinations
  non_self <- combn(rownames(spe_num_ct), 2, simplify = FALSE)
  self <- lapply(rownames(spe_num_ct), function(ct_label) c(ct_label, ct_label))
  pairs <- c(non_self, self)
  
  ## multiply pixel values for each pair of cell types
  affinity_results <- do.call(rbind, lapply(pairs, function(pair) {
    ## create a 2x2 contingency table of counts
    cont_tab <- table(factor(mat_bin[pair[1],], levels = c(0,1)), factor(mat_bin[pair[2],], levels = c(0,1)))
    X <- cont_tab[2,2]
    mA <- sum(cont_tab[2,])
    mB <- sum(cont_tab[,2])
    N <- sum(cont_tab)
    
    out <- CooccurrenceAffinity::ML.Alpha(X,c(mA,mB,N), lev = CI_lev, pvalType = pval_method)
    
    ## set index for the chosen CI method
    CI_idx <- CI_method_dict[CI_method][[1]]
    
    return(data.frame(
      resolution = res,
      pair = paste(pair, collapse = " & "),
      celltypeA = pair[1],
      celltypeB = pair[2],
      X = X,
      mA = mA,
      mB = mB,
      N = N,
      alpha = out$est,
      ci.min = out[CI_idx][[1]][1], 
      ci.max = out[CI_idx][[1]][2], 
      pval = out$pval)
    )
  }))
}))
saveRDS(cooccurrence_results, file = here("outputs", paste0(dataset_name, "_cooccurrence.RDS")))

# Plot --------------------------------------------------------------------
## pixel-wise autocorrelation for gene expression
df <- readRDS(file = here("outputs", paste0(dataset_name, "_gexp.RDS")))

ggplot(df, aes(x = resolution, y = cor_estimate, col = condition)) +
  geom_point() +
  geom_line() +
  geom_errorbar(data = df, aes(ymin = cor_ci_min, ymax = cor_ci_max)) +
  ylim(c(-1,1)) +
  labs(x = "Rasterization Resolution",
       y = "Pearson's correlation") +
  theme_bw()
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_gexp.pdf")), width = 8, height = 5, dpi = 300)

## cell type colocalization
df <- readRDS(file = here("outputs", paste0(dataset_name, "_cooccurrence.RDS")))
df <- df[df$pair %in% c("celltype1 & celltype2", "celltype1 & celltype3"),]

ggplot(df, aes(x = resolution, y = alpha, col = pair)) +
  geom_point() +
  geom_line() +
  geom_errorbar(data = df, aes(ymin = ci.min, ymax = ci.max)) +
  labs(x = "Rasterization Resolution",
       y = "Alpha MLE") +
  theme_bw()
  