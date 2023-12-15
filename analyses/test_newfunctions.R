setwd("~/Desktop/SEraster")
devtools::document()
devtools::load_all()

setwd("~/Desktop/SEraster-analyses/")

library(here)

dataset_name = "codex_humanSpleen"
method = "CooccurrenceAffinity"

res <- 100
spe_rast <- readRDS(file = here("outputs", paste0(dataset_name, "_rasterized_resolution_", res, "_binarized.RDS")))
ct_label <- "Fol B cells"

plotRaster(spe_rast, assay_name = "pixelval", feature_name = "Fol B cells", name = "cell-type counts", option = "inferno")
plotRaster(spe_rast, assay_name = "re", feature_name = "Fol B cells", name = "RE", option = "inferno")
plotRaster(spe_rast, assay_name = "bin", feature_name = "Fol B cells", factor_levels = c(0,1), name = "binarized", option = "inferno")
