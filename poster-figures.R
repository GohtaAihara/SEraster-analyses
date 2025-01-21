
# MERGINUE SVG analysis ---------------------------------------------------

suppressMessages(library(SpatialExperiment))
suppressMessages(library(SEraster))
suppressMessages(library(Matrix))
suppressMessages(library(ggplot2))

dataset_name <- "merfish_mousePOA"
data("merfish_mousePOA")
class(merfish_mousePOA)

## pos is a matrix specifying a cell's x and y positions
pos <- spatialCoords(merfish_mousePOA)
## counts is a genes-by-cells matrix where the gene expression was normalized by cell volume and scaled by 1000
counts <- assays(merfish_mousePOA)$volnorm
df <- data.frame(pos, total_gexp = colSums(counts))
## plot
ggplot(df, aes(x = x, y = y, color = total_gexp)) +
  coord_fixed() +
  scale_color_viridis_c() +
  geom_point(size = 0.5) +
  labs(title = "Single-cell resolution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        axis.text = element_blank()
       )
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_sc-total_gexp.pdf")), width = 6, heigh = 5, dpi = 300)

## let's rasterize the dataset at 100 um resolution.
res <- 100
rastGexp <- SEraster::rasterizeGeneExpression(merfish_mousePOA, 
                                              assay_name = "volnorm", 
                                              resolution = res, 
                                              square = FALSE, 
                                              fun = "mean")
SEraster::plotRaster(rastGexp, plotTitle = "200 um")
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_200um_gexp.pdf")), width = 6, heigh = 5, dpi = 300)

suppressMessages(library(MERINGUE))
## Get neighbor-relationships
rastPos <- spatialCoords(rastGexp)
w <- MERINGUE::getSpatialNeighbors(rastPos, filterDist = 105)
plotNetwork(rastPos, w)
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_spatial_neighbors.pdf")), width = 6, heigh = 5, dpi = 300)

## rasterized genes-by-cells matrix, with cell pixels at 100 um resolution
mat <- assays(rastGexp)$pixel
I <- MERINGUE::getSpatialPatterns(mat, w)

results.filter <- MERINGUE::filterSpatialPatterns(mat = mat,
                                                  I = I,
                                                  w = w,
                                                  adjustPv = TRUE,
                                                  alpha = 0.05,
                                                  minPercentCells = 0.05,
                                                  verbose = TRUE, 
                                                  details = TRUE)


## sort spatially auto-correlated genes with Moran's I statistic
sortedResults <- results.filter[order(results.filter$observed, decreasing = TRUE), ]
## select 3 SVGs with the highest Moran's I
sgenes <- rownames(sortedResults)
for (sg in head(sgenes, 4)) {
  ## plot gene expression pattern
  plt <- SEraster::plotRaster(rastGexp, 
                              feature_name = sg,
                              name = sg,
                              plotTitle = sg)
  show(plt)
  ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_", sg, "_rasterized_gene.pdf")), width = 6, heigh = 5, dpi = 300)
  
}

SEraster::plotRaster(rastGexp, 
                     feature_name = "Cenpe",
                     name = "Cenpe",
                     plotTitle = "Cenpe")
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_Cenpe_rasterized_gene.pdf")), width = 6, heigh = 5, dpi = 300)



# Spatial Bootstrapping ---------------------------------------------------

suppressMessages(library(SpatialExperiment))
suppressMessages(library(SEraster))
suppressMessages(library(Matrix))
suppressMessages(library(ggplot2))

data("merfish_mousePOA")
class(merfish_mousePOA)

ct <- merfish_mousePOA$celltype; names(ct) <- colnames(merfish_mousePOA)
ct <- as.factor(ct)
head(ct)
length(ct)
length(levels(ct))

## plot
suppressMessages(library(ggplot2))
pos <- spatialCoords(merfish_mousePOA)
df <- data.frame(pos, ct = ct)
ggplot(df, aes(x = x, y = y, color = ct)) +
  coord_fixed() +
  geom_point(size = 0.5) +
  labs(title = "Single-cell resolution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        axis.text = element_blank()
  )
##ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_sc-celltype.pdf")), width = 6, heigh = 5, dpi = 300)

suppressMessages(library(SEraster))
## rasterize at 100um resolution with hexagons
rastCt <- SEraster::rasterizeCellType(merfish_mousePOA,
                                      col_name = "celltype",
                                      resolution = 400,
                                      fun = "sum",
                                      square = FALSE)
## plot
SEraster::plotRaster(rastCt, plotTitle = "400 um", name = "Total cells")
##ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_400um_celltype.pdf")), width = 6, heigh = 5, dpi = 300)

## find the hexagonal pixel with the most number of cells
maxCells <- max(rastCt$num_cell) ## 102
pixel <- colnames(rastCt)[colData(rastCt)$num_cell == maxCells]
pixelIdx <- which(colnames(rastCt) == pixel)
print(paste(pixel, " at index ", pixelIdx, " has ", maxCells, " cells."))

## 54th pixel
cells <- colData(rastCt)$cellID_list[[pixelIdx]]
## double check if the number of cells in this pixel matches maxCell
length(cells) == maxCells
## visualize the cells in the 54th pixel
dfsub <- data.frame(pos[cells,], ct = ct[cells])
ggplot(dfsub, aes(x = x, y = y, color = ct)) +
  coord_fixed() +
  geom_point(size = 2.5) +
  labs(title = "400um") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        axis.text = element_blank())
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_400um_celltype_p26.pdf")), width = 6, heigh = 5, dpi = 300)


ctsub <- ct
ctsub[!(names(ctsub) %in% cells)] <- NA
df <- data.frame(pos, ctsub = ctsub)
ggplot(df, aes(x = x, y = y, color = ctsub)) +
  coord_fixed() +
  geom_point(size = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        axis.text = element_blank())
ggsave(filename = here("plots", dataset_name, paste0(dataset_name, "_400um_celltype_isolated.pdf")), width = 6, heigh = 5, dpi = 300)




