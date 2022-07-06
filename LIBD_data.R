library(ggplot2)
library(grid)
library(SpatialExperiment)
library(spatialLIBD)
library(tiff)
# library(foreach)
# library(doParallel)
library(scran)
library(scuttle)
library(ggspavis)
library(torch)
source("../calculate_adj.R")

## Connect to ExperimentHub
# ehub <- ExperimentHub::ExperimentHub()
#> snapshotDate(): 2022-04-26
## Download the small example sce data
# sce <- fetch_data(type = "sce", eh = ehub)
#> 2022-05-03 10:44:08 loading file /home/biocbuild/.cache/R/BiocFileCache/3f32f6ee28586_sce_sub_for_vignette.Rdata%3Fdl%3D1

## Convert to a SpatialExperiment object
# spe <- sce_to_spe(sce)
# spe
# colData(spe)
# rowData(spe)
# dim(spe)
# subset a specific slice
# spe[,colData(spe)$sample_id == "151673"]
# assay(spe)
# spatialLIBD: layer annotation

# spaData <- spe[,colData(spe)$sample_id == "151673"]
# vis_clus(
#   spe = spe,
#   clustervar = "layer_guess_reordered",
#   sampleid = "151673",
#   colors = libd_layer_colors,
#   ... = " LIBD Layers"
# )
# colData(spaData)
# rowData(spaData)
# dim(spaData)
# assay(spaData)
# spaData
# saveRDS(spaData, file = "spaData.rds")
spaData <- readRDS("spaData.rds")


# url <- "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151673_full_image.tif"
# spaData <- addImg(spaData, 
#               sample_id = "151673", 
#               image_id = "fullres",
#               imageSource = url, 
#               scaleFactor = NA_real_, 
#               load = TRUE)
# 
# imgData(spaData)[2,4] <- 1

# img <- imgRaster(spaData, 
#                  sample_id = "151673", 
#                  image_id = "fullres")

img <- imgRaster(spaData, 
                 sample_id = "151673", 
                 image_id = "lowres")
png("histology.png")
plot(img)
dev.off()

colnames(colData(spaData))

# select captured samples
spaData <- spaData[, colData(spaData)$in_tissue==1]
rownames(spaData) <- lapply(rownames(spaData), toupper)
colnames(spaData) <- lapply(colnames(spaData), toupper)

# set coordinates
scale.fac <- scaleFactors(spaData)
x_array <- colData(spaData)$array_row
y_array <- colData(spaData)$array_col
x_pixel <- spatialCoords(spaData)[,2]*scale.fac
y_pixel <- spatialCoords(spaData)[,1]*scale.fac


# add rgb dimension
# array: rgb*pcol*prow
# after permutation:xp*yp*rgb
img.rgb <- aperm(array(col2rgb(img), dim=c(3,ncol(img),nrow(img))),c(3,2,1))
# img.rgb <- aperm(array(col2rgb(img), dim=c(3,ncol(img),nrow(img))),c(2,3,1))


# test coordinates on the image
# a <- as.integer(20*scale.fac)
# img.new <- img
# for (i in 1:length(x_pixel)){
#   x <- x_pixel[i]
#   y <- y_pixel[i]
#   img.new[as.integer(x-a):as.integer(x+a),as.integer(y-a):as.integer(y+a)]="#000000ff"
# }
# png("histology_map.png")
# plot(img.new)
# dev.off()


# img.rgb.new <- img.rgb
# for (i in 1:length(x_pixel)){
#   x <- x_pixel[i]
#   y <- y_pixel[i]
#   img.rgb.new[as.integer(x-a):as.integer(x+a),as.integer(y-a):as.integer(y+a)]=0
# }
# png("histology_map.png")
# plot(img.new)
# dev.off()


# 1. calculate adjacent matrix
# s parameter determines weight given to histology when calculating Euclidean distance between two spots
# s = 1: histology pixel intensity value has the same scale variance as the x,y coordinates
# higher s: higher weights to histology
s <- 1
# b parameter determines the area of each spot when extracting color intensity
b <- 49*scale.fac
adj <- calculate.adj.matrix(x=x_pixel, 
                            y=y_pixel, 
                            x_pixel=x_pixel, 
                            y_pixel=y_pixel, 
                            image=img.rgb, 
                            beta=b, 
                            alpha=s, 
                            histology=TRUE)
# adj <- calculate.adj.matrix(x=x_pixel,
#                             y=y_pixel,
#                             x_pixel=x_pixel,
#                             y_pixel=y_pixel,
#                             image=img.rgb,
#                             beta=b,
#                             alpha=s,
#                             histology=TRUE)
# write.csv(adj, "adj.csv", row.names = FALSE, col.names = FALSE)

# 2. Spatial domain detection
# 2.1 Expression data processing
# adj <- read.csv("adj.csv", header = TRUE)
#filter genes; PYTHON: scanpy filter: min cells, max cells, min counts, max counts
#exclude genes expressed in fewer than three spots
per.gene <- perFeatureQCMetrics(spaData)
discard <- per.gene$mean*dim(spaData)[2] < 3
spaData <- spaData[!discard,]
#exclude special genes: ERCC, MT
# rowData(spadata) --> gene_biotype, only protein_coding
is.ERCC <- base::grepl("ERCC", rownames(spaData))
is.mito <- base::grepl("MT-", rownames(spaData))
spaData <- spaData[!is.ERCC,]
spaData <- spaData[!is.mito,]
#normalize and take log for UMI (exist in spatialLIBD)
# set.seed(100)
# cluster <- quickCluster(spaData) 
# spaData <- computeSumFactors(spaData, cluster=cluster, min.mean=0.1)
# spaData <- logNormCounts(spaData)

  
source("../util.R")
# 2.2 Set hyper-parameters
# p: Percentage of total expression contributed by neighborhoods; average relative contribution of other spots for one spot across all spots
# l: Parameter to control p
p <- 0.5
# Find the l value given p
l <- search.l(p, adj, start = 0.01, end = 1000, tol = 0.01, max.run = 100)
# n.clusters: number of spatial domains wanted
# res: resolution in the initial Louvain's clustering methods. If the number of clusters is known, use search.res() function to search for suitable resolution
# For 151673 slice, set the number of clusters = 7 since this tissue has 7 layers
n.clusters <- 7
seed=100
# Search for suitable resolution
res <- search.res(spaData, adj, l, n.clusters, 
                  start = 0.9, step = 0.1, tol = 5e-3, lr = 0.05, 
                  max.epochs = 20, seed = seed)
source("../spaGCN.R")
# run clustering
set.seed(seed)
torch_manual_seed(seed)
spa.clf <- spaGCN(spaData, adj, init.spa = TRUE, init = "louvain", res = res, tol = 5e-3, lr = 0.05, max.epochs = 200)
spa.clf <- set_l(spa.clf, l)
spa.clf <- train.spaGCN(spa.clf, seed)
spa.pred <- predict.spaGCN(spa.clf)
spa.y_pred <- spa.pred$y_pred
spa.prob <- spa.pred$prob
spa.clf$spaData$spa.y_pred <- as.factor(spa.y_pred)
# Do cluster refinement(optional)
# shape = "hexagon" for Visium data, "square" for ST data
adj.2d <- calculate.adj.matrix(x=x_array, y=y_array, histology = FALSE)
refined.pred <- refine(sample_id = colnames(spaData), pred = spa.clf$spaData$spa.y_pred, 
                       dis = adj.2d, shape = "hexagon")
spa.clf$spaData$refined.pred <- as.factor(refined.pred)
saveRDS(spa.clf$spaData, file = "results.spa.rds")
spaData <- readRDS("results.spa.rds")

# plot spatial domains
plot_color <- c("#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785")
num_cluster <- length(unique(spaData$spa.y_pred))
colors <- plot_color[1:num_cluster]
png("pred.png")
plotSpots(spaData, annotate = "spa.y_pred",
          palette = colors)
dev.off()
plotDimRed(spaData, type = "PCA",
            annotate = "spa.y_pred", palette = colors)

# plot refined spatial domains
num_cluster <- length(unique(spaData$refined.pred))
colors <- plot_color[1:num_cluster]
png("refined.pred.png")
plotSpots(spaData, annotate = "refined.pred",
          palette = colors)
dev.off()
plotDimRed(spaData, type = "PCA",
            annotate = "refined.pred", palette = colors)
















