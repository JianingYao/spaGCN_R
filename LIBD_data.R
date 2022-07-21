library(ggplot2)
library(grid)
library(SpatialExperiment)
library(spatialLIBD)
# library(tiff)
# library(foreach)
# library(doParallel)
library(scran)
library(scuttle)
library(ggspavis)
library(torch)
library(raster)

c = c(100, 123, 345, 478, 763, 999, 5674, 12321)

for (i in c){
  seed = i

sink("out-100(1).txt")
source("../calculate_adj.R")
source("../util.R")
source("../spaGCN.R")



# # Connect to ExperimentHub
# ehub <- ExperimentHub::ExperimentHub()
# # Download the small example sce data
# sce <- fetch_data(type = "sce", eh = ehub)
# 
# # Convert to a SpatialExperiment object
# spe <- sce_to_spe(sce)
# # subset a specific slice: "151507" "151508" "151509" "151510" "151669" "151670" "151671" "151672" "151673" "151674" "151675" "151676"
# sample.id <- "151673"
# spe[,colData(spe)$sample_id == sample.id]
# 
# spaData <- spe[,colData(spe)$sample_id == sample.id]
# vis_clus(
#   spe = spe,
#   clustervar = "layer_guess_reordered",
#   sampleid = sample.id,
#   colors = libd_layer_colors,
#   ... = " LIBD Layers"
# )

# saveRDS(spaData, file = "spaData-151673.rds")
  sample.id <- "151673"
  spaData <- readRDS(paste0("spaData-",sample.id,".rds"))
  
  
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
  
  # img <- imgRaster(spaData,
  #                  sample_id = sample.id,
  #                  image_id = "lowres")
  # pdf("histology.pdf")
  # plot(img)
  # dev.off()
  
  colnames(colData(spaData))
  
  # select captured samples
  spaData <- spaData[, colData(spaData)$in_tissue == 1]
  rownames(spaData) <- lapply(rownames(spaData), toupper)
  colnames(spaData) <- lapply(colnames(spaData), toupper)
  
  # set coordinates
  scale.fac <- scaleFactors(spaData)
  x_array <- colData(spaData)$array_row
  y_array <- colData(spaData)$array_col
  # x_pixel <- spatialCoords(spaData)[, 2] * scale.fac
  # y_pixel <- spatialCoords(spaData)[, 1] * scale.fac
  x_pixel <- spatialCoords(spaData)[, 2]
  y_pixel <- spatialCoords(spaData)[, 1]
  
  # add rgb dimension
  # array: rgb*pcol*prow
  # after permutation:xp*yp*rgb
  # img.rgb <-
  #   aperm(array(col2rgb(img), dim = c(3, ncol(img), nrow(img))), c(3, 2, 1))
  
  # full resolution
  # img <- brick(x = "histology.tif")
  # ext <- extent(img)
  # img.rgb <- raster::extract(img, ext)
  # img.rgb[, c(1,3)] <- img.rgb[, c(3,1)]
  # img.rgb <- aperm(`dim<-`(t(img.rgb), c(3, dim(img)[1], dim(img)[2])), c(3, 2, 1))
  # saveRDS(img.rgb, file = paste0("img.rgb-",sample.id,".rds"))
  img.rgb <- readRDS(paste0("img.rgb-",sample.id,".rds"))
  
  # test coordinates on the image
  # a <- as.integer(20*scale.fac)
  # img.new <- img
  # for (i in 1:length(x_pixel)){
  #   x <- x_pixel[i]
  #   y <- y_pixel[i]
  #   img.new[as.integer(x-a):as.integer(x+a),as.integer(y-a):as.integer(y+a)]="#000000ff"
  # }
  # pdf("histology_map.pdf")
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
  # b <- 49 * scale.fac
  b <- 49
  # browser()
  # adj <- calculate.adj.matrix(
  #   x = x_pixel,
  #   y = y_pixel,
  #   x_pixel = x_pixel,
  #   y_pixel = y_pixel,
  #   image = img.rgb,
  #   beta = b,
  #   alpha = s,
  #   histology = TRUE
  # )
  adj <- as.matrix(read.csv("adj.csv", header = FALSE))

  # 2. Spatial domain detection
  # 2.1 Expression data processing
  #filter genes; PYTHON: scanpy filter: min cells, max cells, min counts, max counts
  #exclude genes expressed in fewer than three spots; exclude special genes: ERCC, MT
  # rowData(spadata) --> gene_biotype, only protein_coding
  is.spike <- base::grepl("^ERCC", rownames(spaData))
  is.mito <- base::grepl("^mt-", rownames(spaData))
  per.cell <- perCellQCMetrics(spaData)
  per.gene <- perFeatureQCMetrics(spaData)
  libsize.drop <- isOutlier(per.cell$total, nmads = 3, type = "lower", log = TRUE)
  numcells <- nexprs(spaData, byrow = TRUE)
  keep <- numcells >= 3
  spaData <- spaData[!(is.spike | is.mito), ]
  spaData <- spaData[keep, ]
  # spaData <- spaData[rowData(spaData)$gene_biotype == "protein_coding", ]
  #normalize and take log for UMI (exist in spatialLIBD??)
  set.seed(seed)
  cluster <- quickCluster(spaData)
  spaData <- computeSumFactors(spaData, cluster=cluster, min.mean=0.1)
  spaData <- logNormCounts(spaData)
  
  
  # 2.2 Set hyper-parameters
  # p: Percentage of total expression contributed by neighborhoods; average relative contribution of other spots for one spot across all spots
  # l: Parameter to control p
  p <- 0.5
  # Find the l value given p
  l <-
    search.l(
      p,
      adj,
      start = 0.01,
      end = 1000,
      tol = 0.01,
      max.run = 100
    )
  # n.clusters: number of spatial domains wanted
  # res: resolution in the initial Louvain's clustering methods. If the number of clusters is known, use search.res() function to search for suitable resolution
  # For 151673 slice, set the number of clusters = 7 since this tissue has 7 layers
  n.clusters <- 7
  # Search for suitable resolution
  res <- search.res(
    spaData,
    adj,
    l,
    n.clusters,
    start = 1.1,
    step = 0.1,
    tol = 5e-3,
    lr = 0.05,
    max.epochs = 10,
    seed = seed
  )
  # run clustering
  set.seed(seed)
  torch_manual_seed(seed)
  spa.clf <-
    spaGCN(
      spaData,
      adj,
      init.spa = TRUE,
      init = "louvain",
      res = res,
      tol = 5e-3,
      lr = 0.05,
      max.epochs = 15
    )
  spa.clf <- set_l(spa.clf, l)
  spa.clf <- train.spaGCN(spa.clf, seed)
  spa.pred <- predict.spaGCN(spa.clf)
  spa.y_pred <- spa.pred$y_pred
  spa.prob <- spa.pred$prob
  spa.clf$spaData$spa.y_pred <- as.factor(spa.y_pred)
  # Do cluster refinement(optional)
  # shape = "hexagon" for Visium data, "square" for ST data
  # adj.2d <-
  #   calculate.adj.matrix(x = x_array, y = y_array, histology = FALSE)
  adj.2d <- as.matrix(read.csv("adj.2d.csv", header = FALSE))
  
  refined.pred <-
    refine(
      sample_id = colnames(spaData),
      pred = spa.clf$spaData$spa.y_pred,
      dis = adj.2d,
      shape = "hexagon"
    )
  spa.clf$spaData$refined.pred <- as.factor(refined.pred)
  saveRDS(spa.clf$spaData, file = "results.spa.rds")
  spaData <- readRDS("results.spa.rds")
  
  # plot spatial domains
  plot_color <-
    c(
      "#F56867",
      "#FEB915",
      "#C798EE",
      "#59BE86",
      "#7495D3",
      "#D1D1D1",
      "#6D1A9C",
      "#15821E",
      "#3A84E6",
      "#997273",
      "#787878",
      "#DB4C6C",
      "#9E7A7A",
      "#554236",
      "#AF5F3C",
      "#93796C",
      "#F9BD3F",
      "#DAB370",
      "#877F6C",
      "#268785"
    )
  num_cluster <- length(unique(spaData$spa.y_pred))
  colors <- plot_color[1:num_cluster]
  pdf(file = paste0("set_pred_", seed, ".pdf"))
  print(plotSpots(spaData, annotate = "spa.y_pred",
            palette = colors, size = 1))
  dev.off()
  plotDimRed(spaData,
             type = "PCA",
             annotate = "spa.y_pred",
             palette = colors)
  
  # plot refined spatial domains
  num_cluster <- length(unique(spaData$refined.pred))
  colors <- plot_color[1:num_cluster]
  pdf(file = paste0("set_refined.pred_", seed, ".pdf"))
  print(plotSpots(spaData, annotate = "refined.pred",
            palette = colors, size = 1))
  dev.off()
  plotDimRed(spaData,
             type = "PCA",
             annotate = "refined.pred",
             palette = colors)
  sink()
  rm(list=ls())
  gc()
  c = c(100, 123, 345, 478, 763, 999, 5674, 12321)
}




