library(spam)
library(scater)
library(bluster)

spaGCN <- function(spaData, 
                   adj, 
                   num.pcs = 50, 
                   lr = 0.005,
                   max.epochs = 2000,
                   weight.decay = 0,
                   opt = "admin",
                   init.spa = TRUE, 
                   init = "louvain", # louvain or kmeans
                   n.neighbors = 10, # for louvain
                   n.clusters = NULL,# for kmeans
                   res = 0.4,
                   tol = 1e-3,
                   l = NULL,
                   adj.exp = NULL){
  structure(list(spaData = spaData,
                 adj = adj,
                 num.pcs = num.pcs,
                 lr = lr,
                 max.epochs = max.epochs,
                 weight.decay = weight.decay, 
                 opt = opt, 
                 init.spa = init.spa,
                 int = init,
                 n.neighbors = n.neighbors,
                 n.clusters = n.clusters,
                 res = res,
                 tol = tol,
                 l = l, 
                 adj.exp = adj.exp),
                 class = "spaGCN")
}

set_l <- function(spaGCN, l){
  spaGCN$l <- l
  invisible(spaGCN)
}

train.spaGCN <- function(clf) {
  stopifnot((dim(clf$spaData)[2] == dim(clf$adj)[1]) && (dim(clf$adj)[1]== dim(clf$adj)[2]))
  set.seed(123)
  # pca_data <- prcomp(t(logcounts(clf$spaData)), rank = 50, scale = TRUE)
  # reducedDims(clf$spaData) <- list(PCA = pca_data$x)

  clf$spaData <- runPCA(clf$spaData, exprs_values = "logcounts", ncomponents = 50,scale = TRUE)
  embed <- reducedDim(clf$spaData, "PCA")
  if (is.null(clf$l)){
    stop("l should be set before fitting the model")
  }
  adj.exp <- exp(-1*(clf$adj^2)/(2*(clf$l^2)))
  clf$adj.exp <- adj.exp
  # train model
  
}
  
  
