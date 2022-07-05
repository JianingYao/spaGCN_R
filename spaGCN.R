library(spam)
library(scater)
library(bluster)
library(torch)

source("../models.R")


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
                   adj.exp = NULL,
                   model = NULL){
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
                 adj.exp = adj.exp,
                 model = model, 
                 embed = NULL),
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
  clf$embed <- embed
  if (is.null(clf$l)){
    stop("l should be set before fitting the model")
  }
  adj.exp <- exp(-1*(clf$adj^2)/(2*(clf$l^2)))
  clf$adj.exp <- adj.exp
  # train model
  clf$model=simple_GC_DEC(dim(embed)[2], dim(embed)[2])
  clf$model$fit(embed, adj.exp, lr = clf$lr, max_epochs = clf$max.epochs, 
                weight_decay = clf$weight.decay, opt = clf$opt, init_spa = clf$init.spa,
                init = clf$init, n_neighbors = clf$n.neighbors, n_clusters = clf$n.clusters, 
                res = clf$res, tol = clf$tol)
  invisible(clf)
}


predict.spaGCN <- function(clf){
  mylist = clf$model$predict(clf$embed, clf$adj.exp)
  z = mylist$x
  q = mylist$q
  y_pred = as.array(torch_argmax(z, dim = 2)$data()$cpu())
  prob = as.array(q)
  mylist = list("y_pred" = y_pred, "prob" = prob)
  return(mylist)
}
  











