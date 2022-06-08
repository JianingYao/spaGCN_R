
spaGCN.train <- function(spaData, adj, 
                         num.pcs = 50, 
                         lr = 0.005,
                         max.epochs = 2000,
                         weight.decay = 0,
                         init.spa = TRUE, 
                         init = "louvain", # louvain or kmeans
                         n_neighbors = 10, # for louvain
                         n.clusters = NULL,# for kmeans
                         res = 0.4,
                         tol = 1e-3,
                         l = l) {
  stopifnot((dim(spaData)[2] == dim(adj)[1]) && (dim(adj)[1]== dim(adj)[2]))
  prcomp()
  # check sparse matrix
  # transform sparse matrix
  # do pca: whether t(logcounts(spaData))
  # embeddings of pca
  if (is.null(l)){
    stop("l should be set before fitting the model")
  }
  adj.exp <- exp(-1*(adj^2)/(2*(l^2)))
}
  
  
