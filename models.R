source("../layers.R")

library(torch)


simple_GC_DEC <- nn_module(
  initialize = function(nfeat, nhid, alpha = 0.2){
    self$gc = GraphConvolution(nfeat, nhid)
    self$nhid = nhid
    self$alpha= alpha
    self$mu = NULL
  },
  
  forward = function(x, adj){
    x = self$gc(x, adj)
    q = 1/((1+torch_sum((x$unsqueenz(1)-self$mu)^2, dim = 2)/self$alpha) + 1e-08)
    q = q^(self$alpha+1)/2
    q = q/torch_sum(q, dim = 1, keepdim = TRUE)
    mylist = list("x" = x, "q" = q)
    return(mylist)
  },
  
  loss_function = function(target, pred){
    loss = torch_mean(torch_sum(target*torch_log(target/(pred+1e-06)), dim = 1))
    return(loss)
  },
  
  target_distribution = function(q){
    p = q^2/torch_sum(q, dim = 0)
    p = p/torch_sum(p, dim = 1, keepdim = TRUE)
    return(p)
  },
  
  fit = function(X, 
                 adj, 
                 lr = 0.001, 
                 max_epochs = 5000, 
                 update_interval = 3, 
                 trajectory_interval = 50,
                 weight_decay = 5e-04, 
                 opt = "sgd",
                 init = "louvain", 
                 n_neighbors = 10,
                 res = 0.4, 
                 n_clusters = 10, 
                 init_spa = TRUE, 
                 tol = 1e-03){
    # TO REVISE INITIATING TRA
    self$trajectory = NULL
    if (opt == "sgd"){
      optmz = optim_sgd(self$parameters, lr = lr, momentum = 0.9)
    }
    else if (opt == "admin"){
      optmz = optim_adam(self$parameters, lr = lr, weight_decay = weight_decay)
    }
    
    features = self$gc(torch_tensor(X), torch_tensor(adj))
    
    if (init == "kmeans"){
      cat("Initializing cluster centers with kmeans, n_clusters known")
    } 
    else if (init == "louvain"){
      cat("Initializing cluster centers with louvain, resolution = ", res)
      if (init_spa){
        mat = as.matrix(features)
      }
      else {
        mat = X
      }
    }
  }
)


gc <- GraphConvolution(50,50)
features = gc(torch_tensor(embed), torch_tensor(adj.exp))
f = as.matrix(features)
library(bluster)
set.seed(123) # just in case there are ties.
np <- NNGraphParam(k=10, cluster.fun="louvain", cluster.args = list(resolution=1.1))
np
graph.out <- clusterRows(f, np)
# plotUMAP(clf$spaData, colour_by=I(graph.out))
table(graph.out)
colLabels(clf$spaData) <- factor(graph.out)
library(ggspavis)
plotSpots(clf$spaData, annotate = "label", 
          palette = "libd_layer_colors")
plotDimRed(clf$spaData, type = "PCA", 
           annotate = "label", palette = "libd_layer_colors")



set.seed(123)
np <- NNGraphParam(k=10, cluster.fun="louvain", cluster.args = list(resolution=1.1))
np
graph.out <- clusterRows(mat1, np)
# plotUMAP(clf$spaData, colour_by=I(graph.out))
table(graph.out)
colLabels(spaData) <- factor(graph.out)
plotSpots(spaData, annotate = "label", 
          palette = "libd_layer_colors")




set.seed(90)
np <- NNGraphParam(k=10, cluster.fun="louvain", cluster.args = list(resolution=1.1))
np
graph.out <- clusterRows(mat, np)
# plotUMAP(clf$spaData, colour_by=I(graph.out))
table(graph.out)
colLabels(spaData1) <- factor(graph.out)
plotSpots(spaData1, annotate = "label", 
          palette = "libd_layer_colors")

