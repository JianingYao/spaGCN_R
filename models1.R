source("../layers.R")

library(torch)
library(bluster)
library(tidyverse)
library(ggspavis)



simple_GC_DEC <- nn_module(
  initialize = function(nfeat, nhid, alpha = 0.2){
    self$gc = GraphConvolution(nfeat, nhid)
    self$nhid = nhid
    self$alpha= alpha
  },
  
  forward = function(x, adj){
    x = self$gc(x, adj)
    print("x is")
    print(x)
    q = 1/((1+torch_sum((x$unsqueeze(2)-self$mu)^2, dim = 3)/self$alpha) + 1e-08)
    q = q^(self$alpha+1)/2
    q = q/torch_sum(q, dim = 2, keepdim = TRUE)
    mylist = list("x" = x, "q" = q)
    return(mylist)
  }
)  

loss_function = function(target, pred){
  loss = torch_mean(torch_sum(target*torch_log(target/(pred+1e-06)), dim = 2))
  print("loss is")
  print(loss)
  return(loss)
}
  
target_distribution = function(q){
  p = q^2/torch_sum(q, dim = 1)
  p = p/torch_sum(p, dim = 2, keepdim = TRUE)
  return(p)
}
  
fit.simple_GC_DEC = function(X,
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
               tol = 1e-03,
               seed = seed,
               spaData = spaData,
               model = model
              ){
  model$trajectory = NULL
  if (opt == "sgd") {
    optmz = optim_sgd(model$parameters, lr = lr, momentum = 0.9)
  } else if (opt == "admin") {
    optmz = optim_adam(model$parameters, lr = lr, weight_decay = weight_decay)
  }
  print("self$parameters is")
  print(model$parameters)
  features = model$gc(torch_tensor(X), torch_tensor(adj))
  if (init == "kmeans") {
    cat("Initializing cluster centers with kmeans, n_clusters known")
  } else if (init == "louvain") {
    cat("Initializing cluster centers with louvain, resolution = ", res)
    if (init_spa) {
      mat = as.matrix(features)
    } else {
      mat = X
    }
    set.seed(seed)
    np <- NNGraphParam(k = n_neighbors, cluster.fun = "louvain", cluster.args = list(resolution = res))
    y_pred <- clusterRows(mat, np)
    model$n_clusters <- length(unique(y_pred))
  }
    
  y_lastPred = y_pred
  model$mu = NULL
  model$mu = nn_parameter(torch_empty(model$n_clusters, model$nhid))
  print("self$mu = nn_parameter(torch_empty(self$n_clusters, self$nhid)) is")
  print(model$mu)
  X = torch_tensor(X)
  adj = torch_tensor(adj)
  model$trajectory = c(model$trajectory, y_lastPred)
  features = as.data.frame(as.matrix(features))
  Group = as.data.frame(y_pred)
  colnames(Group) = "groups"
  Mergefeature = cbind(features, Group)
  cluster_centers = as.data.frame(Mergefeature %>% group_by(groups) %>% summarize_all(mean))
  cluster_centers$groups = NULL
    
  # print and save initialized clustering
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
  colors <- plot_color[1:model$n_clusters]
  print("self$n_clusters is ")
  print(model$n_clusters)
  colLabels(spaData) <- factor(y_lastPred)
  pdf(file = paste0("0cluster_", seed, ".pdf"))
  plotSpots(spaData, annotate = "label",
              palette = colors, size = 1)
  dev.off()
    
  model$mu$data()$copy_(torch_tensor(as.matrix(cluster_centers)))
  print("self$mu$data()$copy_(torch_tensor(as.matrix(cluster_centers))) is")
  print(model$mu)
  model$train(TRUE)
  for (epoch in 0:(max_epochs-1)){
    if (epoch%%update_interval == 0){
      result = model$forward(X, adj)
      q = result$q
      p = target_distribution(q)$data()
    }
    if (epoch%%10==0){
      cat("Epoch ", epoch)
    }
    optmz$zero_grad()
    result3 = model(X, adj)
    z = result3$x
    q = result3$q
    loss = loss_function(p,q)
    loss$backward()
    optmz$step()
    if (epoch%%trajectory_interval == 0){
      model$trajectory = rbind(model$trajectory, as.array(torch_argmax(q, dim = 2)$data()$cpu()))
    }
      
    # check stop criterion
    y_pred = as.array(torch_argmax(q, dim = 2)$data()$cpu())
    delta.label = sum(y_pred != y_lastPred)/dim(X)[1]
    y_lastPred = y_pred
    if (epoch>0 && ((epoch-1)%%update_interval==0) && delta.label < tol){
      cat("delta label ", delta.label, " < tol ", tol)
      cat("Reach tolerance threshold. Stopping training.")
      cat("Total epoch: ", epoch)
      break
    }
  }
  # return (model)
}

predict.simple_GC_DEC = function(X, adj, model){
  print("PREDICTION TENSORS")
  mylist = model(torch_tensor(X), torch_tensor(adj))
  # return value: forward: x, q
  return (mylist)
}












