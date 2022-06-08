source("../spaGCN.R")

calculate.p <- function(adj, l){
  adj.exp <- exp(-1*(adj^2)/(2*(l^2)))
  p <- mean(apply(adj.exp, 1, sum))-1
  return(p)
}


search.l <- function(p, adj, start=0.01, end=1000, tol=0.01, max.run=100){
  run <- 0
  p.low <- calculate.p(adj, start)
  p.high <- calculate.p(adj, end)
  if (p.low > p+tol) {
    cat("l not found, try smaller start point.")
    return(NULL)
  }
  else if (p.high < p+tol) {
    cat("l not found, try bigger end point.")
    return(NULL)
  }
  else if (abs(p.low-p) <= tol) {
    cat(paste("recommended l = ", start))
    return(start)
  }
  else if (abs(p.high-p) <= tol) {
    cat(paste("recommended l = ", end))
    return(end)
  }
  while ((p.low+tol)<p & p<(p.high-tol)) {
    run <- run + 1
    cat(paste0("Run ",run,": l [",start,", ",end,"], p [",p.low,", ",p.high,"]\n"))
    if (run > max.run){
      cat(paste0("Exact l not found, closest values are:",'\n', "l = ",start, ": p = ",p.low, '\n', "l = ",end,": p = ",p.high))
      return(NULL)
    }
    mid <- (start+end)/2
    p.mid <- calculate.p(adj,mid)
    if (abs(p.mid-p) <= tol){
      cat(paste0("recommended l = ", mid))
      return(mid)
    }
    if (p.mid <= p){
      start <- mid
      p.low <- p.mid
    }
    else {
      end <- mid
      p.high <- p.mid
    }
  }
}



search.res <- function(spaData, adj, l, target.num, 
                       start = 0.4, step = 0.1, tol = 5e-3, lr = 0.05, 
                       max.epochs = 10, seed = 100, max.run = 10){
  set.seed(seed)
  res <- start
  cat(paste0("Start at res = ", res, " step = ", step))
  spaGCN.train(spaData, adj, init.spa = TRUE, init = "louvain", 
               res = res, tol = tol, lr = lr, max.epochs = max.epochs, l = l)
  
}



















