source("../spaGCN.R")
library(torch)

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
                       max.epochs = 10, seed = seed, max.run = 10){
  set.seed(seed)
  torch_manual_seed(seed)
  res <- start
  cat(paste0("Start at res = ", res, " step = ", step))
  clf <- spaGCN(spaData, adj, init.spa = TRUE, init = "louvain",
                res = res, tol = tol, lr = lr, max.epochs = max.epochs)
  clf
  clf <- set_l(clf, l)
  clf <- train.spaGCN(clf, seed)
  pred <- predict.spaGCN(clf)
  y_pred <- pred$y_pred
  old.num <- length(unique(y_pred))
  cat("Res = ", res, "Num of clusters = ", old.num)
  run <- 0
  while (old.num != target.num) {
    set.seed(seed)
    torch_manual_seed(seed)
    old.sign <- ifelse(old.num < target.num, 1, -1)
    clf <- spaGCN(spaData, adj, init.spa = TRUE, init = "louvain",
                  res = res+step*old.sign, tol = tol, lr = lr, max.epochs = max.epochs)
    clf <- set_l(clf, l)
    clf <- train.spaGCN(clf, seed)
    pred <- predict.spaGCN(clf)
    y_pred <- pred$y_pred
    new.num <- length(unique(y_pred))
    cat("Res = ", res+step*old.sign, "Num of clusters = ", new.num)
    if (new.num == target.num){
      res = res+step*old.sign
      cat("Recommended res = ", res)
      return(res)
    } 
    new.sign = ifelse(new.num < target.num, 1, -1)
    if (new.sign == old.sign){
      res = res+step*old.sign
      cat("Res changed to ", res)
      old.num = new.num
    } else {
      step = step/2
      cat("step changed to ", step)
    }
    if (run > max.run){
      cat("Exact resolution not found")
      cat("Recommended res = ", res)
      return(res)
    }
    run <- run + 1
  }
  cat("Recommended res = ", res)
  return(res)
}


refine <- function(sample_id, pred, dis, shape = "hexagon"){
  refined.pred = NULL
  pred = as.data.frame(pred, row.names = sample_id)
  dis_df = as.data.frame(dis, row.names = sample_id)
  colnames(dis_df) = sample_id
  if (shape == "hexagon"){
    num_nbs = 6
  } else if (shape == "square"){
    num_obs = 4
  } else {
    cat("Shape not recognized, shape = 'hexagon' for Visium data, 'square' for ST data")
  }
  for (i in 1:length(sample_id)){
    index = sample_id[i]
    dis_tmp = as.data.frame(cbind(dis_df[index],sample_id))
    dis_tmp = dis_tmp[order(dis_tmp[,1]),]
    nbs = dis_tmp[1:num_nbs+1,]
    nbs_pred = pred[nbs$sample_id,]
    self_pred = pred[index,]
    v_c = as.data.frame(table(nbs_pred))
    if ((v_c[v_c$nbs_pred == self_pred, "Freq"] < num_nbs/2) && (max(v_c$Freq) > num_nbs/2)){
      refined.pred = c(refined.pred, v_c[v_c$Freq == max(v_c$Freq),"nbs_pred"])
    } else {
      refined.pred = c(refined.pred, self_pred)
    }
  }
  return(refined.pred)
}

















