var.n <- function(array) {
  variance <- var(array)*(length(array)-1)/length(array)
  return(variance)
}


sd.n <- function(array) {
  return(sqrt(var.n(array)))
}


# pairwise.distance <- function(X) {
#   # run for loop in parallel
#   n.cores <- parallel::detectCores() - 1
#   my.cluster <- parallel::makeCluster(
#     n.cores,
#     type = "PSOCK"
#   )
#   doParallel::registerDoParallel(cl = my.cluster)
#   n <- nrow(X)
#   adj <- 
#     foreach(i=1:n, .combine = 'rbind') %:% 
#     foreach(j=1:n, .combine = 'c') %dopar% (
#       sqrt(sum((X[i,]-X[j,])^2))
#     )
#   parallel::stopCluster(cl = my.cluster)
#   return(adj)
# }

# make symmetric matrix?
# replace function dist
pairwise.distance <- function(X) {
  n <- nrow(X)
  adj <- matrix(,nrow = n, ncol=n)
  for (i in 1:n) {
    for (j in 1:n){
      adj[i,j] <- sqrt(sum((X[i,]-X[j,])^2))
    }
  }
  return(adj)
}

# extract_color <- function(x_pixel = NULL, y_pixel = NULL, image = NULL, beta = 49) {
#   
# }

calculate.adj.matrix <- function(x, y, x_pixel = NULL, y_pixel = NULL, image = NULL,
                                 beta = 49*scale.fac, alpha = 1, histology = TRUE){
  # x, y, x_pixel, y_pixel are lists
  if (histology == TRUE){
    stopifnot(!is.null(x_pixel) & !is.null(y_pixel) & !is.null(image))
    ######## WHY?
    stopifnot((length(x) == length(x_pixel)) & (length(y) == length(y_pixel)))
    cat("Calculating adj matrix using histology image...")
    # beta to control the range of neighborhood when calculate grey value for one spot
    # alpha to control the color scale
    beta_half <- round(beta/2,3)
    g <- NULL
    for (i in 1:length(x_pixel)){
      max_x <- dim(image)[1]
      max_y <- dim(image)[2]
      nbs <- image[max(0,round(x_pixel[i]-beta_half)):min(max_x, round(x_pixel[i]+beta_half+1)),
                   max(0,round(y_pixel[i]-beta_half)):min(max_y, round(y_pixel[i]+beta_half+1)),]
      # create average rgb for each block
      g <- rbind(g, apply(nbs,3,mean))
    }
    
    var.r <- var.n(g[,1])
    var.g <- var.n(g[,2])
    var.b <- var.n(g[,3])
    print(paste("Var of r, g, b = ", var.r, var.g, var.b))
    # raw score
    z0 <- (g[,1]*var.r+g[,2]*var.g+g[,3]*var.b)/(var.r + var.g + var.b)
    # scaled z score
    z1 <- (z0 - mean(z0))/sd.n(z0)
    z_scale <- max(sd.n(x), sd.n(y))*alpha
    z <- z1*z_scale
    print(paste("Var of x, y, z = ", var.n(x), var.n(y), var.n(z)))
    mat <- as.data.frame(cbind(x, y, z))
  }else {
    cat("Calculating adj matrix using xy only...")
    mat <- as.data.frame(cbind(x, y))
  }
  adj <- as.matrix(dist(mat, method = "euclidean", diag = TRUE, upper = TRUE))
  return(adj)
}





