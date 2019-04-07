rifle_func <- function(x, y, type = 'sir'){
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  sigma <- cov(x)
  # Generate matrix A and B
  if(type == 'sir'){
    y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
    y <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
    nclass <- as.integer(length(unique(y)))
    prior <- sapply(seq_len(nclass), function(i){mean(y == i)})
    mu <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      mu[, i] <- apply(x[y == i, ], 2, mean) - colMeans(x)
    }
    mu <- mu %*% diag(prior) %*% t(mu)
    A <- mu
    B <- sigma
  }
  
  a <- initial.convex(A, B, lambda=2*sqrt(log(nvars)/nobs), K=1, nu=1, trace=FALSE)
  init <- eigen(a$Pi+t(a$Pi))$vectors[,1]
  k <- 10
  final.estimator <- rifle(A, B, init, k, 0.01, 1e-3)
}