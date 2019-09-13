rifle_func <- function(x, y, k = 10, categorical = FALSE, type = 'sir', H = 5){
  # k is always given
  
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  
  if(categorical == FALSE){
    ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
    yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    nclass <- as.integer(length(unique(yclass)))
  }else if(categorical == TRUE){
    y_unique <- unique(y)
    nclass <- H <- length(y_unique)
    yclass <- y
  }
  
  # Generate matrix A and B
  if(type == 'sir'){
    prior <- sapply(seq_len(nclass), function(i){mean(yclass == i)})
    mu <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      mu[, i] <- apply(x[yclass == i, ], 2, mean) - colMeans(x)
    }
    # Get the SIR matrix
    mu <- mu %*% diag(prior) %*% t(mu)
    sigma <- cov(x)
    A <- mu
    B <- sigma
  }
  
  a <- initial.convex(A, B, lambda=2*sqrt(log(nvars)/nobs), K=1, nu=1, trace=FALSE)
  init <- eigen(a$Pi+t(a$Pi))$vectors[,1]
  final.estimator <- rifle(A, B, init, k, 0.01, 1e-3)
  final.estimator
}