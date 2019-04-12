DT_SIR <- function(x, y, k, H = 5){
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
  y <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  nclass <- as.integer(length(unique(y)))
  prior <- sapply(seq_len(nclass), function(i){mean(y == i)})
  mu <- matrix(0, nvars, nclass)
  for (i in 1:nclass){
    mu[, i] <- apply(x[y == i, ], 2, mean) - colMeans(x)
  }
  # Get the SIR matrix
  mu <- mu %*% diag(prior) %*% t(mu)
  highest1 <- order(round(diag(mu),3), decreasing = T)[1:k]
  Sigma <- cov(x)
  vec_s <- solve(Sigma[highest1, highest1]) %*% eigen(mu[highest1, highest1])$vec[,1, drop = FALSE]
  B <- matrix(0, nvars, 1)
  B[highest1,] <- vec_s
  B
}