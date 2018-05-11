# generate sigma, delta and mu from x, y
msda.prep <- function(x, y) {
  # data setup
  x <- as.matrix(x)
  y <- drop(y)
  nclass <- as.integer(length(unique(y)))
  prior <- rep(0, nclass)
  for (k in 1:nclass) {
    prior[k] <- mean(y == k)
  }
  nvars <- as.integer(ncol(x))
  nobs <- nrow(x)
  nres <- length(y)
  if (nres != nobs) 
    stop("x and y have different number of observations")
  # prepare sigma and delta
  mu <- matrix(0, nvars, nclass)
  sigma <- matrix(0, nvars, nvars)
  for (i in 1:nclass) {
    mu[, i] <- apply(x[y == i, ], 2, mean)
    sigma <- sigma + (sum(y == i) - 1) * cov(x[y == i, ])
  }
  sigma <- sigma/(nobs - nclass)
  delta <- sweep(mu[, -1], 1, mu[, 1], "-")
  delta <- t(delta)
  outlist <- list(sigma = sigma, delta = delta, mu = mu, prior = prior)
  outlist
}