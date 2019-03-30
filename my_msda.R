my_msda <- function(x, y, y_class, nlambda = 100, lambda.factor = ifelse((nobs - nclass)<=nvars, 0.2, 1e-03), lambda = NULL, 
                    dfmax = nobs, pmax = min(dfmax*2 + 20, nvars), pf = rep(1, nvars), eps = 1e-04, maxit = 1e+06, 
                    sml = 1e-06, verbose = FALSE, perturb = NULL) {
  ## data setup
  this.call <- match.call()
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  nclass <- as.integer(length(unique(y_class)))
  vnames <- colnames(x)
  prior <- rep(0, nclass)
  for (k in 1:nclass) {
    prior[k] <- mean(y_class == k)
  }
  sigma <- cov(x)
  #######################
  # mu <- matrix(0, nvars, nclass)
  # for (i in 1:nclass){
  #   mu[, i] <- apply(x[y == i, ], 2, mean) - colMeans(x)
  # }
  # # #######################
  mu <- matrix(0, nvars, nclass)
  for (i in 1:nclass){
    y_tmp <- y
    y_tmp[y_class!=i] <- 0
    mu[, i] <- cov(y_tmp, x)
  }
  ########################
  if (!is.null(perturb)) 
    diag(sigma) <- diag(sigma) + perturb
  if (is.null(vnames)) 
    vnames <- paste("V", seq(nvars), sep = "")
  nk <- as.integer(dim(mu)[2])
  ## parameter setup
  if (length(pf) != nvars) 
    stop("The size of penalty factor must be same as the number of input variables")
  maxit <- as.integer(maxit)
  verbose <- as.integer(verbose)
  sml <- as.double(sml)
  pf <- as.double(pf)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) 
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)  #ulam=0 if lambda is missing
  } else {
    # flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) 
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))  #lambda is in descending order
    nlam <- as.integer(length(lambda))
  }
  ## call Fortran core
  fit <- .Fortran("msda", obj = double(nlam), nk, nvars, as.double(sigma), as.double(t(mu)), 
                  pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1), 
                  theta = double(pmax * nk * nlam), itheta = integer(pmax), ntheta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))
  ## output
  outlist <- formatoutput(fit, maxit, pmax, nvars, vnames, nk)
  outlist <- c(outlist, list(x = x, y = y, npasses = fit$npass, jerr = fit$jerr, 
                             sigma = sigma, mu = mu, prior = prior, call = this.call))
  if (is.null(lambda)) 
    outlist$lambda <- lamfix(outlist$lambda)
  class(outlist) <- c("msda")
  outlist
}