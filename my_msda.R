my_msda <- function(x, y, H = 5, nlambda = 100, type = 'sir', lambda.factor = ifelse((nobs - nclass)<=nvars, 0.2, 1e-03), lambda = NULL, 
                    dfmax = nobs, pmax = min(dfmax*2 + 20, nvars), pf = rep(1, nvars), eps = 1e-04, maxit = 1e+06, 
                    sml = 1e-06, verbose = FALSE, perturb = NULL, cut_y = FALSE) {
  this.call <- match.call()
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  vnames <- colnames(x)
  sigma <- cov(x)
  ######################################
  # different mu
  y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
  yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  nclass <- as.integer(length(unique(yclass)))
  
  if(type == 'sir'){
    prior <- sapply(seq_len(nclass), function(i){mean(yclass == i)})
    mu <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      mu[, i] <- apply(x[yclass == i, ], 2, mean) - colMeans(x)
    }
  }else if(type == 'save'){
    
  }else if(type == 'pfc'){
    cut_func <- function(x, lb, ub){
      if(x < lb){
        return(lb)
      } else if(x > ub){
        return(ub)
      } else{
        return(x)
      }
    }
    if(cut_y){
      lb <- quantile(y, 0.2)[[1]]
      ub <- quantile(y, 0.8)[[1]]
      y <- sapply(y, cut_func, lb = lb, ub = ub)
    }
    # Fmat <- cbind(y, (y^2), (y^3))
    Fmat <- cbind(y,y^2, y^3)
    Fmat_c <- scale(Fmat,scale = FALSE)
    # Fmat_c <- scale(Fmat)
    x_c <- scale(x, scale = FALSE)
    invhalf_func <- function(Sigma){
      S <- eigen(Sigma)
      pos <- 1/sqrt(S$val[S$val > 0.001])
      zer <- rep(0,length(S$val < 0.001))
      mid <- diag(c(pos,zer), ncol(Sigma), ncol(Sigma))
      S$vec%*%mid%*%t(S$vec)
    }
    # tmp <- svd(t(Fmat_c) %*% Fmat_c)
    # invhalf <- tmp$u %*% diag(1/sqrt(tmp$d)) %*% t(tmp$u)
    # invhalf1 <- invhalf_func(t(Fmat_c) %*% Fmat_c)
    # invhalf <-  pracma::sqrtm(t(Fmat_c) %*% Fmat_c)$B
    
    # A <- matrix(rnorm(3*3), 3, 3)
    # A <- qr.Q(qr(A))
    # invhalf1 <- A %*% diag(c(0.01,1,10),3,3) %*% t(A)
    
    # A <- svd(t(Fmat_c) %*% Fmat_c)$u
    # invhalf <- A %*% diag(c(1,1,1),3,3) %*% t(A)
    # mu <- (t(x_c) %*% Fmat_c %*% invhalf)
    mu <- (t(x_c) %*% Fmat_c)/sqrt(nobs)
  }else if(type == 'intra'){
    mu <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      y_copy <- y
      y_copy[yclass!=i] <- 0
      mu[, i] <- cov(y_copy, x)
    }
  }
  return(mu)
  ######################################
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
    ulam <- as.double(rev(sort(lambda)))  #lambda is declining
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
                             sigma = sigma, mu = mu, prior = ifelse(type %in% c('sir','save'), prior, NA), call = this.call))
  if (is.null(lambda)) 
    outlist$lambda <- lamfix(outlist$lambda)
  class(outlist) <- c("msda")
  outlist
}