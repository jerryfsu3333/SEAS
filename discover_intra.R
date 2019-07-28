rm(list = ls())
source("/Users/cengjing/Documents/GitHub/ssdr/models.R")
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/my_msda.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/rifle_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/lasso_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/CovSIR.R")
library(msda)
library(MASS)

set.seed(1)
p <- 20
N <- 200
H <- 5
d <- 2
r <- 5
Gamma1 <- rep(1,p)/sqrt(p)
Gamma2 <- rep(c(1,-1), p/2)/sqrt(p)
Gamma <- cbind(Gamma1, Gamma2)
# Gamma1 <- rep(c(1,0,0),p/3)/sqrt(7)
# Gamma2 <- rep(c(0,1,0), p/3)/sqrt(7)
# Gamma3 <- rep(c(0,0,1), p/3)/sqrt(7)
# Gamma <- cbind(Gamma1, Gamma2, Gamma3)

Beta <- matrix(rnorm(d*r, 0, 1), d, r)




output <- lapply(1:100, function(i){
  y <- runif(N, 0, 5)
  y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
  yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  f <- c()
  for (i in 1:5){
    tmp <- y
    tmp[yclass!=i] <- 0
    f <- cbind(f, tmp)
  }
  eps <- mvrnorm(N, rep(0,p), diag(rep(1,p), p, p))
  x <- f %*% t(Beta) %*% t(Gamma) + eps
  
  
  #####################
  # PFC
  sigx <- cov(x)
  Fmat <- cbind(y,y^2,y^3)
  Fmat_c <- scale(Fmat,scale = FALSE)
  x_c <- scale(x, scale = FALSE)
  sigfit <- (t(x_c) %*% Fmat_c %*% solve(t(Fmat_c) %*% Fmat_c) %*% t(Fmat_c) %*% x_c)/N
  sigx_half <- pracma::sqrtm(solve(sigx))$B
  Beta_pfc1 <- sigx_half %*% svd(sigx_half %*% sigfit %*% sigx_half)$u[,1:d]
  Beta_pfc2 <- svd((solve(sigx) %*% t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)/sqrt(N))$u[,1:d]
  Beta_pfc3 <- solve(sigx) %*% svd((t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)/sqrt(N))$u[,1:d]
  # Beta_pfc3 <- solve(sigx) %*% sigfit
  
  # sir_mu <- my_msda(x, y, H=H, type = 'sir', lambda = 1e-5, maxit=1e3)
  # intra_mu <- my_msda(x, y, H=H, type = 'intra', lambda = 1e-5, maxit=1e3)
  # pfc_mu <- my_msda(x, y, H=H, type = 'pfc', lambda = 1e-5, maxit=1e3, cut_y = FALSE)
  # 
  # Beta_sir <- solve(sigx) %*% sir_mu
  # Beta_intra <- solve(sigx) %*% intra_mu
  # Beta_pfc3 <- solve(sigx) %*% pfc_mu
  
  # Beta_sir <- as.matrix(sir_mu$theta[[1]])
  # Beta_intra <- as.matrix(intra_mu$theta[[1]])
  # Beta_pfc <- as.matrix(pfc_mu$theta[[1]])
  
  # ########################
  # SIR
  sigx <- cov(x)
  y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
  yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  nclass <- as.integer(length(unique(yclass)))
  prior <- sapply(seq_len(nclass), function(i){mean(yclass == i)})
  mu <- matrix(0, p, nclass)
  for (i in 1:nclass){
    mu[, i] <- sqrt(prior[i]) * (apply(x[yclass == i, ], 2, mean) - colMeans(x))
  }
  U <- mu %*% t(mu)
  sigx_half <- pracma::sqrtm(solve(sigx))$B
  Beta_sir1 <- sigx_half %*% svd(sigx_half %*% U %*% sigx_half)$u[,1:d]
  Beta_sir2 <- svd(solve(sigx) %*% mu)$u[,1:d]
  Beta_sir3 <- solve(sigx) %*% svd(mu)$u[,1:d]
  
  #########################
  # Intraslice
  
  mu <- matrix(0, p, nclass)
  for (i in 1:nclass){
    y_copy <- y
    y_copy[yclass!=i] <- 0
    mu[, i] <- cov(y_copy, x)
  }

  Beta_intra1 <- svd(solve(sigx) %*% mu)$u[,1:d]
  Beta_intra2 <- solve(sigx) %*% svd(mu)$u[,1:d]
  
  
  dist_pfc1 <- subspace(Gamma, Beta_pfc1)
  dist_pfc2 <- subspace(Gamma, Beta_pfc2)
  dist_pfc3 <- subspace(Gamma, Beta_pfc3)
  dist_sir1 <- subspace(Gamma, Beta_sir1)
  dist_sir2 <- subspace(Gamma, Beta_sir2)
  dist_sir3 <- subspace(Gamma, Beta_sir3)
  dist_intra1 <- subspace(Gamma, Beta_intra1)
  dist_intra2 <- subspace(Gamma, Beta_intra2)
  
  
  
  c(dist_pfc1 = dist_pfc1, dist_pfc2 = dist_pfc2, dist_pfc3 = dist_pfc3, dist_sir1 = dist_sir1, dist_sir2 = dist_sir2, dist_sir3 = dist_sir3, 
    dist_intra1 = dist_intra1, dist_intra2 = dist_intra2)
})

output <- do.call(rbind, output)