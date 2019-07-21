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
p <- 20
N <- 200
H <- 5
d <- 2
r <- 3
Gamma1 <- rep(1,p)/sqrt(p)
Gamma2 <- rep(c(1,-1), p/2)/sqrt(p)
Gamma <- cbind(Gamma1, Gamma2)
# Gamma1 <- rep(c(1,0,0),p/3)/sqrt(7)
# Gamma2 <- rep(c(0,1,0), p/3)/sqrt(7)
# Gamma3 <- rep(c(0,0,1), p/3)/sqrt(7)
# Gamma <- cbind(Gamma1, Gamma2, Gamma3)

Beta <- matrix(rnorm(d*r, 0, 1), d, r)




output <- lapply(1:100, function(i){
  y <- runif(N, 0, 4)
  # f <- cbind(y, y^2, y^3*sqrt(0.01))
  f <- cbind(y, y^2, exp(y))
  # f <- cbind(y, y^2)
  eps <- mvrnorm(N, rep(0,p), diag(rep(1,p), p, p))
  x <- f %*% t(Beta) %*% t(Gamma) + eps
  
  
  sigx <- cov(x)
  Fmat <- cbind(y,y^2,y^3)
  Fmat_c <- scale(f,scale = FALSE)
  x_c <- scale(x, scale = FALSE)
  sigfit <- (t(x_c) %*% Fmat_c %*% solve(t(Fmat_c) %*% Fmat_c) %*% t(Fmat_c) %*% x_c)/N
  sigx_half <- pracma::sqrtm(solve(sigx))$B
  Beta_pfc2 <- sigx_half %*% svd(sigx_half %*% sigfit %*% sigx_half)$u[,1:d]
  
  Beta_pfc3 <- solve(sigx) %*% sigfit
  Beta_pfc4 <- solve(sigx) %*% t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B/sqrt(N)
  
  sir_mu <- my_msda(x, y, H=H, type = 'sir', lambda = 1e-5, maxit=1e3)
  intra_mu <- my_msda(x, y, H=H, type = 'intra', lambda = 1e-5, maxit=1e3)
  pfc_mu <- my_msda(x, y, H=H, type = 'pfc', lambda = 1e-5, maxit=1e3, cut_y = FALSE)

  # Beta_sir <- solve(sigx) %*% sir_mu
  # Beta_intra <- solve(sigx) %*% intra_mu
  # Beta_pfc <- solve(sigx) %*% pfc_mu
  
  Beta_sir <- as.matrix(sir_mu$theta[[1]])
  Beta_intra <- as.matrix(intra_mu$theta[[1]])
  Beta_pfc <- as.matrix(pfc_mu$theta[[1]])
  
  
  dist_sir <- subspace_2(Gamma, svd(Beta_sir)$u[,1:d, drop = FALSE])
  dist_intra <- subspace_2(Gamma, svd(Beta_intra)$u[,1:d, drop = FALSE])
  dist_pfc <- subspace_2(Gamma, svd(Beta_pfc)$u[,1:d, drop = FALSE])
  dist_pfc2 <- subspace_2(Gamma, Beta_pfc2)
  dist_pfc3 <- subspace_2(Gamma, svd(Beta_pfc3)$u[,1:d, drop = FALSE])
  dist_pfc4 <- subspace_2(Gamma, svd(Beta_pfc4)$u[,1:d, drop = FALSE])
  
  c(sir = dist_sir, intra = dist_intra, pfc = dist_pfc, pfc2 = dist_pfc2, pfc3 = dist_pfc3, pfc4 = dist_pfc4)
})

output <- do.call(rbind, output)
