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
r <- 3
# Gamma1 <- rep(1,p)/sqrt(p)
# Gamma <- matrix(Gamma1, ncol = 1)

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
  f <- cbind(y, y^2, y^3)
  # f <- cbind(y, y^2, exp(y))
  # f <- cbind(y, y^2)

  # y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
  # yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  # f <- c()
  # for (i in 1:H){
  #   tmp <- as.double(yclass == i)
  #   f <- cbind(f, tmp)
  # }

  # v <- matrix(exp(y), ncol = 1)

  Delta <- diag(rep(1,p), p, p)
  eps <- mvrnorm(N, rep(0,p), Delta)
  x <- f %*% t(Beta) %*% t(Gamma) + 1 * eps
  
  # PFC
  sigx <- cov(x)
  Fmat <- t(sapply(y, function(x){c(x, x^2, x^3)}))
  Fmat_c <- scale(Fmat, scale = FALSE)
  x_c <- scale(x, scale = FALSE)
  sigfit <- (t(x_c) %*% Fmat_c %*% solve(t(Fmat_c) %*% Fmat_c) %*% t(Fmat_c) %*% x_c)/N
  sigx_invhalf <- pracma::sqrtm(solve(sigx))$B
  Beta_pfc1 <- sigx_invhalf %*% svd(sigx_invhalf %*% sigfit %*% sigx_invhalf)$u[,1:d]
  ## Angle
  # angle <- acos( sum(Gamma * Beta_pfc1) / ( norm(Gamma, '2') * norm(Beta_pfc1, '2') ) )
  ##
  
  Beta_pfc2 <- svd(solve(sigx) %*% sigfit)$u[,1:d]
  Beta_pfc3 <- svd(solve(sigx) %*% sigfit)$v[,1:d]
  # tmp <- matrix(rnorm(p*d), p, d)
  # Beta_pfc3 <- qr.Q(qr(tmp))
  Beta_pfc4 <- svd((solve(sigx) %*% t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)/sqrt(N))$u[,1:d]
  # tmp_u <- svd(pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)$u
  # # tmp_d <-  c(0.01,0.1,1)
  # tmp_d <- c(1,5,5)
  # tmp <- tmp_u %*% diag(tmp_d) %*% t(tmp_u)
  # Beta_pfc4 <- svd((solve(sigx) %*% t(x_c) %*% Fmat_c %*% tmp)/sqrt(N))$u[,1:d]

  Beta_pfc5 <- svd(solve(sigx) %*% t(x_c) %*% Fmat_c/sqrt(N))$u[,1:d]
  Beta_pfc6 <- solve(sigx) %*% svd(sigfit)$u[,1:d]
  Beta_pfc7 <- solve(sigx) %*% svd((t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)/sqrt(N))$u[,1:d]

  dist_pfc1 <- subspace(Gamma, Beta_pfc1)
  dist_pfc2 <- subspace(Gamma, Beta_pfc2)
  dist_pfc3 <- subspace(Gamma, Beta_pfc3)
  dist_pfc4 <- subspace(Gamma, Beta_pfc4)
  dist_pfc5 <- subspace(Gamma, Beta_pfc5)
  dist_pfc6 <- subspace(Gamma, Beta_pfc6)
  dist_pfc7 <- subspace(Gamma, Beta_pfc7)
  dist_pfc12 <- subspace(Beta_pfc1, Beta_pfc2)
  dist_pfc13 <- subspace(Beta_pfc1, Beta_pfc3)
  dist_pfc14 <- subspace(Beta_pfc1, Beta_pfc4)
  dist_pfc15 <- subspace(Beta_pfc1, Beta_pfc5)
  dist_pfc16 <- subspace(Beta_pfc1, Beta_pfc6)
  dist_pfc17 <- subspace(Beta_pfc1, Beta_pfc7)

  c(dist_pfc1 = dist_pfc1, dist_pfc2 = dist_pfc2, dist_pfc3 = dist_pfc3, dist_pfc4 = dist_pfc4, dist_pfc5 = dist_pfc5, dist_pfc6 = dist_pfc6, dist_pfc7 = dist_pfc7,
    dist_pfc12 =dist_pfc12, dist_pfc13 = dist_pfc13, dist_pfc14 =dist_pfc14, dist_pfc15 = dist_pfc15, dist_pfc16 = dist_pfc16, dist_pfc17 = dist_pfc17)

  # ######################
  # # PFC scale-invariant
  # sigx <- cov(x)
  # Fmat <- cbind(100*y,y^2,y^3)
  # Fmat_c <- scale(Fmat,scale = FALSE)
  # x_c <- scale(x, scale = FALSE)
  # sigfit <- (t(x_c) %*% Fmat_c %*% solve(t(Fmat_c) %*% Fmat_c) %*% t(Fmat_c) %*% x_c)/N
  # sigx_invhalf <- pracma::sqrtm(solve(sigx))$B
  # Beta1_pfc1 <- sigx_invhalf %*% svd(sigx_invhalf %*% sigfit %*% sigx_invhalf)$u[,1:d]
  # Beta1_pfc2 <- svd(solve(sigx) %*% sigfit)$u[,1:d]
  # Beta1_pfc3 <- svd(solve(sigx) %*% sigfit)$v[,1:d]
  # Beta1_pfc4 <- svd((solve(sigx) %*% t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)/sqrt(N))$u[,1:d]
  # Beta1_pfc5 <- svd(solve(sigx) %*% t(x_c) %*% Fmat_c)$u[,1:d]
  # Beta1_pfc6 <- solve(sigx) %*% svd(sigfit)$u[,1:d]
  # Beta1_pfc7 <- solve(sigx) %*% svd((t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)/sqrt(N))$u[,1:d]
  # 
  # # mu <- t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B/sqrt(N)
  # # U <- mu%*% t(mu)
  # # Beta_pfczz <- svd(solve(sigx) %*% U)$u[,1:d]
  # # norm(sigfit-U, type = 'F')
  # 
  # dist_pfc1 <- subspace(Beta1_pfc1, Beta_pfc1)
  # dist_pfc2 <- subspace(Beta1_pfc2, Beta_pfc2)
  # dist_pfc3 <- subspace(Beta1_pfc3, Beta_pfc3)
  # dist_pfc4 <- subspace(Beta1_pfc4, Beta_pfc4)
  # dist_pfc5 <- subspace(Beta1_pfc5, Beta_pfc5)
  # dist_pfc6 <- subspace(Beta1_pfc6, Beta_pfc6)
  # dist_pfc7 <- subspace(Beta1_pfc7, Beta_pfc7)

# # ########################
# # SIR
# sigx <- cov(x)
# y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
# yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
# nclass <- as.integer(length(unique(yclass)))
# prior <- sapply(seq_len(nclass), function(i){mean(yclass == i)})
# mu <- matrix(0, p, nclass)
# for (i in 1:nclass){
#   mu[, i] <- sqrt(prior[i]) * (apply(x[yclass == i, ], 2, mean) - colMeans(x))
# }
# U <- mu %*% t(mu)
# sigx_invhalf <- pracma::sqrtm(solve(sigx))$B
# Beta_sir1 <- sigx_invhalf %*% svd(sigx_invhalf %*% U %*% sigx_invhalf)$u[,1:d]
# Beta_sir2 <- svd(solve(sigx) %*% U)$u[,1:d]
# Beta_sir3 <- svd(solve(sigx) %*% U)$v[,1:d]
# # tmp <- matrix(rnorm(p*d), p, d)
# # Beta_sir3 <- qr.Q(qr(tmp))
# Beta_sir4 <- svd(solve(sigx) %*% mu)$u[,1:d]
# Beta_sir5 <- solve(sigx) %*% svd(U)$u[,1:d]
# Beta_sir6 <- solve(sigx) %*% svd(mu)$u[,1:d]
# 
# # subspace_test <- function(A,B){
# #   Pa <- qr.Q(qr(A))
# #   Pb <- qr.Q(qr(B))
# #   tmp <- svd(t(Pa) %*% Pb)$d
# #   d <- dim(A)[2]
# #   dist <- sqrt(sum(1-tmp^2)) / sqrt(d)
# #   return(dist)
# # }
# 
# dist_sir1 <- subspace(Gamma, Beta_sir1)
# dist_sir2 <- subspace(Gamma, Beta_sir2)
# dist_sir3 <- subspace(Gamma, Beta_sir3)
# dist_sir4 <- subspace(Gamma, Beta_sir4)
# dist_sir5 <- subspace(Gamma, Beta_sir5)
# dist_sir6 <- subspace(Gamma, Beta_sir6)
# dist_sir12 <- subspace(Beta_sir1, Beta_sir2)
# dist_sir13 <- subspace(Beta_sir1, Beta_sir3)
# dist_sir14 <- subspace(Beta_sir1, Beta_sir4)
# dist_sir15 <- subspace(Beta_sir1, Beta_sir5)
# dist_sir16 <- subspace(Beta_sir1, Beta_sir6)
# dist_sir23 <- subspace(Beta_sir2, Beta_sir3)
# 
# 
# c(dist_sir1 = dist_sir1, dist_sir2 = dist_sir2, dist_sir3 = dist_sir3, dist_sir4 = dist_sir4, dist_sir5 = dist_sir5, dist_sir6 = dist_sir6,
#   dist_sir12 =dist_sir12, dist_sir13 = dist_sir13, dist_sir14 =dist_sir14, dist_sir15 = dist_sir15, dist_sir16 = dist_sir16, dist_sir23 = dist_sir23)


})

output <- do.call(rbind, output)

print(apply(output, 2, mean))
print(apply(output, 2, function(x){sd(x)/sqrt(100)}))
