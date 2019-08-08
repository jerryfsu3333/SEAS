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
N <- 500
H <- 5
d <- 2
r <- H-1

Gamma1 <- rep(1,p)/sqrt(p)
Gamma2 <- rep(c(1,-1), p/2)/sqrt(p)
Gamma <- cbind(Gamma1, Gamma2)

Beta <- matrix(rnorm(d*r, 0, 1), d, r)

output <- lapply(1:100, function(i){
  y <- runif(N, 0, 4)
  
  y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
  yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  f <- c()
  for (i in 1:(H-1)){
    tmp <- as.double(yclass == i)
    f <- cbind(f, tmp)
  }
  
  Delta <- diag(rep(1,p), p, p)
  eps <- mvrnorm(N, rep(0,p), Delta)
  x <- f %*% t(Beta) %*% t(Gamma) + 0.2 * eps
  
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
  sigx_invhalf <- pracma::sqrtm(solve(sigx))$B
  Beta_sir1 <- sigx_invhalf %*% svd(sigx_invhalf %*% U %*% sigx_invhalf)$u[,1:d]
  Beta_sir2 <- svd(solve(sigx) %*% U)$u[,1:d]
  Beta_sir3 <- svd(solve(sigx) %*% U)$v[,1:d]
  Beta_sir4 <- svd(solve(sigx) %*% mu)$u[,1:d]
  Beta_sir5 <- solve(sigx) %*% svd(U)$u[,1:d]
  Beta_sir6 <- solve(sigx) %*% svd(mu)$u[,1:d]

  dist_sir1 <- subspace(Gamma, Beta_sir1)
  dist_sir2 <- subspace(Gamma, Beta_sir2)
  dist_sir3 <- subspace(Gamma, Beta_sir3)
  dist_sir4 <- subspace(Gamma, Beta_sir4)
  dist_sir5 <- subspace(Gamma, Beta_sir5)
  dist_sir6 <- subspace(Gamma, Beta_sir6)
  dist_sir12 <- subspace(Beta_sir1, Beta_sir2)
  dist_sir13 <- subspace(Beta_sir1, Beta_sir3)
  dist_sir14 <- subspace(Beta_sir1, Beta_sir4)
  dist_sir15 <- subspace(Beta_sir1, Beta_sir5)
  dist_sir16 <- subspace(Beta_sir1, Beta_sir6)
  dist_sir23 <- subspace(Beta_sir2, Beta_sir3)


  c(dist_sir1 = dist_sir1, dist_sir2 = dist_sir2, dist_sir3 = dist_sir3, dist_sir4 = dist_sir4, dist_sir5 = dist_sir5, dist_sir6 = dist_sir6,
    dist_sir12 =dist_sir12, dist_sir13 = dist_sir13, dist_sir14 =dist_sir14, dist_sir15 = dist_sir15, dist_sir16 = dist_sir16, dist_sir23 = dist_sir23)
  
  
})

output <- do.call(rbind, output)

print(apply(output, 2, mean))
print(apply(output, 2, function(x){sd(x)/sqrt(100)}))
