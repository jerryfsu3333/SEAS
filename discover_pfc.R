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

p <- 24
N <- 120
H <- 5

output <- lapply(1:100, function(i){
  
  # # Model 17
  # d <- 2
  # r <- 3
  # Gamma <- matrix(0, p, d)
  # Gamma[1:6,1] <- 1
  # Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Beta <- matrix(rnorm(d*r, 0, 1), d, r)
  # Delta <- diag(rep(1,p), p, p)
  # 
  # y <- runif(N,0,4)
  # f <- cbind(y,y^2,y^3)
  # eps <- mvrnorm(N, rep(0,p), Delta)
  # x <- f %*% t(Beta) %*% t(Gamma) + eps
  # Fmat <-cbind(y,y^2,y^3)
  
  # # ############################
  # # Model 17_2
  # d <- 2
  # r <- 3
  # Gamma <- matrix(0, p, d)
  # Gamma[1:6,1] <- 1
  # Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Beta <- matrix(runif(d*r), d, r)
  # Delta <- diag(1,p,p)
  # 
  # y <- rnorm(N,0,1)
  # f <- cbind(y,y^2,y^3)
  # eps <- Train(N, rep(0,p), Delta)
  # x <- f %*% t(Beta) %*% t(Gamma) + 0.1*eps
  # Fmat <-cbind(y,y^2,y^3)
  # # ######################################
  # # # Model 17_3
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:6,1] <- 1
  # Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- diag(1,p,p)
  # 
  # True_sp <- Gamma
  # 
  # y <- runif(N,0,4)
  # v <- cbind(y + y^2/4, y^3/10)
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + 1*eps
  # Fmat <-cbind(y,y^2,y^3)
  
  # ############################
  # # Model Cook&Forzani
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:5,1] <- c(1,1,-1,-1,0)
  # Gamma[1:5,2] <- c(1,0,1,0,1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # A <- matrix(rnorm(p^2), p, p)
  # Delta <- t(A) %*% A
  # 
  # y <- rnorm(N,0,2)
  # v <- cbind(y, abs(y))
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y,abs(y), y^3)
  # ############################
  # # Model 22
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:5,1] <- c(1,1,-1,-1,0)
  # Gamma[1:5,2] <- c(1,0,1,0,1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- diag(1,p,p)
  # 
  # y <- rnorm(N,0,2)
  # v <- cbind(y, abs(y))
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y,abs(y))
  ############################
  # Model 22_2
  # d <- 2
  # tmp <- matrix(0, p, d)
  # tmp[1:6,1] <- 1
  # tmp[1:6,2] <- c(1,-1,1,-1,1,-1)
  # tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
  # tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
  # Delta <- AR(0.1,p)
  # Gamma <- Delta %*% tmp
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y + 1/5*y^2, 2*abs(y))
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y, y^2, abs(y))
  # ############################
  # # Model 22_3
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:5,1] <- c(1,1,-1,-1,0)
  # Gamma[1:5,2] <- c(1,0,1,0,1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- AR(0.5,p)
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y, abs(y))
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y, abs(y))
  # ############################
  # # Model 22_4
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:5,1] <- c(1,1,-1,-1,0)
  # Gamma[1:5,2] <- c(1,0,1,0,1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- AR(0.5,p)
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y + 1/5*y^2, 2*abs(y))
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y, y^2, abs(y))
  # ############################
  # # Model 22_6
  # d <- 2
  # tmp <- matrix(0, p, d)
  # tmp[1:5,1] <- c(1,1,-1,-1,0)
  # tmp[1:5,2] <- c(1,0,1,0,1)
  # tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
  # tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
  # Delta <- AR(0.5,p)
  # 
  # Gamma <- Delta %*% tmp
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y, abs(y))
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y, abs(y))
  # ###############################
  # # Model 24
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:6,1] <- 1
  # Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- CS(0.5,p)
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y+1/5*abs(y), y^2/2)
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + 1*eps
  # Fmat <- cbind(y, abs(y), y^2)
  # ###############################
  # # Model 24_2
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:6,1] <- 1
  # Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- CS(0.5,p)
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y+1/5*exp(y), y^2/2)
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + 1*eps
  # Fmat <- cbind(y, y^2)
  # ###############################
  # # Model 25
  # d <- 2
  # p0 <- 5
  # Gamma <- matrix(0, p, d)
  # Gamma[1:p0,1] <- 1
  # Gamma[(p0+1):(2*p0),2] <- 1
  # # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # 
  # # A <- matrix(rnorm(p^2), p, p)
  # # Delta <- t(A) %*% A
  # Delta <- diag(c(rep(2,p0), rep(20,p0), rep(2,p-2*p0)))
  # Beta <- diag(c(2/sqrt(20), 1/sqrt(20)))
  # 
  # y <- runif(N,0,3)
  # f <- cbind(y, exp(y))
  # eps <- Train(N, rep(0,p), Delta)
  # x <- f %*% t(Beta) %*% t(Gamma) + eps
  # Fmat <-cbind(y, exp(y))
  # # ##############################
  # # Model 25_2
  # d <- 2
  # p0 <- 5
  # Gamma <- matrix(0, p, d)
  # Gamma[1:p0,1] <- 1
  # Gamma[(p0+1):(2*p0),2] <- 1
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- diag(c(rep(2,p0), rep(20,p0), rep(2,p-2*p0)))
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(2/sqrt(20)*y + 1/sqrt(20)*y^3, 3/sqrt(20)*y^2)
  # eps <- mvrnorm(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + 0.2*eps
  # Fmat <-cbind(y, y^2, y^3)
  # ######################################
  # # Model 26
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:6,1] <- 1
  # Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- diag(1,p,p)
  # Delta[1:6,1:6] <- AR(0.3,6)
  # 
  # y <- runif(N,0,3)
  # v <- cbind(y, y^2/2)
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + 0.5*eps
  # Fmat <-cbind(y, y^2)
  # ######################################
  # # # Model 27
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:4,1] <- 1
  # Gamma[1:4,2] <- c(1,-1,1,-1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- AR(0.5,p)
  # 
  # True_sp <- solve(Delta) %*% Gamma
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y, y^2)
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y, y^2)
  # ######################################
  # # # Model 27_2
  # d <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:4,1] <- 1
  # Gamma[1:4,2] <- c(1,-1,1,-1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- AR(0.5,p)
  # 
  # True_sp <- solve(Delta)%*%Gamma
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y, abs(y))
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y, abs(y))
  # ############################
  # # Model 27_3
  # d <- 2
  # r <- 2
  # Gamma <- matrix(0, p, d)
  # Gamma[1:4,1] <- 1
  # Gamma[1:4,2] <- c(1,-1,1,-1)
  # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Delta <- AR(0.5,p)
  # 
  # True_sp <- solve(Delta) %*% Gamma
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y + 1/5*y^2, 2*abs(y))
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y, y^2, abs(y))
  # ############################
  # # Model 27_4
  # d <- 2
  # tmp <- matrix(0, p, d)
  # tmp[1:4,1] <- 1
  # tmp[1:4,2] <- c(1,-1,1,-1)
  # tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
  # tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
  # Delta <- AR(0.5,p)
  # 
  # Gamma <- Delta %*% tmp
  # True_sp <- tmp
  # 
  # y <- rnorm(N,0,1)
  # v <- cbind(y, y^2)
  # eps <- Train(N, rep(0,p), Delta)
  # x <- v %*% t(Gamma) + eps
  # Fmat <-cbind(y, y^2)
  # ###################################
  
  sigx <- cov(x)
  Fmat_c <- scale(Fmat, scale = FALSE)
  x_c <- scale(x, scale = FALSE)
  sigfit <- (t(x_c) %*% Fmat_c %*% solve(t(Fmat_c) %*% Fmat_c) %*% t(Fmat_c) %*% x_c)/N
  sigx_invhalf <- pracma::sqrtm(solve(sigx))$B
  
  Beta_pfc1 <- sigx_invhalf %*% svd(sigx_invhalf %*% sigfit %*% sigx_invhalf)$u[,1:d]
  Beta_pfc2 <- svd(solve(sigx) %*% sigfit)$u[,1:d]
  Beta_pfc3 <- svd(solve(sigx) %*% sigfit)$v[,1:d]
  Beta_pfc4 <- svd((solve(sigx) %*% t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)/sqrt(N))$u[,1:d]
  # tmp_u <- svd(pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)$u
  # # tmp_d <-  c(0.01,0.1,1)
  # tmp_d <- c(1,5,5)
  # tmp <- tmp_u %*% diag(tmp_d) %*% t(tmp_u)
  # Beta_pfc4 <- svd((solve(sigx) %*% t(x_c) %*% Fmat_c %*% tmp)/sqrt(N))$u[,1:d]

  Beta_pfc5 <- svd(solve(sigx) %*% t(x_c) %*% Fmat_c/sqrt(N))$u[,1:d]
  Beta_pfc6 <- solve(sigx) %*% svd(sigfit)$u[,1:d]
  Beta_pfc7 <- solve(sigx) %*% svd((t(x_c) %*% Fmat_c %*% pracma::sqrtm(solve(t(Fmat_c) %*% Fmat_c))$B)/sqrt(N))$u[,1:d]
  

  True_sp <- solve(Delta) %*% Gamma
  
  dist_pfc1 <- subspace(True_sp, Beta_pfc1)
  dist_pfc2 <- subspace(True_sp, Beta_pfc2)
  dist_pfc3 <- subspace(True_sp, Beta_pfc3)
  dist_pfc4 <- subspace(True_sp, Beta_pfc4)
  dist_pfc5 <- subspace(True_sp, Beta_pfc5)
  dist_pfc6 <- subspace(True_sp, Beta_pfc6)
  dist_pfc7 <- subspace(True_sp, Beta_pfc7)
  dist_pfc12 <- subspace(Beta_pfc1, Beta_pfc2)
  dist_pfc13 <- subspace(Beta_pfc1, Beta_pfc3)
  dist_pfc14 <- subspace(Beta_pfc1, Beta_pfc4)
  dist_pfc15 <- subspace(Beta_pfc1, Beta_pfc5)
  dist_pfc16 <- subspace(Beta_pfc1, Beta_pfc6)
  dist_pfc17 <- subspace(Beta_pfc1, Beta_pfc7)
  
  # SIR
  y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
  yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  nclass <- as.integer(length(unique(yclass)))
  prior <- sapply(seq_len(nclass), function(i){mean(yclass == i)})
  mu <- matrix(0, p, nclass)
  for (i in 1:nclass){
    mu[, i] <- sqrt(prior[i]) * (apply(x[yclass == i, ], 2, mean) - colMeans(x))
  }
  Beta_sir <- svd(solve(sigx) %*% mu)$u[,1:d]
  dist_sir <- subspace(True_sp, Beta_sir)
  

  c(dist_pfc1 = dist_pfc1, dist_pfc2 = dist_pfc2, dist_pfc3 = dist_pfc3, dist_pfc4 = dist_pfc4, dist_pfc5 = dist_pfc5, dist_pfc6 = dist_pfc6, dist_pfc7 = dist_pfc7,
    dist_pfc12 =dist_pfc12, dist_pfc13 = dist_pfc13, dist_pfc14 =dist_pfc14, dist_pfc15 = dist_pfc15, dist_pfc16 = dist_pfc16, dist_pfc17 = dist_pfc17, dist_sir = dist_sir)

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

})

output <- do.call(rbind, output)

print(apply(output, 2, mean))
print(apply(output, 2, function(x){sd(x)/sqrt(100)}))
