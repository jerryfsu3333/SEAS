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

set.seed(1)
p <- 500
N <- 500

# # ###########################################
# # # Model 16
# d <- 2
# r <- 2
# Gamma <- matrix(0, p, d)
# Gamma[1:6,1] <- 1
# Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
# Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
# Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
# Beta <- matrix(runif(d*r), d, r)
# Delta <- 0.1^2*diag(1,p,p)
# 
# True_sp <- Gamma
# 
# y <- runif(N,0,4)
# f <- cbind(y,exp(y))/2
# eps <- Train(N, rep(0,p), Delta)
# x <- f %*% t(Beta) %*% t(Gamma) + eps

# # ###########################################
# # # Model 17
# d <- 2
# r <- 3
# Gamma <- matrix(0, p, d)
# Gamma[1:6,1] <- 1
# Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
# Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
# Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
# Beta <- matrix(rnorm(d*r), d, r)
# Delta <- diag(1,p,p)
# 
# True_sp <- Gamma
# 
# y <- runif(N,0,4)
# f <- cbind(y,y^2,y^3)
# eps <- Train(N, rep(0,p), Delta)
# x <- f %*% t(Beta) %*% t(Gamma) + eps

# ###########################################
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
# True_sp <- Gamma
# 
# y <- rnorm(N,0,1)
# f <- cbind(y,y^2,y^3)
# eps <- Train(N, rep(0,p), Delta)
# x <- f %*% t(Beta) %*% t(Gamma) + 0.1*eps

# # ######################################
# # # Model 17_6
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
# True_sp <- Gamma
# 
# y <- rnorm(N,0,1)
# f <- cbind(y,y^2,y^3)
# eps <- Train(N, rep(0,p), Delta)
# x <- f %*% t(Beta) %*% t(Gamma) + 0.02*eps

#######################################
# # Model 25_2
# d <- 2
# r <- 2
# p0 <- 5
# Gamma <- matrix(0, p, d)
# Gamma[1:p0,1] <- 1
# Gamma[(p0+1):(2*p0),2] <- 1
# Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
# Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
# Delta <- diag(c(rep(2,p0), rep(20,p0), rep(2,p-2*p0)))
# 
# # y <- runif(N,0,3)
# y <- rnorm(N,0,1)
# # v <- cbind(2/sqrt(20)*y, 1/sqrt(20)*exp(2*y))
# v <- cbind(2/sqrt(20)*y + 1/sqrt(20)*y^3, 3/sqrt(20)*y^2)
# eps <- mvrnorm(N, rep(0,p), Delta)
# x <- v %*% t(Gamma) + 0.2*eps
# 
# True_sp <- solve(Delta) %*% Gamma

#######################################
# # Model 22
# d <- 2
# r <- 2
# Gamma <- matrix(0, p, d)
# Gamma[1:5,1] <- c(1,1,-1,-1,0)
# Gamma[1:5,2] <- c(1,0,1,0,1)
# Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
# Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
# Delta <- diag(1,p,p)
# 
# True_sp <- Gamma
# 
# y <- rnorm(N,0,2)
# v <- cbind(y, abs(y))
# eps <- Train(N, rep(0,p), Delta)
# x <- v %*% t(Gamma) + eps

# ######################################
# Model 22_2
# d <- 2
# r <- 2
# tmp <- matrix(0, p, d)
# tmp[1:6,1] <- 1
# tmp[1:6,2] <- c(1,-1,1,-1,1,-1)
# tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
# tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
# 
# Delta <- AR(0.1,p)
# # Delta <- CS(0.5,p)
# 
# Gamma <- Delta %*% tmp
# True_sp <- tmp
# 
# y <- rnorm(N,0,1)
# v <- cbind(y + 1/5*y^2, 2*abs(y))
# eps <- Train(N, rep(0,p), Delta)
# x <- v %*% t(Gamma) + eps

# ######################################
# # Model 22_3
# d <- 2
# r <- 2
# Gamma <- matrix(0, p, d)
# Gamma[1:5,1] <- c(1,1,-1,-1,0)
# Gamma[1:5,2] <- c(1,0,1,0,1)
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

# ######################################
# # Model 22_4
# d <- 2
# r <- 2
# Gamma <- matrix(0, p, d)
# Gamma[1:5,1] <- c(1,1,-1,-1,0)
# Gamma[1:5,2] <- c(1,0,1,0,1)
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
# ######################################
# # # Model 22_6
# d <- 2
# tmp <- matrix(0, p, d)
# tmp[1:5,1] <- c(1,1,-1,-1,0)
# tmp[1:5,2] <- c(1,0,1,0,1)
# tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
# tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
# Delta <- AR(0.5,p)
# Gamma <- Delta %*% tmp
# 
# True_sp <- tmp
# 
# y <- rnorm(N,0,1)
# v <- cbind(y, abs(y))
# eps <- Train(N, rep(0,p), Delta)
# x <- v %*% t(Gamma) + eps
# ######################################
# # Model 24
# d <- 2
# r <- 2
# tmp <- matrix(0, p, d)
# tmp[1:6,1] <- 1
# tmp[1:6,2] <- c(1,-1,1,-1,1,-1)
# tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
# tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
# Delta <- CS(0.5,p)
# 
# Gamma <- Delta %*% tmp
# # Gamma <- tmp
# 
# True_sp <- tmp
# 
# y <- rnorm(N,0,1)
# v <- cbind(y+1/5*abs(y), y^2/2)
# eps <- Train(N, rep(0,p), Delta)
# x <- v %*% t(Gamma) + 1*eps

# ######################################
# # Model 24_2
# d <- 2
# r <- 2
# Gamma <- matrix(0, p, d)
# Gamma[1:6,1] <- 1
# Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
# Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
# Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
# 
# # Delta <- diag(1,p,p)
# Delta <- CS(0.5,p)
# True_sp <- solve(Delta) %*% Gamma
# 
# y <- rnorm(N,0,1)
# v <- cbind(y+1/5*exp(y), y^2/2)
# eps <- Train(N, rep(0,p), Delta)
# x <- v %*% t(Gamma) + 1*eps

# ######################################
# # Model 26
# d <- 2
# r <- 2
# Gamma <- matrix(0, p, d)
# Gamma[1:6,1] <- 1
# Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
# Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
# Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
# Delta <- diag(1,p,p)
# Delta[1:6,1:6] <- AR(0.3,6)
# 
# True_sp <- solve(Delta) %*% Gamma
# 
# # y <- rnorm(N,0,1)
# y <- runif(N,0,3)
# v <- cbind(y, y^2/2)
# eps <- Train(N, rep(0,p), Delta)
# x <- v %*% t(Gamma) + 1*eps

# # ######################################
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

# # ######################################
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

# ######################################
# # Model 27_3
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
# v <- cbind(y + 1/5*y^2, 2*abs(y))
# eps <- Train(N, rep(0,p), Delta)
# x <- v %*% t(Gamma) + eps

# ######################################
# # Model 27_4
d <- 2
tmp <- matrix(0, p, d)
tmp[1:4,1] <- 1
tmp[1:4,2] <- c(1,-1,1,-1)
tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
Delta <- AR(0.5,p)

Gamma <- Delta %*% tmp
True_sp <- tmp

y <- rnorm(N,0,1)
v <- cbind(y, y^2)
eps <- Train(N, rep(0,p), Delta)
x <- v %*% t(Gamma) + eps



cat(dcor(y,x), dcor(y,x%*%True_sp), dcor(y,x%*%True_sp[,1]), dcor(y,x%*%True_sp[,2]), sep = '&')
