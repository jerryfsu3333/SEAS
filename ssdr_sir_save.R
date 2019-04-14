library(msda)
library(MASS)
library(methods)
library(glmnet)
library(rifle)
library(LassoSIR)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/my_msda.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/rifle_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/DT_SIR.R")
source("/Users/cengjing/Documents/GitHub/ssdr/lasso_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/CovSIR.R")

###################################################
# Functions
###################################################

# Simulating train data function
Train <- function(n, mu, Sigma){
  m <- mvrnorm(n, mu, Sigma)
  return(m)
}

# AR function
AR <- function(rho, p){
  m <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      m[i,j] <- rho**(abs(i-j))
    }
  }
  return(m)
}

# CS function
CS <- function(rho, p){
  m <- matrix(rho,p,p)
  diag(m) <- 1
  return(m)
}

# Block CS function

CS_blk <- function(rho,p,s){
  block <- CS(rho,s)
  A <- diag(rep(1,p/s),p/s,p/s)
  m <- kronecker(A,block)
  return(m)
}

CS_blk2 <- function(rho,p,s){
  block <- CS(rho,s)
  m <- diag(rep(1,p),p,p)
  m[1:s,1:s] <- block
  return(m)
}

######################## evaluation ########################

# eval_val <- function(Beta, x, y, slices){
#   y_sliced <- cut(y, breaks = slices, include.lowest = TRUE, labels = FALSE)
#   tmp <- prep(x, y_sliced)
#   sigma <- tmp$sigma
#   mu <- tmp$mu
#   l <- length(Beta)
#   result <- rep(0, l)
#   for (i in 1:l){
#     mat <- as.matrix(Beta[[i]])
#     result[i] <- 0.5 * sum(diag(t(mat) %*% sigma %*% mat - 2 * t(mu) %*% mat), na.rm = TRUE)
#   }
#   return(result)
# }

# eval_val_rmse <- function(Beta, x, y){
#   l <- length(Beta)
#   # seq_len(l) is 1:l, but faster
#   result <- sapply(seq_len(l), function(i){
#     if(is.null(Beta[[i]])){
#       NA
#     }else{
#       mat <- as.matrix(Beta[[i]])
#       xnew <- x %*% mat
#       fit <- lm(y~xnew)
#       rmse <- sqrt(mean((fit$residuals)^2))
#       rmse 
#     }
#   })
#   result
# }

# eval_val_rmse_rank <- function(Beta, x, y, rank, slices = NULL){
#   l <- length(Beta)
#   result <- rep(0, l)
#   for (i in 1:l){
#     r <- rank[i]
#     mat <- as.matrix(Beta[[i]])
#     if(r != 0){
#       mat <- svd(as.matrix(Beta[[i]]))$u[,1:r]
#     }
#     xnew <- x %*% mat
#     fit <- lm(y~xnew)
#     rmse <- sqrt(mean((fit$residuals)^2))
#     result[i] <- rmse
#   }
#   return(result)
# }


# eval_val_rmse_2 <- function(Beta, x, y, slices){
#   y_sliced <- cut(y, breaks = slices, include.lowest = TRUE, labels = FALSE)
#   l <- length(Beta)
#   result <- rep(0, l)
#   for (i in 1:l){
#     mat <- as.matrix(Beta[[i]])
#     xnew <- x %*% mat
#     fit <- lm(y_sliced~xnew)
#     rmse <- sqrt(mean((fit$residuals)^2))
#     result[i] <- rmse
#   }
#   return(result)
# }

# eval_val_cart <- function(Beta, xtrain, ytrain, xval, yval, slices_val){
#   y_val_sliced <- cut(yval, breaks = slices_val, include.lowest = TRUE, labels = FALSE)
#   l <- length(Beta)
#   errors <- rep(0, l)
#   for (i in 1:l){
#     mat <- as.matrix(Beta[[i]])
#     xnew_tr <- xtrain %*% mat
#     xnew_val <- as.data.frame(xval %*% mat)
#     fit <- rpart(ytrain~xnew_tr, method = "class", control=rpart.control(minsplit=20, cp=0.001))
#     pred <- predict(fit, xnew_val, type = "class")
#     errors[i] <- mean(y_val_sliced != pred)
#   }
#   return(errors)
# }

##########################################################################################
#                                    Data structure                                      #
##########################################################################################

# #############  Model 1 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 200  # Sample size
# N_val <- 200  # Sample size of validation dataset
# H <- 10
# 
# # Construct true Beta
# Mu <- rep(0,p)
# # Sigma <- CS_blk(0.5,800,4)
# Sigma <- AR(0.5, p)
# Beta <- matrix(0, p, 2)
# Beta[1:10,1] <- 1
# Beta[1:10,2] <- c(1,-1,1,-1,1,-1,1,-1,1,-1)
# # Beta[11:20,2] <- -1
# # Beta[11:20,2] <- 1
# nz_vec <- 1:10
# # nz_vec <- 1:20
# r <- 2
# 
# model <- function(x, Beta){
#   nobs <- dim(x)[1]
#   y <- sign(x %*% Beta[,1]) * log(abs(x %*% Beta[,2] + 5)) + 0.1 * rnorm(nobs)
#   return(y)
# }

# #############  Model 2 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500  # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 1)
# Beta[1:20,1] <- 1
# Beta <- sqrt(0.5)*Beta/norm(Beta, '2')
# 
# nz_vec <- 1:20
# s <- 20
# r <- 1
# 
# True_sp <- Beta
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta + 0.5 * rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 2.1 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500  # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 1)
# Beta[1:20,1] <- 1
# Beta <- sqrt(0.5)*Beta/norm(Beta, '2')
# 
# nz_vec <- 1:20
# s <- 20
# r <- 1
# 
# True_sp <- Beta
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- (x %*% Beta)^3/2 + rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 2.2 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500  # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 1)
# Beta[1:20,1] <- 1
# Beta <- Beta/norm(Beta, '2')
# 
# nz_vec <- 1:20
# s <- 20
# r <- 1
# 
# True_sp <- Beta
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- exp((x%*%Beta)/10) + rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 3 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:6,1] <- 1
# Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
# nz_vec <- 1:6
# s <- 6
# r <- 2
# 
# True_sp <- Beta
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- (x %*% Beta[,1])/(0.5+(x %*% Beta[,2] + 1.5)^2) + 0.2 * rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)

#############  Model 4 #############
set.seed(1)

p <- 100  # Dimension of observations
N <- 500 # Sample size
N_val <- 500  # Sample size of validation dataset
H <- 5

Mu <- rep(0,p)
Sigma <- AR(0.5, p)

# Construct true Beta
Beta <- matrix(0, p, 2)
Beta[1:6,1] <- 1
Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')

nz_vec <- 1:6
s <- 6
r <- 2
True_sp <- Beta

Data <- function(N){
  x <- Train(N, Mu, Sigma)
  nobs <- dim(x)[1]
  y <- abs((x %*% Beta[,1]) / 4 + 2)^3 * sign(x %*% Beta[,2]) + 0.2 * rnorm(nobs)
  list(x = x, y = y)
}

sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)

# #############  Model 4.1 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:4,1] <- 1
# Beta[5:7,2] <- 1
# Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
# Beta[,2] <- sqrt(0.5)*Beta[,2]/norm(Beta[,2], '2')
# nz_vec <- 1:7
# s <- 7
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- abs((x %*% Beta[,1]) / 4 + 2)^3 * sign(x %*% Beta[,2]) + rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)

# #############  Model 5 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:6,1] <- 1
# Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
# Beta[,1] <- Beta[,1]/norm(Beta[,1], '2')
# Beta[,2] <- Beta[,2]/norm(Beta[,2], '2')
# 
# nz_vec <- 1:6
# s <- 6
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2]) + 0.2 * rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 5.1 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:7,1] <- 1
# Beta[8:12,2] <- 1
# Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
# Beta[,2] <- sqrt(0.5)*Beta[,2]/norm(Beta[,2], '2')
# 
# nz_vec <- 1:12
# s <- 12
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2]) + 0.2 * rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 5.2 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:7,1] <- 1
# Beta[8:12,2] <- 1
# Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
# Beta[,2] <- sqrt(0.5)*Beta[,2]/norm(Beta[,2], '2')
# 
# nz_vec <- 1:12
# s <- 12
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2]) + rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 6 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:6,1] <- 1
# Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
# # Beta[,1] <- Beta[,1]/norm(Beta[,1], '2')
# # Beta[,2] <- Beta[,2]/norm(Beta[,2], '2')
# 
# nz_vec <- 1:6
# s <- 6
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.2 *  rnorm(nobs) )
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 6.1 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:7,1] <- 1
# Beta[8:12,2] <- 1
# Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
# Beta[,2] <- sqrt(0.5)*Beta[,2]/norm(Beta[,2], '2')
# 
# nz_vec <- 1:12
# s <- 12
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.2 *  rnorm(nobs) )
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 6.2 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:7,1] <- 1
# Beta[8:12,2] <- 1
# Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
# Beta[,2] <- sqrt(0.5)*Beta[,2]/norm(Beta[,2], '2')
# 
# nz_vec <- 1:12
# s <- 12
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + rnorm(nobs) )
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 7 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:6,1] <- 1
# Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
# nz_vec <- 1:6
# s <- 6
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * (2 + (x %*% Beta[,2]) / 3)^2 + 0.2 * rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 7.1 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:8,1] <- 1
# Beta[9:12,2] <- 1
# 
# nz_vec <- 1:12
# s <- 12
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * (2 + (x %*% Beta[,2]) / 3)^2 + rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 7.2 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# 
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:8,1] <- 1
# Beta[9:12,2] <- 1
# Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
# Beta[,2] <- sqrt(0.5)*Beta[,2]/norm(Beta[,2], '2')
# 
# nz_vec <- 1:12
# s <- 12
# r <- 2
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * (2 + (x %*% Beta[,2]) / 3)^2 + 0.2*rnorm(nobs)
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

# #############  Model 8 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# # Construct true Beta
# Beta <- matrix(0, p, 1)
# Beta[1:6,1] <- 1
# nz_vec <- 1:6
# s <- 6
# r <- 1
# True_sp <- Beta
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- 1 * (x %*% Beta[,1])^2 + 1 * rnorm(nobs)
#   list(x = x, y = y)
# }

# #############  Model 9 #############
# set.seed(1)
# 
# p <- 100  # Dimension of observations
# N <- 500 # Sample size
# N_val <- 500  # Sample size of validation dataset
# d <- 2
# r <- 4
# 
# Delta <- AR(0.5, p)
# Gamma <- matrix(0, p, 2)
# Gamma[1:6,1] <- 1
# Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
# Beta <- matrix(rnorm(d*r), d, r)
# nz_vec <- 1:6
# s <- 6
# True_sp <- solve(Delta) %*% Gamma
# 
# Data <- function(N){
#   # y <- rnorm(N)
#   y <- runif(N,-2,2)
#   p <- dim(Gamma)[1]
#   nobs <- length(y)
#   Fmat <- matrix(0, nobs, 4)
#   Fmat[,1] <- y
#   Fmat[,2] <- y^2
#   Fmat[,3] <- y^3
#   Fmat[,4] <- exp(y)
#   eps <- mvrnorm(nobs, rep(0,p), Delta)
#   x <- Fmat %*% t(Beta) %*% t(Gamma) + eps
#   list(x = x, y = y)
# }
# 
# sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
# intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
# pfc_params <- list(lambda.factor = 0.2, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = FALSE)

# #############################################

times <- 1
output <- sapply(seq_len(times), function(i){
  cat("Time", i, '\n')
  data_train <- Data(N)
  data_val <- Data(N_val)
  ssdrsir_fit <- ssdr_func(data_train$x, data_train$y, data_val$x, data_val$y, H = sir_params$H,
                           type = 'sir',
                           lambda.factor = sir_params$lambda.factor,
                           lam_fac_msda = sir_params$lam_fac_msda,
                           lam_fac_ssdr = sir_params$lam_fac_ssdr)
  ssdrintra_fit <- ssdr_func(data_train$x, data_train$y, data_val$x, data_val$y, H = intra_params$H,
                             type = 'intra',
                             lambda.factor = intra_params$lambda.factor,
                             lam_fac_msda = intra_params$lam_fac_msda,
                             lam_fac_ssdr = intra_params$lam_fac_ssdr,
                             cut_y = intra_params$cut_y)
  ssdrpfc_fit <- ssdr_func(data_train$x, data_train$y, data_val$x, data_val$y, type = 'pfc',
                           lambda.factor = pfc_params$lambda.factor,
                           lam_fac_msda = pfc_params$lam_fac_msda,
                           lam_fac_ssdr = pfc_params$lam_fac_ssdr,
                           cut_y = pfc_params$cut_y)
  LassoSIR_fit <- LassoSIR(data_train$x, data_train$y, H = 5, choosing.d = 'automatic')
  # CovSIR_fit <- CovSIR(data_train$x, data_train$y, Ks = 1, lambdas = seq(0.2,2,by=0.5)*sqrt(log(p)/N))
  # lasso_fit <- lasso_func(data_train$x, data_train$y)[-1,1,drop=FALSE] # the first is zero intercept
  # rifle_fit <- rifle_func(data_train$x, data_train$y, k = s, type = 'sir')
  # # DTSIR_fit <- DT_SIR(data_train$x, data_train$y, s = N/log(N), H = 5)
  # DTSIR_fit <- DT_SIR(data_train$x, data_train$y, k = s, H = 5)
  
  B_ssdrsir <- ssdrsir_fit$mat
  B_ssdrintra <- ssdrintra_fit$mat
  B_ssdrpfc <- ssdrpfc_fit$mat
  B_LassoSIR <- LassoSIR_fit$beta
  # B_CovSIR <- CovSIR_fit$mat
  # B_lasso <- lasso_fit
  # B_rifle <- rifle_fit
  # B_DTSIR <- DTSIR_fit
  
  # calculate C, IC, subspace distance after we obtain estimated matrix from each method.
  if(is.null(B_ssdrsir)){
    C_IC_ssdrsir <- list(C = NA, IC = NA)
    r_ssdrsir <- NA
    dist_ssdrsir <- NA
  }else{
    C_IC_ssdrsir <- C_IC(B_ssdrsir, 1:p, nz_vec)
    r_ssdrsir <- ssdrsir_fit$results$r_ssdr
    dist_ssdrsir <- subspace_2(True_sp, svd(B_ssdrsir)$u[,1:r_ssdrsir, drop = FALSE])
  }

  if(is.null(B_ssdrintra)){
    C_IC_ssdrintra <- list(C = NA, IC = NA)
    r_ssdrintra <- NA
    dist_ssdrintra <- NA
  }else{
    C_IC_ssdrintra <- C_IC(B_ssdrintra, 1:p, nz_vec)
    r_ssdrintra <- ssdrintra_fit$results$r_ssdr
    dist_ssdrintra <- subspace_2(True_sp, svd(B_ssdrintra)$u[,1:r_ssdrintra, drop = FALSE])
  }

  if(is.null(B_ssdrpfc)){
    C_IC_ssdrpfc <- list(C = NA, IC = NA)
    r_ssdrpfc <- NA
    dist_ssdrpfc <- NA
  }else{
    C_IC_ssdrpfc <- C_IC(B_ssdrpfc, 1:p, nz_vec)
    r_ssdrpfc <- ssdrpfc_fit$results$r_ssdr
    dist_ssdrpfc <- subspace_2(True_sp, svd(B_ssdrpfc)$u[,1:r_ssdrpfc, drop = FALSE])
  }

  # C_IC_LassoSIR <- C_IC_cut(B_LassoSIR, 1:p, nz_vec)
  C_IC_LassoSIR <- C_IC(B_LassoSIR, 1:p, nz_vec)
  r_LassoSIR <- LassoSIR_fit$no.dim
  dist_LassoSIR <- subspace_2(True_sp, B_LassoSIR)
  # print(svd(B_LassoSIR)$d)
  
  # C_IC_CovSIR <- C_IC_cut(B_CovSIR, 1:p, nz_vec)
  # r_CovSIR <- CovSIR_fit$r
  # dist_CovSIR <- subspace_2(True_sp, B_CovSIR)
  
  # C_IC_lasso <- C_IC(B_lasso, 1:p, nz_vec)
  # r_lasso <- 1
  # dist_lasso <- subspace_2(True_sp, B_lasso)
  # 
  # C_IC_rifle <- C_IC(B_rifle, 1:p, nz_vec)
  # r_rifle <- 1
  # dist_rifle <- subspace_2(True_sp, B_rifle)
  # 
  # C_IC_DTSIR <- C_IC(B_DTSIR, 1:p, nz_vec)
  # r_DTSIR <- 1
  # dist_DTSIR <- subspace_2(True_sp, B_DTSIR)
  
  # list(C_CovSIR = C_IC_CovSIR$C, IC_CovSIR = C_IC_CovSIR$IC, r_CovSIR = r_CovSIR, dist_CovSIR = dist_CovSIR)
  list(C_ssdrsir = C_IC_ssdrsir$C, IC_ssdrsir = C_IC_ssdrsir$IC, r_ssdrsir = r_ssdrsir, dist_ssdrsir = dist_ssdrsir,
       C_ssdrintra = C_IC_ssdrintra$C, IC_ssdrintra = C_IC_ssdrintra$IC, r_ssdrintra = r_ssdrintra, dist_ssdrintra = dist_ssdrintra,
       C_ssdrpfc = C_IC_ssdrpfc$C, IC_ssdrpfc = C_IC_ssdrpfc$IC, r_ssdrpfc = r_ssdrpfc, dist_ssdrpfc = dist_ssdrpfc,
       C_LassoSIR = C_IC_LassoSIR$C, IC_LassoSIR = C_IC_LassoSIR$IC,  r_LassoSIR = r_LassoSIR, dist_LassoSIR = dist_LassoSIR)
  #      C_lasso = C_IC_lasso$C, IC_lasso = C_IC_lasso$IC,  r_lasso = r_lasso, dist_lasso = dist_lasso,
  #      C_rifle = C_IC_rifle$C, IC_rifle = C_IC_rifle$IC,  r_rifle = r_rifle, dist_rifle = dist_rifle,
  #      C_DTSIR = C_IC_DTSIR$C, IC_DTSIR = C_IC_DTSIR$IC,  r_DTSIR = r_DTSIR, dist_DTSIR = dist_DTSIR)
})

# prof2 <- profvis(a <- replicate(2, run_func()))


# save(output, file = "/Users/cengjing/Desktop/output")


write.table(output, "/Users/cengjing/Desktop/test_ssdr_1")
# write.table(svB, "/Users/cengjing/Desktop/test_ssdr_1")
# write.table(svC, "/Users/cengjing/Desktop/test_ssdr_2")

