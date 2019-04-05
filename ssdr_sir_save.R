library(msda)
library(MASS)
library(methods)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/my_msda.R")

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

# Prediction function
lda_pred <- function(xtrain, ytrain, theta, xtest){
  xfit <- xtrain %*% theta
  l <- lda(xfit, ytrain)
  pred <- predict(l, xtest %*% theta)$class
  return(pred)
}

# subspace function
subspace <- function(A,B){
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  d <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*d))
}

subspace_2 <- function(Beta, B){
  Pa <- qr.Q(qr(Beta))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  d <- dim(Beta)[2]
  result <- norm(Pa - Pa %*% Pb %*% Pa, type = 'F')/sqrt(d)
  return(result)
}

# input the ssdr returned Beta matrix, returns corresponding predictions
predict_ssdr <- function(x_train, y_train, fit, newx){
  mat <- fit$beta
  rank <- fit$rank
  prior <- sapply(1:K, function(x) sum(y_train==x)/length(y_train))
  n.col <- length(mat)
  n.row <- nrow(newx)
  pred <- matrix(0,n.row,n.col)
  for(i in 1:n.col){
    beta <- mat[[i]]
    nz <- sum(beta[,1] != 0)
    if(is.null(beta) || nz == 0){
      # Sometimes ssdr is not converge, so we leave the corresponding matrix null, in this case we use
      # prior to predict. Or if matrix is a zero matrix then we also use prior.
      pred[,i] <- which.max(prior)
    }else{
      r <- rank[i]
      if(r == 0){
        pred[,i] <- which.max(prior)
      }else{
        subset <- svd(beta)$u[,1:r,drop = FALSE]
        pred[,i] <- lda_pred(x_train,y_train,subset,newx)
      }
    }
  }
  return(pred)
}

# rank function

rank_func <- function(B, thrd){
  d <- svd(B)$d
  d[d<thrd] <- 0
  i <- sum(d!=0)
  return(i)
}


# Draw the plot of the ratio of singular values

sv_plot <- function(sv){
  tmp <- rep(0,length(sv)-1)
  for(i in 1:(length(sv)-1)){
    tmp[i] <- (sv[i] + 0.001)/(sv[i+1] + 0.001)
  }
  plot(tmp, main = "Ratio of consecutive singular values", ylab = "ratio")
  plot(sv, main = "Singular values", ylab = "singular values")
}

# Calculate the marginal covariance matrix for x and within-class means
prep <- function(x, y){
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  nclass <- as.integer(length(unique(y)))
  if (nclass != H)
    stop(cat('Class number should be equal to', as.character(H), '!'))
  sigma <- cov(x)
  mu <- matrix(0, nvars, nclass)
  for (i in 1:nclass){
    mu[, i] <- apply(x[y == i, ], 2, mean)
  }
  output <- list(sigma = sigma, mu = mu)
  return(output)
}

# If some whole column is significantly small then we cut the columns to zero
cut_mat <- function(Beta, thrd, rank){
  l <- length(Beta)
  for (i in 1:l){
    if(is.null(Beta[[i]])) next
    mat <- as.matrix(Beta[[i]])
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    r <- rank[i]
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
    vec <- as.vector(mat)
    vec[abs(vec) < thrd] <- 0
    Beta[[i]] <- matrix(vec, nobs, nvars)
    # tmp <- apply(mat, 2, function(x){all(abs(x) < thrd)})
    # mat[,tmp] <- 0
    }
  }
  return(Beta)
}

orth_mat <- function(Beta, rank){
  l <- length(Beta)
  for (i in 1:l){
    mat <- Beta[[i]]
    r <- rank[i]
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    # If rank is equal to zero, then set Beta to zero
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
      Beta[[i]] <- svd(mat)$u[,1:r, drop=FALSE]
    }
  }
  return(Beta)
}

# ssdr function returns corresponding B matrices and other evaluation stuff

########################################  Error code  ##########################################
# case1: jerr < -10000, exit because of non-sparsity from msda, we stop trying more lam2s.
# case2: jerr = 404, exit because we reach the maximum iteration time of ssdr. And we leave matrix NULL
# case3: jerr = 1, succeed.
###############################################################################################


ssdr <- function(sigma, mu, nobs, nvars, lam1, lam2, gam, pf=rep(1, nvars), dfmax=nobs, 
                 pmax=min(dfmax * 2 + 20, nvars), eps=1e-04, maxit=1e+06, sml=1e-06, verbose = FALSE,
                 maxit_outer=1e+3, eps_outer=1e-3){
  
  flmin <- as.double(1)
  nlam <- as.integer(1)
  vnames <- as.character(1:nvars)
  nk <- as.integer(dim(mu)[2])
  pf <- as.double(pf)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  eps <- as.double(eps)
  maxit <- as.integer(maxit)
  sml <- as.double(sml)
  verbose <- as.integer(verbose)
  maxit_outer <- as.integer(maxit_outer)
  eps_outer <- as.double(eps_outer)
  
  sigma0 <- sigma
  mu0 <- mu
  n1 <- length(lam1)
  n2 <- ncol(lam2)
  n3 <- length(gam)
  nparams <- n1*n2*n3
  
  # mat <- as.list(seq_len(nparams))
  # step_final <- as.list(seq_len(nparams))     # To store the iteration times of each run
  # time_final <- as.list(seq_len(nparams))     # To store the running time of each run
  # lam1_list <- as.list(seq_len(nparams))
  # lam2_list <- as.list(seq_len(nparams))
  # gamma_list <- as.list(seq_len(nparams))
  # r_list <- as.list(seq_len(nparams))
  # 
  # sv_list_B <- as.list(seq_len(nparams))
  # sv_list_C <- as.list(seq_len(nparams))
  
  mat <- vector("list", nparams)
  step_final <- vector("list", nparams)     # To store the iteration times of each run
  time_final <- vector("list", nparams)     # To store the running time of each run
  lam1_list <- vector("list", nparams)
  lam2_list <- vector("list", nparams)
  gamma_list <- vector("list", nparams)
  r_list <- vector("list", nparams)
  
  sv_list_B <- vector("list", nparams)
  sv_list_C <- vector("list", nparams)
  
  nlam_ssdr <- 0
  
  for(i in 1:n1){
    ulam <- as.double(lam1[i])
    
    for(j in 1:n3){
      gamma <- gam[j]
      
      for(k in 1:n2){

        lambda2 <- lam2[j,k]
      
        # Maximal interation for outer loop
        sigma <- sigma0 + gamma*diag(rep(1,ncol(sigma0)), ncol(sigma0),ncol(sigma0))
        ##################################
        # SSDR
        ##################################
        # Initialize three matrices
        Bold <- matrix(0,dim(mu0)[1], dim(mu0)[2])
        Cold <- matrix(0,dim(mu0)[1], dim(mu0)[2])
        etaold <- matrix(0,dim(mu0)[1], dim(mu0)[2])
        
        # The MAIN loop of SSDR method
        step_ssdr <- 0    
        
        start_time <- Sys.time()
      
        repeat{
         
          step_ssdr <- step_ssdr + 1
          
          # Update B
          
          mu <- t(mu0 - etaold + gamma * Cold)
          fit <- .Fortran("msda", obj = double(nlam), nk, nvars, as.double(sigma), 
                          as.double(mu), pf, dfmax, pmax, nlam, flmin, ulam, 
                          eps, maxit, sml, verbose, nalam = integer(1), theta = double(pmax * nk * nlam), 
                          itheta = integer(pmax), ntheta = integer(nlam), 
                          alam = double(nlam), npass = integer(1), jerr = integer(1))
          
          # If jerr != 0, msda function returns abnormal results
          if (fit$jerr != 0){
            jerr <- fit$jerr
            break
          }
  
          outlist <- formatoutput(fit, maxit, pmax, nvars, vnames, nk)
          Bnew <- as.matrix(outlist$theta[[1]])
          
          # Update C
          Btemp <- Bnew + 1/gamma * etaold
          r <- svd(Btemp)
          U <- r$u
          V <- r$v
          D <- r$d
          lamtemp <- sapply(D, FUN = function(x) {max(0, x-lambda2/gamma)})
          Cnew <- U %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(V)
          
          # Update mu
          etanew <- etaold + gamma * (Bnew - Cnew)
          
          # Exit condition
          # the success code is 1
          if(max(abs(Bnew - Bold)) < eps_outer){
            jerr <- 1
            break
          }
          if(step_ssdr > maxit_outer){
            jerr <- 404
            break
          }
          
          Bold <- Bnew
          Cold <- Cnew
          etaold <- etanew
          
        }# End of repeat 
        end_time <- Sys.time()  # The time for each repeat
        
        # print(jerr)
        
        # If we get non-sparse matrix for msda, stop here, and leave the rest of matrices, svB, svC, 
        # etc. as NULL
        if(jerr < -10000){
          break
        }
        
        index <- (i-1)*n2*n3 + (j-1)*n2 + k
        nlam_ssdr <- nlam_ssdr + 1
        step_final[[index]] <- step_ssdr
        time_final[[index]] <- difftime(end_time, start_time, units = "secs")
        
        # If jerr == 404, then maximal iteration is reached, we leave the matrix as null
        # if(jerr==404){
        #   mat[index] <- list(NULL)
        #   r_list[[index]] <- NA
        #   
        #   sv_list_B[[index]] <- rep(NA, min(dim(Bnew)))
        #   sv_list_C[[index]] <- rep(NA, min(dim(Cnew)))
        #   
        #   }
        
        # If jerr == 1, then procedure converges, we save the matrix and sv.
        if(jerr==1){
          mat[[index]] <- Bnew
          r_list[[index]] <- rank_func(Bnew, thrd = 1e-3)
          
          # save the singular values of each candidates matrix B and C
          sv_list_B[[index]] <- svd(Bnew)$d
          sv_list_C[[index]] <- svd(Cnew)$d
          
          lam1_list[[index]] <- ulam
          gamma_list[[index]] <- gamma
          lam2_list[[index]] <- lambda2
        }
        
      }# End of lambda2
      
      # If exit because of non-sparsity from msda, we stop trying more lam2s or gammas, step to the larger lambda1
      if(jerr < -10000) break
      
    }# End of gam
    
  }# End of lambda1
  
  return(list(beta = mat, rank=r_list, step = step_final, time_ssdr = time_final, nlam_ssdr = nlam_ssdr, 
              lam1_list = lam1_list, lam2_list = lam2_list, gamma_list = gamma_list, sv_list_B = sv_list_B,
              sv_list_C = sv_list_C))
  
}

######################## evaluation ########################

eval_val <- function(Beta, x, y, slices){
  y_sliced <- cut(y, breaks = slices, include.lowest = TRUE, labels = FALSE)
  tmp <- prep(x, y_sliced)
  sigma <- tmp$sigma
  mu <- tmp$mu
  l <- length(Beta)
  result <- rep(0, l)
  for (i in 1:l){
    mat <- as.matrix(Beta[[i]])
    result[i] <- 0.5 * sum(diag(t(mat) %*% sigma %*% mat - 2 * t(mu) %*% mat), na.rm = TRUE)
  }
  return(result)
}

eval_val_rmse <- function(Beta, x, y, slices = NULL){
  l <- length(Beta)
  # seq_len(l) is 1:l, but faster
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      # print(c(dim(x), dim(mat)))
      xnew <- x %*% mat
      fit <- lm(y~xnew)
      rmse <- sqrt(mean((fit$residuals)^2))
      rmse 
    }
  })
  result
}

eval_val_rmse_rank <- function(Beta, x, y, rank, slices = NULL){
  l <- length(Beta)
  result <- rep(0, l)
  for (i in 1:l){
    r <- rank[i]
    mat <- as.matrix(Beta[[i]])
    if(r != 0){
      mat <- svd(as.matrix(Beta[[i]]))$u[,1:r]
    }
    xnew <- x %*% mat
    fit <- lm(y~xnew)
    rmse <- sqrt(mean((fit$residuals)^2))
    result[i] <- rmse
  }
  return(result)
}


eval_val_rmse_2 <- function(Beta, x, y, slices){
  y_sliced <- cut(y, breaks = slices, include.lowest = TRUE, labels = FALSE)
  l <- length(Beta)
  result <- rep(0, l)
  for (i in 1:l){
    mat <- as.matrix(Beta[[i]])
    xnew <- x %*% mat
    fit <- lm(y_sliced~xnew)
    rmse <- sqrt(mean((fit$residuals)^2))
    result[i] <- rmse
  }
  return(result)
}

eval_val_cart <- function(Beta, xtrain, ytrain, xval, yval, slices_val){
  y_val_sliced <- cut(yval, breaks = slices_val, include.lowest = TRUE, labels = FALSE)
  l <- length(Beta)
  errors <- rep(0, l)
  for (i in 1:l){
    mat <- as.matrix(Beta[[i]])
    xnew_tr <- xtrain %*% mat
    xnew_val <- as.data.frame(xval %*% mat)
    fit <- rpart(ytrain~xnew_tr, method = "class", control=rpart.control(minsplit=20, cp=0.001))
    pred <- predict(fit, xnew_val, type = "class")
    errors[i] <- mean(y_val_sliced != pred)
  }
  return(errors)
}

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
# p <- 1000  # Dimension of observations
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
# nz_vec <- 1:20
# r <- 1
# 
# params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)
# True_sp <- Beta
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta + 0.5 * rnorm(nobs)
#   list(x = x, y = y)
# }

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
# r <- 2
# 
# params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)
# True_sp <- Beta
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- (x %*% Beta[,1])/(0.5+(x %*% Beta[,2] + 1.5)^2) + 0.2 * rnorm(nobs)
#   list(x = x, y = y)
# }

# #############  Model 4 #############
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
# r <- 2
# True_sp <- Beta
# 
# params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- abs((x %*% Beta[,1]) / 4 + 2)^3 * sign(x %*% Beta[,2]) + 0.2 * rnorm(nobs)
#   list(x = x, y = y)
# }

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
# nz_vec <- 1:6
# r <- 2
# True_sp <- Beta
# 
# params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)
# 
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2]) + 0.2 * rnorm(nobs)
#   list(x = x, y = y)
# }

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
# nz_vec <- 1:6
# r <- 2
# True_sp <- Beta
# 
# params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)
# Data <- function(N){
#   x <- Train(N, Mu, Sigma)
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.2 *  rnorm(nobs) )
#   list(x = x, y = y)
# }

#############  Model 7 #############
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
nz_vec <- 1:6
r <- 2
True_sp <- Beta

params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)

Data <- function(N){
  x <- Train(N, Mu, Sigma)
  nobs <- dim(x)[1]
  y <- x %*% Beta[,1] * (2 + (x %*% Beta[,2]) / 3)^2 + 0.2 * rnorm(nobs)
  list(x = x, y = y)
}

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
# r <- 1
# True_sp <- Beta
# 
# params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = FALSE)
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
# True_sp <- solve(Delta) %*% Gamma
# 
# params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)
# Data <- function(N){
#   y <- rnorm(N)
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

# #############################################

run_func <- function(time=1, lambda.factor=0.5, lam_fac_msda=0.8, lam_fac_ssdr=0.8, cut_y=TRUE){
  
  cat('Times', as.character(time), '...\n')
  # Generate training, validation and testing dataset respectively
  
  data_train <- Data(N)
  x_train <- data_train$x
  y_train <- data_train$y
  
  # # Slice y
  # y_breaks_tr <- as.numeric(quantile(y_train, probs=seq(0,1, by=1/H), na.rm=TRUE))
  # y_train <- cut(y_train, breaks = y_breaks_tr, include.lowest = TRUE, labels = FALSE)
  
  # validation dataset
  data_val <- Data(N_val)
  x_val <- data_val$x
  y_val <- data_val$y

  # # Slice y_val
  # y_breaks_val <- c(-Inf, as.numeric(y_breaks_tr[2:H]), Inf)
  # y_val <- cut(y_val, breaks = quantile(y_val, probs=seq(0,1, by=1/H), na.rm=TRUE), 
  #                include.lowest = TRUE, labels = FALSE)


  
  ##################################
  # Bayes error
  ##################################

  #### The start of our methods
  start_time_tot <- Sys.time()
  
  ################################################
  # MSDA
  ################################################
  
  nlam_msda <- 10 # the number of lambdas in msda

  start_time <- Sys.time()
  fit_1 <- my_msda(x_train, y_train, nlambda=nlam_msda, maxit=1e3, lambda.factor=lambda.factor, cut_y=cut_y)
  end_time <- Sys.time()
  time_msda <- difftime(end_time, start_time, units = "secs")/nlam_msda
  
  sigma0 <- as.matrix(fit_1$sigma)
  mu0 <- as.matrix(fit_1$mu)
  
  ######################################
  ## Print true estimation
  # print(svd(solve(sigma0) %*% mu0)$d)
  # next
  ######################################
  
  lam_msda <- fit_1$lambda
  Beta_msda <- fit_1$theta
  
  # # Count the number of non-zero
  # nz_msda <- rep(0,length(Beta_msda))
  # for (i in 1:length(Beta_msda)){
  #   mat <- Beta_msda[[i]]
  #   nz_msda[i] <- sum(apply(mat, 1, function(x) any(x!=0)))
  # }
  
  rank_msda <- sapply(seq_len(length(Beta_msda)), function(i){
    mat <- Beta_msda[[i]]
    if(is.null(mat)){
      NA
    }else{
      rank_func(mat, thrd = 1e-3)
    }
  })
  
  # Cut negligible entries to zero
  Beta_msda <- cut_mat(Beta_msda, 1e-3, rank_msda)
  
  # rank_msda <- rep(0,length(Beta_msda))
  # for (i in 1:length(Beta_msda)){
  #   mat <- Beta_msda[[i]]
  #   rank_msda[i] <- rank_func(mat, thrd = 1e-3)
  # }
  
  # validata
  start_time <- Sys.time()
  eval_msda <- eval_val_rmse(Beta_msda, x_val, y_val, y_breaks_val)
  end_time <- Sys.time()
  time_eval_msda <- difftime(end_time, start_time, units = "secs")
  
  # The optimal lambda1
  id_min_msda <- which.min(eval_msda)
  lam1_min_msda <- lam_msda[id_min_msda]
  
  # calculate C, IC, Frobenious distance, rank and subspace distance
  B_msda <- as.matrix(Beta_msda[[id_min_msda]])
  tmp <- apply(B_msda, 1, function(x) any(x!=0))
  C_msda <- sum(which(tmp) %in% nz_vec)/length(nz_vec)
  IC_msda <- sum(which(tmp) %in% setdiff(1:p, nz_vec))/(p - length(nz_vec))
  r_msda <- rank_msda[id_min_msda]
  
  
  ################################################
  # SSDR
  ################################################
  
  # Lambda1 candidates
  # lam1 <- (lam1_min_msda)*seq(1.5,0.6,-0.1)
  # n1 <- length(lam1)
  n1 <- 10
  lam_fac_msda <- lam_fac_msda
  lam1 <- lam1_min_msda*lam_fac_msda^seq(0,(n1-1))
  
  # Gamma candidates
  gamma <- c(10,30,50)
  n3 <- length(gamma)
  
  # Lambda2 candidates
  n2 <- 15   # we select n2 lambda2 for each gamma
  lam_fac_ssdr <- lam_fac_ssdr
  d <- svd(B_msda)$d
  lam2 <- d[1] * matrix(gamma, ncol = 1) %*% matrix(lam_fac_ssdr^seq((n2-1),0), nrow = 1)
  
  # if lam2 just contains one single value 0, then ssdr just degenerated to msda
  if (all(lam2 == 0)){
    
    print("All lambda2 are zero, degenerate to msda")
    results <- c(C_msda, IC_msda, C_msda, IC_msda, r_msda, r_msda, NA, lam1_min_msda, id_min_msda,
                 lam1_min_msda, NA, NA, which(lam1 == lam1_min_msda), NA, NA, NA, time_msda, time_eval_msda,
                 NA, NA, NA)
    svB <- NA
    svC <- NA
    
    return(list(results = results, svB = svB, svC = svC))
    
  }else{
    
    # fit with ssdr
    nobs <- as.integer(dim(x_train)[1])
    nvars <- as.integer(dim(x_train)[2])
    
    fit_2 <- ssdr(sigma0, mu0, nobs, nvars, lam1, lam2, gamma)
    
    Beta_ssdr <- fit_2$beta
    
    # In some cases, all the Beta is null because the Fortran code didn't return a converaged B matrix 
    if (sum(sapply(Beta_ssdr, is.null)) == n2*n3) {
      print("No converged matrix returned")
      results <- c(C_msda, IC_msda, NA, NA, r_msda, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      return(list(results = results, svB = NA, svC = NA))
    }
    
    # ##############
    # nz_ssdr <- c()
    # for (i in 1:length(Beta_ssdr)){
    #   B <- Beta_ssdr[[i]]
    #   if(is.null(B)){
    #     nz_ssdr <- c(nz_ssdr, NA)
    #   }else{
    #     nz_ssdr <- c(nz_ssdr, sum(apply(B,1,function(x){any(x!=0)})))
    #   }
    # }
    
    gamma_list <- fit_2$gamma_list
    lam1_list <- fit_2$lam1_list
    lam2_list <- fit_2$lam2_list
    rank_ssdr <- fit_2$rank
    step <- unlist(fit_2$step)
    time_ssdr <- unlist(fit_2$time_ssdr)
    
    sv_list_B <- fit_2$sv_list_B
    sv_list_C <- fit_2$sv_list_C
    
    
    
    # Cut negligible columns to zero
    Beta_ssdr <- cut_mat(Beta_ssdr, 1e-3, rank_ssdr)
    
    # rank_ssdr <- rep(0,length(Beta_ssdr))
    # for (i in 1:length(Beta_ssdr)){
    #   mat <- Beta_ssdr[[i]]
    #   rank_ssdr[i] <- rank_func(mat, thrd = 1e-3)
    # }
    
    # validate
    start_time <- Sys.time()
    eval_ssdr <- eval_val_rmse(Beta_ssdr, x_val, y_val, y_breaks_val)
    end_time <- Sys.time()
    time_eval_ssdr <- difftime(end_time, start_time, units = "secs")
    
    # The optimal lambda1 and lambda2 
    #########################
    

    id_min_ssdr <- which.min(eval_ssdr)
    lam1_min_ssdr <- lam1_list[[id_min_ssdr]]
    lam2_min_ssdr <- lam2_list[[id_min_ssdr]]
    gamma_min_ssdr <- gamma_list[[id_min_ssdr]]
    id_lam1 <- which(lam1_min_ssdr == lam1)
    id_lam2 <- which(lam2_min_ssdr == lam2, arr.ind = TRUE)[2]
    id_gamma <- which(gamma_min_ssdr == gamma)
    
    #########################
    
    B_ssdr <- Beta_ssdr[[id_min_ssdr]]
    
    if(is.null(B_ssdr)){
      print("Optimal matrix is a null matrix")
      
      results <- c(C_msda, IC_msda, NA, NA, r_msda, NA, NA, lam1_min_msda,
                   id_min_msda, lam1_min_ssdr, lam2_min_ssdr, gamma_min_ssdr, id_lam1, id_lam2, id_gamma, mean(step), 
                   time_msda, time_eval_msda, mean(time_ssdr), time_eval_ssdr, NA)
      svB <- NA
      svC <- NA
      return(list(results = results, svB = svB, svC = svC))
      
    }else{
      # Calculate C, IC, Frobinious distance, subspace distance
      tmp <- apply(B_ssdr, 1, function(x) any(x!=0))
      C_ssdr <- sum(which(tmp) %in% nz_vec)/length(nz_vec)
      IC_ssdr <- sum(which(tmp) %in% setdiff(1:p, nz_vec))/(p - length(nz_vec))
      r_ssdr <- rank_ssdr[[id_min_ssdr]]
      
      # sub_ssdr <- subspace_2(Beta, svd(B_ssdr)$u[,1:r_ssdr, drop = FALSE])
      sub_ssdr <- subspace_2(True_sp, svd(B_ssdr)$u[,1:r_ssdr, drop = FALSE])
      
      # save the singular values of each optimal matrix B and C
      svB <- sv_list_B[[id_min_ssdr]]
      svC <- sv_list_C[[id_min_ssdr]]
      
      # record total time iff we got converged matrix
      end_time_tot <- Sys.time()
      time_total <- difftime(end_time_tot, start_time_tot, units = "secs")
      
      results <- c(C_msda, IC_msda, C_ssdr, IC_ssdr, r_msda, r_ssdr, sub_ssdr, lam1_min_msda,
                   id_min_msda, lam1_min_ssdr, lam2_min_ssdr, gamma_min_ssdr, id_lam1, id_lam2, id_gamma, mean(step), 
                   time_msda, time_eval_msda, mean(time_ssdr), time_eval_ssdr, time_total)
      
      return(list(results = results, svB = svB, svC = svC))
    }
  
  }

}

# Use apply function to avoid for loop

times <- 5
output <- sapply(seq_len(times), run_func,
                 lambda.factor = params$lambda.factor, lam_fac_msda = params$lam_fac_msda, 
                 lam_fac_ssdr = params$lam_fac_ssdr, cut_y = params$cut_y)

# The first row of output is results, second one is svB, third one is svC. Use do.call to bind them
results <- do.call(rbind, output[1,])
svB <- do.call(rbind, output[2,])
svC <- do.call(rbind, output[3,])

# prof2 <- profvis(a <- replicate(2, run_func()))

results <- as.data.frame(results)
colnames(results) <- c("C_msda", "IC_msda", "C_ssdr", "IC_ssdr", "r_msda", "r_ssdr", "sub_ssdr", "lam1_min_msda",
                       "id_msda", "lam1_min_ssdr", "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", 
                       "step", "time_msda", "teval_msda", "time_ssdr", "teval_ssdr", "time_total")

# write.table(results, "/Users/cengjing/Desktop/test_ssdr_1")
# write.table(svB, "/Users/cengjing/Desktop/test_ssdr_1")
# write.table(svC, "/Users/cengjing/Desktop/test_ssdr_2")

