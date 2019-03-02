library(msda)
library(MASS)
library(methods)
library(rpart)
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

# Cut negligible entries to zero
cut_mat <- function(Beta, thrd, rank){
  l <- length(Beta)
  for (i in 1:l){
    mat <- Beta[[i]]
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    r <- rank[i]
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
      vec <- as.vector(mat)
      vec[abs(vec) < thrd] <- 0
      Beta[[i]] <- matrix(vec, nobs, nvars) 
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

ssdr <- function(lam1,lam2,gam){
  n1 <- length(lam1)
  n2 <- ncol(lam2)
  n3 <- length(gam)
  mat <- list()
  step_final <- c()     # To store the iteration times of each run
  time_final <- c()     # To store the running time of each run
  # jerr_list <- c()
  nlam_ssdr <- 0
  lam1_list <- c()
  lam2_list <- c()
  gamma_list <- c()
  r_list <- c()
  
  for(i in 1:n1){
    ulam <- as.double(lam1[i])
    
    for(k in 1:n3){
      gamma <- gam[k]
      
      for(j in 1:n2){
        lambda2 <- lam2[k,j]
      
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
        diff_B <- c()   
        
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
          lamtemp <- sapply(D, FUN = function(x) max(0, x-lambda2/gamma))
          Cnew <- U %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(V)
          
          # Update mu
          etanew <- etaold + gamma * (Bnew - Cnew)
          
          diff_B <- c(diff_B, norm(Bnew-Bold, type = "F"))
          
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
        
        # If we get non-sparse matrix for msda, stop here
        if(jerr < -10000){
          # jerr_list <- c(jerr_list, jerr)
          break
        }
        
        # If not, we save the matrix
        nlam_ssdr <- nlam_ssdr + 1
        step_final <- c(step_final, step_ssdr)
        time_final <- c(time_final, difftime(end_time, start_time, units = "secs"))
        # jerr_list <- c(jerr_list, jerr)
        
        # If jerr == 404, then maximal iteration is reached, we leave the matrix as null
        if(jerr==404){mat <- c(mat, list())}
        # If jerr == 1, then procedure converges.
        if(jerr==1){
          mat <- c(mat, list(Bnew))
        }
  
        lam1_list <- c(lam1_list, ulam)
        gamma_list <- c(gamma_list, gamma)
        lam2_list <- c(lam2_list, lambda2)
      }# End of lambda2
      
      # If exit because of non-sparsity from msda, we stop trying more lam2s.
      if(jerr < -10000) break
      
    }# End of gam
    
  }# End of lambda1
  
  # Record the rank
  for(i in 1:length(mat)){
     if(is.null(mat[[i]])){
        r_list <- c(r_list, NA)
     }else{
        r_list <- c(r_list, rank_func(mat[[i]], thrd = 1e-3))
     }
  }
  
  return(list(beta = mat, rank=r_list, step = step_final, time_ssdr = time_final, nlam_ssdr = nlam_ssdr, 
              lam1_list = lam1_list, lam2_list = lam2_list, gamma_list = gamma_list))
  
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
  result <- rep(0, l)
  for (i in 1:l){
    mat <- as.matrix(Beta[[i]])
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

#############  Model 2 #############
set.seed(1)

p <- 100  # Dimension of observations
N <- 1000  # Sample size
N_val <- 1000  # Sample size of validation dataset
H <- 5

Mu <- rep(0,p)
# Sigma <- CS_blk(0.5,800,4)
Sigma <- AR(0.5, p)
# Sigma <- diag(rep(1,p),p,p)
# Construct true Beta
Beta <- matrix(0, p, 1)
Beta[1:20,1] <- 1
# Beta[,1] <- Beta[,1]/norm(Beta[,1], type = '2')
nz_vec <- 1:20
r <- 1

model <- function(x, Beta){
  nobs <- dim(x)[1]
  y <- x %*% Beta + 0.5 * rnorm(nobs)
  return(y)
}

# #############  Model 3 #############
# set.seed(1)
# 
# p <- 800  # Dimension of observations
# N <- 1000 # Sample size
# N_val <- 1000  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# # Sigma <- diag(rep(1,p),p,p)
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# # Beta[1:10,1] <- 1
# # Beta[1:10,2] <- c(1,-1,1,-1,1,-1,1,-1,1,-1)
# Beta[1:6,1] <- 1
# Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
# nz_vec <- 1:6
# r <- 2
# 
# model <- function(x, Beta){
#   nobs <- dim(x)[1]
#   y <- (x %*% Beta[,1])/(0.5+(x %*% Beta[,2] + 1.5)^2) + 0.2 * rnorm(nobs)
#   return(y)
# }

# #############  Model 4 #############
# set.seed(1)
# 
# p <- 80  # Dimension of observations
# N <- 1000 # Sample size
# N_val <- 1000  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# # Sigma <- diag(rep(1,p),p,p)
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:6,1] <- 1
# # Beta[11:20,2] <- 1
# Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
## Beta[,1] <- Beta[,1]/norm(Beta[,1], type = '2')
## Beta[,2] <- Beta[,2]/norm(Beta[,2], type = '2')
# # nz_vec <- 1:20
# nz_vec <- 1:6
# r <- 2
# 
# model <- function(x, Beta){
#   nobs <- dim(x)[1]
#   y <- abs((x %*% Beta[,1]) / 4 + 2)^3 * sign(x %*% Beta[,2]) + 0.02 * rnorm(nobs)
#   return(y)
# }

# #############  Model 5 #############
# set.seed(1)
# 
# p <- 20  # Dimension of observations
# N <- 1000 # Sample size
# N_val <- 1000  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# # Sigma <- diag(rep(1,p),p,p)
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:6,1] <- 1
# # Beta[11:20,2] <- 1
# Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
# # Beta[,1] <- Beta[,1]/norm(Beta[,1], type = '2')
# # Beta[,2] <- Beta[,2]/norm(Beta[,2], type = '2')
# # nz_vec <- 1:20
# nz_vec <- 1:6
# r <- 2
# 
# model <- function(x, Beta){
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2]) + 0.2 * rnorm(nobs)
#   return(y)
# }

# #############  Model 6 #############
# set.seed(1)
# 
# p <- 20  # Dimension of observations
# N <- 1000 # Sample size
# N_val <- 1000  # Sample size of validation dataset
# H <- 5
# 
# Mu <- rep(0,p)
# Sigma <- AR(0.5, p)
# # Sigma <- diag(rep(1,p),p,p)
# # Construct true Beta
# Beta <- matrix(0, p, 2)
# Beta[1:6,1] <- 1
# # Beta[11:20,2] <- 1
# Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
# # nz_vec <- 1:20
# nz_vec <- 1:6
# r <- 2
# 
# model <- function(x, Beta){
#   nobs <- dim(x)[1]
#   y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.2 *  rnorm(nobs) )
#   return(y)
# }

#############  Model 7 #############
set.seed(1)

p <- 500  # Dimension of observations
N <- 1000 # Sample size
N_val <- 1000  # Sample size of validation dataset
H <- 5

Mu <- rep(0,p)
Sigma <- AR(0.5, p)
# Sigma <- diag(rep(1,p),p,p)
# Construct true Beta
Beta <- matrix(0, p, 2)
Beta[1:6,1] <- 1
# Beta[11:20,2] <- 1
Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
# nz_vec <- 1:20
nz_vec <- 1:6
r <- 2

model <- function(x, Beta){
  nobs <- dim(x)[1]
  y <- x %*% Beta[,1] * (2 + (x %*% Beta[,2]) / 3)^2 + 0.2 * rnorm(nobs)
  return(y)
}

# #############################################

times <- 10 # Simulation times
results <- matrix(0, times, 18)

nlam_ssdr <- c()

# sv_msda_list <- c()
# sv_ssdr_list <- c()

for(t in 1:times){
  
  cat("Time", as.character(t),'...... :)', '\n')

  # Generate training, validation and testing dataset respectively
  
  x_train <- Train(N, Mu, Sigma)
  y_train <- model(x_train, Beta)
  # Slice y
  # break points of y
  y_breaks_tr <- as.numeric(quantile(y_train, probs=seq(0,1, by=1/H), na.rm=TRUE))
  y_train <- cut(y_train, breaks = y_breaks_tr, include.lowest = TRUE, labels = FALSE)
  
  # validation dataset
  x_val <- Train(N_val, Mu, Sigma)
  y_val <- model(x_val, Beta)
  # break points of y_val
  y_breaks_val <- c(-Inf, as.numeric(y_breaks_tr[2:H]), Inf)
  # y_val <- cut(y_val, breaks = quantile(y_val, probs=seq(0,1, by=1/H), na.rm=TRUE), 
  #                include.lowest = TRUE, labels = FALSE)
  # 
  # x_test <- Train(N_test, Mu, Sigma)
  # y_test <- sign(x_test %*% Beta[,1]) * log(abs(x_test %*% Beta[,2] + 5)) + 0.2 * rnorm(N_test)

  
  ##################################
  # Bayes error
  ##################################

  #### The start of our methods
  start_time <- Sys.time()
  
  ################################################
  # MSDA
  ################################################
  
  nlam_msda <- 10 # the number of lambdas in msda

  fit_1 <- my_msda(x_train, y_train, nlambda = nlam_msda, maxit = 1e3, lambda.factor = 0.5)
  
  sigma0 <- as.matrix(fit_1$sigma)
  mu0 <- as.matrix(fit_1$mu)
  
  Beta_msda <- fit_1$theta
  
  # Count the number of non-zero
  nz_msda <- rep(0,length(Beta_msda))
  for (i in 1:length(Beta_msda)){
    mat <- Beta_msda[[i]]
    nz_msda[i] <- sum(apply(mat, 1, function(x) any(x!=0)))
  }
  
  rank_msda <- rep(0,length(Beta_msda))
  for (i in 1:length(Beta_msda)){
    mat <- Beta_msda[[i]]
    rank_msda[i] <- rank_func(mat, thrd = 0.001)
  }
  
  # Cut negligible entries to zero
  # Beta_msda <- cut_mat(Beta_msda, 1e-6, rank_msda)
  lam_msda <- fit_1$lambda
  
  # Beta_msda <- orth_mat(Beta_msda, rank_msda)
  Beta_msda <- cut_mat(Beta_msda, 1e-3, rank_msda)
  
  # validata
  
  # eval_msda <- eval_val(Beta_msda, x_val, y_val, y_breaks_val)
  eval_msda <- eval_val_rmse(Beta_msda, x_val, y_val, y_breaks_val)
  # eval_msda <- eval_val_rmse_2(Beta_msda, x_val, y_val, y_breaks_val)
  # eval_msda <- eval_val_cart(Beta_msda, x_train, y_train, x_val, y_val, y_breaks_val)
  
  # eval_true <- eval_val_rmse_2(list(Beta), x_val, y_val, y_breaks_val)
  
  # The optimal lambda1
  id_min_msda <- which.min(eval_msda)
  lam1_min_msda <- lam_msda[id_min_msda]
  
  # calculate C, IC, Frobenious distance, rank and subspace distance
  B_msda <- as.matrix(Beta_msda[[id_min_msda]])
  tmp <- apply(B_msda, 1, function(x) any(x!=0))
  C_msda <- sum(which(tmp) %in% nz_vec)/length(nz_vec)
  IC_msda <- sum(which(tmp) %in% setdiff(1:p, nz_vec))/(p - length(nz_vec))
  # Fnorm_msda <- norm(Beta - B_msda, type = 'F')
  r_msda <- rank_func(B_msda, thrd = 1e-3)
  
  # if(r_msda==0){
  #   sub_msda <- NA
  # }else{
  #   sub_msda <- subspace(svd(Beta)$u[,1:r, drop=FALSE], svd(B_msda)$u[,1:r_msda, drop=FALSE])
  # }
  
  ########### Save sv ##############
  # sv_msda_list <- rbind(sv_msda_list, svd(B_msda)$d)
  #################################
  
  
  # # Prediction error
  # pred_msda <- predict(fit_1, x_test)[,id_min_msda]
  # e_msda <- 1 - sum(pred_msda == y_test)/length(y_test)
  
  # Draw the singular values plot
  # sv_plot(svd(B_msda)$d)

  ################################################
  # SSDR
  ################################################
  
  flmin <- as.double(1)
  nlam <- as.integer(1)
  nk <- as.integer(dim(mu0)[2])
  nobs <- as.integer(dim(x_train)[1])
  nvars <- as.integer(dim(x_train)[2])
  pf <- as.double(rep(1, nvars))
  dfmax <- as.integer(nobs)
  # dfmax <- as.integer(nvars)
  pmax <- as.integer(min(dfmax * 2 + 20, nvars))
  # pmax <- as.integer(nvars)
  eps <- as.double(1e-04)
  maxit <- as.integer(1e+06)
  sml <- as.double(1e-06)
  verbose <- as.integer(FALSE)
  maxit_outer <- as.integer(1e+3) 
  eps_outer <- as.double(1e-3)
  vnames <- as.character(1:p)
  
  lam_fac_ssdr <- 0.5
  # We may need to shrink lam1 a little bit
  lam1 <- (lam1_min_msda)*seq(1.5,0.6,-0.1)
  n1 <- length(lam1)
  
  gamma <- c(10,30,50,80)
  # gamma <- 50
  # gamma <- c(10,20,30)
  n3 <- length(gamma)

  # Construct lambda2 candidates
  n2 <- 10   # we select n2 lambda2 for each gamma
  d <- svd(B_msda)$d
  # lam2 <- d[1]*gamma*lam_fac_ssdr^seq((n2-1),0)
  lam2 <- d[1] * matrix(gamma, ncol = 1) %*% matrix(lam_fac_ssdr^seq((n2-1),0), nrow = 1)
  
  # if lam2 just contains one single value 0, then ssdr just degenerated to msda
  if (all(lam2 == 0)){
    C_ssdr <- C_msda
    IC_ssdr <- IC_msda
    # e_ssdr <- e_msda
    r_ssdr <- r_msda
    sub_ssdr <- NA
    # Fnorm_ssdr <- Fnorm_msda
    lam1_min_ssdr <- lam1_min_msda
    lam2_min_ssdr <- NA
    ###
    id_lam1 <- which(lam1 == lam1_min_msda)
    id_lam2 <- NA
    ###
    gamma_min_ssdr <- NA
    step <- NA
    time_ssdr <- NA
  }else{
    
    # fit with ssdr
    fit_2 <- ssdr(lam1, lam2, gamma)
    
    Beta_ssdr <- fit_2$beta
    
    # In some cases, all the Beta is null because the Fortran code didn't return a converaged B matrix 
    if (sum(sapply(Beta_ssdr, is.null)) == n2*n3) {
      results[t,] <- c(C_msda, IC_msda, NA, NA, r_msda, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    ##############
    nz_ssdr <- c()
    for (i in 1:length(Beta_ssdr)){
      B <- Beta_ssdr[[i]]
      if(is.null(B)){
        nz_ssdr <- c(nz_ssdr, NA)
      }else{
        nz_ssdr <- c(nz_ssdr, sum(apply(B,1,function(x){any(x!=0)})))
      }
    }
    
    gamma_list <- fit_2$gamma_list
    lam1_list <- fit_2$lam1_list
    lam2_list <- fit_2$lam2_list
    rank_list <- fit_2$rank
    nlam_ssdr <- c(nlam_ssdr, fit_2$nlam_ssdr)
    step <- fit_2$step
    time_ssdr <- fit_2$time_ssdr
    
    # Cut negligible entries to zero
    # Beta_ssdr <- cut_mat(Beta_ssdr, 1e-6, rank_list)
    
    # Beta_ssdr <- orth_mat(Beta_ssdr, rank_list)
    Beta_ssdr <- cut_mat(Beta_ssdr, 1e-3, rank_list)
    
    # validate
    # eval_ssdr <- eval_val(Beta_ssdr, x_val, y_val, y_breaks_val)
    eval_ssdr <- eval_val_rmse(Beta_ssdr, x_val, y_val, y_breaks_val)
    # eval_ssdr <- eval_val_rmse_2(Beta_ssdr, x_val, y_val, y_breaks_val)
    
    # eval_true <- eval_val_rmse_2(list(Beta), x_val, y_val, y_breaks_val)
    # print(c(eval_true <= min(eval_msda), eval_true <= min(eval_ssdr)))
    
    # The optimal lambda1 and lambda2 
    #########################
    

    id_min_ssdr <- which.min(eval_ssdr)
    lam1_min_ssdr <- lam1_list[id_min_ssdr]
    lam2_min_ssdr <- lam2_list[id_min_ssdr]
    gamma_min_ssdr <- gamma_list[id_min_ssdr]
    # gamma_min_ssdr <- 10
    id_lam1 <- which(lam1_min_ssdr == lam1)
    id_lam2 <- which(lam2_min_ssdr == lam2, arr.ind = TRUE)[2]
    id_gamma <- which(gamma_min_ssdr == gamma)
    # id_gamma <- 1
    
    #########################
    
    B_ssdr <- Beta_ssdr[[id_min_ssdr]]
    
    if(is.null(B_ssdr)){
      C_ssdr <- NA
      IC_ssdr <- NA
      r_ssdr <- NA
      # sub_ssdr <- NA
      # Fnorm_ssdr <- NA
    }else{

      # Calculate C, IC, Frobinious distance, subspace distance
      tmp <- apply(B_ssdr, 1, function(x) any(x!=0))
      C_ssdr <- sum(which(tmp) %in% nz_vec)/length(nz_vec)
      IC_ssdr <- sum(which(tmp) %in% setdiff(1:p, nz_vec))/(p - length(nz_vec))
      # Fnorm_ssdr <- norm(Beta - B_ssdr, type = 'F')
      r_ssdr <- fit_2$rank[id_min_ssdr]
      
      sub_ssdr <- subspace_2(Beta, svd(B_ssdr)$u[,1:r_ssdr, drop = FALSE])
      # if(r_ssdr==0){
      #   sub_ssdr <- NA
      # }else{
      #   sub_ssdr <- subspace(svd(Beta)$u[,1:r, drop=FALSE], svd(B_ssdr)$u[,1:r_ssdr, drop=FALSE])
      # }
      
      ########### Save sv ##############
      # sv_ssdr_list <- rbind(sv_ssdr_list, svd(B_ssdr)$d)
      #################################
      
    }
    
    # B_ssdr_final <- list(beta=list(B_ssdr), rank=r_ssdr)
    # pred_ssdr <- predict_ssdr(x_train, y_train, B_ssdr_final, x_test)
    # e_ssdr <- 1 - sum(pred_ssdr == y_test)/length(y_test)
    
    # Draw the singular values plot
    # sv_plot(svd(B_ssdr)$d)
  }
  
  #####################################################################
  
  end_time <- Sys.time()
  time_total <- difftime(end_time, start_time, units = "secs")
  # store the prediction errors
  results[t,] <- c(C_msda, IC_msda, C_ssdr, IC_ssdr, r_msda, r_ssdr, sub_ssdr, lam1_min_msda,
                   id_min_msda, lam1_min_ssdr, lam2_min_ssdr, gamma_min_ssdr, id_lam1, id_lam2, id_gamma, mean(step), 
                   mean(time_ssdr), time_total)
}

results <- as.data.frame(results)
colnames(results) <- c("C_msda", "IC_msda", "C_ssdr", "IC_ssdr", "r_msda", "r_ssdr", "sub_ssdr", "lam1_min_msda",
                       "id_msda", "lam1_min_ssdr", 
                       "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", "step", "time_ssdr", "time_total")
write.table(results, "/Users/cengjing/Desktop/test_ssdr_1")