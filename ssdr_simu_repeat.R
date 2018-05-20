library(msda)
library(MASS)
library(methods)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")

p <- 800  #Dimension of observations
K <- 3    # The number of class
Nperclass <- 30  # The number of training observations in each class
Nperclass_test <- 300   # The number of testing data in each class


# Simulating train data function
Train <- function(n, mu, Sigma){
  m <- mvrnorm(n, mu[,1], Sigma)
  for (i in 2:ncol(mu)){
    m <- rbind(m, mvrnorm(n, mu[,i], Sigma))
  }
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
  u <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*u))
}

# ssdr function returns corresponding B matrices, training dataset and prior
ssdr <- function(lam1,lam2,gam){
  n1 <- length(lam1)
  n2 <- length(lam2)
  n3 <- length(gam)
  # prior <- sapply(1:K, function(x) sum(y==x)/length(y))
  mat <- vector(mode = "list", length = n1*n2*n3)
  
  # exe_times <- rep(0, n1*n2*n3)
  exe_times <- matrix(0, n1*n2*n3, 7)
  for(i in 1:n1){
    lambda1 <- lam1[i]
    ulam <- as.double(lambda1)
    
    for(j in 1:n2){
      lambda2 <- lam2[j]
      
      for(k in 1:n3){
        gamma <- gam[k]
        # Maximal interation for outer loop
        sigma <- sigma0 + gamma*diag(rep(1,ncol(sigma0)), ncol(sigma0),ncol(sigma0))
        ##################################
        # SSDR
        ##################################
        # Initialize three matrices
        Bold <- matrix(0,dim(delta0)[2], dim(delta0)[1])
        Cold <- matrix(0,dim(delta0)[2], dim(delta0)[1])
        muold <- matrix(0,dim(delta0)[2], dim(delta0)[1])
        
        # The MAIN loop of SSDR method
        step_ssdr <- 0
        
        start_time_all <- Sys.time()

        repeat{
          temp <- c()
          step_ssdr <- step_ssdr + 1
          
          # Update B
          start_time <- Sys.time()
          
          delta <- delta0 - t(muold) + gamma * t(Cold)
          fit <- .Fortran("msda", obj = double(nlam), nk, nvars, as.double(sigma), 
                          as.double(delta), pf, dfmax, pmax, nlam, flmin, ulam, 
                          eps, maxit, sml, verbose, nalam = integer(1), theta = double(pmax * nk * nlam), 
                          itheta = integer(pmax), ntheta = integer(nlam), 
                          alam = double(nlam), npass = integer(1), jerr = integer(1))
          
          
          if (fit$jerr != 0){
            jerr <- fit$jerr
            break
          }
          
          end_time <- Sys.time()
          time_1 <- difftime(end_time, start_time, units = "secs")
          
          start_time <- Sys.time()

          outlist <- formatoutput(fit, maxit, pmax, nvars, vnames, nk)
          Bnew <- as.matrix(outlist$theta[[1]])
          
          end_time <- Sys.time()
          time_2 <- difftime(end_time, start_time, units = "secs")
          
          # Update C
          start_time <- Sys.time()
          
          Btemp <- Bnew + 1/gamma * muold
          r <- svd(Btemp)
          U <- r$u
          V <- r$v
          D <- r$d
          lamtemp <- sapply(D, FUN = function(x) max(0, x-lambda2/gamma))
          Cnew <- U %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(V)
          
          end_time <- Sys.time()
          time_3 <- difftime(end_time, start_time, units = "secs")
          
          # Update mu
          start_time <- Sys.time()

          munew <- muold + gamma * (Bnew - Cnew)
          
          
          
          # Exit condition
          if(max(abs(Bnew - Bold)) < eps_outer){
            jerr <- 1
            temp <- rbind(temp, matrix(c(time_1, time_2, time_3, 0, 0), nrow = 1))
            break
          }
          if(step_ssdr > maxit_outer){
            jerr <- -2
            temp <- rbind(temp, matrix(c(time_1, time_2, time_3, 0, 0), nrow = 1))
            break
          }
          
          end_time <- Sys.time()
          time_4 <- difftime(end_time, start_time, units = "secs")
          
          start_time <- Sys.time()
          
          Bold <- Bnew
          Cold <- Cnew
          muold <- munew
          
          end_time <- Sys.time()
          time_5 <- difftime(end_time, start_time, units = "secs")
          temp <- rbind(temp, matrix(c(time_1, time_2, time_3, time_4, time_5), nrow = 1))
        }# End of repeat 
        
        temp <- colMeans(temp)
        end_time_all <- Sys.time()
        time_all <- difftime(end_time_all, start_time_all, units = "secs")
        
        exe_times[(i-1)*n2*n3+(j-1)*n3+k,] <- c(temp, time_all, step_ssdr)
        # exe_times[(i-1)*n2*n3+(j-1)*n3+k] <- (end_time - start_time)

        
        # If jerr == 1, then procedure converges. And if not, we leave the matrix NULL.
        if(jerr==1){
          mat[[(i-1)*n2*n3+(j-1)*n3+k]] <- Bnew
        }
        
      }
    }
  }
  return(list(Beta = mat, exe_time = colMeans(exe_times)))
}

# input the ssdr object, returns corresponding predictions
predict_ssdr <- function(x_train, y_train, mat, newx){
  prior <- sapply(1:K, function(x) sum(y_train==x)/length(y_train))
  n.col <- length(mat)
  n.row <- nrow(newx)
  pred <- matrix(0,n.row,n.col)
  for(i in 1:n.col){
    beta <- mat[[i]]
    nz <- sum(beta[,1] != 0)
    if(is.null(beta) || nz == 0){
      pred[,i] <- which.max(prior)
    }else{
      subset <- svd(beta)$u[,1,drop = FALSE]    # since we fix rank at 1
      pred[,i] <- lda_pred(x_train,y_train,subset,newx)
    }
  }
  return(pred)
}

#############################################
#             Simulate data             #
#############################################

# ##############    Situation 1    ############
# nz <- 6
# Sigma <- AR(0.5, p)
# Theta <- matrix(0, p, K)
# for (i in 1:3){
#   Theta[c(2*i-1,2*i),i] <- 1.6
# }
# for (i in 4:K){
#   # random <- runif(2, -1, 1)
#   # Theta[,i] <- ((i-3)+random[1])*(Theta[,2] - Theta[,1]) + ((i-3)+random[2])*(Theta[,3] - Theta[,1]) + Theta[,1]
#   Theta[,i] <- (1+(i-4)/2)*(Theta[,2] - Theta[,1]) + (1+(i-4)/2)*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# #############################################
# 
##############    Situation 1.1   for K = 3 ############
nz <- 6
Sigma <- AR(0.5, p)
Theta <- matrix(0, p, K)
for (i in 1:2){
  Theta[((i-1)*nz/2+1):(i*nz/2),i] <- 1.6
}
for (i in 3:K){
  Theta[,i] <- (i/2)*(Theta[,2] - Theta[,1]) + Theta[,1]
}

Mu <- Sigma%*%Theta

Beta <- matrix(0, p, K-1)
for(i in 1:(K-1)){
  Beta[,i] <- Theta[,i+1] - Theta[,1]
}
#############################################

# ##############    Situation 1.2   for K = 10   ############
# nz <- 20
# Sigma <- AR(0.5, p)
# Theta <- matrix(0, p, K)
# for (i in 1:2){
#   Theta[c(3*i-2,3*i-1,3),i] <- 1.6
# }
# for (i in 3:K){
#   Theta[,i] <- (1+2*(i-2)/3)*(Theta[,2] - Theta[,1]) + Theta[,1]
# }
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# #############################################


# ###############    Situation 2    ############
# nz <- 6
# Sigma <- kronecker(diag(5), CS(0.5, p/5))
# Theta <- matrix(0, p, K)
# for (i in 1:3){
#   Theta[c(2*i-1,2*i),i] <- 2.5
# }
# for (i in 4:K){
#   random <- sample(1:10, 2)
#   Theta[,i] <- random[1]*(Theta[,2] - Theta[,1]) + random[2]*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# #############################################

# ###############    Situation 3    ############
# nz <- 8   # The number of non-variables
# Sigma <- CS(0.5, p)
# Theta <- matrix(0, p, K)
# 
# for (i in 1:3){
#   random <- runif(nz, -1/8, 1/8)
#   Theta[1:K,i] <- i + random
# }
# 
# for (i in 4:K){
#   random <- runif(2, -1/8, 1/8)
#   Theta[,i] <- random[1]*(Theta[,2] - Theta[,1]) + random[2]*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# #############################################

# ###############    Situation 4    ############
# nz <- 8   # The number of non-variables
# Sigma <- CS(0.8, p)
# Theta <- matrix(0, p, K)
# 
# for (i in 1:3){
#   random <- runif(nz, -1/8, 1/8)
#   Theta[1:K,i] <- i + random
# }
# 
# for (i in 4:K){
#   random <- runif(2, -1/8, 1/8)
#   Theta[,i] <- random[1]*(Theta[,2] - Theta[,1]) + random[2]*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# #############################################

# ###############    Situation 5    ############
# nz <- 8   # The number of non-variables
# Sigma <- AR(0.5, p)
# Theta <- matrix(0, p, K)
# 
# Theta[,1] <- 0
# Theta[1:8,2] <- 1.2
# Theta[1:4,3] <- -1.2; Theta[5:8,3] <- 1.2
# 
# for (i in 4:K){
#   random <- runif(2, -1/8, 1/8)
#   Theta[,i] <- random[1]*(Theta[,2] - Theta[,1]) + random[2]*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# #############################################

# ###############    Situation 6    ############
# nz <- 8   # The number of non-variables
# Sigma <- AR(0.8, p)
# Theta <- matrix(0, p, K)
# 
# Theta[,1] <- 0
# Theta[1:8,2] <- 1.2
# Theta[1:4,3] <- -1.2; Theta[5:8,3] <- 1.2
# 
# for (i in 4:K){
#   random <- runif(2, -1/8, 1/8)
#   Theta[,i] <- random[1]*(Theta[,2] - Theta[,1]) + random[2]*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# #############################################

times <- 1
results <- matrix(0,times,14)
# used to store the execution time of each parts in ssdr function
ssdr_results <- matrix(0,times,7)
for(t in 1:times){

  # Create training, validation and testing dataset respectively
  start_time <- Sys.time()
  
  x_train <- Train(Nperclass, Mu, Sigma)
  y_train <- rep(1:K, each = Nperclass)

  x_val <- Train(Nperclass, Mu, Sigma)
  y_val <- rep(1:K, each = Nperclass)

  x_test <- Train(Nperclass_test, Mu, Sigma)
  y_test <- rep(1:K, each = Nperclass_test)
  
  end_time <- Sys.time()
  time_1 <- difftime(end_time, start_time, units = "secs")
  
  ##################################
  # Bayes error
  ##################################
  start_time <- Sys.time()
  
  B <- Beta
  pi <- rep(1/K,K)
  # Use true parameters to predict the testing data
  b_er <- function(x){
    tmp <- diag(t(replicate(K,x) - 1/2*Mu) %*% Theta) + log(pi)
    return(which.max(tmp))
  }
  pred_bayes <- apply(x_test, 1, b_er)
  e_bayes <- 1 - sum(pred_bayes == y_test)/length(y_test)
  
  end_time <- Sys.time()
  time_2 <- difftime(end_time, start_time, units = "secs")
  
  # the initial estimate
  start_time <- Sys.time()
  
  tmp <- msda.prep(x_train,y_train)
  sigma0 <- as.matrix(tmp$sigma)
  delta0 <- as.matrix(tmp$delta)
  mu0 <- as.matrix(tmp$mu)
  
  end_time <- Sys.time()
  time_3 <- difftime(end_time, start_time, units = "secs")
  ################################################
  # MSDA
  ################################################
  
  nlam_msda <- 10 # the number of lambdas in msda
  
  start_time <- Sys.time()

  fit_1 <- msda(x_train, y_train, nlambda = nlam_msda, maxit = 1e3)
  lam_msda <- fit_1$lambda
  pred_msda_val <- predict(fit_1, x_val)
  e_msda_val <- rep(0, ncol(pred_msda_val))

  end_time <- Sys.time()
  time_4 <- difftime(end_time, start_time, units = "secs")
  
  start_time <- Sys.time()
  
  for (i in 1:ncol(pred_msda_val)){
    pred <- pred_msda_val[,i]
    e_msda_val[i] <- 1 - sum(pred == y_val)/length(y_val)
  }
  
  id_min_msda <- which.min(e_msda_val)
  # calculate C and IC
  B_msda <- as.matrix(fit_1$theta[[id_min_msda]])
  tmp <- apply(B_msda, 1, function(x) any(x!=0))
  C_msda <- sum(which(tmp) %in% 1:nz)
  IC_msda <- sum(tmp) - C_msda
  #######################
  pred_msda <- predict(fit_1, x_test)[,id_min_msda]
  e_msda <- 1 - sum(pred_msda == y_test)/length(y_test)
  
  end_time <- Sys.time()
  time_5 <- difftime(end_time, start_time, units = "secs")
  ################################################
  # SSDR
  ################################################
  
  flmin <- as.double(1)
  nlam <- as.integer(1)
  nk <- as.integer(dim(delta0)[1])
  nobs <- as.integer(dim(x_train)[1])
  nvars <- as.integer(dim(x_train)[2])
  pf <- as.double(rep(1, nvars))
  dfmax <- as.integer(nobs)
  pmax <- as.integer(min(dfmax * 2 + 20, nvars))
  eps <- as.double(1e-04)
  maxit <- as.integer(1e+06)
  sml <- as.double(1e-06)
  verbose <- as.integer(FALSE)
  maxit_outer <- as.integer(1e+3) 
  eps_outer <- as.double(1e-4)
  vnames <- as.character(1:p)
  
  lam1 <- lam_msda
  lam2 <- seq(0.8,1.2,0.1)
  gamma <- c(10,20,30)
  n1 <- length(lam1)
  n2 <- length(lam2)
  n3 <- length(gamma)
  e_ssdr_val <- rep(0,n1*n2*n3)
  
  start_time <- Sys.time()
  fit_2 <- ssdr(lam1, lam2, gamma)
  pred_ssdr_val <- predict_ssdr(x_train, y_train, fit_2$Beta, x_val)
  time_6 <- fit_2$exe_time
  
  end_time <- Sys.time()
  time_7 <- difftime(end_time, start_time, units = "secs")
  
  # prediction error for validation set
  start_time <- Sys.time()
  
  for (i in 1:n1){
    for (j in 1:n2){
      for (k in 1:n3){
        pos <- (i-1)*n2*n3+(j-1)*n3+k
        pred <- pred_ssdr_val[,pos]
        e_ssdr_val[pos] <- 1 - sum(pred == y_val)/length(y_val)
      }
    }
  }
  
  id_min_ssdr <- which.min(e_ssdr_val)
  
  # calculate C and IC
  B_ssdr <- fit_2$Beta[[id_min_ssdr]]
  tmp <- apply(B_ssdr, 1, function(x) any(x!=0))
  C_ssdr <- sum(which(tmp) %in% 1:nz)
  IC_ssdr <- sum(tmp) - C_ssdr
  #######################
  
  pred_ssdr <- predict_ssdr(x_train, y_train, list(B_ssdr), x_test)
  e_ssdr <- 1 - sum(pred_ssdr == y_test)/length(y_test)
  
  end_time <- Sys.time()
  time_8 <- difftime(end_time, start_time, units = "secs")
  # store the prediction errors
  results[t,] <- c(C_msda, IC_msda, C_ssdr, IC_ssdr, e_bayes, e_msda, e_ssdr, time_1, time_2, time_3, time_4, time_5, time_7, time_8)
  ssdr_results[t,] <- time_6
}

# med_result <- apply(results,2, median)
# med_result <- as.data.frame(matrix(med_result, 1))
# colnames(med_result) <- c("C_msda", "IC_msda", "C_ssdr", "IC_ssdr", "error_bayes", "error_msda", "error_ssdr", "time_1", "time_2", "time_3", 
# "time_4", "time_5", "time_6", "time_7", "time_8")
results <- as.data.frame(results)
colnames(results) <- c("C_msda", "IC_msda", "C_ssdr", "IC_ssdr", "error_bayes", "error_msda", "error_ssdr", "time_1", "time_2", "time_3",
                       "time_4", "time_5", "time_7", "time_8")
ssdr_results <- as.data.frame(ssdr_results)
write.table(format(results, digits=4), "/Users/cengjing/Desktop/test")
write.table(format(ssdr_results, digits=4), "/Users/cengjing/Desktop/test")
