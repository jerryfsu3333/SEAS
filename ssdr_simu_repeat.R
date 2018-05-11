library(msda)
library(MASS)
library(methods)
source("/Users/cengjing/Documents/DIS/Record/Code/msda_prep.R")
source("/Users/cengjing/Documents/DIS/Record/Code/utility.R")

p <- 800  #Dimension of observations
Nperclass <- 75  # The number of training observations in each class
K <- 3   # The number of class
Nperclass_test <- 150   # The number of testing data in each class


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

##############    Situation 1.1  K=3,4,5  ############
nz <- 6
Sigma <- AR(0.5, p)
Theta <- matrix(0, p, K)
for (i in 1:2){
  Theta[c(3*i-2,3*i-1,3*i),i] <- 1.6
}
#for (i in 3:K){
#  Theta[,i] <- (1+(i-2)/2)*(Theta[,2] - Theta[,1]) + Theta[,1]
#
#}
Theta[,3] <- 1.1*(Theta[,2] - Theta[,1]) + Theta[,1]
Mu <- Sigma%*%Theta

Beta <- matrix(0, p, K-1)
for(i in 1:(K-1)){
  Beta[,i] <- Theta[,i+1] - Theta[,1]
}
############################################# 

###############    Situation 1.2 K=10   ############
#nz <- 6
#Sigma <- AR(0.5, p) 
#Theta <- matrix(0, p, K) 
#for (i in 1:2){ 
#  Theta[c(3*i-2,3*i-1,3*i),i] <- 1.6
#}
#for (i in 3:K){
#  Theta[,i] <- (1+2*(i-2)/3)*(Theta[,2] - Theta[,1]) + Theta[,1]
#}
#Mu <- Sigma%*%Theta
#
#Beta <- matrix(0, p, K-1)
#for(i in 1:(K-1)){
#  Beta[,i] <- Theta[,i+1] - Theta[,1]
#}
#############################################

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

lam2 <- seq(0.5,1,0.1)
times <- 200
#The array is used for storing C and IC
array_tmp <- array(0, c(times, 14, length(lam2)))


for (t in 1:times){
  
  x <- Train(Nperclass, Mu, Sigma)
  y <- rep(1:K, each = Nperclass)
  
  x_test <- Train(Nperclass_test, Mu, Sigma)
  y_test <- rep(1:K, each = Nperclass_test)
  
  
  ##################################
  # Bayes error
  ##################################
  B <- Beta
  pi <- rep(1/K,K)
  # Use true parameters to predict the testing data
  b_er <- function(x){
    tmp <- diag(t(replicate(K,x) - 1/2*Mu) %*% Theta) + log(pi)
    return(which.max(tmp))
  }
  pred_bayes <- apply(x_test, 1, b_er)
  error_bayes <- 1 - sum(pred_bayes == y_test)/length(y_test)
  
  tmp_bayes <- apply(B, 1, function(x) any(x!=0))
  C_bayes <- sum(which(tmp_bayes) %in% 1:nz)
  IC_bayes <- sum(tmp_bayes) - C_bayes
  
  
  tmp <- msda.prep(x,y)
  sigma0 <- as.matrix(tmp$sigma)
  delta0 <- as.matrix(tmp$delta)
  mu0 <- as.matrix(tmp$mu)
  
  flag <- 0

  # Try different tuning parameters
  for(lambda1 in 2){
    
    # Initialized some parameters
    flmin <- as.double(1)
    ulam <- as.double(lambda1)
    nlam <- as.integer(1)
    nk <- as.integer(dim(delta0)[1])
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
    pf <- as.double(rep(1, nvars))
    dfmax <- as.integer(nobs)
    pmax <- as.integer(min(dfmax * 2 + 20, nvars))
    eps <- as.double(1e-04)
    maxit <- as.integer(1e+06)
    sml <- as.double(1e-06)
    verbose <- as.integer(FALSE)
    
    ##################################
    # MSDA
    ##################################
    fit_1 <- msda(x,y, lambda = lambda1)
    step_msda <- fit_1$npasses
    msda_mat <- as.matrix(fit_1$theta[[1]])    # The B matrix in msda method
    
    # To obtain the C and IC in MSDA and SSDR methods
    tmp_msda <- apply(msda_mat, 1, function(x) any(x!=0))
    C_msda <- sum(which(tmp_msda) %in% 1:nz)      # Correctly discovered non-zero variables
    IC_msda <- sum(tmp_msda) - C_msda
    
    # We use the whole matrix to predict
    pred_msda <- lda_pred(x,y,msda_mat,x_test)
    error_msda <- 1 - sum(pred_msda == y_test)/length(y_test)     # Prediction error

    # We use left-singular vector to estimate the subspace
    subset_msda <- svd(msda_mat)$u[,1, drop=FALSE]      # Reduced rank subspace estimation
    sub_msda <- subspace(subset_msda, Beta[,1,drop=FALSE])     # Subspace distance
    
    # We use the first column to estimate the subspace
    subset_msda_col <- msda_mat[,1,drop=FALSE]
    sub_msda_col <- subspace(subset_msda_col, Beta[,1,drop=FALSE])     # Subspace distance

    for(lambda2 in lam2){
      for(gamma in 30){

        flag <- flag + 1
        # Maximal interation for outer loop
        maxit_outer <- as.integer(1e+3) 
        eps_outer <- as.double(1e-4)
        vnames <- as.character(1:p)
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
        repeat{
          
          step_ssdr <- step_ssdr + 1
          
          # Update B
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
          
          outlist <- formatoutput(fit, maxit, pmax, nvars, vnames, nk)
          Bnew <- as.matrix(outlist$theta[[1]])
          
          # Update C
          Btemp <- Bnew + 1/gamma * muold
          r <- svd(Btemp)
          U <- r$u
          V <- r$v
          D <- r$d
          lamtemp <- sapply(D, FUN = function(x) max(0, x-lambda2/gamma))
          Cnew <- U %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(V)
          
          # Update mu
          munew <- muold + gamma * (Bnew - Cnew)
         
          # Exit condition
          if(max(abs(Bnew - Bold)) < eps_outer){
            jerr <- 1
            break
          }
          if(step_ssdr > maxit_outer){
            jerr <- -2
            break
          }
        
          Bold <- Bnew
          Cold <- Cnew
          muold <- munew
           
        }
        

        tmp <- apply(Bnew, 1, function(x) any(x!=0))
        C <- sum(which(tmp) %in% 1:nz)      # Correctly discovered non-zero variables
        IC <- sum(tmp) - C
        
        # We use left-singular vector to estimate the subspace
        subset_ssdr <- svd(Bnew)$u[,1,drop = FALSE]
        pred_ssdr <- lda_pred(x,y,subset_ssdr,x_test)
        error_ssdr <- 1 - sum(pred_ssdr == y_test)/length(y_test)   # Prediction error
        sub_ssdr <- subspace(subset_ssdr, Beta[,1,drop=FALSE])     # Subspace distance
        
        # We use the first column to estimate the subspace
        subset_ssdr_col <- Bnew[,1,drop=FALSE]
        pred_ssdr_col <- lda_pred(x,y,subset_ssdr_col,x_test)
        error_ssdr_col <- 1 - sum(pred_ssdr_col == y_test)/length(y_test)   # Prediction error
        sub_ssdr_col <- subspace(subset_ssdr_col, Beta[,1,drop=FALSE]) 

        array_tmp[t,,flag] <- c(C_msda, IC_msda, C, IC, error_bayes, error_msda, error_ssdr, error_ssdr_col, sub_msda, sub_msda_col, sub_ssdr, sub_ssdr_col, step_msda, step_ssdr)

      }
    }
  }

}

results <- as.data.frame(t(colMeans(array_tmp)), row.names = as.character(lam2))
colnames(results) <- c("C_msda", "IC_msda", "C_ssdr", "IC_ssdr", "error_bayes", "error_msda", "error_ssdr", "error_ssdr_col", "sub_msda", "sub_msda_col", "sub_ssdr", "sub_ssdr_col", "step_msda", "step_ssdr")
 
write.table(format(results, digits=4), "/Users/cengjing/Desktop/test")
