library(msda)
library(MASS)
library(methods)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")

p <- 800  #Dimension of observations
Nperclass <- 75  # The number of training observations in each class
K <- 3    # The number of class
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

# ssdr function returns corresponding B matrices, training dataset and prior
ssdr <- function(x,y,lam1,lam2,gam){
  n1 <- length(lam1)
  n2 <- length(lam2)
  n3 <- length(gam)
  prior <- sapply(1:K, function(x) sum(y==x)/length(y))
  mat <- vector(mode = "list", length = n1*n2*n3)
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
        
        }# End of repeat
        # If jerr == 1, then procedure converges. And if not, we leave the matrix NULL.
        if(jerr==1){
          mat[[(i-1)*n2*n3+(j-1)*n3+k]] <- Bnew
        }

      }
    }
  }
  outlist <- list(x = x, y = y, Beta = mat, prior = prior)
  return(outlist)
}

# input the ssdr object, returns corresponding predictions
predict_ssdr <- function(obj, newx){
  x_train <- obj$x
  y_train <- obj$y
  mat <- obj$Beta
  prior <- obj$prior
  n.col <- length(mat)
  n.row <- nrow(newx)
  pred <- matrix(0,n.row,n.col)
  for(i in 1:n.col){
    beta <- mat[[i]]
    nz <- sum(beta[,1] != 0)
    if(is.null(beta) || nz == 0){
      pred[,i] <- which.max(prior)
    }else{
      subset <- svd(beta)$u[,1,drop = FALSE]
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
##############    Situation 1.1   for K = 3,4,5 ############
nz <- 6
Sigma <- AR(0.5, p)
Theta <- matrix(0, p, K)
for (i in 1:2){
  Theta[c(3*i-2,3*i-1,3*i),i] <- 1.6
}
# for (i in 3:K){
#   Theta[,i] <- (1+(i-2)/2)*(Theta[,2] - Theta[,1]) + Theta[,1]
# }
Theta[,3] <- 1.5*(Theta[,2] - Theta[,1]) + Theta[,1]
Mu <- Sigma%*%Theta

Beta <- matrix(0, p, K-1)
for(i in 1:(K-1)){
  Beta[,i] <- Theta[,i+1] - Theta[,1]
}
#############################################

# ##############    Situation 1.2   for K = 10   ############
# nz <- 6
# Sigma <- AR(0.5, p)
# Theta <- matrix(0, p, K)
# for (i in 1:2){
#   Theta[c(3*i-2,3*i-1,3*i),i] <- 1.6
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

# Create training, validation and testing dataset respectively
  
x_train <- Train(Nperclass, Mu, Sigma)
y_train <- rep(1:K, each = Nperclass)

x_val <- Train(Nperclass, Mu, Sigma)
y_val <- rep(1:K, each = Nperclass)

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


tmp <- msda.prep(x_train,y_train)
sigma0 <- as.matrix(tmp$sigma)
delta0 <- as.matrix(tmp$delta)
mu0 <- as.matrix(tmp$mu)

################################################
# MSDA
################################################
nlam_msda <- 10 # the number of lambdas in msda
e_msda_val <- rep(0, nlam_msda)
fit_1 <- msda(x_train, y_train, nlambda = nlam_msda)
lam_msda <- fit_1$lambda
pred_msda_val <- predict(fit_1, x_val)

for (i in 1:nlam_msda){
  pred <- pred_msda_val[,i]
  e_msda_val[i] <- 1 - sum(pred == y_val)/length(y_val)
}

id_min_msda <- which.min(e_msda_val)
pred_msda <- predict(fit_1, x_test)[,id_min_msda]
e_msda <- 1 - sum(pred_msda == y_test)/length(y_test)

################################################
# SSDR
################################################

flmin <- as.double(1)
# ulam <- as.double(lambda1)
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
# sigma <- sigma0 + gamma*diag(rep(1,ncol(sigma0)), ncol(sigma0),ncol(sigma0))
  
lam1 <- lam_msda
lam2 <- seq(0.8,1.3,0.1)
gamma <- c(10,20,30)
n1 <- length(lam1)
n2 <- length(lam2)
n3 <- length(gamma)
e_ssdr_val <- rep(0,n1*n2*n3)
fit_2 <- ssdr(x_train, y_train, lam1, lam2, gamma)
pred_ssdr_val <- predict_ssdr(fit_2, x_val)

# prediction error for validation set
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
pred_ssdr <- predict_ssdr(fit_2, x_test)[,id_min_ssdr]
e_ssdr <- 1 - sum(pred_ssdr == y_test)/length(y_test)