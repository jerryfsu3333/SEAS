library(msda)
library(MASS)
library(methods)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")

p <- 800  #Dimension of observations
K <- 21   # The number of class
Nperclass <- 10  # The number of training observations in each class
Nperclass_test <- 500   # The number of testing data in each class


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
  n2 <- length(lam2)
  n3 <- length(gam)
  mat <- vector(mode = "list", length = n2*n3)
  diff_B_final <- vector(mode = "list", length = n2*n3)
  sv2_final <- vector(mode = "list", length = n2*n3)
  step_final <- rep(0, n2*n3)
  time_final <- rep(0, n2*n3)
  
  lambda1 <- lam1
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
      diff_B <- c()
      sv2 <- c()
      
      start_time <- Sys.time()
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
        
        diff_B <- c(diff_B, norm(Bnew-Bold, type = "F"))
        sv2 <- c(sv2, svd(Bnew)$d[2])
        
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
        mat[[(j-1)*n3+k]] <- Bnew
      }
      
      end_time <- Sys.time()
      diff_time <- difftime(end_time, start_time, units = "secs")
      
      diff_B_final[[(j-1)*n3+k]] <- diff_B
      sv2_final[[(j-1)*n3+k]] <- sv2
      step_final[(j-1)*n3+k] <- step_ssdr
      time_final[(j-1)*n3+k] <- diff_time
      
    }
  }
  
  return(list(Beta = mat, diff = diff_B_final, sv2 = sv2_final, step = step_final, time_ssdr = time_final))
  
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
      # If matrix is null or a zero matrix, then use prior to predict
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
results <- matrix(0,times,16)

# lam_min_all <- matrix(0, times, 8)


for(t in 1:times){

  # Create training, validation and testing dataset respectively
  
  start_time <- Sys.time()

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
  e_bayes <- 1 - sum(pred_bayes == y_test)/length(y_test)
  
  
  # the initial estimate
  
  tmp <- msda.prep(x_train,y_train)
  sigma0 <- as.matrix(tmp$sigma)
  delta0 <- as.matrix(tmp$delta)
  mu0 <- as.matrix(tmp$mu)
  
  ################################################
  # MSDA
  ################################################
  
  nlam_msda <- 10 # the number of lambdas in msda

  fit_1 <- msda(x_train, y_train, nlambda = nlam_msda, maxit = 1e3, lambda.factor = 0.5)
  mat_msda <- fit_1$theta
  lam_msda <- fit_1$lambda
  pred_msda_val <- predict(fit_1, x_val)
  e_msda_val <- rep(0, ncol(pred_msda_val))
  
  for (i in 1:ncol(pred_msda_val)){
    pred <- pred_msda_val[,i]
    e_msda_val[i] <- 1 - sum(pred == y_val)/length(y_val)
  }
  
  # The optimal lambda1
  id_min_msda <- which.min(e_msda_val)
  lam1_min_msda <- lam_msda[id_min_msda]
  # calculate C and IC
  B_msda <- as.matrix(fit_1$theta[[id_min_msda]])
  tmp <- apply(B_msda, 1, function(x) any(x!=0))
  C_msda <- sum(which(tmp) %in% 1:nz)
  IC_msda <- sum(tmp) - C_msda
  #######################
  pred_msda <- predict(fit_1, x_test)[,id_min_msda]
  e_msda <- 1 - sum(pred_msda == y_test)/length(y_test)

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
  eps_outer <- as.double(1e-2)
  vnames <- as.character(1:p)
  lam_fac_ssdr <- 0.5
  n2 <- 5
  
  
  # lam1 <- seq(lam1_min_msda-1, lam1_min_msda+1, by = 0.2)
  lam1 <- lam_msda
  gamma <- c(10,20,30)
  n1 <- length(lam1)
  # n2 <- 5
  n3 <- length(gamma)


  Beta_ssdr <- list()
  diff_B <- list()
  sv2 <- list()
  step <- c()
  new_lam1 <- c()
  new_lam2 <- c()
  time_ssdr <- c()
  # For each lambda1, we use the singular value of B as the candidate lambda2
  for (i in 1:n1){
    B <- as.matrix(mat_msda[[i]])
    d <- svd(B)$d
    # lam2 <- matrix(seq(d[length(d)], d[1], (d[1]-d[length(d)])/(n2-1)), nrow = 1)
    lam2 <- matrix(d[1]*lam_fac_ssdr^seq(0,(n2-1)), nrow = 1)
    # if lam2 just contains one single value 0, we drop this lam1
    if (all(lam2 == 0)) next

    fit_2 <- ssdr(lam1[i], lam2, gamma)

    # We store lambda1 and lambda2 here, since some of the lambda1 are dropped
    # because the singular values are zero
    new_lam1 <- c(new_lam1, lam1[i])
    new_lam2 <- rbind(new_lam2, lam2)

    Beta_ssdr <- c(Beta_ssdr, fit_2$Beta)
    diff_B <- c(diff_B, fit_2$diff)
    sv2 <- c(sv2, fit_2$sv2)
    step <- c(step, fit_2$step)
    time_ssdr <- c(time_ssdr, fit_2$time_ssdr)
  }

  # ############    Draw the plot for each tuning parameter  ############
  # 
  # for (i in 1:new_n1){
  #   for (j in 1:n2){
  #     for (k in 1:n3){
  #       pos <- (i-1)*n2*n3+(j-1)*n3+k
  #       tmp <- diff_B[[pos]]
  #       sv <- sv2[[pos]]
  #       tmp2 <- step[pos]
  #       plot(tmp, xlab = "iteration", ylab = "difference of B", main = paste("i=", i, ", j=", j,", k=", k, ", step=", tmp2),
  #            xlim = c(0,100), ylim = c(0,0.05))
  #       plot(sv, xlab = "iteration", ylab = "2nd singular value of B", main = paste("i=", i, ", j=", j,", k=", k, ", step=", tmp2))
  #     }
  #   }
  # }
  # 
  # #####################################################################
  
  new_n1 <- length(new_lam1)
  pred_ssdr_val <- predict_ssdr(x_train, y_train, Beta_ssdr, x_val)
  
  
  # prediction error for validation set
  
  e_ssdr_val <- rep(0,new_n1*n2*n3)
  for (i in 1:new_n1){
    for (j in 1:n2){
      for (k in 1:n3){
        pos <- (i-1)*n2*n3+(j-1)*n3+k
        pred <- pred_ssdr_val[,pos]
        e_ssdr_val[pos] <- 1 - sum(pred == y_val)/length(y_val)
      }
    }
  }
  
  # We find the optimal lambda1 and lambda2 here
  id_min_ssdr <- which.min(e_ssdr_val)
  id_lam1 <- ceiling(id_min_ssdr/(n2*n3))
  id_lam2 <- ceiling((id_min_ssdr-(id_lam1-1)*n2*n3)/n3)
  lam1_min_ssdr <- new_lam1[id_lam1]
  lam2_min_ssdr <- new_lam2[id_lam1, id_lam2]

  
  B_ssdr <- Beta_ssdr[[id_min_ssdr]]
  pred_ssdr <- predict_ssdr(x_train, y_train, list(B_ssdr), x_test)
  e_ssdr <- 1 - sum(pred_ssdr == y_test)/length(y_test)
  
  #####################################################################
  # store the information about lambdas in msda and ssdr
  if (!(lam1_min_msda %in% new_lam1)){
    tmp <- c(lam1_min_msda, lam1_min_ssdr, lam2_min_ssdr, id_lam1, id_lam2, 
                         NA, e_ssdr)
  }else{
    temp_id <- which(lam1_min_msda == new_lam1)
    temp_id2 <- which.min(e_ssdr_val[((temp_id-1)*n2*n3+1):((temp_id)*n2*n3)])
    B_ssdr_2 <- Beta_ssdr[[(temp_id-1)*n2*n3 + temp_id2]]
    
    # Calculate C and IC
    tmp <- apply(B_ssdr_2, 1, function(x) any(x!=0))
    C_ssdr <- sum(which(tmp) %in% 1:nz)
    IC_ssdr <- sum(tmp) - C_ssdr
    
    pred_ssdr_2 <- predict_ssdr(x_train, y_train, list(B_ssdr_2), x_test)
    e_ssdr_msda <- 1 - sum(pred_ssdr_2 == y_test)/length(y_test)
    
    tmp <- c(lam1_min_msda, lam1_min_ssdr, lam2_min_ssdr, id_lam1, id_lam2,
                         e_ssdr_msda, e_ssdr, mean(step), mean(time_ssdr))
  }
  
  #####################################################################
  
  end_time <- Sys.time()
  time_total <- difftime(end_time, start_time, units = "secs")
  # store the prediction errors
  results[t,] <- c(C_msda, IC_msda, C_ssdr, IC_ssdr, e_bayes, e_msda, tmp, time_total)
}

results <- as.data.frame(results)
colnames(results) <- c("C_msda", "IC_msda", "C_ssdr", "IC_ssdr", "error_bayes", "error_msda", 
                       "lam1_msda", "lam1_ssdr", "lam2_ssdr", "id_lam1_ssdr", "id_lam2_ssdr",
                       "error1", "error2", "step", "time_ssdr", "time_total")

write.table(results, "/Users/cengjing/Desktop/test_ssdr_1")