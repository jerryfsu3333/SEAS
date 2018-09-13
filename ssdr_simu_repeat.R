library(msda)
library(MASS)
library(methods)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")

p <- 800  #Dimension of observations
K <- 21   # The number of class
Nperclass <- 10  # The number of training observations in each class
Nperclass_test <- 10   # The number of testing data in each class

###################################################
# Functions
###################################################

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

# input the ssdr returned Beta matrix, returns corresponding predictions
predict_ssdr <- function(x_train, y_train, mat, newx, r){
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
      # How to decide the rank of estimated B???
      subset <- svd(beta)$u[,1:r,drop = FALSE]
      pred[,i] <- lda_pred(x_train,y_train,subset,newx)
    }
  }
  return(pred)
}

# rank function
# is this criteria reasonable???
# rank_func <- function(mat, thrd){
#   d <- svd(mat)$d
#   for (i in 1:(length(d)-1)){
#     if (d[i]/d[i+1] > thrd) break
#   }
#   return(i)
# }

rank_func <- function(mat, thrd){
  d <- svd(mat)$d
  for (i in 1:(length(d)-1)){
    if ((d[i]+0.001)/(d[i+1]+.001) > thrd) break
  }
  return(i)
}


# Draw the plot of the ratio of singular values

# sv_plot <- function(sv){
#   tmp <- rep(0,length(sv)-1)
#   for(i in 1:(length(sv)-1)){
#     tmp[i] <- sv[i]/sv[i+1]
#   }
#   plot(tmp, main = "Ratio of consecutive singular values", ylab = "ratio")
#   plot(sv, main = "Singular values", ylab = "singular values")
# }

sv_plot <- function(sv){
  tmp <- rep(0,length(sv)-1)
  for(i in 1:(length(sv)-1)){
    tmp[i] <- (sv[i] + 0.001)/(sv[i+1] + 0.001)
  }
  plot(tmp, main = "Ratio of consecutive singular values", ylab = "ratio")
  plot(sv, main = "Singular values", ylab = "singular values")
}

# ssdr function returns corresponding B matrices and other evaluation stuff
ssdr <- function(lam1,lam2,gam){
  # lam1 is selected as the optimal one from msda, so we just have one lambda 1 here
  n2 <- length(lam2)
  n3 <- length(gam)
  mat <- vector(mode = "list", length = n2*n3)    # To store the converged matrix B of each run into a list
  diff_B_final <- vector(mode = "list", length = n2*n3)     # To store the difference of consecutive B sequence of each run into a list
  step_final <- c()     # To store the iteration times of each run
  time_final <- c()     # To store the running time of each run
  
  lambda1 <- lam1
  ulam <- as.double(lambda1)
  jerr_list <- matrix(0, n2, n3)
  
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
        
        # Exit condition
        # the success code is 1
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
        diff_B_final[[(j-1)*n3+k]] <- diff_B 
        step_final <- c(step_final, step_ssdr)
      }
      
      end_time <- Sys.time()
      time_final <- c(time_final, difftime(end_time, start_time, units = "secs"))
      jerr_list[j,k] <- jerr
    
    }
  }
  
  return(list(Beta = mat, diff = diff_B_final, step = step_final, time_ssdr = time_final, jerr = jerr_list))
  
}



#############################################
#             Data structure                #
#############################################

# ##############    Situation 1.1    ############
# nz <- 8
# r <- 3
# Sigma <- AR(0.8, p)
# Theta <- matrix(0, p, K)
# for (i in 1:4){
#   Theta[((i-1)*nz/4+1):(i*nz/4),i] <- 1.3
# }
# for (i in 5:K){
#   j <- i/2
#   u <- runif(3, j-0.5, j+0.5)
#   Theta[,i] <- u[1]*(Theta[,2] - Theta[,1]) + u[2]*(Theta[,3] - Theta[,1]) +  u[3]*(Theta[,4] - Theta[,1]) + Theta[,1]
# }
# 
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# #############################################

# ##############    Situation 1.2     ############
# nz <- 6
# r <- 2
# Sigma <- AR(0.5, p)
# Theta <- matrix(0, p, K)
# for (i in 1:3){
#   Theta[((i-1)*nz/3+1):(i*nz/3),i] <- 1.6
# }
# for (i in 4:K){
#   j <- i/2
#   Theta[,i] <- j*(Theta[,2] - Theta[,1]) + j*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# 
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# 
# sv_plot(svd(Beta)$d)
# #############################################


# # ###############    Situation 2    ############
# nz <- 6
# r <- 2
# Sigma <- kronecker(diag(5), CS(0.5, p/5))
# Theta <- matrix(0, p, K)
# for (i in 1:3){
#   Theta[((i-1)*nz/3+1):(i*nz/3),i] <- 1.6
# }
# for (i in 4:K){
#   j <- i/2
#   u <- runif(2, j-0.5, j+0.5)
#   Theta[,i] <- u[1]*(Theta[,2] - Theta[,1]) + u[2]*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# 
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# # #############################################



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


# ###############    Situation 4  rank = 2  ############
# nz <- 6   # The number of non-variables
# r <- 2
# Sigma <- AR(0.3, p)
# # Sigma <- CS(0.8, p)
# Theta <- matrix(0, p, K)
# 
# Theta[,1] <- 0
# Theta[1:6,2] <- 1
# Theta[1:3,3] <- -1; Theta[4:6,3] <- 1
# 
# for (i in 4:K){
#   j <- i/2
#   u <- runif(2, j-0.5, j+0.5)
#   Theta[,i] <- u[1]*(Theta[,2] - Theta[,1]) + u[2]*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# 
# sv_plot(svd(Beta)$d)
# rank_func(Beta, 1e2)
# #############################################

# ###############    Situation 10   rank = 3 ############
nz <- 8   # The number of non-variables
r <- 3
# Sigma <- AR(0.5, p)
Sigma <- CS(0.5, p)
Theta <- matrix(0, p, K)

Theta[,1] <- 0
Theta[1:8,2] <- 1
Theta[1:4,3] <- -1; Theta[5:8,3] <- 1
Theta[c(1,3,5,7), 4] <- -1; Theta[c(2,4,6,8), 4] <- 1

for (i in 5:K){
  j <- i/2
  u <- runif(3, j-0.5, j+0.5)
  Theta[,i] <- u[1]*(Theta[,2] - Theta[,1]) + u[2]*(Theta[,3] - Theta[,1]) + u[3]*(Theta[,4] - Theta[,1]) + Theta[,1]
}
Mu <- Sigma%*%Theta

Beta <- matrix(0, p, K-1)
for(i in 1:(K-1)){
  Beta[,i] <- Theta[,i+1] - Theta[,1]
}

sv_plot(svd(Beta)$d)
rank_func(Beta, 1e3)
# #############################################

# ###############    Situation   rank = 2  ############
# nz <- 6   # The number of non-variables
# r <- 2
# Sigma <- CS(0.5, p)
# Theta <- matrix(0, p, K)
# 
# for (i in 1:3){
#   random <- runif(nz, -0.5, 0.5)
#   Theta[1:nz,i] <- i/2 + random
# }
# 
# for (i in 4:K){
#   j <- i/2
#   u <- runif(2, j-0.5, j+0.5)
#   Theta[,i] <- u[1]*(Theta[,2] - Theta[,1]) + u[2]*(Theta[,3] - Theta[,1]) + Theta[,1]
# }
# 
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# 
# sv_plot(svd(Beta)$d)
# #############################################

# ###############    Situation 11  rank = 3  ############
# nz <- 8   # The number of non-variables
# r <- 3
# # Sigma <- AR(0.5, p)
# Sigma <- CS(0.8, p)
# 
# Theta <- matrix(0, p, K)
# 
# for (i in 1:(r+1)){
#   random <- runif(nz, -0.5, 0.5)
#   Theta[1:nz,i] <- i/2 + random
# }
# 
# for (i in (r+2):K){
#   j <- (i-1)/3
#   u <- runif(r, j-0.5, j+0.5)
#   Theta[,i] <- (Theta[,2:(r+1)] - Theta[,1] %*% matrix(rep(1,r), nrow = 1)) %*% u + Theta[,1]
#   # Theta[,i] <- u[1]*(Theta[,2] - Theta[,1]) + u[2]*(Theta[,3] - Theta[,1]) + u[3]*(Theta[,4] - Theta[,1]) + Theta[,1]
# }
# 
# Mu <- Sigma%*%Theta
# 
# Beta <- matrix(0, p, K-1)
# for(i in 1:(K-1)){
#   Beta[,i] <- Theta[,i+1] - Theta[,1]
# }
# 
# sv_plot(svd(Beta)$d)
# rank_func(Beta, 1e3)
# #############################################

#############################################
# Test
# nz <- 6
# r <- 2
# Gamma <- rbind(matrix(runif(nz*r), nz, r), matrix(0, nrow = p-nz, ncol = r))
# b <- matrix(runif(r*(K-1)), r, K-1)
# orth_gamma <- qr.Q(qr(Gamma))
# Beta <- orth_gamma %*% b
# 
# Theta <- matrix(0, p, K)
# Theta[,2:K] <- Beta
# Sigma <- CS(0.1, p)
# 
# Mu <- Sigma%*%Theta
# 
# sv_plot(svd(Beta)$d)

# #############################################

times <- 1
results <- matrix(0,times,15)

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
  Beta_msda <- fit_1$theta
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
  B_msda <- as.matrix(Beta_msda[[id_min_msda]])
  tmp <- apply(B_msda, 1, function(x) any(x!=0))
  C_msda <- sum(which(tmp) %in% 1:nz)
  IC_msda <- sum(tmp) - C_msda
  
  # rank of B_msda matrix and distance between subspace
  r_msda <- rank_func(B_msda, thrd = 1e2)
  sub_msda <- subspace(B[,1:r, drop=FALSE], B_msda[,1:r,drop=FALSE])
  # test1 <- subspace(svd(B)$u[,1:r, drop=FALSE], svd(B_msda)$u[,1:r, drop=FALSE])
  #######################
  
  pred_msda <- predict(fit_1, x_test)[,id_min_msda]
  e_msda <- 1 - sum(pred_msda == y_test)/length(y_test)
  
  # Draw the singular values plot
  sv_plot(svd(B_msda)$d)

  ################################################
  # SSDR
  ################################################
  
  flmin <- as.double(1)
  nlam <- as.integer(1)
  nk <- as.integer(dim(delta0)[1])
  nobs <- as.integer(dim(x_train)[1])
  nvars <- as.integer(dim(x_train)[2])
  pf <- as.double(rep(1, nvars))
  # dfmax <- as.integer(nobs)
  dfmax <- as.integer(nvars)
  # pmax <- as.integer(min(dfmax * 2 + 20, nvars))
  pmax <- as.integer(nvars)
  eps <- as.double(1e-04)
  maxit <- as.integer(1e+06)
  sml <- as.double(1e-06)
  verbose <- as.integer(FALSE)
  maxit_outer <- as.integer(1e+3) 
  eps_outer <- as.double(1e-2)
  vnames <- as.character(1:p)
  lam_fac_ssdr <- 0.5
  
  lam1_min_ssdr <- lam1_min_msda
  n2 <- 5   # Choose n2 = 5, we select 5 lambda2
  gamma <- c(40,50,60)
  # gamma <- c(1,2,3)
  # gamma <- c(3,5,8)
  n3 <- length(gamma)

  # Construct lambda2 candidates
  temp <- as.matrix(Beta_msda[[id_min_msda]])
  d <- svd(temp)$d
  lam2 <- d[1]*lam_fac_ssdr^seq(0,(n2-1))
  
  # if lam2 just contains one single value 0, then ssdr just degenerated to msda
  if (all(lam2 == 0)){
    C_ssdr <- C_msda
    IC_ssdr <- IC_msda
    e_ssdr <- e_msda
    lam1_min_ssdr <- lam1_min_msda
    lam2_min_ssdr <- NA
    gamma_min_ssdr <- NA
    step <- NA
    time_ssdr <- NA
  }else{
    
    fit_2 <- ssdr(lam1_min_ssdr, lam2, gamma)
    jerr <- fit_2$jerr
    Beta_ssdr <- fit_2$Beta
    # In some cases, all the Beta is null because the Fortran code didn't return a converaged B matrix 
    if (sum(sapply(Beta_ssdr, is.null)) == n2*n3) {
      results[t,] <- c(C_msda, IC_msda, NA, NA, e_bayes, e_msda, NA, r_msda, NA, sub_msda, 
                       NA, NA, NA, NA, NA)
      next
    }
    
    step <- fit_2$step
    time_ssdr <- fit_2$time_ssdr
    
    # validation
    pred_ssdr_val <- predict_ssdr(x_train, y_train, Beta_ssdr, x_val, r)
    
    # prediction error for validation set
    
    e_ssdr_val <- rep(0,n2*n3)
    
    for (j in 1:n2){
      for (k in 1:n3){
        pos <- (j-1)*n3+k
        pred <- pred_ssdr_val[,pos]
        e_ssdr_val[pos] <- 1 - sum(pred == y_val)/length(y_val)
      }
    }
    
    # We find the optimal lambda1 and lambda2 here
    id_min_ssdr <- which.min(e_ssdr_val)
    # id_lam1 <- ceiling(id_min_ssdr/(n2*n3))
    id_lam2 <- ceiling(id_min_ssdr/n3)
    id_gamma <- id_min_ssdr - (id_lam2 - 1)*n3
    # lam1_min_ssdr <- lam1_min_msda
    lam2_min_ssdr <- lam2[id_lam2]
    gamma_min_ssdr <- gamma[id_gamma]
  
    
    B_ssdr <- Beta_ssdr[[id_min_ssdr]] 
    
    if(is.null(B_ssdr)){
      C_ssdr <- NA
      IC_ssdr <- NA
      r_ssdr <- NA
      sub_ssdr <- NA
    }else{
    # Calculate C and IC
    tmp <- apply(B_ssdr, 1, function(x) any(x!=0))
    C_ssdr <- sum(which(tmp) %in% 1:nz)
    IC_ssdr <- sum(tmp) - C_ssdr
    r_ssdr <- rank_func(B_ssdr, thrd = 1e3)
    sub_ssdr <- subspace(B[,1:r, drop=FALSE], B_ssdr[,1:r,drop=FALSE])
    # test2 <- subspace(svd(B)$u[,1:r, drop=FALSE], svd(B_ssdr)$u[,1:r, drop=FALSE])
    }
    
    pred_ssdr <- predict_ssdr(x_train, y_train, list(B_ssdr), x_test, r)
    e_ssdr <- 1 - sum(pred_ssdr == y_test)/length(y_test)
    
    # Draw the singular values plot
    sv_plot(svd(B_ssdr)$d)
  }
  
  #####################################################################
  
  end_time <- Sys.time()
  time_total <- difftime(end_time, start_time, units = "secs")
  # store the prediction errors
  results[t,] <- c(C_msda, IC_msda, C_ssdr, IC_ssdr, e_bayes, e_msda, e_ssdr, r_msda, r_ssdr, sub_msda, sub_ssdr, gamma_min_ssdr, mean(step), mean(time_ssdr), time_total)
}

results <- as.data.frame(results)
colnames(results) <- c("C_msda", "IC_msda", "C_ssdr", "IC_ssdr", "error_bayes", "error_msda", "error_ssdr",
                       "r_msda", "r_ssdr","sub_msda", "sub_ssdr","gamma_min_ssdr", "step", "time_ssdr", "time_total")


# write.table(results, "/Users/cengjing/Desktop/test_ssdr_1")