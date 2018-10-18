library(msda)
library(MASS)
library(methods)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")

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
predict_ssdr <- function(x_train, y_train, fit, newx){
  prior <- fit$prior
  mat <- fit$beta
  rank <- fit$rank
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
      subset <- svd(beta)$u[,1:r,drop = FALSE]
      pred[,i] <- lda_pred(x_train,y_train,subset,newx)
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

#########
# cv.msda
#########

# ssdr function returns corresponding B matrices and other evaluation stuff
ssdr <- function(x, y, lam1=NA, lam2=NA, gam=NA, grid=TRUE, lam_fac_ssdr=0.5, pf=rep(1, nvars), dfmax=nobs, 
                 pmax = min(dfmax * 2 + 20, nvars), eps=1e-04, maxit=1e+06, sml=1e-6, verbose=FALSE,
                 maxit_outer=1e+3, eps_outer=1e-3){
  
  y <- factor(y)
  x <- as.matrix(x)
  tmp <- msda.prep(x,y)
  prior <- tmp$prior
  sigma0 <- as.matrix(tmp$sigma)
  delta0 <- as.matrix(tmp$delta)
  mu0 <- as.matrix(tmp$mu)
  
  nlam <- as.integer(1)
  nk <- as.integer(dim(delta0)[1])
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  pf <- as.double(pf)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  eps <- as.double(eps)
  maxit <- as.integer(maxit)
  sml <- as.double(sml)
  verbose <- as.integer(verbose)
  maxit_outer <- as.integer(maxit_outer) 
  eps_outer <- as.double(eps_outer)
  vnames <- as.character(1:nvars)
  
  if(is.na(gam)){
    gamma <- 10
    n3 <- length(gamma)
    }
  
  if((is.na(lam1))||(is.na(lam2))){
    cvfit_msda <- cv.msda(x,y,nlambda = 10, maxit = 1e3, lambda.factor = 0.5)
    lam1_msda <- cvfit_msda$lambda.min
    fit_msda <- cvfit_msda$msda.fit
    id_msda <- which(fit_msda$lambda==lam1_msda)
    if(is.na(lam1)){
      lam1 <- (lam1_msda)*seq(1.5,0.6,-0.1)
      n1 <- length(lam1)
    }
    if(is.na(lam2)){
      B_msda <- as.matrix(fit_msda$theta[[id_msda]])
      n2 <- 10
      lam_fac_ssdr <- 0.5
      d <- svd(B_msda)$d
      lam2 <- d[1]*gamma*lam_fac_ssdr^seq((n2-1),0)
    }
    tmp1 <- lam1
    tmp2 <- lam2
    lam1 <- rep(tmp1,each=length(tmp2))
    lam2 <- rep(tmp2,time=length(tmp1))
  }else{
    if(grid){
      tmp1 <- lam1
      tmp2 <- lam2
      lam1 <- rep(tmp1,each=length(tmp2))
      lam2 <- rep(tmp2,time=length(tmp1))
    }else{
      if(length(lam1)!=length(lam2)){
        stop('The length of lambda1 should be equal to the length of lambda2')
      }
    }
  }
  
  # mat <- vector(mode = "list", length = n1*n2*n3)    # To store the converged matrix B of each run into a list
  # diff_B_final <- vector(mode = "list", length = n1*n2*n3)     # To store the difference of consecutive B sequence of each run into a list
  mat <- list()
  step_final <- c()     # To store the iteration times of each run
  time_final <- c()     # To store the running time of each run
  # jerr_list <- c()
  nlam_ssdr <- 0
  lam1_list <- c()
  lam2_list <- c()
  r_list <- c()
  uniq_lam1 <- unique(lam1)
  n1 <- length(uniq_lam1)
  
  for(i in 1:n1){
    ulam <- as.double(uniq_lam1[i])
    
    lam2_tmp <- lam2[lam1==ulam]
    n2 <- length(lam2_tmp)
    
    for(j in 1:n2){
      lambda2 <- lam2_tmp[j]
      
      # for(k in 1:n3){
        # gamma <- gam[k]
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
          fit <- .Fortran("msda", obj = double(1), nk, nvars, as.double(sigma), 
                          as.double(delta), pf, dfmax, pmax, nlam, flmin=as.double(1), ulam, 
                          eps, maxit, sml, verbose, nalam = integer(1), theta = double(pmax * nk * nlam), 
                          itheta = integer(pmax), ntheta = integer(nlam), 
                          alam = double(nlam), npass = integer(1), jerr = integer(1))
          
          # If jerr = 0, msda function returns normal results
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
            jerr <- 404
            break
          }
          
          Bold <- Bnew
          Cold <- Cnew
          muold <- munew
          
        }# End of repeat 
        end_time <- Sys.time()
        
        if(jerr < -10000) break
        
        lam1_list <- c(lam1_list, ulam)
        lam2_list <- c(lam2_list, lambda2)
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
        
    }# End of lambda2
    
  }# End of lambda1
  
  for(i in 1:length(mat)){
    r_list <- c(r_list, rank_func(mat[[i]], thrd = 1e-3))
  }
  
  return(list(beta = mat, rank = r_list, step = step_final, time_ssdr = time_final, nlam_ssdr = nlam_ssdr, 
              lam1_list = lam1_list, lam2_list = lam2_list, prior = prior))
  
}

# msda_best <- function(x_train, y_train, x_val, y_val){
#   nlam_msda <- 10 # the number of lambdas in msda
#   fit_1 <- msda(x_train, y_train, nlambda = nlam_msda, maxit = 1e3, lambda.factor = 0.5)
#   Beta_msda <- fit_1$theta
#   lam_msda <- fit_1$lambda
#   pred_msda_val <- predict(fit_1, x_val)
#   e_msda_val <- rep(0, ncol(pred_msda_val))
#   
#   for (i in 1:ncol(pred_msda_val)){
#     pred <- pred_msda_val[,i]
#     e_msda_val[i] <- 1 - sum(pred == y_val)/length(y_val)
#   }
#   
#   # The optimal lambda1
#   id_min_msda <- which.min(e_msda_val)
#   lam1_min_msda <- lam_msda[id_min_msda]
#   B_msda <- as.matrix(Beta_msda[[id_min_msda]])
# 
#   obj <- list(lam_msda=lam1_min_msda, mat_msda=B_msda)
#   return(obj)
# }



#################  cv.ssdr     ##############################

cv.ssdr <- function(x,y){
  ord <- order(y)
  x <- x[ord,]
  y <- y[ord]
  fit <- ssdr(x,y)  # fit the whole data first so we can get lam1 and lam2
  lam1 <- fit$lam1_list
  lam2 <- fit$lam2_list
  nfolds <- 6
  foldid <- rep(0,nrow(x))
  for (i in 1:21){
    foldid[((i-1)*6+1):(i*6)] <- sample(1:6)
  }
  cv_mat <- matrix(NA,nfolds,length(lam1))  # all the possible lam1 and lam2 candidates on column
  
  for (i in 1:nfolds){
    x_train <- x[foldid!=i,]
    y_train <- y[foldid!=i]
    x_val <- x[foldid==i,]
    y_val <- y[foldid==i]
   
    # fit all the possible models on training dataset and make predictions
    fit_fold <- ssdr(x_train, y_train, lam1 = lam1, lam2 = lam2, grid = FALSE)
    pred <- predict_ssdr(x_train, y_train, fit_fold, x_val)
    
    lam1_fold <- fit_fold$lam1_list
    lam2_fold <- fit_fold$lam2_list
    uniq_lam1 <- unique(lam1_fold)
    # put the prediction error on the corresponding position
    pos <- 0
    for(j in 1:length(uniq_lam1)){
      id_fold <- which(lam1_fold == uniq_lam1[j])
      lam2_tmp <- lam2_fold[id_fold]
      id_all <- which(lam1 == uniq_lam1[j])
      lam2_all <- lam2[id_all]
      cv_mat[i, pos+match(lam2_tmp,lam2_all)] <- colMeans(pred[,id_fold]!=y_val)
      pos <- pos + length(lam2_all)
    }
  }
  
  cv_mat[is.na(cv_mat)] <- 0
  cvm <- colMeans(cv_mat)
  cvsd <- sqrt(colMeans(scale(cv_mat, cvm, FALSE)^2)/(nfolds - 1))
  
  lambda1.min <- lam1[which.min(cvm)]
  lambda2.min <- lam2[which.min(cvm)]
  
  obj <- list(ssdr.fit = fit, lam1 = lam1, lam2 = lam2, lam1.min = lambda1.min, lam2.min = lambda2.min, 
              cvm=cvm, cvsd=cvsd)
}

  
############  Real data  ################
data <- read.table('/Users/cengjing/Documents/DIS/Record/Oct 12/DB_220_IDLH_2mins.txt')
data[,1] <- rep(1:21,each=7)
x <- as.matrix(data[,-1])
y <- data[,1]

# Divide the whole dataset into 7 folds
nfold <- 7
foldid <- rep(0,nrow(data))
for (i in 1:21){
  foldid[((i-1)*7+1):(i*7)] <- sample(1:7)
}

error_rate <- rep(0, nfold) #prediction error for each fold
for (i in 1:nfold){
  x_train <- x[foldid!=i,]
  y_train <- y[foldid!=i]
  x_test <- x[foldid==i,]
  y_test <- y[foldid==i]
  
  # choose the best model using cross-validation
  fit_cv <- cv.ssdr(x_train,y_train)
  lam1.min <- fit_cv$lam1.min
  lam2.min <- fit_cv$lam2.min
  id <- which((fit_cv$lam1 == lam1.min) & (fit_cv$lam2 == lam2.min))
  beta_final <- (fit_cv$ssdr.fit)$beta[[id]]
  r <- rank_func(beta_final, thrd = 1e-3)
  
  # predict on testing data with fitted model
  fit_pred <- list(beta=list(beta_final), rank=r)
  pred <- predict_ssdr(x_train, y_train, fit_pred, x_test)
  error_rate[i] = mean(pred != y_test)
}

##########################################################################################
#                                    Data structure                                      #
##########################################################################################

# ###############    Situation 10   rank = 3 ############
# nz <- 8   # The number of non-variables
# r <- 3
# Sigma <- AR(0.5, p)
# # Sigma <- CS(0.5, p)
# Theta <- matrix(0, p, K)
# 
# Theta[,1] <- 0
# Theta[1:8,2] <- 1
# Theta[1:4,3] <- -1; Theta[5:8,3] <- 1
# Theta[c(1,3,5,7), 4] <- -1; Theta[c(2,4,6,8), 4] <- 1
# 
# for (i in 5:K){
#   j <- i/2
#   u <- runif(3, j-0.5, j+0.5)
#   Theta[,i] <- u[1]*(Theta[,2] - Theta[,1]) + u[2]*(Theta[,3] - Theta[,1]) + u[3]*(Theta[,4] - Theta[,1]) + Theta[,1]
# }
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

##########################################################################################
#                                    Data structure                                      #
##########################################################################################

# set.seed(123)
# 
# p <- 800  #Dimension of observations
# K <- 21   # The number of class
# Nperclass <- 100  # The number of training observations in each class
# Nperclass_test <- 100   # The number of testing data in each class
# nz <- 8
# r <- 3
# 
# # Construct true Beta 
# Gamma <- rbind(matrix(runif(nz*r), nz, r), matrix(0, nrow = p-nz, ncol = r))
# orth_gamma <- qr.Q(qr(Gamma))
# eta <- matrix(runif((K-1)*r),(K-1),r)
# orth_eta <- qr.Q(qr(eta))
# Alpha <- diag(rep(20,r), r, r)
# Beta <- orth_gamma %*% Alpha %*% t(orth_eta)
# 
# 
# Theta <- matrix(0, p, K)
# Theta[,2:K] <- Beta
# Sigma <- AR(0.5,p)
# # Sigma <- CS(0.5,p)
# Mu <- Sigma%*%Theta
# 
# # sv_plot(svd(Beta)$d)
# 
# # #############################################
# 
# set.seed(Sys.time())
# 
# times <- 1 # Simulation times
# results <- matrix(0,times,18)
# 
# sv_msda_list <- c()
# sv_ssdr_list <- c()
# # nz_list <- c()
# 
# for(t in 1:times){
# 
#   # Generate training, validation and testing dataset respectively
#   x_train <- Train(Nperclass, Mu, Sigma)
#   y_train <- rep(1:K, each = Nperclass)
# 
#   x_val <- Train(Nperclass, Mu, Sigma)
#   y_val <- rep(1:K, each = Nperclass)
# 
#   x_test <- Train(Nperclass_test, Mu, Sigma)
#   y_test <- rep(1:K, each = Nperclass_test)
# 
#   
#   ##################################
#   # Bayes error
#   ##################################
# 
#   pi <- rep(1/K,K)
#   # Use true parameters to predict the testing data
#   b_er <- function(x){
#     tmp <- diag(t(replicate(K,x) - 1/2*Mu) %*% Theta) + log(pi)
#     return(which.max(tmp))
#   }
#   pred_bayes <- apply(x_test, 1, b_er)
#   e_bayes <- 1 - sum(pred_bayes == y_test)/length(y_test)
# 
#   #### The start of our methods
#   start_time <- Sys.time()
#   
#   # the initial estimate
#   
#   tmp <- msda.prep(x_train,y_train)
#   sigma0 <- as.matrix(tmp$sigma)
#   delta0 <- as.matrix(tmp$delta)
#   mu0 <- as.matrix(tmp$mu)
#   
#   ################################################
#   # MSDA
#   ################################################
#   
#   nlam_msda <- 10 # the number of lambdas in msda
# 
#   fit_1 <- msda(x_train, y_train, nlambda = nlam_msda, maxit = 1e3, lambda.factor = 0.5)
#   Beta_msda <- fit_1$theta
#   lam_msda <- fit_1$lambda
#   pred_msda_val <- predict(fit_1, x_val)
#   e_msda_val <- rep(0, ncol(pred_msda_val))
#   
#   for (i in 1:ncol(pred_msda_val)){
#     pred <- pred_msda_val[,i]
#     e_msda_val[i] <- 1 - sum(pred == y_val)/length(y_val)
#   }
#   
#   # The optimal lambda1
#   id_min_msda <- which.min(e_msda_val)
#   lam1_min_msda <- lam_msda[id_min_msda]
#   
#   # calculate C and IC
#   B_msda <- as.matrix(Beta_msda[[id_min_msda]])
#   tmp <- apply(B_msda, 1, function(x) any(x!=0))
#   C_msda <- sum(which(tmp) %in% 1:nz)
#   IC_msda <- sum(tmp) - C_msda
#   
#   ########### Save sv ##############
#   sv_msda_list <- rbind(sv_msda_list, svd(B_msda)$d)
#   #################################
#   
#   # rank of B_msda matrix, subspace distance and Frobenius distance
#   r_msda <- rank_func(B_msda, thrd = 1e-3)
#   sub_msda <- subspace(svd(Beta)$u[,1:r, drop=FALSE], svd(B_msda)$u[,1:r, drop=FALSE])
#   Fnorm_msda <- norm(Beta - B_msda, type = 'F')
#   
#   # Prediction error
#   pred_msda <- predict(fit_1, x_test)[,id_min_msda]
#   e_msda <- 1 - sum(pred_msda == y_test)/length(y_test)
#   
#   # Draw the singular values plot
#   # sv_plot(svd(B_msda)$d)
# 
#   ################################################
#   # SSDR
#   ################################################
#   
#   flmin <- as.double(1)
#   nlam <- as.integer(1)
#   nk <- as.integer(dim(delta0)[1])
#   nobs <- as.integer(dim(x_train)[1])
#   nvars <- as.integer(dim(x_train)[2])
#   pf <- as.double(rep(1, nvars))
#   # dfmax <- as.integer(nobs)
#   dfmax <- as.integer(nvars)
#   # pmax <- as.integer(min(dfmax * 2 + 20, nvars))
#   pmax <- as.integer(nvars)
#   eps <- as.double(1e-04)
#   maxit <- as.integer(1e+06)
#   sml <- as.double(1e-06)
#   verbose <- as.integer(FALSE)
#   maxit_outer <- as.integer(1e+3) 
#   eps_outer <- as.double(1e-3)
#   vnames <- as.character(1:p)
#   
#   lam_fac_ssdr <- 0.5
#   # We may need to shrink lam1 a little bit
#   lam1 <- (lam1_min_msda)*seq(0.6,1,0.1)
#   n1 <- length(lam1)
#   
#   gamma <- 10
#   # gamma <- c(10,20,30)
#   n3 <- length(gamma)
# 
#   # Construct lambda2 candidates
#   n2 <- 10   # we select n2 lambda2
#   d <- svd(B_msda)$d
#   lam2 <- d[1]*gamma*lam_fac_ssdr^seq(0,(n2-1))
#   
#   # if lam2 just contains one single value 0, then ssdr just degenerated to msda
#   if (all(lam2 == 0)){
#     C_ssdr <- C_msda
#     IC_ssdr <- IC_msda
#     e_ssdr <- e_msda
#     r_ssdr <- r_msda
#     sub_ssdr <- sub_msda
#     Fnorm_ssdr <- Fnorm_msda
#     lam1_min_ssdr <- lam1_min_msda
#     lam2_min_ssdr <- NA
#     gamma_min_ssdr <- NA
#     step <- NA
#     time_ssdr <- NA
#   }else{
#     
#     fit_2 <- ssdr(lam1, lam2, gamma)
#     # jerr <- rbind(jerr, fit_2$jerr)
#     Beta_ssdr <- fit_2$beta
#     # In some cases, all the Beta is null because the Fortran code didn't return a converaged B matrix 
#     if (sum(sapply(Beta_ssdr, is.null)) == n2*n3) {
#       results[t,] <- c(C_msda, IC_msda, NA, NA, e_bayes, e_msda, NA, r_msda, NA, sub_msda, 
#                        NA, Fnorm_msda, NA, NA, NA, NA, NA, NA)
#       next
#     }
#     
#     # for (i in 1:length(Beta_ssdr)){
#     #   B <- Beta_ssdr[[i]]
#     #   nz_list <- c(nz_list, sum(apply(B,1,function(x){any(x!=0)})))
#     # }
#     
#     
#     step <- fit_2$step
#     time_ssdr <- fit_2$time_ssdr
#     
#     # validation
#     pred_ssdr_val <- predict_ssdr(x_train, y_train, Beta_ssdr, x_val, r)
#     
#     # prediction error for validation set
#     
#     e_ssdr_val <- rep(0,n1*n2*n3)
#     
#     for (j in 1:ncol(pred_ssdr_val)){
#       pred <- pred_ssdr_val[,j]
#       e_ssdr_val[j] <- 1 - sum(pred == y_val)/length(y_val)
#     }
#     
#     # We find the optimal lambda1 and lambda2 here
#     id_min_ssdr <- which.min(e_ssdr_val)
#     id_lam1 <- ceiling(id_min_ssdr/(n2*n3))
#     id_lam2 <- ceiling((id_min_ssdr - (id_lam1-1)*n2*n3)/n3)
#     id_gamma <- id_min_ssdr - (id_lam1-1)*n2*n3 - (id_lam2-1)*n3
#     lam1_min_ssdr <- lam1[id_lam1]
#     lam2_min_ssdr <- lam2[id_lam2]
#     gamma_min_ssdr <- gamma[id_gamma]
#   
#     
#     B_ssdr <- Beta_ssdr[[id_min_ssdr]] 
#     
#     if(is.null(B_ssdr)){
#       C_ssdr <- NA
#       IC_ssdr <- NA
#       r_ssdr <- NA
#       sub_ssdr <- NA
#       Fnorm_ssdr <- NA
#     }else{
#     # Calculate C and IC
#       
#     ########### Save sv ##############
#     sv_ssdr_list <- rbind(sv_ssdr_list, svd(B_ssdr)$d)
#     #################################
#     
#     tmp <- apply(B_ssdr, 1, function(x) any(x!=0))
#     C_ssdr <- sum(which(tmp) %in% 1:nz)
#     IC_ssdr <- sum(tmp) - C_ssdr
#     r_ssdr <- rank_func(B_ssdr, thrd = 1e-3)
#     sub_ssdr <- subspace(svd(Beta)$u[,1:r, drop=FALSE], svd(B_ssdr)$u[,1:r, drop=FALSE])
#     Fnorm_ssdr <- norm(Beta - B_ssdr, type = 'F')
#     }
#     
#     pred_ssdr <- predict_ssdr(x_train, y_train, list(B_ssdr), x_test, r)
#     e_ssdr <- 1 - sum(pred_ssdr == y_test)/length(y_test)
#     
#     # Draw the singular values plot
#     # sv_plot(svd(B_ssdr)$d)
#   }
#   
#   #####################################################################
#   
#   end_time <- Sys.time()
#   time_total <- difftime(end_time, start_time, units = "secs")
#   # store the prediction errors
#   results[t,] <- c(C_msda, IC_msda, C_ssdr, IC_ssdr, e_bayes, e_msda, e_ssdr, r_msda, r_ssdr, sub_msda, 
#                    sub_ssdr, Fnorm_msda, Fnorm_ssdr, lam1_min_ssdr, lam2_min_ssdr, mean(step), 
#                    mean(time_ssdr), time_total)
# }
# 
# results <- as.data.frame(results)
# colnames(results) <- c("C_msda", "IC_msda", "C_ssdr", "IC_ssdr", "error_bayes", "error_msda", "error_ssdr",
#                        "r_msda", "r_ssdr","sub_msda", "sub_ssdr","Fnorm_msda","Fnorm_ssdr","lam1_min_ssdr", 
#                        "lam2_min_ssdr", "step", "time_ssdr", "time_total")
# write.table(results, "/Users/cengjing/Desktop/test_ssdr_1")
# # write.table(sv_msda_list, file = )
# # write.table(sv_ssdr_list, file = )