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
  u <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*u))
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

# ssdr function returns corresponding B matrices and other evaluation stuff
ssdr <- function(lam1,lam2,gam){
  n1 <- length(lam1)
  n2 <- length(lam2)
  n3 <- length(gam)
  # mat <- vector(mode = "list", length = n1*n2*n3)    # To store the converged matrix B of each run into a list
  mat <- list()
  # diff_B_final <- vector(mode = "list", length = n1*n2*n3)     # To store the difference of consecutive B sequence of each run into a list
  step_final <- c()     # To store the iteration times of each run
  time_final <- c()     # To store the running time of each run
  # jerr_list <- c()
  nlam_ssdr <- 0
  lam1_list <- c()
  lam2_list <- c()
  r_list <- c()
  
  for(i in 1:n1){
    ulam <- as.double(lam1[i])
    
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
  
      }# End of Gamma
      # If exit because of non-sparsity from msda, we stop trying more lam2s.
      if(jerr < -10000) break
      lam1_list <- c(lam1_list, ulam)
      lam2_list <- c(lam2_list, lambda2)
    }# End of lambda2
    
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
              lam1_list = lam1_list, lam2_list = lam2_list))
  
}

######################## AIC/BIC ########################

eval <- function(Beta, n, type, d=NULL){
  Beta <- as.matrix(Beta)
  nz <- sum(apply(Beta, 1, function(x) any(x!=0)))
  sample_size <- n
  p <- dim(Beta)[2]
  obj <- 0.5 * sum(diag(t(Beta) %*% sigma0 %*% Beta - 2 * t(mu0) %*% Beta), na.rm = TRUE)
  if(is.null(d)){
    df <- nz * p
  }else{
    df <- nz * d
  }
  # Different information criterion
  if(type=='AIC'){
    IC <- 2/sample_size * df
  }
  if(type=='BIC'){
    IC <- log(sample_size)/sample_size * df
  }
  if(type=='EBIC'){
    IC <- sqrt(log(sample_size)) * log(sample_size)/sample_size * df
  }
  
  e <- obj + IC
  return(c(obj, IC, e))
}


##########################################################################################
#                                    Data structure                                      #
##########################################################################################

# #############  Model 1 #############
# set.seed(123)
# 
# p <- 800  # Dimension of observations
# N <- 200  # Sample size
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
#   y <- sign(x %*% Beta[,1]) * log(abs(x %*% Beta[,2] + 5)) + 0.2 * rnorm(nobs)
#   return(y)
# }

# #############  Model 2 #############
# set.seed(123)
# 
# p <- 800  # Dimension of observations
# N <- 200  # Sample size
# H <- 10
# 
# Mu <- rep(0,p)
# # Sigma <- CS_blk(0.5,800,4)
# Sigma <- AR(0.5, p)
# # Construct true Beta
# Beta <- matrix(0, p, 1)
# Beta[1:20,1] <- 1
# nz_vec <- 1:20
# r <- 1
# 
# model <- function(x, Beta){
#   nobs <- dim(x)[1]
#   y <- x %*% Beta + 0.5 * rnorm(nobs)
#   return(y)
# }

#############  Model 3 #############

p <- 800  # Dimension of observations
N <- 200 # Sample size
H <- 10

Mu <- rep(0,p)
# Sigma <- CS_blk(0.5,800,4)
Sigma <- AR(0.5, p)
# Construct true Beta
Beta <- matrix(0, p, 2)
Beta[1:10,1] <- 1
Beta[11:20,2] <- 1
# Beta[1:10,2] <- c(1,-1,1,-1,1,-1,1,-1,1,-1)
nz_vec <- 1:20
# nz_vec <- 1:10
r <- 2

model <- function(x, Beta){
  nobs <- dim(x)[1]
  y <- (x %*% Beta[,1])/(0.5+(x %*% Beta[,2] + 1.5)^2) + 0.2 * rnorm(nobs)
  return(y)
}

# #############################################

times <- 100 # Simulation times
results <- matrix(0,times,13)

nlam_ssdr <- c()

sv_msda_list <- c()
sv_ssdr_list <- c()

for(t in 1:times){

  # Generate training, validation and testing dataset respectively
  
  x_train <- Train(N, Mu, Sigma)
  y_train <- model(x_train, Beta)
  # Slice y
  y_train <- cut(y_train, breaks = quantile(y_train, probs=seq(0,1, by=1/H), na.rm=TRUE), 
                     include.lowest = TRUE, labels = FALSE)
  
  # x_val <- Train(N, Mu, Sigma)
  # y_val <- sign(x_val %*% Beta[,1]) * log(abs(x_val %*% Beta[,2] + 5)) + 0.2 * rnorm(N)
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
  lam_msda <- fit_1$lambda
  
  # Calculate AIC for each lambda1
  obj_msda <- rep(0,length(Beta_msda))
  IC_msda <- rep(0,length(Beta_msda))
  eval_msda <- rep(0,length(Beta_msda))
  for (i in 1:length(Beta_msda)){
    tmp <- eval(Beta_msda[[i]], n = N, type = 'EBIC')
    obj_msda[i] <- tmp[1]
    IC_msda[i] <- tmp[2]
    eval_msda[i] <- tmp[3]
  }
  
  id_min_msda <- which.min(eval_msda)
  lam1_min_msda <- lam_msda[id_min_msda]
  
  # pred_msda_val <- predict(fit_1, x_val)
  # e_msda_val <- rep(0, ncol(pred_msda_val))
  # 
  # for (i in 1:ncol(pred_msda_val)){
  #   pred <- pred_msda_val[,i]
  #   e_msda_val[i] <- 1 - sum(pred == y_val)/length(y_val)
  # }
  
  # The optimal lambda1
  # id_min_msda <- which.min(e_msda_val)
  # lam1_min_msda <- lam_msda[id_min_msda]
  
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
  sv_msda_list <- rbind(sv_msda_list, svd(B_msda)$d)
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
  
  gamma <- 10
  # gamma <- c(10,20,30)
  n3 <- length(gamma)

  # Construct lambda2 candidates
  n2 <- 10   # we select n2 lambda2
  d <- svd(B_msda)$d
  lam2 <- d[1]*gamma*lam_fac_ssdr^seq((n2-1),0)
  
  # if lam2 just contains one single value 0, then ssdr just degenerated to msda
  if (all(lam2 == 0)){
    C_ssdr <- C_msda
    IC_ssdr <- IC_msda
    # e_ssdr <- e_msda
    r_ssdr <- r_msda
    # sub_ssdr <- sub_msda
    # Fnorm_ssdr <- Fnorm_msda
    lam1_min_ssdr <- lam1_min_msda
    lam2_min_ssdr <- NA
    ###
    id_lam1 <- 1
    id_lam2 <- NA
    ###
    gamma_min_ssdr <- NA
    step <- NA
    time_ssdr <- NA
  }else{
    fit_2 <- ssdr(lam1, lam2, gamma)
    Beta_ssdr <- fit_2$beta
    # In some cases, all the Beta is null because the Fortran code didn't return a converaged B matrix 
    if (sum(sapply(Beta_ssdr, is.null)) == n2*n3) {
      results[t,] <- c(C_msda, IC_msda, NA, NA, r_msda, NA, NA, NA, NA, NA, NA, NA, NA)
      next
    }
    
    # nz_list <- c()
    # for (i in 1:length(Beta_ssdr)){
    #   B <- Beta_ssdr[[i]]
    #   if(is.null(B)){
    #     nz_list <- c(nz_list, NA)
    #   }else{
    #     nz_list <- c(nz_list, sum(apply(B,1,function(x){any(x!=0)})))
    #   }
    # }
    
    lam1_list <- fit_2$lam1_list
    lam2_list <- fit_2$lam2_list
    rank_list <- fit_2$rank
    nlam_ssdr <- c(nlam_ssdr, fit_2$nlam_ssdr)
    step <- fit_2$step
    time_ssdr <- fit_2$time_ssdr
    
    # # validation
    # pred_ssdr_val <- predict_ssdr(x_train, y_train, fit_2, x_val)
    
    # # prediction error for validation set
    # e_ssdr_val <- rep(0, ncol(pred_ssdr_val))
    
    # for (j in 1:ncol(pred_ssdr_val)){
    #   pred <- pred_ssdr_val[,j]
    #   e_ssdr_val[j] <- 1 - sum(pred == y_val)/length(y_val)
    # }
    
    # The optimal lambda1 and lambda2 
    #########################
    
    obj_ssdr <- rep(0,length(Beta_ssdr))
    IC_ssdr <- rep(0,length(Beta_ssdr))
    eval_ssdr <- rep(0,length(Beta_ssdr))
    for (i in 1:length(Beta_ssdr)){
      mat <- Beta_ssdr[[i]]
      if (is.null(mat)){
        obj_ssdr[i] <- NULL
        IC_ssdr[i] <- NULL
        eval_ssdr[i] <- NULL
      }else{
        tmp <- eval(mat, n = N, d = rank_list[i], type = 'EBIC')
        obj_ssdr[i] <- tmp[1]
        IC_ssdr[i] <- tmp[2]
        eval_ssdr[i] <- tmp[3]
      }
    }

    id_min_ssdr <- which.min(eval_ssdr)
    lam1_min_ssdr <- lam1_list[id_min_ssdr]
    lam2_min_ssdr <- lam2_list[id_min_ssdr]
    gamma_min_ssdr <- 10
    id_lam1 <- which(lam1_min_ssdr == lam1)
    id_lam2 <- which(lam2_min_ssdr == lam2)
    id_gamma <- 1
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
      # if(r_ssdr==0){
      #   sub_ssdr <- NA
      # }else{
      #   sub_ssdr <- subspace(svd(Beta)$u[,1:r, drop=FALSE], svd(B_ssdr)$u[,1:r_ssdr, drop=FALSE])
      # }
      
      ########### Save sv ##############
      sv_ssdr_list <- rbind(sv_ssdr_list, svd(B_ssdr)$d)
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
  results[t,] <- c(C_msda, IC_msda, C_ssdr, IC_ssdr, r_msda, r_ssdr, lam1_min_ssdr, lam2_min_ssdr, id_lam1, id_lam2, mean(step), 
                   mean(time_ssdr), time_total)
}

results <- as.data.frame(results)
colnames(results) <- c("C_msda", "IC_msda", "C_ssdr", "IC_ssdr", "r_msda", "r_ssdr", "lam1_min_ssdr", 
                       "lam2_min_ssdr","id1", "id2", "step", "time_ssdr", "time_total")
write.table(results, "/Users/cengjing/Desktop/test_ssdr_1")
# write.table(sv_msda_list, file = )
# write.table(sv_ssdr_list, file = )