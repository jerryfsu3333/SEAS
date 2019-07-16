ssdr.cv <- function(x, y, H=5, type = 'sir', lambda.factor=0.5, 
                      lam_fac_msda=0.8, lam_fac_ssdr=0.8, nlam_msda=10, nlam1=10, nlam2=15, 
                      gamma=c(10,30,50), cut_y=TRUE, nfold = 5){
  
  y <- drop(y)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  # Cross validation with msda to find lambda1_msda
  fit_1 <- msda.cv(x, y, H = H, type = type, nlam = nlam_msda, lambda.factor = lambda.factor, cut_y=cut_y, nfold = nfold, maxit = 1e3)
  sigma0 <- as.matrix(fit_1$sigma)
  mu0 <- as.matrix(fit_1$mu)
  id_min_msda <- fit_1$id
  lam1_min_msda <- fit_1$lambda
  rank_min_msda <- fit_1$rank
  B_msda <- as.matrix(fit_1$Beta)
  
  # Generate tuning parameter candidates
  n1 <- nlam1
  lam_fac_msda <- lam_fac_msda
  lam1 <- lam1_min_msda*lam_fac_msda^seq(0,(n1-1))
  
  gamma <- gamma
  n3 <- length(gamma)
  
  n2 <- nlam2   # we select n2 lambda2 for each gamma
  lam_fac_ssdr <- lam_fac_ssdr
  d <- svd(B_msda)$d
  lam2 <- d[1] * matrix(gamma, ncol = 1) %*% matrix(lam_fac_ssdr^seq((n2-1),0), nrow = 1)
  
  # if lam2 just contains one single value 0, then ssdr just degenerated to msda
  if (all(lam2 == 0)){
    # cat("All lambda2 are zero, degenerate to msda\n")
    # return(list(mat = B_msda, rank = rank_min_msda, cvm = NA, cvsd = NA, id = NA, lam1 = lam1, lam2 = lam2, gamma = gamma,
    #             lam1.min = lam1_min_msda, lam2.min = NA, gamma.min = NA))
    cat("All lambda2 are zero, msda matrix is zero matrix\n")
    return(list(mat = NULL, rank = NA, cvm = NA, cvsd = NA, id = NA, lam1 = lam1, lam2 = lam2, gamma = gamma,
                lam1.min = lam1_min_msda, lam2.min = NA, gamma.min = NA))
  }else{
    # Cross-validation
    if (nfold < 3) stop("nfold must be larger than 3")
    if (nfold > nobs) stop("nfold is larger than the sample size")
    
    fold <- sample(rep(seq(nfold), length = nobs))

    eval_ssdr <- sapply(1:nfold, function(k){
      x_train <- x[which(fold!=k),,drop=FALSE]
      x_val <- x[which(fold==k),,drop=FALSE]
      y_train <- y[which(fold!=k)]
      y_val <- y[which(fold==k)]
      
      nobs_fold <- as.integer(dim(x_train)[1])
      nvars_fold <- as.integer(dim(x_train)[2])
      prep_out <- prep(x_train,y_train,type,H,cut_y)
      sigma0 <- prep_out$sigma
      mu0 <- prep_out$mu
      fit_fold <- ssdr(sigma0, mu0, nobs_fold, nvars_fold, lam1, lam2, gamma)
      
      Beta_fold <- fit_fold$beta
      if (all(sapply(Beta_fold, is.null))) {
        cat("Fold",k,":No converged matrix returned\n")
        return(rep(NA, length(Beta_fold)))
      }
      rank_fold <- fit_fold$rank
      Beta_fold <- cut_mat(Beta_fold, 1e-3, rank_fold)
      return(eval_val_rmse(Beta_fold, x_val, y_val))
    })
    
    # If no matrix is converged in any fold, return NULL matrix
    if(all(is.na(eval_ssdr))){
      cat("No converged matrix returned in the process of cross-validation\n")
      return(list(mat = NULL, rank = NA, cvm = NA, cvsd = NA, id = NA, lam1 = lam1, lam2 = lam2, gamma = gamma,
                         lam1.min = lam1_min_msda, lam2.min = NA, gamma.min = NA))
    }
    # Calculate cv mean and cv std
    cvm <- apply(eval_ssdr, 1, mean, na.rm=TRUE)
    cvsd <- sqrt(colMeans(scale(t(eval_ssdr), cvm, FALSE)^2, na.rm = TRUE)/(nfold-1))
    
    # Find the optimal lam1, lam2 and gamma
    id_min_ssdr <- which.min(cvm)
    id_lam1 <- ceiling(id_min_ssdr/(n2*n3))
    id_gamma <- ceiling((id_min_ssdr-(id_lam1-1)*(n2*n3))/n2)
    id_lam2 <- id_min_ssdr-(id_lam1-1)*(n2*n3)-(id_gamma-1)*n2
    lam1_min_ssdr <- lam1[id_lam1]
    gamma_min_ssdr <- gamma[id_gamma]
    lam2_min_ssdr <- lam2[id_gamma,id_lam2]
    
    # Refit with the optimal parameters
    prep_full <- prep(x,y,type,H,cut_y)
    sigma0 <- prep_full$sigma
    mu0 <- prep_full$mu
    fit_full <- ssdr(sigma0, mu0, nobs, nvars, lam1_min_ssdr, matrix(lam2_min_ssdr,1,1), gamma_min_ssdr)
    
    Beta_ssdr <- fit_full$beta
    rank_ssdr <- fit_full$rank
    Beta_ssdr <- cut_mat(Beta_ssdr, 1e-3, rank_ssdr)
    
    r_ssdr <- rank_ssdr[[1]]
    B_ssdr <- Beta_ssdr[[1]]
    
    id <- data.frame(id_lam1 = id_lam1, id_lam2 = id_lam2, id_gamma = id_gamma)
    
    if(is.null(B_ssdr)){
      cat("Optimal matrix is a null matrix\n")
      return(list(mat = NULL, rank = NA, cvm = cvm, cvsd = cvsd, id = id, lam1 = lam1, lam2 = lam2, gamma = gamma,
                  lam1.min = lam1_min_ssdr, lam2.min = lam2_min_ssdr, gamma.min = gamma_min_ssdr))
    }else{
      return(list(mat = B_ssdr, rank = r_ssdr, cvm = cvm, cvsd = cvsd, id = id, lam1 = lam1, lam2 = lam2, gamma = gamma,
                  lam1.min = lam1_min_ssdr, lam2.min = lam2_min_ssdr, gamma.min = gamma_min_ssdr))
    }
    
  }
}


msda.cv <- function(x,y,H,type,nlam,lambda.factor,cut_y=FALSE,nfold=5, maxit=1e3){
  
  # Fit full data, obtain the msda lambda candidates
  fit <- msda_func(x, y, H, type, nlam = nlam, lambda.factor = lambda.factor, cut_y = cut_y, maxit = maxit)
  sigma0 <- fit$sigma
  mu0 <- fit$mu
  lam_msda <- fit$lambda
  rank_msda <- fit$rank
  Beta_msda <- fit$Beta
  
  # Cross-validation
  fold <- sample(rep(1:nfold,ceiling(length(y)/nfold))[1:length(y)])
  eval_msda <- sapply(1:nfold, function(k){
    x_train <- x[which(fold!=k),]		
    x_val <- x[which(fold==k),]		
    y_train <- y[which(fold!=k)]
    y_val <- y[which(fold==k)]
    
    fit_fold <- msda_func(x_train, y_train, H, type, lambda = lam_msda, nlam = nlam, lambda.factor = lambda.factor,cut_y = cut_y, maxit = maxit)
    Beta_fold <- fit_fold$Beta
    # return evaluation of each fold
    eval_val_rmse(Beta_fold, x_val, y_val)
  })
  
  eval <- apply(eval_msda, 1, mean)
  
  # The optimal lambda1
  id_min_msda <- which.min(eval)
  lam1_min_msda <- lam_msda[id_min_msda]
  rank_min_msda <- rank_msda[id_min_msda]
  
  # calculate C, IC, Frobenious distance, rank and subspace distance
  B_msda <- as.matrix(Beta_msda[[id_min_msda]])
  
  list(id = id_min_msda, lambda = lam1_min_msda, rank = rank_min_msda, Beta = B_msda, sigma = sigma0, mu = mu0)
}


msda_func <-function(x,y,H,type,nlam,lambda.factor,lambda = NULL, cut_y=FALSE, maxit=1e3){
  
  fit <- my_msda(x, y, H=H, type = type, lambda=lambda, nlambda=nlam, maxit=maxit, lambda.factor=lambda.factor, cut_y=cut_y)
  sigma0 <- as.matrix(fit$sigma)
  mu0 <- as.matrix(fit$mu)
  lam_msda <- fit$lambda
  Beta_msda <- fit$theta
  
  rank_msda <- sapply(seq_len(length(Beta_msda)), function(i){
    mat <- Beta_msda[[i]]
    if(is.null(mat)){
      NA
    }else{
      rank_func(mat, thrd = 1e-3)
    }
  })
  
  Beta_msda <- cut_mat(Beta_msda, 1e-3, rank_msda)
  list(lambda = lam_msda, Beta = Beta_msda, sigma = sigma0, mu = mu0, rank = rank_msda)
}


# ssdr algorithm function

########################################  Error code  ##########################################
# case1: jerr < -10000, exit because of non-sparsity from msda, we stop trying more lam2s.
# case2: jerr = 404, exit because we reach the maximum iteration time of ssdr. And we leave matrix NULL
# case3: jerr = 1, succeed.
###############################################################################################

ssdr <- function(sigma, mu, nobs, nvars, lam1, lam2, gam, pf=rep(1, nvars), dfmax=nobs, pmax=min(dfmax * 2 + 20, nvars), eps=1e-04, maxit=1e+06, sml=1e-06, verbose = FALSE, maxit_outer=1e+3, eps_outer=1e-3){
  
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
  
  mat <- vector("list", nparams)
  step_final <- vector("list", nparams)     # To store the iteration times of each run
  time_final <- vector("list", nparams)     # To store the running time of each run
  lam1_list <- vector("list", nparams)
  lam2_list <- vector("list", nparams)
  gamma_list <- vector("list", nparams)
  r_list <- vector("list", nparams)
  
  sv_list_B <- vector("list", nparams)
  sv_list_C <- vector("list", nparams)
  
  # The number of converged matrices
  nlam_ssdr <- 0
  
  for(i in 1:n1){
    ulam <- as.double(lam1[i])
    
    for(j in 1:n3){
      gamma <- gam[j]
      
      for(k in 1:n2){
        
        lambda2 <- lam2[j,k]
        
        # Maximal interation for outer loop
        sigma <- sigma0 + gamma*diag(rep(1,ncol(sigma0)), ncol(sigma0),ncol(sigma0))

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
        
        # If we get non-sparse matrix for msda, stop here, and leave the rest of matrices, svB, svC, 
        # etc. as NULL
        index <- (i-1)*n2*n3 + (j-1)*n2 + k
        if(jerr < -10000){
          cat('Iteration',index,':lam1 is too small\n')
          break
        }
        # If jerr == 404, then maximal iteration is reached, we leave the matrix as null
        # If jerr == 1, then procedure converges, we save the matrix and sv.
        if(jerr==1){
          nlam_ssdr <- nlam_ssdr + 1
          step_final[[index]] <- step_ssdr
          time_final[[index]] <- difftime(end_time, start_time, units = "secs")
          
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
      
      # If exit because of non-sparsity from msda, we stop trying more lam2s or gammas
      if(jerr < -10000) break
      
    }# End of gam
    
  }# End of lambda1
  
  return(list(beta = mat, rank=r_list, step = mean(unlist(step_final)), time_ssdr = mean(unlist(time_final)), nlam_ssdr = nlam_ssdr, lam1_list = lam1_list, lam2_list = lam2_list, gamma_list = gamma_list, sv_list_B = sv_list_B, sv_list_C = sv_list_C))
  
}
