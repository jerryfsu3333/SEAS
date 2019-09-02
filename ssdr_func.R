# The complete ssdr function, consisting of discovering the tunining parameter candidates.

ssdr_func <- function(x_train, y_train, x_val, y_val, H=5, type = 'sir', lambda.factor=0.5, nlam_msda=10,
                      lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10),
                      gamma=c(10,30,50), cut_y=TRUE){
  
  #### The start of our methods
  start_time_tot <- Sys.time()
  
  ################################################
  # MSDA
  ################################################
  
  start_time <- Sys.time()
  fit_1 <- my_msda(x_train, y_train, H=H, type = type, nlambda=nlam_msda, maxit=1e3, lambda.factor=lambda.factor, cut_y=cut_y)
  end_time <- Sys.time()
  time_msda <- difftime(end_time, start_time, units = "secs")/nlam_msda
  
  sigma0 <- as.matrix(fit_1$sigma)
  mu0 <- as.matrix(fit_1$mu)
  
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
  eval_msda <- eval_val_rmse(Beta_msda, x_val, y_val)
  end_time <- Sys.time()
  time_eval_msda <- difftime(end_time, start_time, units = "secs")
  
  # The optimal lambda1
  id_min_msda <- which.min(eval_msda)
  lam1_min_msda <- lam_msda[id_min_msda]
  
  # #####
  plot(1:length(eval_msda), eval_msda)
  points(id_min_msda, eval_msda[id_min_msda], col = 'red')
  # ####
  
  # calculate C, IC, Frobenious distance, rank and subspace distance
  B_msda <- as.matrix(Beta_msda[[id_min_msda]])
  # tmp <- apply(B_msda, 1, function(x) any(x!=0))
  # C_msda <- sum(which(tmp) %in% nz_vec)/length(nz_vec)
  # IC_msda <- sum(which(tmp) %in% setdiff(1:p, nz_vec))/(p - length(nz_vec))
  # r_msda <- rank_msda[id_min_msda]
  
  
  ################################################
  # SSDR
  ################################################
  
  # Lambda1 candidates
  lam1 <- (lam1_min_msda)*lam1_fac
  n1 <- length(lam1)
  
  # Gamma candidates
  gamma <- gamma
  n3 <- length(gamma)
  
  # Lambda2 candidates
  d <- svd(B_msda)$d
  lam2 <- d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
  n2 <- dim(lam2)[2]
  
  # if lam2 just contains one single value 0, then msda matrix is zero matrix
  if (all(lam2 == 0)){
    
    cat("All lambda2 are zero, msda matrix is zero matrix\n")
    results <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    results <- as.data.frame(t(results))
    colnames(results) <- c("r_ssdr", "lam1_min_msda","id_msda", "lam1_min_ssdr", "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", "step", "time_msda", "teval_msda", "time_ssdr", "teval_ssdr", "time_total")
    
    return(list(mat = NULL, results = results, eval = NA, svB = NULL, svC = NULL))
    
  }else{
    
    # fit with ssdr
    nobs <- as.integer(dim(x_train)[1])
    nvars <- as.integer(dim(x_train)[2])
    
    fit_2 <- ssdr(sigma0, mu0, nobs, nvars, lam1, lam2, gamma)
    
    Beta_ssdr <- fit_2$beta
    
    # In some cases, all the Beta is null because the Fortran code didn't return a converaged B matrix 
    if (all(sapply(Beta_ssdr, is.null))) {
      print("No converged matrix returned")
      results <- c(NA,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      results <- as.data.frame(t(results))
      colnames(results) <- c("r_ssdr", "lam1_min_msda","id_msda", "lam1_min_ssdr", "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", "step", "time_msda", "teval_msda", "time_ssdr", "teval_ssdr", "time_total")
      
      return(list(mat = NULL, results = results, eval = NA, svB = NULL, svC = NULL))
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
    rank_ssdr_B <- fit_2$rank_B
    step <- fit_2$step
    time_ssdr <- fit_2$time_ssdr
    # rank_ssdr_C <- fit_2$rank_C
    
    sv_list_B <- fit_2$sv_list_B
    sv_list_C <- fit_2$sv_list_C
    
    # Cut negligible columns to zero
    Beta_ssdr <- cut_mat(Beta_ssdr, 1e-3, rank_ssdr_B)
    
    # Recalculate the rank after the cut
    rank_ssdr_B <- vector("list", length(Beta_ssdr))
    for(i in 1:length(Beta_ssdr)){
      if(!is.null(Beta_ssdr[[i]])){
      rank_ssdr_B[[i]] <- rank_func(Beta_ssdr[[i]], thrd = 1e-3)
      }
    }
    
    # validate
    start_time <- Sys.time()
    eval_ssdr <- eval_val_dc(Beta_ssdr, x_val, y_val, d = rank_ssdr_B)
    
    end_time <- Sys.time()
    time_eval_ssdr <- difftime(end_time, start_time, units = "secs")
    
    
    ############################
    ind <- which(sapply(Beta_ssdr, is.null))
    rank_ssdr_B[ind] <- 0
    eval_ssdr[ind] <- max(eval_ssdr, na.rm = TRUE)
    plot(1:length(eval_ssdr), eval_ssdr)
    points(which(rank_ssdr_B > 2), eval_ssdr[rank_ssdr_B > 2], col = 'green')
    points(which(rank_ssdr_B == 2), eval_ssdr[rank_ssdr_B == 2], col = 'red')
    points(which(rank_ssdr_B == 1), eval_ssdr[rank_ssdr_B == 1], col = 'blue')
    points(ind, eval_ssdr[ind], pch=4)

    for(i in 1:(n1-1)){
      abline(v = n2*n3*i+1, lty = 'dashed')
    }
    ##########################
    
    # The optimal lambda1 and lambda2 
    id_min_ssdr <- which.min(eval_ssdr)
    lam1_min_ssdr <- lam1_list[[id_min_ssdr]]
    lam2_min_ssdr <- lam2_list[[id_min_ssdr]]
    gamma_min_ssdr <- gamma_list[[id_min_ssdr]]
    id_lam1 <- which(lam1_min_ssdr == lam1)
    id_lam2 <- which(lam2_min_ssdr == lam2, arr.ind = TRUE)[2]
    id_gamma <- which(gamma_min_ssdr == gamma)
    
    # The optimal ssdr
    B_ssdr <- Beta_ssdr[[id_min_ssdr]]
    
    if(is.null(B_ssdr)){
      print("Optimal matrix is a null matrix")
      
      results <- c(NA, lam1_min_msda, id_min_msda, lam1_min_ssdr, lam2_min_ssdr, gamma_min_ssdr, id_lam1, id_lam2, id_gamma, mean(unlist(step)), time_msda, time_eval_msda, mean(unlist(time_ssdr)), time_eval_ssdr, NA)
      results <- as.data.frame(t(results))
      colnames(results) <- c("r_ssdr", "lam1_min_msda","id_msda", "lam1_min_ssdr", "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", "step", "time_msda", "teval_msda", "time_ssdr", "teval_ssdr", "time_total")

      return(list(mat = NULL, results = results, eval = eval_ssdr, svB = NULL, svC = NULL))
      
    }else{
      # Calculate C, IC, Frobinious distance, subspace distance
      r_ssdr <- rank_ssdr_B[[id_min_ssdr]]
      
      # r_ssdr_C <- rank_ssdr_C[[id_min_ssdr]]
      # save the singular values of each optimal matrix B and C
      svB <- sv_list_B[[id_min_ssdr]]
      svC <- sv_list_C[[id_min_ssdr]]
      
      # record total time iff we got converged matrix
      end_time_tot <- Sys.time()
      time_total <- difftime(end_time_tot, start_time_tot, units = "secs")
      
      results <- c(r_ssdr, lam1_min_msda, id_min_msda, lam1_min_ssdr, lam2_min_ssdr, gamma_min_ssdr, id_lam1, id_lam2, id_gamma, mean(unlist(step)), time_msda, time_eval_msda, mean(unlist(time_ssdr)), time_eval_ssdr, time_total)
      
      results <- as.data.frame(t(results))
      colnames(results) <- c("r_ssdr", "lam1_min_msda","id_msda", "lam1_min_ssdr", "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", "step", "time_msda", "teval_msda", "time_ssdr", "teval_ssdr", "time_total")
      
      return(list(mat = B_ssdr, results = results, eval = eval_ssdr, svB = svB, svC = svC))
    }
    
  }
}


ssdr.cv <- function(x, y, H=5, categorical=FALSE, type = 'sir', lambda.factor=0.5, nlam_msda=10, 
                    lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10),
                    gamma=c(10,30,50), cut_y=TRUE, nfold = 5){
  col.names <- colnames(x)
  x <- as.matrix(x)
  y <- drop(y)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  # Cross validation with msda to find lambda1_msda
  fit_1 <- msda.cv(x, y, H=H, categorical=categorical, type=type, nlam=nlam_msda, lambda.factor=lambda.factor, cut_y=cut_y, nfold=nfold, maxit=1e3)
  sigma0 <- as.matrix(fit_1$sigma)
  mu0 <- as.matrix(fit_1$mu)
  id_min_msda <- fit_1$id
  lam1_min_msda <- fit_1$lambda
  rank_min_msda <- fit_1$rank
  B_msda <- as.matrix(fit_1$Beta)
  
  # Generate tuning parameter candidates
  lam1 <- (lam1_min_msda)*lam1_fac
  n1 <- length(lam1)
  
  gamma <- gamma
  n3 <- length(gamma)
  
  d <- svd(B_msda)$d
  lam2 <- d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
  n2 <- dim(lam2)[2]
  
  # if lam2 just contains one single value 0, then ssdr just degenerated to msda
  if (all(lam2 == 0)){
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
      prep_fold <- prep(x_train, y_train, H=H, categorical = categorical, type = type, cut_y=cut_y)
      sigma_fold <- prep_fold$sigma
      mu_fold <- prep_fold$mu
      fit_fold <- ssdr(sigma_fold, mu_fold, nobs_fold, nvars_fold, lam1, lam2, gamma)
      
      Beta_fold <- fit_fold$beta
      if (all(sapply(Beta_fold, is.null))) {
        cat("Fold",k,":No converged matrix returned\n")
        return(rep(NA, length(Beta_fold)))
      }
      rank_fold <- fit_fold$rank_B
      Beta_fold <- cut_mat(Beta_fold, 1e-3, rank_fold)
      
      # Recalculate the rank after the cut
      rank_fold<- vector("list", length(Beta_fold))
      for(i in 1:length(Beta_fold)){
        if(!is.null(Beta_fold[[i]])){
          rank_fold[[i]] <- rank_func(Beta_fold[[i]], thrd = 1e-3)
        }
      }
      eval <- eval_val_dc(Beta_fold, x_val, y_val, d = rank_fold)
      
      ############################
      ############################
      ind <- which(sapply(Beta_fold, is.null))
      rank_fold_copy <- rank_fold
      eval_copy <- eval
      rank_fold_copy[ind] <- 0
      eval_copy[ind] <- max(eval_copy, na.rm = TRUE)
      plot(1:length(eval_copy), eval_copy)
      points(which(rank_fold_copy > 2), eval_copy[rank_fold_copy > 2], col = 'green')
      points(which(rank_fold_copy == 2), eval_copy[rank_fold_copy == 2], col = 'red')
      points(which(rank_fold_copy == 1), eval_copy[rank_fold_copy == 1], col = 'blue')
      points(ind, eval_copy[ind], pch=4)
      
      for(i in 1:(n1-1)){
        abline(v = n2*n3*i+1, lty = 'dashed')
      }
      ##########################
      ############################
      eval
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
    
    # prep_full <- prep(x,y,type,H,cut_y)
    # sigma0 <- prep_full$sigma
    # mu0 <- prep_full$mu
    fit_full <- ssdr(sigma0, mu0, nobs, nvars, lam1_min_ssdr, matrix(lam2_min_ssdr,1,1), gamma_min_ssdr)
    
    Beta_ssdr <- fit_full$beta
    rank_ssdr <- fit_full$rank_B
    Beta_ssdr <- cut_mat(Beta_ssdr, 1e-3, rank_ssdr)
    
    rank_ssdr <- vector("list", 1)
    if(!is.null(Beta_ssdr[[1]])){
      rank_ssdr[[1]] <- rank_func(Beta_ssdr[[1]], thrd = 1e-3)
    }
    
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


msda.cv <- function(x, y, H, categorical, type, nlam, lambda.factor, cut_y=FALSE, nfold=5, maxit=1e3){
  
  # Fit full data, obtain the msda lambda candidates
  fit <- msda_func(x, y, H=H, categorical=categorical, type=type, nlam=nlam, lambda.factor=lambda.factor, cut_y=cut_y, maxit=maxit)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  sigma0 <- as.matrix(fit$sigma)
  mu0 <- as.matrix(fit$mu)
  
  lam_msda <- fit$lambda
  rank_msda <- fit$rank
  Beta_msda <- fit$Beta
  
  # Cross-validation
  fold <-   sample(rep(seq(nfold), length = nobs))
  eval_msda <- sapply(1:nfold, function(k){
    x_train <- x[which(fold!=k),,drop=FALSE]		
    x_val <- x[which(fold==k),,drop=FALSE]		
    y_train <- y[which(fold!=k)]
    y_val <- y[which(fold==k)]
    
    # matrix is already cut inside msda_func
    fit_fold <- msda_func(x_train, y_train, H=H, categorical = categorical, type=type, lambda = lam_msda, nlam = nlam, lambda.factor = lambda.factor,cut_y = cut_y, maxit = maxit)
    Beta_fold <- fit_fold$Beta
    # return evaluation of each fold
    eval_val_rmse(Beta_fold, x_val, y_val)
  })
  
  eval <- apply(eval_msda, 1, mean)
  
  # The optimal lambda1
  id_min_msda <- which.min(eval)
  lam1_min_msda <- lam_msda[id_min_msda]
  rank_min_msda <- rank_msda[id_min_msda]
  
  # #####
  plot(1:length(eval), eval)
  points(id_min_msda, eval[id_min_msda], col = 'red')
  # ####
  
  # calculate C, IC, Frobenious distance, rank and subspace distance
  B_msda <- as.matrix(Beta_msda[[id_min_msda]])
  
  list(id = id_min_msda, lambda = lam1_min_msda, rank = rank_min_msda, Beta = B_msda, sigma = sigma0, mu = mu0)
}


msda_func <-function(x, y, H, categorical, type, nlam, lambda.factor, lambda=NULL, cut_y=FALSE, maxit=1e3){
  
  fit <- my_msda(x, y, H=H, categorical, type= type, lambda=lambda, nlambda=nlam, maxit=maxit, lambda.factor=lambda.factor, cut_y=cut_y)
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
  r_list_B <- vector("list", nparams)
  # r_list_C <- vector("list", nparams)
  
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
          print('lam1 are too small!')
          break
        }
        
        # If jerr == 404, then maximal iteration is reached, we leave the matrix as null
        
        # If jerr == 1, then procedure converges, we save the matrix and sv.
        if(jerr==1){
          
          index <- (i-1)*n2*n3 + (j-1)*n2 + k
          nlam_ssdr <- nlam_ssdr + 1
          step_final[[index]] <- step_ssdr
          time_final[[index]] <- difftime(end_time, start_time, units = "secs")
          
          mat[[index]] <- Bnew
          r_list_B[[index]] <- rank_func(Bnew, thrd = 1e-3)
          
          # tol_rank <- max(dim(Cnew)) * .Machine$double.eps
          # r_list_C[[index]] <- rank_func2(Cnew, thrd = tol_rank)
          
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
  
  return(list(beta = mat, rank_B = r_list_B, step = step_final, time_ssdr = time_final, nlam_ssdr = nlam_ssdr, lam1_list = lam1_list, lam2_list = lam2_list, gamma_list = gamma_list, sv_list_B = sv_list_B, sv_list_C = sv_list_C))
  
}