# The complete ssdr function, consisting of discovering the tunining parameter candidates.

ssdr_func <- function(x_train, y_train, x_val, y_val, type = 'sir', lambda.factor=0.5, 
                      lam_fac_msda=0.8, lam_fac_ssdr=0.8, nlam_msda=10, nlam1=10, nlam2=15, 
                      gamma=c(10,30,50), cut_y=TRUE){
  
  #### The start of our methods
  start_time_tot <- Sys.time()
  
  ################################################
  # MSDA
  ################################################
  
  start_time <- Sys.time()
  fit_1 <- my_msda(x_train, y_train, nlambda=nlam_msda, maxit=1e3, lambda.factor=lambda.factor, cut_y=cut_y, type = type)
  end_time <- Sys.time()
  time_msda <- difftime(end_time, start_time, units = "secs")/nlam_msda
  
  sigma0 <- as.matrix(fit_1$sigma)
  mu0 <- as.matrix(fit_1$mu)
  
  ######################################
  ## Print true estimation
  # print(svd(solve(sigma0) %*% mu0)$d)
  # next
  ######################################
  
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
  # lam1 <- (lam1_min_msda)*seq(1.5,0.4,-0.1)
  # n1 <- length(lam1)
  n1 <- nlam1
  lam_fac_msda <- lam_fac_msda
  lam1 <- lam1_min_msda*lam_fac_msda^seq(0,(n1-1))
  
  # Gamma candidates
  gamma <- gamma
  n3 <- length(gamma)
  
  # Lambda2 candidates
  n2 <- nlam2   # we select n2 lambda2 for each gamma
  lam_fac_ssdr <- lam_fac_ssdr
  d <- svd(B_msda)$d
  lam2 <- d[1] * matrix(gamma, ncol = 1) %*% matrix(lam_fac_ssdr^seq((n2-1),0), nrow = 1)
  
  # if lam2 just contains one single value 0, then ssdr just degenerated to msda
  if (all(lam2 == 0)){
    
    print("All lambda2 are zero, degenerate to msda")
    results <- c(r_msda, lam1_min_msda, id_min_msda, lam1_min_msda, NA, NA, which(lam1 == lam1_min_msda), NA, NA, NA, time_msda, time_eval_msda, NA, NA, NA)
    results <- as.data.frame(t(results))
    colnames(results) <- c("r_ssdr", "lam1_min_msda","id_msda", "lam1_min_ssdr", "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", "step", "time_msda", "teval_msda", "time_ssdr", "teval_ssdr", "time_total")

    return(list(mat = B_msda, results = results, svB = NA, svC = NA))
    
  }else{
    
    # fit with ssdr
    nobs <- as.integer(dim(x_train)[1])
    nvars <- as.integer(dim(x_train)[2])
    
    fit_2 <- ssdr(sigma0, mu0, nobs, nvars, lam1, lam2, gamma)
    
    Beta_ssdr <- fit_2$beta
    
    # In some cases, all the Beta is null because the Fortran code didn't return a converaged B matrix 
    if (sum(sapply(Beta_ssdr, is.null)) == n2*n3) {
      print("No converged matrix returned")
      results <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      results <- as.data.frame(t(results))
      colnames(results) <- c("r_ssdr", "lam1_min_msda","id_msda", "lam1_min_ssdr", "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", "step", "time_msda", "teval_msda", "time_ssdr", "teval_ssdr", "time_total")
      
      return(list(mat = NULL, results = results, svB = NA, svC = NA))
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
    rank_ssdr <- fit_2$rank
    step <- unlist(fit_2$step)
    time_ssdr <- unlist(fit_2$time_ssdr)
    
    sv_list_B <- fit_2$sv_list_B
    sv_list_C <- fit_2$sv_list_C
    
    # Cut negligible columns to zero
    Beta_ssdr <- cut_mat(Beta_ssdr, 1e-3, rank_ssdr)
    
    # rank_ssdr <- rep(0,length(Beta_ssdr))
    # for (i in 1:length(Beta_ssdr)){
    #   mat <- Beta_ssdr[[i]]
    #   rank_ssdr[i] <- rank_func(mat, thrd = 1e-3)
    # }
    
    # validate
    start_time <- Sys.time()
    eval_ssdr <- eval_val_rmse(Beta_ssdr, x_val, y_val)
    end_time <- Sys.time()
    time_eval_ssdr <- difftime(end_time, start_time, units = "secs")
    
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
      
      results <- c(NA, lam1_min_msda, id_min_msda, lam1_min_ssdr, lam2_min_ssdr, gamma_min_ssdr, id_lam1, id_lam2, id_gamma, mean(step), time_msda, time_eval_msda, mean(time_ssdr), time_eval_ssdr, NA)
      results <- as.data.frame(t(results))
      colnames(results) <- c("r_ssdr", "lam1_min_msda","id_msda", "lam1_min_ssdr", "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", "step", "time_msda", "teval_msda", "time_ssdr", "teval_ssdr", "time_total")

      return(list(mat = NULL, results = results, svB = NA, svC = NA))
      
    }else{
      # Calculate C, IC, Frobinious distance, subspace distance
      r_ssdr <- rank_ssdr[[id_min_ssdr]]
      
      # save the singular values of each optimal matrix B and C
      svB <- sv_list_B[[id_min_ssdr]]
      svC <- sv_list_C[[id_min_ssdr]]
      
      # record total time iff we got converged matrix
      end_time_tot <- Sys.time()
      time_total <- difftime(end_time_tot, start_time_tot, units = "secs")
      
      results <- c(r_ssdr, lam1_min_msda, id_min_msda, lam1_min_ssdr, lam2_min_ssdr, gamma_min_ssdr, id_lam1, id_lam2, id_gamma, mean(step), time_msda, time_eval_msda, mean(time_ssdr), time_eval_ssdr, time_total)
      
      results <- as.data.frame(t(results))
      colnames(results) <- c("r_ssdr", "lam1_min_msda","id_msda", "lam1_min_ssdr", "lam2_min_ssdr", "gam_min_ssdr", "id1", "id2", "id_gam", "step", "time_msda", "teval_msda", "time_ssdr", "teval_ssdr", "time_total")
      
      return(list(mat = B_ssdr, results = results, svB = svB, svC = svC))
    }
    
  }
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
          r_list[[index]] <- rank_func(Bnew, thrd = 1e-3)
          
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
  
  return(list(beta = mat, rank=r_list, step = step_final, time_ssdr = time_final, nlam_ssdr = nlam_ssdr, lam1_list = lam1_list, lam2_list = lam2_list, gamma_list = gamma_list, sv_list_B = sv_list_B, sv_list_C = sv_list_C))
  
}