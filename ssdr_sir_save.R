rm(list = ls())
library(parallel)
library(msda)
library(MASS)
library(methods)
library(glmnet)
library(rifle)
library(LassoSIR)
library(energy)
source("/Users/cengjing/Documents/GitHub/ssdr/models.R")
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/my_msda.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/rifle_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/lasso_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/CovSIR.R")

# #############################################
# Set random seed for parallel computing
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

p <- 500
N <- 500
N_val <- 500

model <- Model1(p)
Data <- model$Data
sir_params <- model$sir_params
intra_params <- model$intra_params
pfc_params <- model$pfc_params
nz_vec <- model$nz_vec
True_sp <- model$True_sp

times <- 1
  
# output <- mclapply(seq_len(times), function(i){
output <- lapply(seq_len(times), function(i){
    
  cat("Time", i, '\n')
  
  data_train <- Data(N)
  data_val <- Data(N_val)
  
  start_time <- Sys.time()
  ssdrsir_fit <- ssdr_func(data_train$x, data_train$y, data_val$x, data_val$y, H = sir_params$H,
                           type = 'sir',
                           lambda.factor = sir_params$lambda.factor,
                           lam1_fac = sir_params$lam1_fac,
                           lam2_fac = sir_params$lam2_fac)
  end_time <- Sys.time()
  time_sir <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  ssdrintra_fit <- ssdr_func(data_train$x, data_train$y, data_val$x, data_val$y, H = intra_params$H,
                             type = 'intra',
                             lambda.factor = intra_params$lambda.factor,
                             lam1_fac = intra_params$lam1_fac,
                             lam2_fac = intra_params$lam2_fac)
  end_time <- Sys.time()
  time_intra <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  ssdrpfc_fit <- ssdr_func(data_train$x, data_train$y, data_val$x, data_val$y, type = 'pfc',
                           lambda.factor = pfc_params$lambda.factor,
                           lam1_fac = pfc_params$lam1_fac,
                           lam2_fac = pfc_params$lam2_fac,
                           cut_y = pfc_params$cut_y)
  end_time <- Sys.time()
  time_pfc <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  LassoSIR_fit <- LassoSIR(data_train$x, data_train$y, H = 5, choosing.d = 'automatic')
  end_time <- Sys.time()
  time_lassosir <- difftime(end_time, start_time, units = "secs")
  
  # start_time <- Sys.time()
  # lasso_fit <- lasso_func(data_train$x, data_train$y)[-1,1,drop=FALSE] # the first is zero intercept
  # end_time <- Sys.time()
  # time_lasso <- difftime(end_time, start_time, units = "secs")
  # 
  # start_time <- Sys.time()
  # rifle_fit <- rifle_func(data_train$x, data_train$y, k = length(nz_vec), type = 'sir') # For Rifle, use true sparsity as k and H = 5.
  # end_time <- Sys.time()
  # time_rifle <- difftime(end_time, start_time, units = "secs")
  
  B_ssdrsir <- ssdrsir_fit$mat
  B_ssdrintra <- ssdrintra_fit$mat
  B_ssdrpfc <- ssdrpfc_fit$mat
  B_LassoSIR <- LassoSIR_fit$beta
  # B_lasso <- lasso_fit
  # B_rifle <- rifle_fit
  
  # plot(ssdrsir_fit$eval)
  # plot(ssdrintra_fit$eval)
  # plot(ssdrpfc_fit$eval)
  
  cat(c(ssdrsir_fit$results$id1, ssdrsir_fit$results$id2, ssdrsir_fit$results$id_gam, '\n'))
  cat(c(ssdrintra_fit$results$id1, ssdrintra_fit$results$id2, ssdrintra_fit$results$id_gam, '\n'))
  cat(c(ssdrpfc_fit$results$id1, ssdrpfc_fit$results$id2, ssdrpfc_fit$results$id_gam, '\n'))
  
  # cat(ssdrsir_fit$svB)
  # cat(ssdrsir_fit$svC)
  # cat(ssdrintra_fit$svB)
  # cat(ssdrintra_fit$svC)
  # cat(ssdrpfc_fit$svB)
  # cat(ssdrpfc_fit$svC)

  # calculate C, IC, subspace distance after we obtain estimated matrix from each method.
  if(is.null(B_ssdrsir)){
    C_IC_ssdrsir <- list(C = NA, IC = NA)
    r_ssdrsir <- NA
    r_ssdrsir_C <- NA
    dist_ssdrsir <- NA
    distord_ssdrsir <- NA
  }else{
    C_IC_ssdrsir <- C_IC(B_ssdrsir, 1:p, nz_vec)
    r_ssdrsir <- ssdrsir_fit$results$r_ssdr
    r_ssdrsir_C <- ssdrsir_fit$results$r_ssdr_C
    dist_ssdrsir <- subspace_2(True_sp, svd(B_ssdrsir)$u[,1:r_ssdrsir, drop = FALSE])
    distord_ssdrsir <- subspace(True_sp, svd(B_ssdrsir)$u[,1:r_ssdrsir, drop = FALSE])
  }

  if(is.null(B_ssdrintra)){
    C_IC_ssdrintra <- list(C = NA, IC = NA)
    r_ssdrintra <- NA
    r_ssdrintra_C <- NA
    dist_ssdrintra <- NA
    distord_ssdrintra <- NA
  }else{
    C_IC_ssdrintra <- C_IC(B_ssdrintra, 1:p, nz_vec)
    r_ssdrintra <- ssdrintra_fit$results$r_ssdr
    r_ssdrintra_C <- ssdrintra_fit$results$r_ssdr_C
    dist_ssdrintra <- subspace_2(True_sp, svd(B_ssdrintra)$u[,1:r_ssdrintra, drop = FALSE])
    distord_ssdrintra <- subspace(True_sp, svd(B_ssdrintra)$u[,1:r_ssdrintra, drop = FALSE])
  }

  if(is.null(B_ssdrpfc)){
    C_IC_ssdrpfc <- list(C = NA, IC = NA)
    r_ssdrpfc <- NA
    r_ssdrpfc_C <- NA
    dist_ssdrpfc <- NA
    distord_ssdrpfc <- NA
  }else{
    C_IC_ssdrpfc <- C_IC(B_ssdrpfc, 1:p, nz_vec)
    r_ssdrpfc <- ssdrpfc_fit$results$r_ssdr
    r_ssdrpfc_C <- ssdrpfc_fit$results$r_ssdr_C
    dist_ssdrpfc <- subspace_2(True_sp, svd(B_ssdrpfc)$u[,1:r_ssdrpfc, drop = FALSE])
    distord_ssdrpfc <- subspace(True_sp, svd(B_ssdrpfc)$u[,1:r_ssdrpfc, drop = FALSE])
  }

  C_IC_LassoSIR <- C_IC(B_LassoSIR, 1:p, nz_vec)
  r_LassoSIR <- LassoSIR_fit$no.dim
  dist_LassoSIR <- subspace_2(True_sp, B_LassoSIR)
  distord_LassoSIR <- subspace(True_sp, B_LassoSIR)

  # C_IC_lasso <- C_IC(B_lasso, 1:p, nz_vec)
  # r_lasso <- 1
  # dist_lasso <- subspace_2(True_sp, B_lasso)
  # distord_lasso <- subspace(True_sp, B_lasso)
  # 
  # C_IC_rifle <- C_IC(B_rifle, 1:p, nz_vec)
  # r_rifle <- 1
  # dist_rifle <- subspace_2(True_sp, B_rifle)
  # distord_rifle <- subspace(True_sp, B_rifle)
  
  c(C_ssdrsir = C_IC_ssdrsir$C, IC_ssdrsir = C_IC_ssdrsir$IC, r_ssdrsir = r_ssdrsir, dist_ssdrsir = dist_ssdrsir, distord_ssdrsir = distord_ssdrsir, time_sir=time_sir,
    C_ssdrintra = C_IC_ssdrintra$C, IC_ssdrintra = C_IC_ssdrintra$IC, r_ssdrintra = r_ssdrintra, dist_ssdrintra = dist_ssdrintra, distord_ssdrintra = distord_ssdrintra, time_intra = time_intra,
    C_ssdrpfc = C_IC_ssdrpfc$C, IC_ssdrpfc = C_IC_ssdrpfc$IC, r_ssdrpfc = r_ssdrpfc, dist_ssdrpfc = dist_ssdrpfc, distord_ssdrpfc = distord_ssdrpfc, time_pfc = time_pfc,
    C_LassoSIR = C_IC_LassoSIR$C, IC_LassoSIR = C_IC_LassoSIR$IC, r_LassoSIR = r_LassoSIR, dist_LassoSIR = dist_LassoSIR, distord_LassoSIR = distord_LassoSIR, time_lassosir = time_lassosir)
  #   C_lasso = C_IC_lasso$C, IC_lasso = C_IC_lasso$IC, r_lasso = r_lasso, dist_lasso = dist_lasso, distord_lasso = distord_lasso, time_lasso = time_lasso,
  #   C_rifle = C_IC_rifle$C, IC_rifle = C_IC_rifle$IC, r_rifle = r_rifle, dist_rifle = dist_rifle, distord_rifle = distord_rifle, time_rifle = time_rifle
  # )
  
  # start_time <- Sys.time()
  # CovSIR_fit <- CovSIR(data_train$x, data_train$y, Ks = 1:3, lambdas = seq(0.2,2,by=0.5)*sqrt(log(p)/N))
  # end_time <- Sys.time()
  # time_covsir <- difftime(end_time, start_time, units = "secs")
  # 
  # B_CovSIR <- CovSIR_fit$mat
  # 
  # C_IC_CovSIR <- C_IC_cut(B_CovSIR, 1:p, nz_vec)
  # r_CovSIR <- CovSIR_fit$r
  # dist_CovSIR <- subspace_2(True_sp, B_CovSIR)
  # distord_CovSIR <- subspace(True_sp, B_CovSIR)
  # 
  # c(C_CovSIR = C_IC_CovSIR$C, IC_CovSIR = C_IC_CovSIR$IC, r_CovSIR = r_CovSIR, dist_CovSIR = dist_CovSIR, distord_CovSIR = distord_CovSIR, time_covsir = time_covsir)
})

output <- do.call(rbind, output)
write.table(output, "/Users/cengjing/Desktop/test3")
# write.table(svB, "/Users/cengjing/Desktop/test_ssdr_1")
# write.table(svC, "/Users/cengjing/Desktop/test_ssdr_2")

