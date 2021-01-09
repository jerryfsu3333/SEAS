rm(list = ls())
library(parallel)
library(msda)
library(MASS)
library(methods)
library(glmnet)
library(rifle)
library(LassoSIR)
library(energy)
source("R/models.R")
source("R/utility.R")
source("R/ssdr_utility.R")
source("R/ssdr_func.R")
source("R/rifle_func.R")
source("R/lasso_func.R")
source("R/CovSIR.R")

# #############################################
# Set random seed for parallel computing
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

p <- 500
N <- 500

model <- Model1(p)
Data <- model$Data
sir_params <- model$sir_params
intra_params <- model$intra_params
pfc_params <- model$pfc_params
nz_vec <- model$nz_vec
True_sp <- model$True_sp

times <- 8

# output <- mclapply(seq_len(times), function(i){
output <- lapply(seq_len(times), function(i){
  
  cat("Time", i, '\n')
  
  data_train <- Data(N)
  
  start_time <- Sys.time()
  ssdrsir_fit <- ssdr.cv(data_train$x, data_train$y, H = sir_params$H,
                         type = 'sir',  nfolds = 5,
                         lambda.factor = sir_params$lambda.factor,
                         lam1_fac = sir_params$lam1_fac,
                         lam2_fac = sir_params$lam2_fac)
  end_time <- Sys.time()
  time_sir <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  ssdrintra_fit <- ssdr.cv(data_train$x, data_train$y, H = intra_params$H,
                         type = 'intra',  nfolds = 5,
                         lambda.factor = intra_params$lambda.factor,
                         lam1_fac = intra_params$lam1_fac,
                         lam2_fac = intra_params$lam2_fac)
  end_time <- Sys.time()
  time_intra <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  ssdrpfc_fit <- ssdr.cv(data_train$x, data_train$y, type = 'pfc', nfolds = 5,
                         lambda.factor = pfc_params$lambda.factor,
                         lam1_fac = pfc_params$lam1_fac,
                         lam2_fac = pfc_params$lam2_fac,
                         cut_y = pfc_params$cut_y)
  end_time <- Sys.time()
  time_pfc <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  LassoSIR_fit <- LassoSIR(data_train$x, data_train$y, H = 5, nfolds = 5, choosing.d = 'automatic')
  end_time <- Sys.time()
  time_lassosir <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  lasso_fit <- lasso_func(data_train$x, data_train$y, nfolds=5)[-1,1,drop=FALSE] # the first is zero intercept
  end_time <- Sys.time()
  time_lasso <- difftime(end_time, start_time, units = "secs")

  start_time <- Sys.time()
  rifle_fit <- rifle_func(data_train$x, data_train$y, k = length(nz_vec), type = 'sir') # For Rifle, use true sparsity as k and H = 5.
  end_time <- Sys.time()
  time_rifle <- difftime(end_time, start_time, units = "secs")

  B_ssdrsir <- ssdrsir_fit$Beta
  B_ssdrintra <- ssdrintra_fit$Beta
  B_ssdrpfc <- ssdrpfc_fit$Beta
  B_LassoSIR <- LassoSIR_fit$beta
  B_lasso <- lasso_fit
  B_rifle <- rifle_fit

  # calculate C, IC, subspace distance after we obtain estimated matrix from each method.
  if(is.null(B_ssdrsir)){
    C_IC_ssdrsir <- list(C = NA, IC = NA)
    r_ssdrsir <- NA
    dist_ssdrsir <- NA
    distord_ssdrsir <- NA
  }else{
    C_IC_ssdrsir <- C_IC(B_ssdrsir, 1:p, nz_vec)
    r_ssdrsir <- ssdrsir_fit$rank
    dist_ssdrsir <- subspace_2(True_sp, B_ssdrsir)
    distord_ssdrsir <- subspace(True_sp, B_ssdrsir)
  }

  if(is.null(B_ssdrintra)){
    C_IC_ssdrintra <- list(C = NA, IC = NA)
    r_ssdrintra <- NA
    dist_ssdrintra <- NA
    distord_ssdrintra <- NA
  }else{
    C_IC_ssdrintra <- C_IC(B_ssdrintra, 1:p, nz_vec)
    r_ssdrintra <- ssdrintra_fit$rank
    dist_ssdrintra <- subspace_2(True_sp, B_ssdrintra)
    distord_ssdrintra <- subspace(True_sp, B_ssdrintra)
  }

  if(is.null(B_ssdrpfc)){
    C_IC_ssdrpfc <- list(C = NA, IC = NA)
    r_ssdrpfc <- NA
    dist_ssdrpfc <- NA
    distord_ssdrpfc <- NA
  }else{
    C_IC_ssdrpfc <- C_IC(B_ssdrpfc, 1:p, nz_vec)
    r_ssdrpfc <- ssdrpfc_fit$rank
    dist_ssdrpfc <- subspace_2(True_sp, B_ssdrpfc)
    distord_ssdrpfc <- subspace(True_sp, B_ssdrpfc)
  }

  C_IC_LassoSIR <- C_IC(B_LassoSIR, 1:p, nz_vec)
  r_LassoSIR <- LassoSIR_fit$no.dim
  dist_LassoSIR <- subspace_2(True_sp, B_LassoSIR)
  distord_LassoSIR <- subspace(True_sp, B_LassoSIR)

  C_IC_lasso <- C_IC(B_lasso, 1:p, nz_vec)
  r_lasso <- 1
  dist_lasso <- subspace_2(True_sp, B_lasso)
  distord_lasso <- subspace(True_sp, B_lasso)

  C_IC_rifle <- C_IC(B_rifle, 1:p, nz_vec)
  r_rifle <- 1
  dist_rifle <- subspace_2(True_sp, B_rifle)
  distord_rifle <- subspace(True_sp, B_rifle)

  c(C_ssdrsir = C_IC_ssdrsir$C, IC_ssdrsir = C_IC_ssdrsir$IC, r_ssdrsir = r_ssdrsir, dist_ssdrsir = dist_ssdrsir, distord_ssdrsir = distord_ssdrsir, time_sir=time_sir,
    C_ssdrintra = C_IC_ssdrintra$C, IC_ssdrintra = C_IC_ssdrintra$IC, r_ssdrintra = r_ssdrintra, dist_ssdrintra = dist_ssdrintra, distord_ssdrintra = distord_ssdrintra, time_intra = time_intra,
    C_ssdrpfc = C_IC_ssdrpfc$C, IC_ssdrpfc = C_IC_ssdrpfc$IC, r_ssdrpfc = r_ssdrpfc, dist_ssdrpfc = dist_ssdrpfc, distord_ssdrpfc = distord_ssdrpfc, time_pfc = time_pfc,
    C_LassoSIR = C_IC_LassoSIR$C, IC_LassoSIR = C_IC_LassoSIR$IC, r_LassoSIR = r_LassoSIR, dist_LassoSIR = dist_LassoSIR, distord_LassoSIR = distord_LassoSIR, time_lassosir = time_lassosir,
    C_lasso = C_IC_lasso$C, IC_lasso = C_IC_lasso$IC, r_lasso = r_lasso, dist_lasso = dist_lasso, distord_lasso = distord_lasso, time_lasso = time_lasso,
    C_rifle = C_IC_rifle$C, IC_rifle = C_IC_rifle$IC, r_rifle = r_rifle, dist_rifle = dist_rifle, distord_rifle = distord_rifle, time_rifle = time_rifle
  )
  
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
# }, mc.cores = 4)

output <- do.call(rbind, output)
# write.table(output, "/Users/cengjing/Desktop/test")

