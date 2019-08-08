rm(list = ls())
library(parallel)
library(msda)
library(MASS)
library(methods)
library(glmnet)
library(rifle)
library(LassoSIR)
source("/Users/cengjing/Documents/GitHub/ssdr/models.R")
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/my_msda.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/rifle_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/lasso_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/CovSIR.R")

######################## evaluation ########################

# eval_val <- function(Beta, x, y, slices){
#   y_sliced <- cut(y, breaks = slices, include.lowest = TRUE, labels = FALSE)
#   tmp <- prep(x, y_sliced)
#   sigma <- tmp$sigma
#   mu <- tmp$mu
#   l <- length(Beta)
#   result <- rep(0, l)
#   for (i in 1:l){
#     mat <- as.matrix(Beta[[i]])
#     result[i] <- 0.5 * sum(diag(t(mat) %*% sigma %*% mat - 2 * t(mu) %*% mat), na.rm = TRUE)
#   }
#   return(result)
# }

# eval_val_rmse <- function(Beta, x, y){
#   l <- length(Beta)
#   # seq_len(l) is 1:l, but faster
#   result <- sapply(seq_len(l), function(i){
#     if(is.null(Beta[[i]])){
#       NA
#     }else{
#       mat <- as.matrix(Beta[[i]])
#       xnew <- x %*% mat
#       fit <- lm(y~xnew)
#       rmse <- sqrt(mean((fit$residuals)^2))
#       rmse 
#     }
#   })
#   result
# }

# eval_val_rmse_rank <- function(Beta, x, y, rank, slices = NULL){
#   l <- length(Beta)
#   result <- rep(0, l)
#   for (i in 1:l){
#     r <- rank[i]
#     mat <- as.matrix(Beta[[i]])
#     if(r != 0){
#       mat <- svd(as.matrix(Beta[[i]]))$u[,1:r]
#     }
#     xnew <- x %*% mat
#     fit <- lm(y~xnew)
#     rmse <- sqrt(mean((fit$residuals)^2))
#     result[i] <- rmse
#   }
#   return(result)
# }


# eval_val_rmse_2 <- function(Beta, x, y, slices){
#   y_sliced <- cut(y, breaks = slices, include.lowest = TRUE, labels = FALSE)
#   l <- length(Beta)
#   result <- rep(0, l)
#   for (i in 1:l){
#     mat <- as.matrix(Beta[[i]])
#     xnew <- x %*% mat
#     fit <- lm(y_sliced~xnew)
#     rmse <- sqrt(mean((fit$residuals)^2))
#     result[i] <- rmse
#   }
#   return(result)
# }

# eval_val_cart <- function(Beta, xtrain, ytrain, xval, yval, slices_val){
#   y_val_sliced <- cut(yval, breaks = slices_val, include.lowest = TRUE, labels = FALSE)
#   l <- length(Beta)
#   errors <- rep(0, l)
#   for (i in 1:l){
#     mat <- as.matrix(Beta[[i]])
#     xnew_tr <- xtrain %*% mat
#     xnew_val <- as.data.frame(xval %*% mat)
#     fit <- rpart(ytrain~xnew_tr, method = "class", control=rpart.control(minsplit=20, cp=0.001))
#     pred <- predict(fit, xnew_val, type = "class")
#     errors[i] <- mean(y_val_sliced != pred)
#   }
#   return(errors)
# }

# #############################################
# Set random seed for parallel computing
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

p <- 500
N <- 500
N_val <- 500

model <- Model22(p)
Data <- model$Data
sir_params <- model$sir_params
intra_params <- model$intra_params
pfc_params <- model$pfc_params
nz_vec <- model$nz_vec
True_sp <- model$True_sp

times <- 5
  
output <- mclapply(seq_len(times), function(i){
# output <- lapply(seq_len(times), function(i){
    
  cat("Time", i, '\n')
  
  data_train <- Data(N)
  data_val <- Data(N_val)
  
  start_time <- Sys.time()
  ssdrsir_fit <- ssdr_func(data_train$x, data_train$y, data_val$x, data_val$y, H = sir_params$H,
                           type = 'sir',
                           lambda.factor = sir_params$lambda.factor,
                           lam_fac_msda = sir_params$lam_fac_msda,
                           lam_fac_ssdr = sir_params$lam_fac_ssdr)
  end_time <- Sys.time()
  time_sir <- difftime(end_time, start_time, units = "secs")
  
  start_time <- Sys.time()
  ssdrintra_fit <- ssdr_func(data_train$x, data_train$y, data_val$x, data_val$y, H = intra_params$H,
                             type = 'intra',
                             lambda.factor = intra_params$lambda.factor,
                             lam_fac_msda = intra_params$lam_fac_msda,
                             lam_fac_ssdr = intra_params$lam_fac_ssdr)
  end_time <- Sys.time()
  time_intra <- difftime(end_time, start_time, units = "secs")
  
  start_time <- Sys.time()
  ssdrpfc_fit <- ssdr_func(data_train$x, data_train$y, data_val$x, data_val$y, type = 'pfc',
                           lambda.factor = pfc_params$lambda.factor,
                           lam_fac_msda = pfc_params$lam_fac_msda,
                           lam_fac_ssdr = pfc_params$lam_fac_ssdr,
                           cut_y = pfc_params$cut_y)
  end_time <- Sys.time()
  time_pfc <- difftime(end_time, start_time, units = "secs")
  
  start_time <- Sys.time()
  LassoSIR_fit <- LassoSIR(data_train$x, data_train$y, H = 5, choosing.d = 'automatic')
  end_time <- Sys.time()
  time_lassosir <- difftime(end_time, start_time, units = "secs")
  
  # start_time <- Sys.time()
  # CovSIR_fit <- CovSIR(data_train$x, data_train$y, Ks = 1:3, lambdas = seq(0.2,2,by=0.5)*sqrt(log(p)/N))
  # end_time <- Sys.time()
  # time_covsir <- difftime(end_time, start_time, units = "secs")
  
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
  # B_CovSIR <- CovSIR_fit$mat
  # B_lasso <- lasso_fit
  # B_rifle <- rifle_fit
  
  # plot(ssdrsir_fit$eval)
  # plot(ssdrintra_fit$eval)
  # plot(ssdrpfc_fit$eval)
  
  # cat(c(ssdrsir_fit$results$id1, ssdrsir_fit$results$id2, ssdrsir_fit$results$id_gam, '\n'))
  # cat(c(ssdrintra_fit$results$id1, ssdrintra_fit$results$id2, ssdrintra_fit$results$id_gam, '\n'))
  # cat(c(ssdrpfc_fit$results$id1, ssdrpfc_fit$results$id2, ssdrpfc_fit$results$id_gam, '\n'))

  # calculate C, IC, subspace distance after we obtain estimated matrix from each method.
  if(is.null(B_ssdrsir)){
    C_IC_ssdrsir <- list(C = NA, IC = NA)
    r_ssdrsir <- NA
    dist_ssdrsir <- NA
  }else{
    C_IC_ssdrsir <- C_IC(B_ssdrsir, 1:p, nz_vec)
    r_ssdrsir <- ssdrsir_fit$results$r_ssdr
    dist_ssdrsir <- subspace_2(True_sp, svd(B_ssdrsir)$u[,1:r_ssdrsir, drop = FALSE])
  }

  if(is.null(B_ssdrintra)){
    C_IC_ssdrintra <- list(C = NA, IC = NA)
    r_ssdrintra <- NA
    dist_ssdrintra <- NA
  }else{
    C_IC_ssdrintra <- C_IC(B_ssdrintra, 1:p, nz_vec)
    r_ssdrintra <- ssdrintra_fit$results$r_ssdr
    dist_ssdrintra <- subspace_2(True_sp, svd(B_ssdrintra)$u[,1:r_ssdrintra, drop = FALSE])
  }

  if(is.null(B_ssdrpfc)){
    C_IC_ssdrpfc <- list(C = NA, IC = NA)
    r_ssdrpfc <- NA
    dist_ssdrpfc <- NA
  }else{
    C_IC_ssdrpfc <- C_IC(B_ssdrpfc, 1:p, nz_vec)
    r_ssdrpfc <- ssdrpfc_fit$results$r_ssdr
    dist_ssdrpfc <- subspace_2(True_sp, svd(B_ssdrpfc)$u[,1:r_ssdrpfc, drop = FALSE])
  }

  # C_IC_LassoSIR <- C_IC_cut(B_LassoSIR, 1:p, nz_vec)
  C_IC_LassoSIR <- C_IC(B_LassoSIR, 1:p, nz_vec)
  r_LassoSIR <- LassoSIR_fit$no.dim
  dist_LassoSIR <- subspace_2(True_sp, B_LassoSIR)
  # print(svd(B_LassoSIR)$d)
  
  # C_IC_CovSIR <- C_IC_cut(B_CovSIR, 1:p, nz_vec)
  # r_CovSIR <- CovSIR_fit$r
  # dist_CovSIR <- subspace_2(True_sp, B_CovSIR)
  # 
  # C_IC_lasso <- C_IC(B_lasso, 1:p, nz_vec)
  # r_lasso <- 1
  # dist_lasso <- subspace_2(True_sp, B_lasso)
  # 
  # C_IC_rifle <- C_IC(B_rifle, 1:p, nz_vec)
  # r_rifle <- 1
  # dist_rifle <- subspace_2(True_sp, B_rifle)
  
  # list(C_CovSIR = C_IC_CovSIR$C, IC_CovSIR = C_IC_CovSIR$IC, r_CovSIR = r_CovSIR, dist_CovSIR = dist_CovSIR, time_covsir = time_covsir)
  c(C_ssdrsir = C_IC_ssdrsir$C, IC_ssdrsir = C_IC_ssdrsir$IC, r_ssdrsir = r_ssdrsir, dist_ssdrsir = dist_ssdrsir, time_sir=time_sir,
       C_ssdrintra = C_IC_ssdrintra$C, IC_ssdrintra = C_IC_ssdrintra$IC, r_ssdrintra = r_ssdrintra, dist_ssdrintra = dist_ssdrintra, time_intra = time_intra,
       C_ssdrpfc = C_IC_ssdrpfc$C, IC_ssdrpfc = C_IC_ssdrpfc$IC, r_ssdrpfc = r_ssdrpfc, dist_ssdrpfc = dist_ssdrpfc, time_pfc = time_pfc,
       C_LassoSIR = C_IC_LassoSIR$C, IC_LassoSIR = C_IC_LassoSIR$IC,  r_LassoSIR = r_LassoSIR, dist_LassoSIR = dist_LassoSIR, time_lassosir = time_lassosir)
       # C_CovSIR = C_IC_CovSIR$C, IC_CovSIR = C_IC_CovSIR$IC, r_CovSIR = r_CovSIR, dist_CovSIR = dist_CovSIR,
       # C_lasso = C_IC_lasso$C, IC_lasso = C_IC_lasso$IC,  r_lasso = r_lasso, dist_lasso = dist_lasso, time_lasso = time_lasso,
       # C_rifle = C_IC_rifle$C, IC_rifle = C_IC_rifle$IC,  r_rifle = r_rifle, dist_rifle = dist_rifle, time_rifle = time_rifle)
})

output <- do.call(rbind, output)
write.table(output, "/Users/cengjing/Desktop/test3")
# write.table(svB, "/Users/cengjing/Desktop/test_ssdr_1")
# write.table(svC, "/Users/cengjing/Desktop/test_ssdr_2")

