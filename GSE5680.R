rm(list = ls())
library(parallel)
library(msda)
library(MASS)
library(methods)
library(glmnet)
library(energy)
library(caret)
library(LassoSIR)
library(R.matlab)
library(rifle)
library(ggplot2)
library(dplyr)
library(reshape2)
library(np)
library(e1071)
library(randomForest)
source("/Users/cengjing/Documents/GitHub/ssdr/models.R")
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/rifle_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/lasso_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/CovSIR.R")
load('~/Documents/GitHub/ssdr/Real_dataset/GSE5680x_fix')
load('~/Documents/GitHub/ssdr/Real_dataset/GSE5680y_fix')

###### Prediction ########
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
times <- 1
train_size <- 0.7

err_func <- function(x, y, newx, newy){

  m <- apply(x, 2, mean)
  se <- apply(x, 2, sd)
  x <- (x-m)/se
  newx <- (newx-m)/se

  x <- data.frame(x)
  newx <- data.frame(newx)
  colnames(x) <- paste0('X', seq_len(ncol(x)))
  colnames(newx) <- paste0('X', seq_len(ncol(newx)))

  data <- cbind(y = y, x)


  model <- npreg(as.formula(paste(c(names(data)[1], paste(names(data)[-1], collapse = '+')), collapse = '~')),
                 data = data, regtype = 'lc', ckernel='gaussian')
  prediction <- predict(model, newdata = newx)
  ker_err <- mean((prediction - newy)^2)

  model <- randomForest(y~., data = data, sample_size = nrow(x)*0.6)
  prediction <- predict(model, newdata = newx)
  rf_err <- mean((prediction - newy)^2)

  model <- svm(y~., data = data, scale=FALSE, type = 'eps-regression', kernel = 'linear')
  prediction <- predict(model, newdata = newx)
  svm_err <- mean((prediction - newy)^2)

  c(ker_err = ker_err, rf_err = rf_err, svm_err = svm_err)
}

err_list <- lapply(seq_len(times), function(i){
  # for(i in seq_len(times)){
  cat(c('Time', i, '\n'))
  train_index <- drop(createDataPartition(y, p = train_size, list = FALSE))
  valid_index <- setdiff(seq_len(length(y)), train_index)
  train_x <- x[train_index,, drop=FALSE]
  train_y <- y[train_index]
  test_x <- x[valid_index,, drop=FALSE]
  test_y <- y[valid_index]

  # train_x <- log(x[train_index,])
  # train_y <- log(y[train_index])
  # test_x <- log(x[valid_index,])
  # test_y <- log(y[valid_index])
  #
  # # Screen variables
  # dist_cor <- sapply(seq_len(dim(train_x)[2]), function(i){
  #   dcor(exp(train_y), exp(train_x[,i]))
  # })

  # Screen variables
  dist_cor <- sapply(seq_len(dim(train_x)[2]), function(i){
    dcor(train_y, train_x[,i])
  })
  ord <- order(dist_cor, decreasing = TRUE)[1:1500]

  train_x <- scale(log(train_x))
  train_y <- scale(log(train_y))
  test_x <- scale(log(test_x))
  test_y <- scale(log(test_y))

  train_x <- train_x[,ord]
  test_x <- test_x[,ord]

  fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.8, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), gamma = c(1,2,3), nfolds = 5, type = 'sir')
  if(!is.numeric(fit_sir$Beta)){
    warning('A NULL matrix is returned in time ', i, ' (SIR).')
    err_sir <- NA
  }else{
    # d_sir <- fit_sir$rank
    directions_sir <- fit_sir$Beta
    new_train <- train_x %*% directions_sir
    new_test <- test_x %*% directions_sir
    err_sir <- err_func(new_train, train_y, new_test, test_y)
  }

  
  fit_intra <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.8, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), gamma = c(1,2,3), nfolds = 5, type = 'intra')
  if(!is.numeric(fit_intra$Beta)){
    warning('A NULL matrix is returned in time ', i, ' (intra).')
    err_intra <- NA
  }else{
    d_intra <- fit_intra$rank
    directions_intra <- fit_intra$Beta
    new_train <- train_x %*% directions_intra
    new_test <- test_x %*% directions_intra
    err_intra <- err_func(new_train, train_y, new_test, test_y)
  }
  
  fit_pfc <- ssdr.cv(train_x, train_y, lam1_fac = seq(2.0,0.8, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), gamma = c(1,2,3), type = 'pfc', nfolds = 5, cut_y = TRUE, maxit_outer = 1e+4)
  if(!is.numeric(fit_pfc$Beta)){
    warning('A NULL matrix is returned in time ', i, ' (intra).')
    err_pfc <- NA
  }else{
    d_pfc <- fit_pfc$rank
    directions_pfc <- fit_pfc$Beta
    new_train <- train_x %*% directions_pfc
    new_test <- test_x %*% directions_pfc
    err_pfc <- err_func(new_train, train_y, new_test, test_y)
  }



  fit_lassosir <- LassoSIR(train_x, train_y, H = 5, nfolds = 5, choosing.d = 'automatic')
  if(!is.numeric(fit_lassosir$beta)){
    warning('A NULL matrix is returned in time ', i, ' (lassosir).')
    err_lassosir <- NA
  }else{
    directions_lassosir <- fit_lassosir$beta
    new_train <- train_x %*% directions_lassosir
    new_test <- test_x %*% directions_lassosir
    err_lassosir <- err_func(new_train, train_y, new_test, test_y)
  }



  # fit_lasso <- lasso_func(train_x, train_y)[-1,1,drop=FALSE]
  # directions_lasso <- unname(fit_lasso)
  # new_train <- train_x %*% directions_lasso
  # new_test <- test_x %*% directions_lasso
  # err_lasso <- err_func(new_train, train_y, new_test, test_y)
  #
  #
  #
  # fit_rifle <- rifle_func(train_x, train_y, k=15, type = 'sir')
  # directions_rifle <- fit_rifle
  # new_train <- train_x %*% directions_rifle
  # new_test <- test_x %*% directions_rifle
  # err_rifle <- err_func(new_train, train_y, new_test, test_y)

  # list(err_sir = err_sir, err_intra = err_intra, err_pfc = err_pfc, err_lassosir = err_lassosir, err_lasso = err_lasso, err_rifle = err_rifle)
  list(err_sir = err_sir, err_intra = err_intra, err_pfc = err_pfc, err_lassosir = err_lassosir)

})
# }

save(err_list, file = "/Users/cengjing/Desktop/test3")


# ##### Estimation consistency ##########
# output_func <- function(x, y){
# 
#   # x <- log(x)
#   # y <- log(y)
#   # x <- scale(x)
#   # y <- scale(y)
# 
#   # Screen variables
#   # dist_cor <- sapply(seq_len(dim(x)[2]), function(i){
#   #   dcor(exp(y), exp(x[,i]))
#   # })
#   # dist_cor <- sapply(seq_len(dim(x)[2]), function(i){
#   #   dcor(y, x[,i])
#   # })
# 
#   dist_cor <- sapply(seq_len(dim(x)[2]), function(i){
#     dcor(y, x[,i])
#   })
#   ord <- order(dist_cor, decreasing = TRUE)[1:1500]
# 
#   x <- scale(log(x))
#   y <- scale(log(y))
#   x <- x[,ord]
# 
#   fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.8, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), gamma = c(1,2,3), nfolds = 5, type = 'sir')
#   if(!is.numeric(fit_sir$Beta)){
#     warning('A NULL matrix is returned (intra).')
#     d_sir <- NA
#     directions_sir <- NA
#     nz_sir <- NA
#     ord_sir <- NA
#     s_sir <- NA
#   }else{
#     d_sir <- fit_sir$rank
#     directions_sir <- fit_sir$Beta
#     nz_sir <- nz_func(directions_sir)
#     ord_sir <- ord[nz_sir]
#     s_sir <- length(nz_sir)
#   }
# 
#   fit_intra <- ssdr.cv(x, y, lam1_fac = seq(2,0.8, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), gamma = c(1,2,3), nfolds = 5, type = 'intra')
#   if(!is.numeric(fit_intra$Beta)){
#     warning('A NULL matrix is returned (intra).')
#     d_intra <- NA
#     directions_intra <- NA
#     nz_intra <- NA
#     ord_intra <- NA
#     s_intra <- NA
#   }else{
#     d_intra <- fit_intra$rank
#     directions_intra <- fit_intra$Beta
#     nz_intra <- nz_func(directions_intra)
#     ord_intra <- ord[nz_intra]
#     s_intra <- length(nz_intra)
#   }
# 
# 
#   fit_pfc <- ssdr.cv(x, y, lam1_fac = seq(2.0,0.8, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), gamma = c(1,2,3), nfolds = 5, type = 'pfc', cut_y = TRUE, maxit_outer = 1e+4)
#   if(!is.numeric(fit_pfc$Beta)){
#     warning('A NULL matrix is returned (intra).')
#     d_pfc <- NA
#     directions_pfc <- NA
#     nz_pfc <- NA
#     ord_pfc <- NA
#     s_pfc <- NA
#   }else{
#     d_pfc <- fit_pfc$rank
#     directions_pfc <- fit_pfc$Beta
#     nz_pfc <- nz_func(directions_pfc)
#     ord_pfc <- ord[nz_pfc]
#     s_pfc <- length(nz_pfc)
#   }
# 
#   fit_lassosir <- LassoSIR(x, y, H = 5, nfolds = 5, choosing.d = 'automatic')
#   if(!is.numeric(fit_lassosir$beta)){
#     warning('A NULL matrix is returned (lassosir).')
#     d_lassosir <- NA
#     directions_lassosir <- NA
#     nz_lassosir <- NA
#     ord_lassosir <- NA
#     s_lassosir <- NA
#   }else{
#     d_lassosir <- fit_lassosir$no.dim
#     directions_lassosir <- fit_lassosir$beta
#     nz_lassosir <- nz_func(directions_lassosir)
#     ord_lassosir <- ord[nz_lassosir]
#     s_lassosir <- length(nz_lassosir)
#   }
# 
#   # fit_lasso <- lasso_func(x, y, nfolds = 5)[-1,1,drop=FALSE]
#   # directions_lasso <- unname(fit_lasso)
#   # nz_lasso <- nz_func(directions_lasso)
#   # if(length(nz_lasso)==1 && is.na(nz_lasso[1])){
#   #   s_lasso <- NA
#   # }else{
#   #   s_lasso <- length(nz_lasso)
#   # }
#   #
#   # fit_rifle <- rifle_func(x, y, k=15, type = 'sir')
#   # directions_rifle <- fit_rifle
#   # nz_rifle <- nz_func(directions_rifle)
#   # if(length(nz_rifle)==1 && is.na(nz_rifle[1])){
#   #   s_rifle <- NA
#   # }else{
#   #   s_rifle <- length(nz_rifle)
#   # }
# 
#   # output <- list(rank = c(d_sir, d_intra, d_pfc, d_lassosir, 1, 1), s = c(s_sir, s_intra, s_pfc, s_lassosir, s_lasso, s_rifle),
#   #                nz = list(nz_sir, nz_intra, nz_pfc, nz_lassosir, nz_lasso, nz_rifle),
#   #                directions = list(directions_sir, directions_intra, directions_pfc, directions_lassosir, directions_lasso, directions_rifle))
#   output <- list(rank = c(d_sir, d_intra, d_pfc, d_lassosir), s = c(s_sir, s_intra, s_pfc, s_lassosir), nz = list(nz_sir, nz_intra, nz_pfc, nz_lassosir),
#                  ord = ord, ord_est = list(ord_sir, ord_intra, ord_pfc, ord_lassosir),
#                  directions = list(directions_sir, directions_intra, directions_pfc, directions_lassosir))
#   output
# }
# 
# RNGkind("L'Ecuyer-CMRG")
# set.seed(1)
# 
# # Full dataset
# true_output <- output_func(x,y)
# 
# # Bootstrap samples
# times <- 100
# samples <- createResample(y, times = times)
# 
# output <- mclapply(seq_len(times), function(i){
# # output <- lapply(seq_len(times), function(i){
#   cat('Time', i, '\n')
#   index <- samples[[i]]
#   boot_x <- x[index,]
#   boot_y <- y[index]
#   boot_output <- output_func(boot_x, boot_y)
# 
#   directions <- boot_output$directions
#   dist <- sapply(1:length(directions), function(i){
#     if(!is.numeric(true_output$directions[[i]]) | !is.numeric(directions[[i]])){
#       NA
#     }else{
#       subspace(true_output$directions[[i]], directions[[i]])
#     }
#   })
#   
#   list(rank=boot_output$rank, s = boot_output$s, dist = unname(dist), nz = boot_output$nz, ord = boot_output$ord, 
#        ord_est = boot_output$ord_est)
# })
# 
# save(true_output, file = '')
# save(output, file = '')

