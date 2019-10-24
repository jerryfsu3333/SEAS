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
setwd("~/Documents/GitHub/ssdr/R/")
source("models.R")
source("utility.R")
source("ssdr_utility.R")
source("ssdr_func.R")
source("rifle_func.R")
source("lasso_func.R")
source("CovSIR.R")
dat<-scan("~/Documents/GitHub/ssdr/data/blood.txt")
dat <- matrix(dat,ncol=71)
dat <- t(dat)
# people with malaria
dat <- dat[23:71,]
y <- dat[,2059]
x <- dat[,-2059]
y <- log(y)
x <- log(x)


###### Prediction ########
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
# times <- 16
# train_size <- 0.8

err_func <- function(x, y, newx, newy){

  m <- apply(x, 2, mean)
  se <- apply(x, 2, sd)
  x <- sweep(sweep(x, 2, m), 2, se, '/')
  newx <- sweep(sweep(newx, 2, m), 2, se, '/')

  x <- data.frame(x)
  newx <- data.frame(newx)
  colnames(x) <- paste0('X', seq_len(ncol(x)))
  colnames(newx) <- paste0('X', seq_len(ncol(newx)))

  data <- cbind(y = y, x)

  model <- lm(y~., data = data)
  prediction <- predict(model, newdata = newx)
  lm_err <- mean((prediction - newy)^2)
  
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

  c(lm_err = lm_err, ker_err = ker_err, rf_err = rf_err, svm_err = svm_err)
}


err_list <- mclapply(seq_len(NROW(x)), function(i){
# err_list <- mclapply(seq_len(times), function(i){
# err_list <- lapply(seq_len(times), function(i){
  cat(c('Time', i, '\n'))
  # index <- sample(1:length(y), train_size*length(y), replace = FALSE)
  index <- -i
  
  train_x <- x[index,, drop=FALSE]
  train_y <- y[index]
  test_x <- x[-index,, drop=FALSE]
  test_y <- y[-index]
  
  # Screen variables
  dist_cor <- sapply(seq_len(dim(train_x)[2]), function(i){
    dcor(exp(train_y), exp(train_x[,i]))
  })
  ord <- order(dist_cor, decreasing = TRUE)[1:200]

  train_x <- train_x[,ord, drop=FALSE]
  test_x <- test_x[,ord, drop=FALSE]
  

  fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                     gamma = c(1,2), nfolds = 5, type = 'sir', pmax = 400)
  # fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.3,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  # gamma = c(1,2), nfolds = 5, type = 'sir', plot = TRUE, pmax = 400)
  # fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.3,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  #                    gamma = c(1,2), nfolds = 5, type = 'sir', pmax = 400)

  if(!is.numeric(fit_sir$Beta)){
    print('A NULL matrix is returned (sir).')
    err_sir <- NA
  }else if(fit_sir$rank==0){
    print('rank is 0 (sir).')
    err_sir <- NA
  }else{
    directions_sir <- fit_sir$Beta
    new_train <- train_x %*% directions_sir
    new_test <- test_x %*% directions_sir
    err_sir <- err_func(new_train, train_y, new_test, test_y)
  }


  fit_intra <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                       gamma = c(1,2), nfolds = 5, type = 'intra', pmax = 400)
  # fit_intra <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.3,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  #                      gamma = c(1,2), nfolds = 5, type = 'intra', plot = TRUE, pmax = 400)
  # fit_intra <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.3,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  #                      gamma = c(1,2), nfolds = 5, type = 'intra', pmax = 400)
  
  if(!is.numeric(fit_intra$Beta)){
    print('A NULL matrix is returned (intra).')
    err_intra <- NA
  }else if(fit_intra$rank==0){
    print('rank is 0 (intra).')
    err_intra <- NA
  }else{
    directions_intra <- fit_intra$Beta
    new_train <- train_x %*% directions_intra
    new_test <- test_x %*% directions_intra
    err_intra <- err_func(new_train, train_y, new_test, test_y)
  }

  
  fit_pfc <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.3,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                     gamma = 0.5, type = 'pfc', nfolds = 5, cut_y = TRUE, maxit_outer = 1e+4, pmax = 400)
  # fit_pfc <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.3,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  # gamma = 0.5, type = 'pfc', nfolds = 5, cut_y = TRUE, maxit_outer = 1e+4, plot = TRUE, pmax = 400)
  # fit_pfc <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.3,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  #                    gamma = 0.5, type = 'pfc', nfolds = 5, cut_y = TRUE, maxit_outer = 1e+4, pmax = 400)
  
  if(!is.numeric(fit_pfc$Beta)){
    print('A NULL matrix is returned (intra).')
    err_pfc <- NA
  }else if(fit_pfc$rank==0){
    print('rank is 0 (pfc).')
    err_pfc <- NA
  }else{
    directions_pfc <- fit_pfc$Beta
    new_train <- train_x %*% directions_pfc
    new_test <- test_x %*% directions_pfc
    err_pfc <- err_func(new_train, train_y, new_test, test_y)
  }



  fit_lassosir <- LassoSIR(train_x, train_y, H = 5, nfolds = 5, choosing.d = 'automatic')
  if(!is.numeric(fit_lassosir$beta)){
    print('A NULL matrix is returned (lassosir).')
    err_lassosir <- NA
  }else{
    directions_lassosir <- fit_lassosir$beta
    new_train <- train_x %*% directions_lassosir
    new_test <- test_x %*% directions_lassosir
    err_lassosir <- err_func(new_train, train_y, new_test, test_y)
  }

  list(err_sir = err_sir, err_intra = err_intra, err_pfc = err_pfc, err_lassosir = err_lassosir)

}, mc.cores = 16)
# })

save(err_list, file = "")
