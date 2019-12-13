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
setwd("~/Documents/GitHub/seas/R/")
source("models.R")
source("utility.R")
source("ssdr_utility.R")
source("ssdr_func.R")
source("rifle_func.R")
source("lasso_func.R")
source("CovSIR.R")

data <- readMat('../data/NIR.mat')$data

hist_plot <- function(x, y, title){
  if(!is.factor(y)){y <- factor(y)}
  if(!is.null(dim(x))){x <- drop(x)}
  df <- data.frame(component = x, class = y)
  means <- sapply(unique(df$class), function(i){
    mean(df$component[df$class==i])
  })
  means_df <- data.frame(means = means, class=unique(df$class))
  g <- ggplot(df, aes(x=component, colour=class, fill=class)) +
    geom_histogram(aes(y=..density..), bins = 50, position = 'identity', alpha=0.5) +
    geom_density(alpha=0.3) +
    geom_vline(data = means_df, aes(xintercept=means, color=class), linetype='dashed')+
    theme(
      plot.title = element_text(size=16),
      axis.title.x = element_text(size=16),
      axis.title.y = element_text(size=16)
      # legend.position = 'none'
    )+
    labs(title = title)
  g
}
  

# Pork (y=1) only
data <- data[data[,1] == 1,]
y <- data[,4]
x <- data[,-c(1,4,5)]
y <- log(y)
x <- log(x)


######## Prediction ##########
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
times <- 100
train_size <- 0.8

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

# err_list <- lapply(seq_len(NROW(x)), function(i){
err_list <- mclapply(seq_len(times), function(i){
# err_list <- lapply(seq_len(times), function(i){
  cat(c('Time', i, '\n'))
  index <- sample(1:length(y), train_size*length(y), replace = FALSE)
  
  ## For LOO
  # index <- -i
  
  train_x <- x[index,, drop=FALSE]
  train_y <- y[index]
  test_x <- x[-index,, drop=FALSE]
  test_y <- y[-index]
  
  m <- apply(train_x, 2, mean)
  se <- apply(train_x, 2, sd)
  train_x <- sweep(sweep(train_x, 2, m), 2, se, '/')
  test_x <- sweep(sweep(test_x, 2, m), 2, se, '/')

  m <- mean(train_y)
  se <- sd(train_y)
  train_y <- (train_y-m)/se
  test_y <- (test_y-m)/se
  
  # fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.7, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  #                    gamma = c(0.1), nfolds = 3, type = 'sir', plot=TRUE)
  fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.7, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                     gamma = c(0.1), nfolds = 3, type = 'sir')
  
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


  # fit_intra <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.7, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  #                      gamma = c(0.1), nfolds = 3, type = 'intra', plot=TRUE)
  fit_intra <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.7, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                       gamma = c(0.1), nfolds = 3, type = 'intra')
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
  
  
  # fit_pfc <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.2,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  #                    gamma = c(0.1), nfolds = 3, type = 'pfc', cut_y = TRUE, plot = TRUE, maxit_outer = 1e+4)
  fit_pfc <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.2,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                     gamma = c(0.01), nfolds = 3, type = 'pfc', cut_y = FALSE, maxit_outer = 1e+4)
  if(!is.numeric(fit_pfc$Beta)){
    print('A NULL matrix is returned (pfc).')
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



  fit_lassosir <- LassoSIR(train_x, train_y, H = 5, nfolds = 3, choosing.d = 'automatic')
  if(!is.numeric(fit_lassosir$beta)){
    print('A NULL matrix is returned (lassosir).')
    err_lassosir <- NA
  }else if(fit_lassosir$no.dim==0){
    print('rank is 0 (lassosir).')
    err_lassosir <- NA
  }else{
    directions_lassosir <- fit_lassosir$beta
    new_train <- train_x %*% directions_lassosir
    new_test <- test_x %*% directions_lassosir
    err_lassosir <- err_func(new_train, train_y, new_test, test_y)
  }


  list(err_sir = err_sir, err_intra = err_intra, err_pfc = err_pfc, err_lassosir = err_lassosir)

}, mc.cores=16)
# })

save(err_list, file = "")
