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
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/rifle_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/lasso_func.R")
data <- readMat('~/Documents/GitHub/ssdr/Real_dataset/NIR.mat')$data

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
    theme(legend.position = 'none') +
    labs(title = title)
  g
}

# Pork (y=1) only
data <- data[data[,1] == 1,]
y <- data[,2]
# x <- data[,-c(1,2)]
x <- data[,-c(1,2,3,4,5)]


############## Visualization ####################
# # fit <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), type = 'sir')
# # fit <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), type = 'intra')
# # fit <- ssdr.cv(x, y, lam1_fac = seq(1.5,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), type = 'pfc', cut_y = FALSE, maxit_outer = 1e+4)
# d <- fit$rank
# directions <- svd(fit$Beta)$u[,1:d, drop=FALSE]
# nz_func(directions)
# x_new <- as.matrix(x) %*% directions
# # hist_plot(x_new, y, 'SSDR-SIR')
# # hist_plot(x_new, y, 'SSDR-intra')
# # hist_plot(x_new, y, 'SSDR-PFC')

# # #### LassoSIR ####
# LassoSIR_fit <- LassoSIR(x, y, H = 5, nfolds = 5, choosing.d = 'automatic')
# d_LassoSIR <- LassoSIR_fit$no.dim
# directions_Lassosir <- LassoSIR_fit$beta
# nz_func(directions_Lassosir)
# x_new_Lassosir <- as.matrix(x) %*% directions_Lassosir
# # df <- data.frame(Component=x_new_Lassosir, class=factor(y))
# hist_plot(x_new_Lassosir, y, 'Lasso-SIR')


#### CovSIR ####
# N <- dim(x)[1]
# p <- dim(x)[2]
# CovSIR_fit <- CovSIR(x, y, Ks = 1:3, lambdas = seq(0.2,2,by=0.5)*sqrt(log(p)/N))
# d_CovSIR <- CovSIR_fit$r
# directions_Covsir <- CovSIR_fit$mat
# sum(directions_Covsir != 0)
# x_new_Covsir <- as.matrix(x) %*% directions_Covsir
# plot(x_new_Covsir[,1], col=y, xlab = 'Index', ylab = 'Component 1')
# abline(h=0, lty = 'dashed')

# #### Lasso ####
# lasso_fit <- lasso_func(x, y, nfolds = 5)[-1,1,drop=FALSE] # the first is zero intercept
# directions_lasso <- lasso_fit
# nz_func(directions_lasso)
# x_new_lasso <- as.matrix(x) %*% directions_lasso
# # hist_plot(x_new_lasso, y, 'Lasso')

# #### Rifle ####
# rifle_fit <- rifle_func(x, y, k=15, type = 'sir')
# directions_rifle <- rifle_fit
# nz_func(directions_rifle)
# x_new_rifle <- as.matrix(x) %*% directions_rifle
# hist_plot(x_new_rifle, y, 'Rifle-SIR')

######## Prediction ##########
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
times <- 5
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

err_list <- mclapply(seq_len(times), function(i){
# err_list <- lapply(seq_len(times), function(i){
  cat(c('Time', i, '\n'))
  index <- sample(1:length(y), train_size*length(y), replace = FALSE)
  # train_index <- drop(createDataPartition(y, p = train_size, list = FALSE))
  # valid_index <- setdiff(seq_len(length(y)), train_index)
  train_x <- x[index,, drop=FALSE]
  train_y <- y[index]
  test_x <- x[-index,, drop=FALSE]
  test_y <- y[-index]
  # train_index <- drop(createDataPartition(y, p = train_size, list = FALSE))
  # valid_index <- setdiff(seq_len(length(y)), train_index)
  # train_x <- x[train_index,]
  # train_y <- y[train_index]
  # test_x <- x[valid_index,]
  # test_y <- y[valid_index]

  train_x <- log(train_x)
  train_y <- log(train_y)
  test_x <- log(test_x)
  test_y <- log(test_y)

  # train_x <- scale(train_x)
  # train_y <- scale(train_y)
  # test_x <- scale(test_x)
  # test_y <- scale(test_y)
  
  fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.7, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), 
                     gamma = c(0.01), nfolds = 3, type = 'sir', plot=TRUE, pmax=100)
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
  
  
  fit_intra <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.7, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                       gamma = c(0.01), nfolds = 3, type = 'intra', plot=TRUE, pmax=100)
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

  
  
  fit_pfc <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.2,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                     gamma = c(0.01), nfolds = 3, type = 'pfc', cut_y = TRUE, plot = TRUE, maxit_outer = 1e+4, pmax=100)
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



  fit_lassosir <- LassoSIR(train_x, train_y, H = 5, nfolds = 5, choosing.d = 'automatic')
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

save(err_list, file = "")
