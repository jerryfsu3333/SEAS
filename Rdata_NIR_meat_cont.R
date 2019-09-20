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
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/my_msda.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/CovSIR.R")
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

# import dataset
y <- data[,2]
# x <- data[,-c(1,3)]
x <- data[,-(1:5)]
y <- scale(y)
x <- scale(x)


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

err_list <- lapply(seq_len(times), function(i){
# for(i in seq_len(times)){
  cat(c('Time', i, '\n'))
  train_index <- drop(createDataPartition(y, p = train_size, list = FALSE))
  valid_index <- setdiff(seq_len(length(y)), train_index)
  train_x <- x[train_index,]
  train_y <- y[train_index]
  test_x <- x[valid_index,]
  test_y <- y[valid_index]

  fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.8, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), type = 'sir')
  if(!is.numeric(fit_sir$Beta)){
    warning('A NULL matrix is returned in time ', i, ' (SIR).')
    err_sir <- NA
  }else{
    d_sir <- fit_sir$rank
    directions_sir <- svd(fit_sir$Beta)$u[,1:d_sir, drop=FALSE]
    new_train <- train_x %*% directions_sir
    new_test <- test_x %*% directions_sir
    err_sir <- err_func(new_train, train_y, new_test, test_y)
  }


  fit_intra <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.8, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), type = 'intra')
  if(!is.numeric(fit_intra$Beta)){
    warning('A NULL matrix is returned in time ', i, ' (intra).')
    err_intra <- NA
  }else{
    d_intra <- fit_intra$rank
    directions_intra <- svd(fit_intra$Beta)$u[,1:d_intra, drop=FALSE]
    new_train <- train_x %*% directions_intra
    new_test <- test_x %*% directions_intra
    err_intra <- err_func(new_train, train_y, new_test, test_y)
  }

  
  fit_pfc <- ssdr.cv(train_x, train_y, lam1_fac = seq(1.5,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), type = 'pfc', cut_y = FALSE, maxit_outer = 1e+4)
  if(!is.numeric(fit_pfc$Beta)){
    warning('A NULL matrix is returned in time ', i, ' (intra).')
    err_pfc <- NA
  }else{
    d_pfc <- fit_pfc$rank
    directions_pfc <- svd(fit_pfc$Beta)$u[,1:d_pfc, drop=FALSE]
    new_train <- train_x %*% directions_pfc
    new_test <- test_x %*% directions_pfc
    err_pfc <- err_func(new_train, train_y, new_test, test_y)
  }

  

  fit_lassosir <- LassoSIR(train_x, train_y, H = 5, nfolds = 5, choosing.d = 'automatic')
  directions_lassosir <- fit_lassosir$beta
  new_train <- train_x %*% directions_lassosir
  new_test <- test_x %*% directions_lassosir
  err_lassosir <- err_func(new_train, train_y, new_test, test_y)



  fit_lasso <- lasso_func(train_x, train_y)[-1,1,drop=FALSE]
  directions_lasso <- unname(fit_lasso)
  new_train <- train_x %*% directions_lasso
  new_test <- test_x %*% directions_lasso
  err_lasso <- err_func(new_train, train_y, new_test, test_y)



  fit_rifle <- rifle_func(train_x, train_y, k=15, type = 'sir')
  directions_rifle <- fit_rifle
  new_train <- train_x %*% directions_rifle
  new_test <- test_x %*% directions_rifle
  err_rifle <- err_func(new_train, train_y, new_test, test_y)

  list(err_sir = err_sir, err_intra = err_intra, err_pfc = err_pfc, err_lassosir = err_lassosir, err_lasso = err_lasso, err_rifle = err_rifle)

})
# }

save(err_list, file = "/Users/cengjing/Desktop/test3")


######## Estimation consistency ##########
# output_func <- function(x, y){
# 
#   fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), type = 'sir')
#   d_sir <- fit_sir$rank
#   directions_sir <- svd(fit_sir$Beta)$u[,1:d_sir, drop=FALSE]
#   nz_sir <- nz_func(directions_sir)
#   if(length(nz_sir)==1 && is.na(nz_sir[1])){
#     s_sir <- NA
#   }else{
#     s_sir <- length(nz_sir)
#   }
# 
# 
#   fit_intra <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), type = 'intra')
#   d_intra <- fit_intra$rank
#   directions_intra <- svd(fit_intra$Beta)$u[,1:d_intra, drop=FALSE]
#   nz_intra <- nz_func(directions_intra)
#   if(length(nz_intra)==1 && is.na(nz_intra[1])){
#     s_intra <- NA
#   }else{
#     s_intra <- length(nz_intra)
#   }
# 
# 
#   fit_pfc <- ssdr.cv(x, y, lam1_fac = seq(1.5,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), type = 'pfc', cut_y = TRUE, maxit_outer = 1e+4)
#   d_pfc <- fit_pfc$rank
#   directions_pfc <- svd(fit_pfc$Beta)$u[,1:d_pfc, drop=FALSE]
#   nz_pfc <- nz_func(directions_pfc)
#   if(length(nz_pfc)==1 && is.na(nz_pfc[1])){
#     s_pfc <- NA
#   }else{
#     s_pfc <- length(nz_pfc)
#   }
# 
#   fit_lassosir <- LassoSIR(x, y, H = 5, nfolds = 5, choosing.d = 'automatic')
#   d_lassosir <- fit_lassosir$no.dim
#   directions_lassosir <- fit_lassosir$beta
#   nz_lassosir <- nz_func(directions_lassosir)
#   if(length(nz_lassosir)==1 && is.na(nz_lassosir[1])){
#     s_lassosir <- NA
#   }else{
#     s_lassosir <- length(nz_lassosir)
#   }
# 
#   fit_lasso <- lasso_func(x, y, nfolds = 5)[-1,1,drop=FALSE]
#   directions_lasso <- unname(fit_lasso)
#   nz_lasso <- nz_func(directions_lasso)
#   if(length(nz_lasso)==1 && is.na(nz_lasso[1])){
#     s_lasso <- NA
#   }else{
#     s_lasso <- length(nz_lasso)
#   }
# 
#   fit_rifle <- rifle_func(x, y, k=15, type = 'sir')
#   directions_rifle <- fit_rifle
#   nz_rifle <- nz_func(directions_rifle)
#   if(length(nz_rifle)==1 && is.na(nz_rifle[1])){
#     s_rifle <- NA
#   }else{
#     s_rifle <- length(nz_rifle)
#   }
# 
#   output <- list(rank = c(d_sir, d_intra, d_pfc, d_lassosir, 1, 1), s = c(s_sir, s_intra, s_pfc, s_lassosir, s_lasso, s_rifle),
#                  nz = list(nz_sir, nz_intra, nz_pfc, nz_lassosir, nz_lasso, nz_rifle),
#                  directions = list(directions_sir, directions_intra, directions_pfc, directions_lassosir, directions_lasso, directions_rifle))
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
# output <- lapply(seq_len(times), function(i){
# # output <- mclapply(seq_len(times), function(i){
#   cat('Time', i, '\n')
#   index <- samples[[i]]
#   boot_x <- x[index,]
#   boot_y <- y[index]
#   boot_output <- output_func(boot_x, boot_y)
# 
#   directions <- boot_output$directions
#   dist <- sapply(1:length(directions), function(i){
#     subspace(true_output$directions[[i]], directions[[i]])
#   })
# 
#   list(rank=boot_output$rank, s = boot_output$s, dist = unname(dist), nz = boot_output$nz)
# })
# 
# boot_rank <- do.call(rbind, lapply(output, '[[', 1))
# boot_s <- do.call(rbind, lapply(output, '[[', 2))
# boot_dist <- do.call(rbind, lapply(output, '[[', 3))
# boot_nz <- lapply(output, '[[', 4)
# 
# # combine bootstrap statistics with true statistics
# boot_rank <- rbind(true_output$rank, boot_rank)
# boot_s <- rbind(true_output$s, boot_s)
# boot_nz <- c(list(true_output$nz), boot_nz)
# 
# 
# # output <- do.call(rbind, prediction)
# write.table(boot_rank, "/Users/cengjing/Desktop/test1")
# write.table(boot_s, "/Users/cengjing/Desktop/test2")
# write.table(boot_dist, "/Users/cengjing/Desktop/test3")
# save(boot_nz, file = '/Users/cengjing/Desktop/test4')

