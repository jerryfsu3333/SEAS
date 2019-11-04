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

data <- readMat('../data/NIR.mat')$data
y <- factor(data[,1])
# x <- data[,-1]
x <- data[,-c(1,5)]
x <- log(x)


# ############# Visualization ####################
# hist_plot <- function(x, y, title){
#   if(!is.factor(y)){y <- factor(y)}
#   if(!is.null(dim(x))){x <- drop(x)}
#   df <- data.frame(component = x, class = y)
#   means <- sapply(unique(df$class), function(i){
#     mean(df$component[df$class==i])
#   })
#   means_df <- data.frame(means = means, class=unique(df$class))
#   g <- ggplot(df, aes(x=component, colour=class, fill=class)) +
#     geom_histogram(aes(y=..density..), bins = 50, position = 'identity', alpha=0.5) +
#     geom_density(alpha=0.3) +
#     geom_vline(data = means_df, aes(xintercept=means, color=class), linetype='dashed')+
#     theme(
#       plot.title = element_text(size=16),
#       axis.title.x = element_text(size=16),
#       axis.title.y = element_text(size=16),
#       legend.position = 'none'
#     )+
#     labs(title = title)
#   g
# }
# 
# x <- scale(x)
# fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), categorical=TRUE, type = 'sir')
# directions_sir <- fit_sir$Beta
# x_new <- as.matrix(x) %*% directions_sir
# hist_plot(x_new, y, 'SSDR-SIR')
# 
# #### LassoSIR ####
# fit_lassosir <- LassoSIR(x, y, H = 5, nfolds = 5, choosing.d = 'automatic')
# directions_Lassosir <- fit_lassosir$beta
# x_new_Lassosir <- as.matrix(x) %*% directions_Lassosir
# hist_plot(x_new_Lassosir, y, 'Lasso-SIR')
# 
# 
# #### CovSIR ####
# # N <- dim(x)[1]
# # p <- dim(x)[2]
# # CovSIR_fit <- CovSIR(x, y, Ks = 1:3, lambdas = seq(0.2,2,by=0.5)*sqrt(log(p)/N))
# # d_CovSIR <- CovSIR_fit$r
# # directions_Covsir <- CovSIR_fit$mat
# # sum(directions_Covsir != 0)
# # x_new_Covsir <- as.matrix(x) %*% directions_Covsir
# # plot(x_new_Covsir[,1], col=y, xlab = 'Index', ylab = 'Component 1')
# # abline(h=0, lty = 'dashed')
# 
# ### Lasso ####
# fit_lasso <- lasso_func(x, y, family='binomial', type.measure='class')[-1,1,drop=FALSE] # the first is zero intercept
# directions_lasso <- fit_lasso
# x_new_lasso <- as.matrix(x) %*% directions_lasso
# hist_plot(x_new_lasso, y, 'Lasso')
# 
# ### Rifle ####
# fit_rifle <- rifle_func(x, y, categorical = TRUE, k=15, type = 'sir')
# directions_rifle <- fit_rifle
# x_new_rifle <- as.matrix(x) %*% directions_rifle
# hist_plot(x_new_rifle, y, 'Rifle-SIR')

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
  
  model <- glm(y~., data = data, family = binomial())
  prediction <- ifelse(predict(model, newdata=newx, type='response') > 0.5, 2, 1)
  log_err <- mean(prediction != newy)
  
  
  model <- svm(y ~., data = data, scale = FALSE, kernel = 'linear', cost = 1)
  prediction <- predict(model, newdata = newx)
  svm_err <- mean(prediction != newy)
  
  model <- lda(y~., data = data)
  prediction <- predict(model, newdata = newx)$class
  lda_err <- mean(prediction != newy)
  
  model <- randomForest(y~., data = data)
  prediction <- predict(model, newdata = newx)
  rf_err <- mean(prediction != newy)
  
  c(log_err = log_err, svm_err = svm_err, lda_err = lda_err, rf_err = rf_err)
}


# err_list <- mclapply(seq_len(NROW(x)), function(i){
err_list <- mclapply(seq_len(times), function(i){
# err_list <- lapply(seq_len(times), function(i){
  cat(c('Time', i, '\n'))
  ## stratified K-fold, 20% testing dataset
  class <- unique(y)
  index <- c()
  for(k in class){
    index <- c(index, sample(which(y==k), train_size*length(y[y==k]), replace = FALSE))
  }
  
  ## For LOO
  # index <- -i
  
  
  train_x <- x[index,,drop=FALSE]
  train_y <- y[index]
  test_x <- x[-index,,drop=FALSE]
  test_y <- y[-index]
  
  m <- apply(train_x, 2, mean)
  se <- apply(train_x, 2, sd)
  train_x <- sweep(sweep(train_x, 2, m), 2, se, '/')
  test_x <- sweep(sweep(test_x, 2, m), 2, se, '/')

  # fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
  # categorical=TRUE, plot = TRUE, type = 'sir')
  fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                     categorical=TRUE, type = 'sir')
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


  fit_lasso <- lasso_func(train_x, train_y, family='binomial', nfolds=5, type.measure='class')[-1,1,drop=FALSE]
  if(!is.numeric(fit_lasso)){
    print('A NULL matrix is returned (lasso).')
    err_lasso <- NA
  }else{
    directions_lasso <- fit_lasso
    new_train <- train_x %*% directions_lasso
    new_test <- test_x %*% directions_lasso
    err_lasso <- err_func(new_train, train_y, new_test, test_y)
  }

  
  fit_rifle <- rifle_func(train_x, train_y, categorical = TRUE, k=15, type = 'sir')
  if(!is.numeric(fit_rifle)){
    print('A NULL matrix is returned (rifle).')
    err_rifle <- NA
  }else{
    directions_rifle <- fit_rifle
    new_train <- train_x %*% directions_rifle
    new_test <- test_x %*% directions_rifle
    err_rifle <- err_func(new_train, train_y, new_test, test_y)
  }
  
  list(err_sir = err_sir, err_lassosir = err_lassosir, err_lasso = err_lasso, err_rifle = err_rifle)
}, mc.cores=16)
# })
  
save(err_list, file="")


