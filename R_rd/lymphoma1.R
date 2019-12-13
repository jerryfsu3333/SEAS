rm(list = ls())
library(parallel)
library(MASS)
library(methods)
library(glmnet)
library(energy)
library(caret)
library(R.matlab)
library(ggplot2)
library(dplyr)
library(reshape2)
#### Packages for methods ####
library(msda)
library(rifle)
library(LassoSIR)
library(np)
library(nnet)
library(e1071)
library(randomForest)
library(sparseLDA)
library(penalizedLDA)
####
setwd("~/Documents/GitHub/seas/R/")
source("models.R")
source("utility.R")
source("ssdr_utility.R")
source("ssdr_func.R")
source("rifle_func.R")
source("lasso_func.R")
source("CovSIR.R")

### Lymphoma data description ###
# The lymphoma dataset consists of 42 samples of diffuse large B-cell lymphoma (DLBCL), 9 samples of follicular lymphoma (FL), and 11 samples of chronic lymphocytic leukemia (CLL). DBLCL, FL, and CLL classes are coded in 0, 1, and 2, respectively, in y vector. Matrix x is gene expression data and arrays were normalized, imputed, log transformed, and standardized to zero mean and unit variance across genes 

library(spls)
data("lymphoma")

x <- lymphoma$x
y <- lymphoma$y
y <- y + 1

y_dummy <- sapply(unique(y), function(i){
  as.numeric(y == i)
})


######## Prediction ##########
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
times <- 5
train_size <- 0.8

# Screen variables
dist_cor <- sapply(seq_len(dim(x)[2]), function(i){
  dcor(y_dummy, x[,i])
})
ord <- order(dist_cor, decreasing = TRUE)[1:200]
x <- x[,ord, drop=FALSE]


err_func <- function(train_x, train_y, test_x, test_y){
  
  # m <- apply(x, 2, mean)
  # se <- apply(x, 2, sd)
  # x <- sweep(sweep(x, 2, m), 2, se, '/')
  # newx <- sweep(sweep(newx, 2, m), 2, se, '/')
  
  train_x <- data.frame(train_x)
  test_x <- data.frame(test_x)
  colnames(train_x) <- paste0('X', seq_len(ncol(train_x)))
  colnames(test_x) <- paste0('X', seq_len(ncol(test_x)))
  
  data <- cbind(y = factor(train_y), train_x)
  
  model <- multinom(y~., data = data)
  prediction <- predict(model, newdata=test_x, type='class')
  log_err <- mean(prediction != test_y)
  
  
  model <- svm(y ~., data = data, scale = FALSE, kernel = 'linear', cost = 1, probability = TRUE)
  prediction <- predict(model, newdata = test_x, probability = TRUE)
  svm_err <- mean(prediction != test_y)
  
  model <- lda(y~., data = data)
  prediction <- predict(model, newdata = test_x)$class
  lda_err <- mean(prediction != test_y)
  
  model <- randomForest(y~., data = data)
  prediction <- predict(model, newdata = test_x)
  rf_err <- mean(prediction != test_y)
  
  c(log_err = log_err, svm_err = svm_err, lda_err = lda_err, rf_err = rf_err)
}



# err_list <- mclapply(seq_len(NROW(x)), function(i){
# err_list <- mclapply(seq_len(times), function(i){
err_list <- lapply(seq_len(times), function(i){
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
  
  fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), categorical=TRUE, nfolds=5, plot = TRUE, type = 'sir')
  # fit_sir <- ssdr.cv(train_x, train_y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), categorical=TRUE, nfolds=5, type = 'sir')
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

  obj.cv <- cv.msda(train_x, train_y, nfolds = 5)
  lambda.min<-obj.cv$lambda.min
  id.min<-which(obj.cv$lambda==lambda.min)
  fit_msda <- as.matrix(obj.cv$msda.fit$theta[[id.min]])
  if(!is.numeric(fit_msda)){
    print('A NULL matrix is returned (lassosir).')
    err_msda <- NA
  }else{
    directions_msda <- fit_msda
    new_train <- train_x %*% directions_msda
    new_test <- test_x %*% directions_msda
    err_msda <- err_func(new_train, train_y, new_test, test_y)
  }

  fit_sos <- sda(train_x, factor(train_y))
  if(!is.numeric(fit_sos)){
    print('A NULL matrix is returned (sos).')
    err_sos <- NA
  }else{
    directions_sos <- fit_sos
    new_train <- train_x %*% directions_sos
    new_test <- test_x %*% directions_sos
    err_sos <- err_func(new_train, train_y, new_test, test_y)
  }
  
  cv.out <- PenalizedLDA.cv(train_x, train_y, lambdas=c(1e-4,1e-3,1e-2,.1,1,10), nfold = 5)
  fit_fisher <- PenalizedLDA(train_x, train_y, lambda=cv.out$bestlambda, K=cv.out$bestK)$discrim
  if(!is.numeric(fit_fisher)){
    print('A NULL matrix is returned (rifle).')
    err_fisher <- NA
  }else{
    directions_fisher <- fit_fisher
    new_train <- train_x %*% directions_fisher
    new_test <- test_x %*% directions_fisher
    err_fisher <- err_func(new_train, train_y, new_test, test_y)
  }
  
  list(err_sir = err_sir, err_lassosir = err_lassosir, err_rifle = err_rifle, err_msda = err_msda, err_fisher = err_fisher)
# }, mc.cores=16)
})

save(err_list, file="")

