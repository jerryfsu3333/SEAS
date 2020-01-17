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

load("../data/Ramaswamy.RData")

train <- Ramaswamy$Train
test <- Ramaswamy$Test

train_x <- train$x
train_y <- train$y

test_x <- test$x
test_y <- test$y

y_dummy <- sapply(unique(train_y), function(i){
  as.numeric(train_y == i)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(1)
# times <- 5
# train_size <- 0.8

### Variable screening (first 200 variables) ###
dist_cor <- sapply(seq_len(dim(train_x)[2]), function(i){
  dcor(y_dummy, train_x[,i])
})
ord <- order(dist_cor, decreasing = TRUE)[1:4000]
train_x <- train_x[,ord, drop=FALSE]
test_x <- test_x[,ord, drop=FALSE]

### Check the distribution of preditors from training and testing dataset ###
for (i in 1:5){
  data <- data.frame(x = c(train_x[,i], test_x[,i]), class = c(rep(1, dim(train_x)[1]), rep(2, dim(test_x)[1])))
  g <- ggplot(data, aes(x = x, color = factor(class), fill = factor(class))) +
    geom_density(alpha=0.3)
  print(g)
}

######## Prediction ##########

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
# err_list <- lapply(seq_len(times), function(i){
  # cat(c('Time', i, '\n'))
  
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
  d_sir <- fit_sir$rank
  nz_sir <- nz_func(directions_sir)
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
  d_lassosir <- fit_lassosir$no.dim
  nz_lassosir <- nz_func(directions_lassosir)
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
  d_rifle <- 1
  nz_rifle <- nz_func(directions_rifle)
}
  
obj.cv <- cv.msda(train_x, train_y, nfolds = 5)
lambda.min<-obj.cv$lambda.min
id.min<-which(obj.cv$lambda==lambda.min)
prediction <- predict(obj.cv$msda.fit, test_x)[,id.min]
err_msda <- mean(prediction != test_y)
  
cv.out <- PenalizedLDA.cv(train_x, train_y, lambdas=c(1e-4,1e-3,1e-2,.1,1,10), nfold = 5)
tmp <- PenalizedLDA(train_x, train_y, xte = test_x, lambda=cv.out$bestlambda, K=cv.out$bestK)$ypred
prediction <- tmp[,dim(tmp)[2]]
err_fda <- mean(prediction != test_y)
  
err_list <- list(err_sir = err_sir, err_lassosir = err_lassosir, err_rifle = err_rifle, err_msda = err_msda, err_fda = err_fda, d = c(d_sir, d_lassosir, d_rifle), nz = c(nz_sir, nz_lassosir, nz_rifle))
# }, mc.cores=16)
# })

save(err_list, file="")
