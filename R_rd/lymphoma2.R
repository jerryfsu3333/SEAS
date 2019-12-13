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
####s
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


RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# Screen variables
dist_cor <- sapply(seq_len(dim(x)[2]), function(i){
  dcor(y_dummy, x[,i])
})
ord <- order(dist_cor, decreasing = TRUE)[1:200]
x <- x[,ord, drop=FALSE]



######## Estimation consistency ##########
output_func <- function(x, y){
  
  
  # fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), categorical=TRUE, plot = TRUE, nfolds=5, type = 'sir')
  fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), categorical=TRUE, nfolds=5, type = 'sir')
  if(!is.numeric(fit_sir$Beta)){
    print('A NULL matrix is returned (sir).')
    d_sir <- NA
    directions_sir <- NA
    nz_sir <- NA
    s_sir <- NA
  }else if(fit_sir$rank==0){
    print('rank is 0 (sir).')
    d_sir <- NA
    directions_sir <- NA
    nz_sir <- NA
    s_sir <- NA
  }else{
    d_sir <- fit_sir$rank
    directions_sir <- fit_sir$Beta
    nz_sir <- nz_func(directions_sir)
    s_sir <- length(nz_sir)
  }


  fit_lassosir <- LassoSIR(x, y, H = 5, nfolds = 5, choosing.d = 'automatic')
  if(!is.numeric(fit_lassosir$beta)){
    print('A NULL matrix is returned (lassosir).')
    d_lassosir <- NA
    directions_lassosir <- NA
    nz_lassosir <- NA
    s_lassosir <- NA
  }else if(fit_lassosir$no.dim==0){
    print('rank is 0 (lassosir).')
    d_lassosir <- NA
    directions_lassosir <- NA
    nz_lassosir <- NA
    s_lassosir <- NA
  }else{
    d_lassosir <- fit_lassosir$no.dim
    directions_lassosir <- fit_lassosir$beta
    nz_lassosir <- nz_func(directions_lassosir)
    s_lassosir <- length(nz_lassosir)
  }

  # fit_lasso <- lasso_func(x, y, family='binomial', nfolds=5, type.measure='class')[-1,1,drop=FALSE]
  # if(!is.numeric(fit_lasso)){
  #   print('A NULL matrix is returned (lassosir).')
  #   d_lasso <- NA
  #   directions_lasso <- NA
  #   nz_lasso <- NA
  #   s_lasso <- NA
  # }else{
  #   d_lasso <- 1
  #   directions_lasso <- fit_lasso
  #   nz_lasso <- nz_func(directions_lasso)
  #   s_lasso <- length(nz_lasso)
  # }


  fit_rifle <- rifle_func(x, y, categorical = TRUE, k=15, type = 'sir')
  if(!is.numeric(fit_rifle)){
    print('A NULL matrix is returned (lassosir).')
    d_rifle <- NA
    directions_rifle <- NA
    nz_rifle <- NA
    s_rifle <- NA
  }else{
    d_rifle <- 1
    directions_rifle <- fit_rifle
    nz_rifle <- nz_func(directions_rifle)
    s_rifle <- length(nz_rifle)
  }

  obj.cv <- cv.msda(x, y, nfolds = 5)
  lambda.min<-obj.cv$lambda.min
  id.min<-which(obj.cv$lambda==lambda.min)
  fit_msda <- as.matrix(obj.cv$msda.fit$theta[[id.min]])
  if(!is.numeric(fit_msda)){
    print('A NULL matrix is returned (msda).')
    d_msda <- NA
    directions_msda <- NA
    nz_msda <- NA
    s_msda <- NA
  }else if(dim(fit_msda)[2]==0){
    print('rank is 0 (msda).')
    d_msda <- NA
    directions_msda <- NA
    nz_msda <- NA
    s_msda <- NA
  }else{
    d_msda <- dim(fit_msda)[2]
    directions_msda <- fit_msda
    nz_msda <- nz_func(directions_msda)
    s_msda <- length(nz_msda)
  }
  
  cv.out <- PenalizedLDA.cv(x, y, lambdas=c(1e-4,1e-3,1e-2,.1,1,10,100), nfold = 5)
  fit_fisher <- PenalizedLDA(x, y, lambda=cv.out$bestlambda, K=cv.out$bestK)$discrim
  if(!is.numeric(fit_fisher)){
    print('A NULL matrix is returned (fisher).')
    d_fisher <- NA
    directions_fisher <- NA
    nz_fisher <- NA
    s_fisher <- NA
  }else if(dim(fit_fisher)[2]==0){
    print('rank is 0 (fisher).')
    d_fisher <- NA
    directions_fisher <- NA
    nz_fisher <- NA
    s_fisher <- NA
  }else{
    d_fisher <- dim(fit_fisher)[2]
    directions_fisher <- fit_fisher
    nz_fisher <- nz_func(directions_fisher)
    s_fisher <- length(nz_fisher)
  }
  
  output <- list(rank = c(d_sir, d_lassosir, d_rifle, d_msda, d_fisher), s = c(s_sir, s_lassosir, s_rifle, s_msda, s_fisher),
                 nz = list(nz_sir, nz_lassosir, nz_rifle, nz_msda, nz_fisher),
                 directions = list(directions_sir, directions_lassosir, directions_rifle, directions_msda, directions_fisher))
  output
}

# Full dataset
true_output <- output_func(x,y)

# Bootstrap samples
times <- 100
# output <- mclapply(seq_len(times), function(i){
output <- lapply(seq_len(times), function(i){
  cat('Time', i, '\n')
  index <- sample(1:length(y), length(y), replace = TRUE)
  boot_x <- x[index,]
  boot_y <- y[index]
  boot_output <- output_func(boot_x, boot_y)
  
  directions <- boot_output$directions
  dist <- sapply(1:length(directions), function(i){
    if(!is.numeric(true_output$directions[[i]]) | !is.numeric(directions[[i]])){
      NA
    }else{
      subspace(true_output$directions[[i]], directions[[i]])
    }
  })
  list(rank=boot_output$rank, s = boot_output$s, nz = boot_output$nz, dist = unname(dist))
})
# }, mc.cores=16)

save(ord, file='')
save(true_output, file = '')
save(output, file = '')
