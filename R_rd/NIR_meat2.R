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
# source("msda_prep.R")
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
# x <- scale(x)


######## Estimation consistency ##########
output_func <- function(x, y){

  x <- scale(x)
  
  # fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), 
                     # categorical=TRUE, plot = TRUE, type = 'sir')
  fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                     categorical=TRUE, type = 'sir')
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
  

  fit_lasso <- lasso_func(x, y, family='binomial', type.measure='class')[-1,1,drop=FALSE]
  if(!is.numeric(fit_lasso)){
    print('A NULL matrix is returned (lassosir).')
    d_lasso <- NA
    directions_lasso <- NA
    nz_lasso <- NA
    s_lasso <- NA
  }else{
    d_lasso <- 1
    directions_lasso <- fit_lasso
    nz_lasso <- nz_func(directions_lasso)
    s_lasso <- length(nz_lasso)
  }
  

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
  

  output <- list(rank = c(d_sir, d_lassosir, d_lasso, d_rifle), s = c(s_sir, s_lassosir, s_lasso, s_rifle),
                 nz = list(nz_sir, nz_lassosir, nz_lasso, nz_rifle),
                 directions = list(directions_sir, directions_lassosir, directions_lasso, directions_rifle))
  output
}

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# Full dataset
true_output <- output_func(x,y)

# Bootstrap samples
times <- 100
output <- mclapply(seq_len(times), function(i){
# output <- lapply(seq_len(times), function(i){
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
# })
}, mc.cores=16)

save(true_output, file = '')
save(output, file = '')
