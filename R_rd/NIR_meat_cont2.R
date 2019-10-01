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

# Pork (y=1) only
data <- data[data[,1] == 1,]
y <- data[,2]
# x <- data[,-c(1,2)]
x <- data[,-c(1,2,3,4,5)]


####### Estimation consistency ##########
output_func <- function(x, y){
  
  # y <- scale(y)
  # x <- scale(x)
  
  y <- log(y)
  x <- log(x)
  
  fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.7, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                     gamma = c(0.01), nfolds = 3, type = 'sir', plot=TRUE, pmax=100)
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
  
  
  fit_intra <- ssdr.cv(x, y, lam1_fac = seq(2,0.7, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                       gamma = c(0.01), nfolds = 3, type = 'intra', plot=TRUE, pmax=100)
  if(!is.numeric(fit_intra$Beta)){
    print('A NULL matrix is returned (intra).')
    d_intra <- NA
    directions_intra <- NA
    nz_intra <- NA
    s_intra <- NA
  }else if(fit_intra$rank==0){
    print('rank is 0 (intra).')
    d_intra <- NA
    directions_intra <- NA
    nz_intra <- NA
    s_intra <- NA
  }else{
    d_intra <- fit_intra$rank
    directions_intra <- fit_intra$Beta
    nz_intra <- nz_func(directions_intra)
    s_intra <- length(nz_intra)
  }
  
  
  fit_pfc <- ssdr.cv(x, y, lam1_fac = seq(1.2,0.5, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10),
                     gamma = c(0.01), nfolds = 3, type = 'pfc', cut_y = TRUE, plot = TRUE, maxit_outer = 1e+4, pmax=100)
  if(!is.numeric(fit_pfc$Beta)){
    print('A NULL matrix is returned (pfc).')
    d_pfc <- NA
    directions_pfc <- NA
    nz_pfc <- NA
    s_pfc <- NA
  }else if(fit_pfc$rank==0){
    print('rank is 0 (pfc).')
    d_pfc <- NA
    directions_pfc <- NA
    nz_pfc <- NA
    s_pfc <- NA
  }else{
    d_pfc <- fit_pfc$rank
    directions_pfc <- fit_pfc$Beta
    nz_pfc <- nz_func(directions_pfc)
    s_pfc <- length(nz_pfc)
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
  
  # fit_lasso <- lasso_func(x, y, nfolds = 5)[-1,1,drop=FALSE]
  # directions_lasso <- unname(fit_lasso)
  # nz_lasso <- nz_func(directions_lasso)
  # if(length(nz_lasso)==1 && is.na(nz_lasso[1])){
  #   s_lasso <- NA
  # }else{
  #   s_lasso <- length(nz_lasso)
  # }
  #
  # fit_rifle <- rifle_func(x, y, k=15, type = 'sir')
  # directions_rifle <- fit_rifle
  # nz_rifle <- nz_func(directions_rifle)
  # if(length(nz_rifle)==1 && is.na(nz_rifle[1])){
  #   s_rifle <- NA
  # }else{
  #   s_rifle <- length(nz_rifle)
  # }
  
  # output <- list(rank = c(d_sir, d_intra, d_pfc, d_lassosir, 1, 1), s = c(s_sir, s_intra, s_pfc, s_lassosir, s_lasso, s_rifle),
  #                nz = list(nz_sir, nz_intra, nz_pfc, nz_lassosir, nz_lasso, nz_rifle),
  #                directions = list(directions_sir, directions_intra, directions_pfc, directions_lassosir, directions_lasso, directions_rifle))
  output <- list(rank = c(d_sir, d_intra, d_pfc, d_lassosir), s = c(s_sir, s_intra, s_pfc, s_lassosir),
                 nz = list(nz_sir, nz_intra, nz_pfc, nz_lassosir),
                 directions = list(directions_sir, directions_intra, directions_pfc, directions_lassosir))
  output
}

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# Full dataset
true_output <- output_func(x,y)

# Bootstrap samples
times <- 100
# samples <- createResample(y, times = times)

output <- mclapply(seq_len(times), function(i){
  # output <- lapply(seq_len(times), function(i){
  cat('Time', i, '\n')
  index <- sample(1:length(y), length(y), replace = TRUE)
  # index <- samples[[i]]
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
  
  list(rank=boot_output$rank, s = boot_output$s, dist = unname(dist), nz = boot_output$nz)
})

# save(true_output, file = '')
save(output, file = '')