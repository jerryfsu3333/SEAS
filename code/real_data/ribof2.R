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

data <- read.csv("../data/riboflavinv100.csv")
sample_name <- names(data)[-1]
## Remove the column1 - the name of each gene
gene_name <- data[,1]
data <- data[,-1]
data <- t(data)

## The first column "q_RIBFLV" is the response
y <- unname(data[,1])
x <- unname(data[,-1])


####### Estimation consistency ##########
output_func <- function(x, y){
  
  y <- scale(y)
  x <- scale(x)
  
  # fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.3, length.out = 10), lam2_fac = seq(0.1,4, length.out = 10),
  #                    gamma = c(0.1), nfolds = 3, type = 'sir', plot=TRUE)
  fit_sir <- ssdr.cv(x, y, lam1_fac = seq(2,0.3, length.out = 10), lam2_fac = seq(0.1,4, length.out = 10),
                     gamma = c(0.1), nfolds = 3, type = 'sir')
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
  
  
  # fit_intra <- ssdr.cv(x, y, lam1_fac = seq(2,0.3, length.out = 10), lam2_fac = seq(0.1,2, length.out = 10),
  #                      gamma = c(0.1), nfolds = 3, type = 'intra', plot=TRUE)
  fit_intra <- ssdr.cv(x, y, lam1_fac = seq(2,0.3, length.out = 10), lam2_fac = seq(0.1,2, length.out = 10),
                       gamma = c(0.1), nfolds = 3, type = 'intra')
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
  
  # fit_pfc <- ssdr.cv(x, y, lam1_fac = seq(1.5,0.3, length.out = 10), lam2_fac = seq(0.1,2, length.out = 10),
  #                    gamma = c(0.1), nfolds = 3, type = 'pfc', cut_y = TRUE, plot = TRUE, maxit_outer = 1e+4)
  fit_pfc <- ssdr.cv(x, y, lam1_fac = seq(1.5,0.3, length.out = 10), lam2_fac = seq(0.1,2, length.out = 10),
                     gamma = c(0.01), nfolds = 3, type = 'pfc', cut_y = FALSE, maxit_outer = 1e+4)
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
  
  fit_lassosir <- LassoSIR(x, y, H = 5, nfolds = 3, choosing.d = 'automatic')
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