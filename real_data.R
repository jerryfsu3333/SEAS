rm(list = ls())
library(msda)
library(MASS)
library(methods)
library(glmnet)
library(energy)
library(caret)
library(LassoSIR)
library(mlbench)
library(R.matlab)
library(nnet)
library(rifle)
# library(plot3D)
# library(plot3Drgl)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/my_msda.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/CovSIR.R")
source("/Users/cengjing/Documents/GitHub/ssdr/rifle_func.R")
source("/Users/cengjing/Documents/GitHub/ssdr/lasso_func.R")
# ####################  Wine data #############################
# data <- read.csv('/Users/cengjing/Documents/GitHub/ssdr/wine.data', header = FALSE, col.names =
#                 c('class', 'Alc', 'Malic', 'Ash', 'Alka', 'Mag', 'Tot_ph', 'Fla', 'NonFla', 'Proant', 'Col', 'Hue', 'OD', 'Proline'))
# 
# y <- data[,1]
# x <- as.matrix(data[,-1])
# x <- scale(x)
# 
# # fit_lda <- MASS::lda(x,y)
# # pred <- predict(fit_lda, x)$class
# # pred <- as.numeric(levels(pred)[pred])
# # sum(y == pred)/length(y)
# 
# fit <- ssdr.cv(x, y, lam1_fac = seq(1,0.01, length.out = 10), categorical=TRUE, type = 'sir')
# d <- fit$rank
# directions <- svd(fit$mat)$u[,1:d, drop=FALSE]
# x_new <- as.matrix(x) %*% directions
# plot(x_new[,1], x_new[,2], col=y, xlab = 'Component 1', ylab = 'Component 2')
# #
# # fit_lda <- MASS::lda(x_new,y)
# # pred <- predict(fit_lda, x_new)$class
# # pred <- as.numeric(levels(pred)[pred])
# # sum(y == pred)/length(y)
# # 
# # 
# # LassoSIR_fit <- LassoSIR(as.matrix(x), y, categorical = TRUE, choosing.d = 'manual')
# # d <- LassoSIR_fit$no.dim
# # direction <- LassoSIR_fit$beta
# # x_new <- as.matrix(x) %*% direction
# # # x_new <- scale(x_new, center = FALSE)
# # plot(x_new[,1], x_new[,2], col=y)
# # 
# # fit_lda <- MASS::lda(x_new,y)
# # pred <- predict(fit_lda, x_new)$class
# # pred <- as.numeric(levels(pred)[pred])
# # sum(y == pred)/length(y)

# ####################  Boston housing data #############################
# data("BostonHousing2")
# data <- BostonHousing2[,c(5,7:19)]
# data <- data[data$crim <= 3.2,]
# y <- data$medv
# x <- data[,-1]
# names <- colnames(x)
# x[,4] <- as.numeric(levels(x[,4])[x[,4]])
# x <- scale(x)
# y <- drop(y)
# 
# fit <- ssdr.cv(x, y, lam1_fac = seq(1,0.01, length.out = 15), categorical=FALSE, type = 'sir')
# d <- fit$rank
# directions <- svd(fit$mat)$u[,1:d, drop=FALSE]
# x_new <- as.matrix(x) %*% directions
# plot(x_new[,1], y)

####################  NIR meat data #############################
set.seed(1)
data <- readMat('~/Documents/GitHub/ssdr/Real_dataset/NIR.mat')$data
y <- data[,1]
x <- data[,-1]
x <- scale(x)


LassoSIR_fit <- LassoSIR(x, y, H = 5, nfolds = 5, choosing.d = 'automatic')
d_LassoSIR <- LassoSIR_fit$no.dim
directions_Lassosir <- LassoSIR_fit$beta
sum(directions_Lassosir != 0)
x_new_Lassosir <- as.matrix(x) %*% directions_Lassosir
plot(x_new_Lassosir[,1], col=y, xlab = 'Index', ylab = 'Component 1')
abline(h=0, lty = 'dashed')

# N <- dim(x)[1]
# p <- dim(x)[2]
# CovSIR_fit <- CovSIR(x, y, Ks = 1:3, lambdas = seq(0.2,2,by=0.5)*sqrt(log(p)/N))
# d_CovSIR <- CovSIR_fit$r
# directions_Covsir <- CovSIR_fit$mat
# sum(directions_Covsir != 0)
# x_new_Covsir <- as.matrix(x) %*% directions_Covsir
# plot(x_new_Covsir[,1], col=y, xlab = 'Index', ylab = 'Component 1')
# abline(h=0, lty = 'dashed')

# lasso_fit <- lasso_func(x, y, family='binomial', type.measure='class')[-1,1,drop=FALSE] # the first is zero intercept
# directions_lasso <- lasso_fit
# sum(directions_lasso != 0)
# x_new_lasso <- as.matrix(x) %*% directions_lasso
# plot(x_new_lasso[,1], col=y, xlab = 'Index', ylab = 'Component 1')
# abline(h=0, lty = 'dashed')

# rifle_fit <- rifle_func(x, y, categorical = TRUE, k=15, type = 'sir')
# directions_rifle <- rifle_fit
# sum(directions_rifle != 0)
# x_new_rifle <- as.matrix(x) %*% directions_rifle
# plot(x_new_rifle[,1], col=y, xlab = 'Index', ylab = 'Component 1')
# abline(h=0, lty = 'dashed')

# fit <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), categorical=TRUE, type = 'intra')
# # fit <- ssdr.cv(x, y, lam1_fac = seq(1.5,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), categorical=TRUE,type = 'pfc', cut_y = FALSE, maxit_outer = 1e+4)
# d <- fit$rank
# directions <- svd(fit$Beta)$u[,1:d, drop=FALSE]
# sum(directions != 0)
# x_new <- as.matrix(x) %*% directions
# plot(x_new[,1], col=y, xlab = 'Index', ylab = 'Component 1')
# abline(h=0, lty = 'dashed')

# # ####################  NIR corn data #############################
# data <- readMat('~/Documents/GitHub/ssdr/Real_dataset/corn.mat')
# data <- cbind(data$propvals$data, data$m5spec$data, data$mp5spec$data, data$mp6spec$data)
# y <- data[,1]
# x <- data[,-1]
# x <- scale(x)
# 
# fit <- ssdr.cv(x, y, lam1_fac = seq(2,0.2, length.out = 10), lam2_fac = seq(0.001,0.2, length.out = 10), categorical=TRUE, type = 'sir')
# d <- fit$rank
# directions <- svd(fit$Beta)$u[,1:d, drop=FALSE]
# x_new <- as.matrix(x) %*% directions
# plot(x_new[,1], col=y, xlab = 'Index', ylab = 'Component 1')

# # ####################  GEO data #############################
# library("GEOquery")
# dat<-getGEO("GDS3709")
# 
# x<-Table(dat)
# x<-x[,-c(1,2)]
# x<-t(x)
# x<-as.numeric(x)
# x<-matrix(x,nrow=79)
# x<-scale(log(x))
# y<-c(rep(1,19),rep(2,20),rep(3,20),rep(4,20))
# K<-4
# 
# # Select the first 1000 variables via distance correlation
# dist_cor <- sapply(seq_len(dim(x)[2]), function(i){
#   dcor(y,x[,i])
# })
# 
# ord <- order(dist_cor, decreasing = TRUE)[1:1500]
# x <- x[,ord]

############################################
set.seed(1)

# save(x, y, file = '~/Documents/GitHub/ssdr/Real_dataset/GDS3709_sel')
load('~/Documents/GitHub/ssdr/Real_dataset/GDS3709_sel')

# fit <- ssdr.cv(x, y, lam1_fac = seq(1.5,0.9, length.out = 15), lam2_fac = seq(0.001,0.1, length.out = 10), categorical=TRUE, type = 'sir')
# fit <- ssdr.cv(x, y, lam1_fac = seq(2,1.2, length.out = 15), lam2_fac = seq(0.001,0.1, length.out = 10), categorical=TRUE, type = 'intra')
fit <- ssdr.cv(x, y, lam1_fac = seq(1.5,0.9, length.out = 15), lam2_fac = seq(1e-4,1e-2, length.out = 10), categorical=TRUE, type = 'pfc', maxit_outer = 1e+4)
d <- fit$rank
directions <- svd(fit$Beta)$u[,1:d, drop=FALSE]

# Sparsity
sparsity <- sapply(seq_len(dim(directions)[1]), function(i){
  any(directions[i,] != 0)
})
sum(sparsity)

x_new <- as.matrix(x) %*% directions
plot(x_new[,1], x_new[,2], col=y, xlab = 'Component 1', ylab = 'Component 2', ylim = c())
legend('topright', legend = 1:4, col = 1:4, pch = 1, cex = 0.6, pt.cex = 1)
# plot(x_new[,1], x_new[,3], col=y, xlab = 'Component 1', ylab = 'Component 3')
# plot(x_new[,2], x_new[,3], col=y, xlab = 'Component 2', ylab = 'Component 3')
# scatter3Drgl(x_new[,1], x_new[,2], x_new[,3], col = y)
# plotrgl()

# fit_lda <- MASS::lda(x_new,y)
# pred <- predict(fit_lda, x_new)$class
# pred <- as.numeric(levels(pred)[pred])
# sum(y == pred)/length(y)

fit_lda <- multinom(y~x_new)
pred <- predict(fit_lda, x_new)
pred <- as.numeric(levels(pred)[pred])
sum(y == pred)/length(y)

