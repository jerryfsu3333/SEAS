# ####################  GEO data #############################
library("GEOquery")
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
# 
# ###########################################
# set.seed(1)
# 
# save(x, y, file = '~/Documents/GitHub/ssdr/Real_dataset/GDS3709_sel')
# load('~/Documents/GitHub/ssdr/Real_dataset/GDS3709_sel')
# 
# fit <- ssdr.cv(x, y, lam1_fac = seq(1.5,0.9, length.out = 15), lam2_fac = seq(0.001,0.1, length.out = 10), categorical=TRUE, type = 'sir')
# fit <- ssdr.cv(x, y, lam1_fac = seq(2,1.2, length.out = 15), lam2_fac = seq(0.001,0.1, length.out = 10), categorical=TRUE, type = 'intra')
# fit <- ssdr.cv(x, y, lam1_fac = seq(1.5,0.9, length.out = 15), lam2_fac = seq(1e-4,1e-2, length.out = 10), categorical=TRUE, type = 'pfc', maxit_outer = 1e+4)
# d <- fit$rank
# directions <- svd(fit$Beta)$u[,1:d, drop=FALSE]
# 
# # Sparsity
# sparsity <- sapply(seq_len(dim(directions)[1]), function(i){
#   any(directions[i,] != 0)
# })
# sum(sparsity)
# 
# x_new <- as.matrix(x) %*% directions
# plot(x_new[,1], x_new[,2], col=y, xlab = 'Component 1', ylab = 'Component 2', ylim = c())
# legend('topright', legend = 1:4, col = 1:4, pch = 1, cex = 0.6, pt.cex = 1)
# # plot(x_new[,1], x_new[,3], col=y, xlab = 'Component 1', ylab = 'Component 3')
# # plot(x_new[,2], x_new[,3], col=y, xlab = 'Component 2', ylab = 'Component 3')
# # scatter3Drgl(x_new[,1], x_new[,2], x_new[,3], col = y)
# # plotrgl()
# 
# # fit_lda <- MASS::lda(x_new,y)
# # pred <- predict(fit_lda, x_new)$class
# # pred <- as.numeric(levels(pred)[pred])
# # sum(y == pred)/length(y)
# 
# fit <- multinom(y~x_new)
# pred <- predict(fit, x_new)
# pred <- as.numeric(levels(pred)[pred])
# sum(y == pred)/length(y)
# 
# ### LassoSIR ####
# LassoSIR_fit <- LassoSIR(x, y, H = 5, nfolds = 5, choosing.d = 'automatic', categorical = TRUE)
# d_LassoSIR <- LassoSIR_fit$no.dim
# directions_Lassosir <- LassoSIR_fit$beta
# sum(directions_Lassosir != 0)
# x_new_Lassosir <- as.matrix(x) %*% directions_Lassosir
# plot(x_new_Lassosir[,1], col=y, xlab = 'Index', ylab = 'Component 1')
# 
# fit_lassosir <- multinom(y~x_new_Lassosir)
# pred <- predict(fit_lassosir, x_new_Lassosir)
# pred <- as.numeric(levels(pred)[pred])
# sum(y == pred)/length(y)
# 
# # #### CovSIR ####
# # N <- dim(x)[1]
# # p <- dim(x)[2]
# # CovSIR_fit <- CovSIR(x, y, Ks = 1:3, lambdas = seq(0.2,2,by=0.5)*sqrt(log(p)/N))
# # d_CovSIR <- CovSIR_fit$r
# # directions_Covsir <- CovSIR_fit$mat
# # sum(directions_Covsir != 0)
# # x_new_Covsir <- as.matrix(x) %*% directions_Covsir
# # plot(x_new_Covsir[,1], col=y, xlab = 'Index', ylab = 'Component 1')

######################################
dat<-getGEO("GDS3916")

