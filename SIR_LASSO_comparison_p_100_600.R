library(MASS)
library(lars)
library(fastclime)

half <- function(Sigma){
  S <- eigen(Sigma)
  pos <- sqrt(S$val[S$val > 0.0001])
  zer <- rep(0,length(S$val < 0.0001))
  mid <- diag(c(pos,zer), ncol(Sigma), ncol(Sigma))
  S$vec%*%mid%*%t(S$vec)
}

H <- 10

# SIR for linear regression
# this function returns the SIR estimate of the covariance matrix
SIR <- function(Y,X,H,SS.inv.half){
  
  n <- length(Y)
  p <- ncol(X)
  
  Y.ord <- order(Y)
  
  m <- floor(n/H)
  
  start.ind <- cumsum(c(1,(rep(m,H) + c(rep(1, n%%m), rep(0, H - n%%m)))[-H]))
  end.ind <- cumsum(rep(m,H) + c(rep(1, n%%m), rep(0, H - n%%m)))
  
  slice.means <- t(sapply(1:H, function(i) colMeans(X[Y.ord[start.ind[i]:end.ind[i]],]%*%SS.inv.half)))
  
  cov(slice.means)
  
}

SIR.LL <- function(Y,X,H,SS.inv.half){
  
  n <- length(Y)
  p <- ncol(X)
  
  Y.ord <- order(Y)
  
  m <- floor(n/H)
  
  start.ind <- cumsum(c(1,(rep(m,H) + c(rep(1, n%%m), rep(0, H - n%%m)))[-H]))
  end.ind <- cumsum(rep(m,H) + c(rep(1, n%%m), rep(0, H - n%%m)))
  
  slice.means <- t(sapply(1:H, function(i) colMeans(X[Y.ord[start.ind[i]:end.ind[i]],]%*%SS.inv.half)))
  
  C <- cov(slice.means)
  V <- eigen(C)$vec[,1]
  
  X.est <- matrix(0,p,n)
  for (j in 1:n){
    X.est[,j] <- slice.means[which((j >= start.ind)&(j <= end.ind)),]
  }
  
  #browser()
  
  t(X.est[,order(Y.ord)])%*%t(t(V)%*%SS.inv.half)
}


err1.final <- list()
err2.final <- list()
err3.final <- list()
err4.final <- list()

signs1.final <- list()
signs2.final <- list()
signs3.final <- list()
signs4.final <- list()

err1.LASSO.final <- list()
err2.LASSO.final <- list()
err3.LASSO.final <- list()
err4.LASSO.final <- list()

signs1.LASSO.final <- list()
signs2.LASSO.final <- list()
signs3.LASSO.final <- list()
signs4.LASSO.final <- list()

err1.LL.final <- list()
err2.LL.final <- list()
err3.LL.final <- list()
err4.LL.final <- list()

signs1.LL.final <- list()
signs2.LL.final <- list()
signs3.LL.final <- list()
signs4.LL.final <- list()

for (j in 1:2){
  p <- c(100,600)[j]#,600,1200)[j]
  
  err1.null <- NULL
  err2.null <- NULL
  err3.null <- NULL
  err4.null <- NULL
  
  signs1.null <- NULL
  signs2.null <- NULL
  signs3.null <- NULL
  signs4.null <- NULL
  
  err1.LASSO.null <- NULL
  err2.LASSO.null <- NULL
  err3.LASSO.null <- NULL
  err4.LASSO.null <- NULL
  
  signs1.LASSO.null <- NULL
  signs2.LASSO.null <- NULL
  signs3.LASSO.null <- NULL
  signs4.LASSO.null <- NULL
  
  err1.LL.null <- NULL
  err2.LL.null <- NULL
  err3.LL.null <- NULL
  err4.LL.null <- NULL
  
  signs1.LL.null <- NULL
  signs2.LL.null <- NULL
  signs3.LL.null <- NULL
  signs4.LL.null <- NULL
  
  for (theta in c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)){
    
    # s <- round(log(p))
    s <- round(sqrt(p))
    n <- round(theta*s*log(p-s))
    n <- H*floor(n/H)
    
    # type Sigma Toeplitz
    Sigma <- matrix(0, p, p)
    Sigma <- .5**abs(row(Sigma)-col(Sigma))
    
    beta <- 1/sqrt(s)*c(rep(1,s),rep(0,p-s))
    beta[1] <- -beta[1]
    beta <- beta/sqrt(beta%*%Sigma%*%beta)
    
    # errors for 4 different models 
    err1 <- NULL
    err2 <- NULL
    err3 <- NULL
    err4 <- NULL
    
    signs1 <- NULL
    signs2 <- NULL
    signs3 <- NULL
    signs4 <- NULL
    
    err1LASSO <- NULL
    err2LASSO <- NULL
    err3LASSO <- NULL
    err4LASSO <- NULL
    
    signs1LASSO <- NULL
    signs2LASSO <- NULL
    signs3LASSO <- NULL
    signs4LASSO <- NULL
    
    err1LL <- NULL
    err2LL <- NULL
    err3LL <- NULL
    err4LL <- NULL
    
    signs1LL <- NULL
    signs2LL <- NULL
    signs3LL <- NULL
    signs4LL <- NULL
    
    for (i in 1:100){
      
      #X <- matrix(rnorm(n*p), nrow = n)
      X <- mvrnorm(n, rep(0,p), Sigma)
      
      SS.inv <- fastclime(half(cov(X)), lambda.min = 0.08, nlambda = 50)
      SS.h.inv <- (fastclime.lambda(SS.inv$lambdamtx, SS.inv$icovlist, 0.14)[[1]])
      #SS.h.inv <- (SS.h.inv + t(SS.h.inv))/2
      
      Y1 <- (X%*%beta) + sin((X%*%beta)) + rnorm(n)
      Y2 <- 2*atan(X%*%beta) + rnorm(n)
      Y3 <- (X%*%beta)**3 + rnorm(n)
      Y4 <- sinh(X%*%beta) + rnorm(n)
      
      Cov1 <- SIR(Y1, X, H = 10, SS.h.inv)
      Cov2 <- SIR(Y2, X, H = 10, SS.h.inv)
      Cov3 <- SIR(Y3, X, H = 10, SS.h.inv)
      Cov4 <- SIR(Y4, X, H = 10, SS.h.inv)
      
      highest1 <- order(round(diag(Cov1),3), decreasing = T)[1:s]
      highest2 <- order(round(diag(Cov2),3), decreasing = T)[1:s]
      highest3 <- order(round(diag(Cov3),3), decreasing = T)[1:s]
      highest4 <- order(round(diag(Cov4),3), decreasing = T)[1:s]
      
      err1 <- c(err1, length(setdiff(1:s,highest1[1:s])))
      err2 <- c(err2, length(setdiff(1:s,highest2[1:s])))
      err3 <- c(err3, length(setdiff(1:s,highest3[1:s])))
      err4 <- c(err4, length(setdiff(1:s,highest4[1:s])))
      
      pvec1 <- eigen(Cov1[highest1, highest1], symmetric = T)$vec[,1]
      pvec2 <- eigen(Cov1[highest2, highest2], symmetric = T)$vec[,1]
      pvec3 <- eigen(Cov1[highest3, highest3], symmetric = T)$vec[,1]
      pvec4 <- eigen(Cov1[highest4, highest4], symmetric = T)$vec[,1]
      
      signs1 <- c(signs1, all(sign(pvec1)*sign(pvec1[1]) == sign(beta[highest1])*sign(beta[highest1][1])))
      signs2 <- c(signs2, all(sign(pvec2)*sign(pvec2[1]) == sign(beta[highest2])*sign(beta[highest2][1])))
      signs3 <- c(signs3, all(sign(pvec3)*sign(pvec3[1]) == sign(beta[highest3])*sign(beta[highest3][1])))
      signs4 <- c(signs4, all(sign(pvec4)*sign(pvec4[1]) == sign(beta[highest4])*sign(beta[highest4][1])))
      
      ## LASSO
      coefs_m1 <- lars(X,Y1,type = "lasso", normalize = FALSE, intercept = FALSE, max.steps = s+1)$beta[s + 1,]
      coefs_m2 <- lars(X,Y2,type = "lasso", normalize = FALSE, intercept = FALSE, max.steps = s+1)$beta[s + 1,]
      coefs_m3 <- lars(X,Y3,type = "lasso", normalize = FALSE, intercept = FALSE, max.steps = s+1)$beta[s + 1,]
      coefs_m4 <- lars(X,Y4,type = "lasso", normalize = FALSE, intercept = FALSE, max.steps = s+1)$beta[s + 1,]
      
      err1LASSO <- c(err1LASSO, length(setdiff(1:s,order(round(abs(coefs_m1),3), decreasing = T)[1:s])))
      err2LASSO <- c(err2LASSO, length(setdiff(1:s,order(round(abs(coefs_m2),3), decreasing = T)[1:s])))
      err3LASSO <- c(err3LASSO, length(setdiff(1:s,order(round(abs(coefs_m3),3), decreasing = T)[1:s])))
      err4LASSO <- c(err4LASSO, length(setdiff(1:s,order(round(abs(coefs_m4),3), decreasing = T)[1:s])))
      
      signs1LASSO <- c(signs1LASSO, (coefs_m1[1] < 0)*all(coefs_m1[-1] >= 0))
      signs2LASSO <- c(signs2LASSO, (coefs_m2[1] < 0)*all(coefs_m2[-1] >= 0))
      signs3LASSO <- c(signs3LASSO, (coefs_m3[1] < 0)*all(coefs_m3[-1] >= 0))
      signs4LASSO <- c(signs4LASSO, (coefs_m4[1] < 0)*all(coefs_m4[-1] >= 0))
      
      ## LEXIN LI
      #v1 <- eigen(Cov1)$vec[,1]
      #v2 <- eigen(Cov2)$vec[,1]
      #v3 <- eigen(Cov3)$vec[,1]
      #v4 <- eigen(Cov4)$vec[,1]
      
      #Y1.tilde <- X%*%v1
      #Y2.tilde <- X%*%v2
      #Y3.tilde <- X%*%v3
      #Y4.tilde <- X%*%v4
      Y1.tilde <- SIR.LL(Y1, X, H = 10, SS.h.inv)#X%*%v1
      Y2.tilde <- SIR.LL(Y2, X, H = 10, SS.h.inv)#X%*%v2
      Y3.tilde <- SIR.LL(Y3, X, H = 10, SS.h.inv)#X%*%v3
      Y4.tilde <- SIR.LL(Y4, X, H = 10, SS.h.inv)#X%*%v4
      
      coefs_m1_ll <- lars(X,Y1.tilde,type = "lasso", normalize = FALSE, intercept = FALSE, max.steps = s+1)$beta[s + 1,]
      coefs_m2_ll <- lars(X,Y2.tilde,type = "lasso", normalize = FALSE, intercept = FALSE, max.steps = s+1)$beta[s + 1,]
      coefs_m3_ll <- lars(X,Y3.tilde,type = "lasso", normalize = FALSE, intercept = FALSE, max.steps = s+1)$beta[s + 1,]
      coefs_m4_ll <- lars(X,Y4.tilde,type = "lasso", normalize = FALSE, intercept = FALSE, max.steps = s+1)$beta[s + 1,]
      
      err1LL <- c(err1LL, length(setdiff(1:s,order(round(abs(coefs_m1_ll),3), decreasing = T)[1:s])))
      err2LL <- c(err2LL, length(setdiff(1:s,order(round(abs(coefs_m2_ll),3), decreasing = T)[1:s])))
      err3LL <- c(err3LL, length(setdiff(1:s,order(round(abs(coefs_m3_ll),3), decreasing = T)[1:s])))
      err4LL <- c(err4LL, length(setdiff(1:s,order(round(abs(coefs_m4_ll),3), decreasing = T)[1:s])))
      
      signs1LL <- c(signs1LL, (coefs_m1[1] < 0)*all(coefs_m1[-1] >= 0))
      signs2LL <- c(signs2LL, (coefs_m2[1] < 0)*all(coefs_m2[-1] >= 0))
      signs3LL <- c(signs3LL, (coefs_m3[1] < 0)*all(coefs_m3[-1] >= 0))
      signs4LL <- c(signs4LL, (coefs_m4[1] < 0)*all(coefs_m4[-1] >= 0))
    }
    
    err1.null <- cbind(err1.null, err1)
    err2.null <- cbind(err2.null, err2)
    err3.null <- cbind(err3.null, err3)
    err4.null <- cbind(err4.null, err4)
    
    signs1.null <- cbind(signs1.null, signs1)
    signs2.null <- cbind(signs2.null, signs2)
    signs3.null <- cbind(signs3.null, signs3)
    signs4.null <- cbind(signs4.null, signs4)
    
    err1.LASSO.null <- cbind(err1.LASSO.null, err1LASSO)
    err2.LASSO.null <- cbind(err2.LASSO.null, err2LASSO)
    err3.LASSO.null <- cbind(err3.LASSO.null, err3LASSO)
    err4.LASSO.null <- cbind(err4.LASSO.null, err4LASSO)
    
    signs1.LASSO.null <- cbind(signs1.LASSO.null, signs1LASSO)
    signs2.LASSO.null <- cbind(signs2.LASSO.null, signs2LASSO)
    signs3.LASSO.null <- cbind(signs3.LASSO.null, signs3LASSO)
    signs4.LASSO.null <- cbind(signs4.LASSO.null, signs4LASSO)
    
    err1.LL.null <- cbind(err1.LL.null, err1LL)
    err2.LL.null <- cbind(err2.LL.null, err2LL)
    err3.LL.null <- cbind(err3.LL.null, err3LL)
    err4.LL.null <- cbind(err4.LL.null, err4LL)
    
    signs1.LL.null <- cbind(signs1.LL.null, signs1LL)
    signs2.LL.null <- cbind(signs2.LL.null, signs2LL)
    signs3.LL.null <- cbind(signs3.LL.null, signs3LL)
    signs4.LL.null <- cbind(signs4.LL.null, signs4LL)
  }
  
  err1.final[[j]] <- err1.null
  err2.final[[j]] <- err2.null
  err3.final[[j]] <- err3.null
  err4.final[[j]] <- err4.null
  
  signs1.final[[j]] <- signs1.null
  signs2.final[[j]] <- signs2.null
  signs3.final[[j]] <- signs3.null
  signs4.final[[j]] <- signs4.null
  
  err1.LASSO.final[[j]] <- err1.LASSO.null
  err2.LASSO.final[[j]] <- err2.LASSO.null
  err3.LASSO.final[[j]] <- err3.LASSO.null
  err4.LASSO.final[[j]] <- err4.LASSO.null
  
  signs1.LASSO.final[[j]] <- signs1.LASSO.null
  signs2.LASSO.final[[j]] <- signs2.LASSO.null
  signs3.LASSO.final[[j]] <- signs3.LASSO.null
  signs4.LASSO.final[[j]] <- signs4.LASSO.null
  
  err1.LL.final[[j]] <- err1.LL.null
  err2.LL.final[[j]] <- err2.LL.null
  err3.LL.final[[j]] <- err3.LL.null
  err4.LL.final[[j]] <- err4.LL.null
  
  signs1.LL.final[[j]] <- signs1.LL.null
  signs2.LL.final[[j]] <- signs2.LL.null
  signs3.LL.final[[j]] <- signs3.LL.null
  signs4.LL.final[[j]] <- signs4.LL.null
}

L1 <- list(err1.final, signs1.final)
L2 <- list(err2.final, signs2.final)
L3 <- list(err3.final, signs3.final)
L4 <- list(err4.final, signs4.final)

save(L1, file="sqrterr1SIR_p_100_600.RData")
save(L2, file="sqrterr2SIR_p_100_600.RData")
save(L3, file="sqrterr3SIR_p_100_600.RData")
save(L4, file="sqrterr4SIR_p_100_600.RData")

LAS1 <- list(err1.LASSO.final, signs1.LASSO.final)
LAS2 <- list(err2.LASSO.final, signs2.LASSO.final)
LAS3 <- list(err3.LASSO.final, signs3.LASSO.final)
LAS4 <- list(err4.LASSO.final, signs4.LASSO.final)

save(LAS1, file="sqrterr1LASSO_p_100_600.RData")
save(LAS2, file="sqrterr2LASSO_p_100_600.RData")
save(LAS3, file="sqrterr3LASSO_p_100_600.RData")
save(LAS4, file="sqrterr4LASSO_p_100_600.RData")

LL1 <- list(err1.LL.final, signs1.LL.final)
LL2 <- list(err2.LL.final, signs2.LL.final)
LL3 <- list(err3.LL.final, signs3.LL.final)
LL4 <- list(err4.LL.final, signs4.LL.final)

save(LL1, file="sqrterr1LL_p_100_600.RData")
save(LL2, file="sqrterr2LL_p_100_600.RData")
save(LL3, file="sqrterr3LL_p_100_600.RData")
save(LL4, file="sqrterr4LL_p_100_600.RData")
