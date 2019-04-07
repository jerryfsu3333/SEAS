library(MASS)

# SIR for linear regression
# this function returns the SIR estimate of the covariance matrix
SIR <- function(x,y,H){
  
  n <- length(y)
  p <- ncol(x)
  
  y.ord <- order(y)
  
  m <- floor(n/H)
  
  start.ind <- cumsum(c(1,(rep(m,H) + c(rep(1, n%%m), rep(0, H - n%%m)))[-H]))
  end.ind <- cumsum(rep(m,H) + c(rep(1, n%%m), rep(0, H - n%%m)))
    
  slice.means <- t(sapply(1:H, function(i) colMeans(x[y.ord[start.ind[i]:end.ind[i]],])))
  
  cov(slice.means)
}

# err1.final <- list()
# err2.final <- list()
# err3.final <- list()
# err4.final <- list()
# 
# signs1.final <- list()
# signs2.final <- list()
# signs3.final <- list()
# signs4.final <- list()

# for (j in 1:5){
#   p <- c(100,200,300,600,1200)[j]
#   
#   err1.null <- NULL
#   err2.null <- NULL
#   err3.null <- NULL
#   err4.null <- NULL
#   
#   signs1.null <- NULL
#   signs2.null <- NULL
#   signs3.null <- NULL
#   signs4.null <- NULL
  
  # for (theta in c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)){
    # s <- round(sqrt(p))
    # n <- round(theta*s*log(p-s))
    # 
    # #beta <- 1/sqrt(s)*c(rep(1,s-1),-1,rep(0,p-s))
    # beta <- c(runif(s,.5,1), rep(0,p-s))
    # beta[s] <- -beta[s]
    # beta <- beta/sqrt(sum(beta**2))
    
    # errors for 4 different models 
    # err1 <- NULL
    # err2 <- NULL
    # err3 <- NULL
    # err4 <- NULL
    # 
    # signs1 <- NULL
    # signs2 <- NULL
    # signs3 <- NULL
    # signs4 <- NULL
    
    # for (i in 1:100){
    #         
    #   X <- matrix(rnorm(n*p), nrow = n)
    #   
    #   Y1 <- (X%*%beta) + sin((X%*%beta)) + rnorm(n)
    #   Y2 <- 2*atan(X%*%beta) + rnorm(n)
    #   Y3 <- (X%*%beta)**3 + rnorm(n)
    #   Y4 <- sinh(X%*%beta) + rnorm(n)
DT_SIR <- function(x, y, s, H = 10){
  
      Cov1 <- SIR(x, y, H)
      # Cov2 <- SIR(Y2, X, H = 10)
      # Cov3 <- SIR(Y3, X, H = 10)
      # Cov4 <- SIR(Y4, X, H = 10)
      
      highest1 <- order(round(diag(Cov1),3), decreasing = T)[1:s]
      # highest2 <- order(round(diag(Cov2),3), decreasing = T)[1:s]
      # highest3 <- order(round(diag(Cov3),3), decreasing = T)[1:s]
      # highest4 <- order(round(diag(Cov4),3), decreasing = T)[1:s]
      
      err1 <- length(setdiff(1:s,highest1[1:s]))
      # err2 <- c(err2, length(setdiff(1:s,highest2[1:s])))
      # err3 <- c(err3, length(setdiff(1:s,highest3[1:s])))
      # err4 <- c(err4, length(setdiff(1:s,highest4[1:s])))
      
      pvec1 <- eigen(Cov1[highest1, highest1], symmetric = T)$vec[,1, drop = FALSE]
      # pvec2 <- eigen(Cov1[highest2, highest2], symmetric = T)$vec[,1]
      # pvec3 <- eigen(Cov1[highest3, highest3], symmetric = T)$vec[,1]
      # pvec4 <- eigen(Cov1[highest4, highest4], symmetric = T)$vec[,1]
      B <- matrix(0, p, dim(pvec1)[2])
      B[highest1,] <- pvec1
      B
      # signs1 <- all(sign(pvec1)*sign(pvec1[1]) == sign(beta[highest1])*sign(beta[highest1][1]))
}
      # signs2 <- c(signs2, all(sign(pvec2)*sign(pvec2[1]) == sign(beta[highest2])*sign(beta[highest2][1])))
      # signs3 <- c(signs3, all(sign(pvec3)*sign(pvec3[1]) == sign(beta[highest3])*sign(beta[highest3][1])))
      # signs4 <- c(signs4, all(sign(pvec4)*sign(pvec4[1]) == sign(beta[highest4])*sign(beta[highest4][1])))
    # }
    
  #   err1.null <- cbind(err1.null, err1)
  #   err2.null <- cbind(err2.null, err2)
  #   err3.null <- cbind(err3.null, err3)
  #   err4.null <- cbind(err4.null, err4)
  #   
  #   signs1.null <- cbind(signs1.null, signs1)
  #   signs2.null <- cbind(signs2.null, signs2)
  #   signs3.null <- cbind(signs3.null, signs3)
  #   signs4.null <- cbind(signs4.null, signs4)
  #   
  # }
  
  # err1.final[[j]] <- err1.null
  # err2.final[[j]] <- err2.null
  # err3.final[[j]] <- err3.null
  # err4.final[[j]] <- err4.null
  # 
  # signs1.final[[j]] <- signs1.null
  # signs2.final[[j]] <- signs2.null
  # signs3.final[[j]] <- signs3.null
  # signs4.final[[j]] <- signs4.null
# }

# L1 <- list(err1.final, signs1.final)
# L2 <- list(err2.final, signs2.final)
# L3 <- list(err3.final, signs3.final)
# L4 <- list(err4.final, signs4.final)

# save(L1, file="sqrterr1.RData")
# save(L2, file="sqrterr2.RData")
# save(L3, file="sqrterr3.RData")
# save(L4, file="sqrterr4.RData")
