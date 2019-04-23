#############################################
# Oct 25, 2018
# Preliminary Code for Convex Sliced Inverse Regression
# L-ADMM for sparse Sliced Inverse Regression
# Covxy is the conditional covariance
# Covx is marginal covariance of x
# Lambda is the sparsity tuning parameter
# K constraint the rank of the subspace
# Covxy can be constructed using calsigmafit function
#############################################

CovSIR <- function(x, y, Ks = 1:3, lambdas = seq(0.2,3,by=0.5)*sqrt(log(p)/n), nfold = 5, nslice = 5){
  # nfold <- 5
  # nslice <- 5
  n <- dim(x)[1]
  p <- dim(x)[2]
  nk <- length(Ks)
  nlam <- length(lambdas)
  a <- ssir.cv(x,y,Ks,lambdas,nfold,nslice)
  calmean <- function(mat){apply(mat,2,mean)}
  calse <- function(mat){apply(mat,2,sd)/nfold}
  temp1 <- lapply(a,calmean)
  temp1se <- lapply(a,calse)
  cverror <- NULL
  
  for(k in 1:nk){
    tempcv <- min(temp1[[k]])
    tempcv2 <- temp1[[k]][tail(which(temp1[[k]]<=tempcv),1)]
    cverror <- c(cverror,tempcv2)
  }
  
  #	temp2 <- unlist(lapply(temp1,min))
  chosenK_ind <- which(cverror==min(cverror))[1]
  chosenK <- Ks[chosenK_ind]
  # print(chosenK)
  
  templambda <- min(temp1[[chosenK_ind]])
  lambda_ind <- tail(which(temp1[[chosenK_ind]]<=templambda),1)
  chosenlambda <- lambdas[lambda_ind]
  # chosenlambda <- chosenlambda[1]
  # print(lambda_ind)
  
  Sigmax <- cov(x)
  Sigmafit <- calsigmafit(y,x,nslice=5)
  res <- ssir(Sigmafit,Sigmax,chosenlambda,chosenK,epsilon=1e-04,maxiter=1000,trace=FALSE)
  temp <- eigen(round((res$Pi+t(res$Pi))/2,2))
  est <- temp$vectors[,1:chosenK]%*%t(temp$vectors[,1:chosenK])
  mat <- temp$vectors[,1:chosenK]
  
  return(list(mat=mat, mat_full = est, r=chosenK))
}


ssir <- function(covxy,covx,lambda,K,nu=1,epsilon=1e-3,maxiter=1000,trace=FALSE,init=FALSE,initPi=NULL,initH=NULL,initGamma=NULL){
  
  p <- nrow(covx)
  eigencovx <- eigen(covx)
  sqcovx <- eigencovx$vectors%*%sqrt(diag(pmax(eigencovx$values,0)))%*%t(eigencovx$vectors)	
  
  tau <- 4*nu*eigencovx$values[1]^2	
  criteria <- 1e10
  i <- 1
  
  
  # Initialize parameters
  H <- Pi <- oldPi <-  diag(1,p,p)
  Gamma <- matrix(0,p,p)
  
  if(init==TRUE){
    H <- initH
    Pi <- initPi
    Gamma <- initGamma
  }
  
  # While loop for the iterations
  while(criteria > epsilon && i <= maxiter){
    Pi <- updatePi(covx,sqcovx,covxy,H,Gamma,nu,lambda,Pi,tau)
    
    H <- updateH(sqcovx,Gamma,nu,Pi,K)
    Gamma <- Gamma + sqcovx%*%Pi%*%sqcovx-H	
    criteria <- sqrt(sum((Pi-oldPi)^2))
    oldPi <- Pi
    i <- i+1
    if(trace==TRUE)
    {
      print(i)
      print(criteria)
    }
    
  }
  
  return(list(Pi=Pi,H=H,Gamma=Gamma,iteration=i,convergence=criteria))
  
}



######################################################
# Update Pi
######################################################
updatePi <- function(covx,sqcovx,covxy,H,Gamma,nu,lambda,Pi,tau){
  
  A <- Pi + 1/tau*covxy-nu/tau*covx%*%Pi%*%covx+nu/tau*sqcovx%*%(H-Gamma)%*%sqcovx
  B <- lambda/tau
  return(Soft(A,B))
}


######################################################
# Update H
######################################################
updateH <- function(sqcovx,Gamma,nu,Pi,K){
  
  temp <- Gamma + sqcovx%*%Pi%*%sqcovx
  temp <- (temp+t(temp))/2
  svdtemp <- eigen(temp)
  d <- svdtemp$values
  p <- length(d)
  
  if(sum(pmin(1,pmax(d,0)))<=K){
    dfinal <- pmin(1,pmax(d,0))
    return(svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors))
  }
  
  fr <- function(x){
    sum(pmin(1,pmax(d-x,0)))
  }
  # Vincent Vu Fantope Projection
  knots <- unique(c((d-1),d))
  knots <- sort(knots,decreasing=TRUE)
  temp <- which(sapply(knots,fr)<=K)
  lentemp <- tail(temp,1)
  a=knots[lentemp]
  b=knots[lentemp+1]
  fa <- sum(pmin(pmax(d-a,0),1))
  fb <- sum(pmin(pmax(d-b,0),1))
  theta <- a+ (b-a)*(K-fa)/(fb-fa)
  dfinal <- pmin(1,pmax(d-theta,0))
  res <- svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors)
  return(res)
}	



######################################################
# Soft-thresholding Operator
######################################################
Soft <- function(a,b){
  if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
  return(sign(a)*pmax(0,abs(a)-b))
}

############################################
# A function to calculate the conditional covariance
# for sliced inverse regression
# If y is categorical, nslice is # of unique values in y
# If y is continuous, then the conditional covariance is constructed with Sigma_x - T where T is constructed based on Equation (7) 
# X is first mean centered
############################################
calsigmafit <- function(y,X,nslice,standardize=TRUE){
  
  n <- length(y)
  
  if(standardize==TRUE){ 	
    X <- scale(X,T,F)
  }
  
  # Test if this is integer or continuous
  if(sum(y-round(y))==0){
    nslice <- length(unique(y))
    quany <- sort(unique(y))
  }
  
  # For continuous response y
  else if(sum(y-round(y))!=0){
    quany <- quantile(y,seq(1/nslice,1,length.out=nslice))
  }
  
  indexy <- vector("list",nslice)
  
  # Assign indices into indexy
  indexy[[1]] <- which(y<=quany[1])
  for(k in 2:nslice){
    indexy[[k]] <- which(y >quany[k-1] & y<=quany[k])
  }
  
  nindexy <- lapply(indexy,length)
  
  # f matrix 
  f <- matrix(0,n,nslice)
  for(k1 in 1:(nslice-1)){
    for(k2 in 1:nslice){
      if(k1==k2){
        f[indexy[[k1]],k2] <- 1 - nindexy[[k2]]/n
      }
      if(k1!=k2){
        f[indexy[[k1]],k2] <- -nindexy[[k2]]/n			
      }
    }
  }
  for(k in 1:nslice){
    f[indexy[[nslice]],k] <- -nindexy[[k]]/n
  }
  
  bigF <- f%*%solve(t(f)%*%f)%*%t(f)
  Sigmafit <- t(X)%*%bigF%*%X/(n)
  return(Sigmafit)
}


############################################
# Prediction using PFC under X|y Normality Assumption 
# See Cook Statistical Science Paper 2007 and 2008
# Insert output of PFC object and desired rank
# Xnew is a matrix
############################################
predictpfc <- function(pfcobject,K,y,X,Xnew){
  
  Pi <- pfcobject$Pi
  Pi <- (t(Pi)+Pi)/2
  temp <- eigen(Pi)
  if(max(temp$values)<0.01){
    return(rep(mean(y),nrow(Xnew)))
  }
  temp <- temp$vectors[,1:K]
  RhatX <- X%*%temp
  
  Xnew <- as.list(data.frame(t(Xnew)))
  
  predicty <- function(x){ 		
    temp2 <- x%*%temp	
    residual <- t(t(RhatX)-as.vector(temp2))
    weights <- exp(-0.5*apply((residual)^2,1,sum))
    weights <- weights/sum(weights)
    return(sum(weights*y))
  }
  yhat <- unlist(lapply(Xnew,predicty))
  return(yhat)	
} 

############################################
# Apr 4, 2017
# Contain functions for fitting sliced inverse regression 
############################################
############################################
# Cross-validation for choosing K and lambda
# nfold is the number of fold  
# K is a vector of possible values of the rank
# lambda is a vector of possible values of lambda
############################################
ssir.cv <- function(X,y,Ks,lambdas,nfold,nslice){
  nk <- length(Ks)
  nlam <- length(lambdas)
  p <- ncol(X)
  # Index for training vs test set
  fold <- rep(NA,length(y))
  fold <- sample(rep(1:nfold,length(y)/nfold))
  cv.error <- vector("list",nk)
  for(K in 1:nk){
    cv.error[[K]] <- matrix(NA,nrow=nfold,ncol=nlam)	
  }
  for(j in 1:nfold){
    print(paste("cross validation for dataset ",j,sep=""))
    tmp <- 1
    ytrain <- y[which(fold!=j)]
    ytest <- y[which(fold==j)]		
    Xtrain <- X[which(fold!=j),]		
    Xtest <- X[which(fold==j),]		
    Sigmax <- cov(Xtrain)
    Sigmafit <- calsigmafit(ytrain,Xtrain,nslice=nslice)
    
    for(k in 1:nk){
      K = Ks[k]
      initPi = diag(1,p,p)
      initH = diag(1,p,p)
      initGamma = diag(0,p,p)
      
      # tmp <- 1
      for(i in 1:nlam){
        lambda = lambdas[i]
        res <- ssir(Sigmafit,Sigmax,lambda,K,epsilon=5e-04,maxiter=1000,init=TRUE,initPi=initPi,initH=initH,initGamma=initGamma,trace=FALSE)
        initPi = res$Pi
        initH = res$H
        initGamma = res$Gamma
        yhat <- predictpfc(res,K,ytrain,Xtrain,Xtest)
        # Test if this is integer or continuous
        #if(sum(y-round(y))==0){
        #	yhat <- sign(yhat)
        #	cv.error[[K]][j,tmp] <- 	1-sum(diag(table(yhat,ytest)))/sum(table(yhat,ytest))
        #}
        #else if(sum(y-round(y))!=0){
        cv.error[[k]][j,i] <- sum((ytest-yhat)^2)
        #}
        # tmp <- tmp + 1
      }
    }
  }
  return(cv.error)
}
  

############################################
# Simulation for linear regression
############################################
linearcase <- function(n,p,lambdas,seed){
  
  set.seed(seed)
  beta <- rep(0,p)
  beta[1:3] <- 1/sqrt(3)
  
  # Generate AR-1 type covariance for X
  Sigma <- matrix(0.5,p,p)
  tmpmat <- matrix(0,p,p)
  for(i in 1:(p-1)){
    tmpmat[i,i:p] <- c(0:(length(i:p)-1))
  }
  tmpmat = tmpmat+t(tmpmat)
  Sigma <- Sigma^tmpmat	
  
  X <- mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
  y <- X[,1]+X[,2]+X[,3]+2*rnorm(n,0,1)
  nfold <- 5
  nslice <- 5
  
  ########## Cv to select lambda #############
  Ks <- 1:3
  nfold <- 5
  nslice <- 5
  lambdas <- lambdas
  a <- ssir.cv(X,y,Ks,lambdas,nfold,nslice)
  calmean <- function(mat){ apply(mat,2,mean)}
  calse <- function(mat){ apply(mat,2,sd)/nfold}
  temp1 <- lapply(a,calmean)
  temp1se <- lapply(a,calse)
  cverror <- NULL
  
  for(k in Ks){
    
    tempcv <- min(temp1[[k]])#+temp1se[[k]][which(temp1[[k]]==min(temp1[[k]]))][1]
    
    tempcv2 <- temp1[[k]][tail(which(temp1[[k]]<=tempcv),1)]
    
    cverror <- c(cverror,tempcv2)
  }
  
  #	temp2 <- unlist(lapply(temp1,min))
  chosenK <- which(cverror==min(cverror))
  chosenK <- chosenK[1]
  print(chosenK)
  
  templambda <- min(temp1[[chosenK]])#+temp1se[[chosenK]][which(temp1[[chosenK]]==min(temp1[[chosenK]]))][1]
  
  chosenlambda <- lambdas[tail(which(temp1[[chosenK]]<=templambda),1)]
  chosenlambda <- chosenlambda[1]
  print(chosenlambda)
  
  ######### Fit it using CV selected tuning parameters
  Sigmax <- cov(X)
  Sigmafit <- calsigmafit(y,X,nslice=5)
  res <- ssir(Sigmafit,Sigmax,chosenlambda,chosenK,epsilon=1e-04,maxiter=1000,trace=FALSE)
  temp <- eigen(round((res$Pi+t(res$Pi))/2,2))
  est <- temp$vectors[,1:chosenK]%*%t(temp$vectors[,1:chosenK])
  
  cvspec <- sum(abs(diag(est))<1e-5 & beta< 1e-5)/sum(beta==0)
  cvsens <- sum(abs(diag(est))>1e-5 & abs(beta)>1e-5)/sum(beta!=0)
  
  if(chosenK==1){
    abscor <-abs(cor(temp$vectors[,1],beta))
    
  }
  if(chosenK>1){
    abscor <- apply(abs(cor(temp$vectors[,1:chosenK],beta)),2,max)
  }
  
  
  
  return(list(cvsens=cvsens,cvspec=cvspec,K=chosenK,abscor=abscor,lambda=chosenlambda))
}




############################################
expcase <- function(n,p,lambdas,seed){
  
  set.seed(seed)
  beta <- rep(0,p)
  beta[1:3] <- 1/sqrt(3)
  
  # Generate AR-1 type covariance for X
  Sigma <- matrix(0.5,p,p)
  tmpmat <- matrix(0,p,p)
  for(i in 1:(p-1)){
    tmpmat[i,i:p] <- c(0:(length(i:p)-1))
  }
  tmpmat = tmpmat+t(tmpmat)
  Sigma <- Sigma^tmpmat	
  
  X <- mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
  #	y <- X[,1]+X[,2]+X[,3]+2*rnorm(n,0,1)
  y <- 1+exp(X%*%beta)+rnorm(n,0,1)
  nfold <- 5
  nslice <- 5
  
  ########## Cv to select lambda #############
  Ks <- 1:3
  nfold <- 5
  nslice <- 5
  lambdas <- lambdas
  #lambdas <- seq(0.02,0.5,by=0.02)
  #lambdas <- lambdas*sqrt(log(p)/n)
  a <- ssir.cv(X,y,Ks,lambdas,nfold,nslice)
  calmean <- function(mat){ apply(mat,2,mean)}
  calse <- function(mat){ apply(mat,2,sd)/nfold}
  temp1 <- lapply(a,calmean)
  temp1se <- lapply(a,calse)
  cverror <- NULL
  
  for(k in Ks){
    
    tempcv <- min(temp1[[k]])#+temp1se[[k]][which(temp1[[k]]==min(temp1[[k]]))][1]
    
    tempcv2 <- temp1[[k]][tail(which(temp1[[k]]<=tempcv),1)]
    
    cverror <- c(cverror,tempcv2)
  }
  
  #	temp2 <- unlist(lapply(temp1,min))
  chosenK <- which(cverror==min(cverror))
  chosenK <- chosenK[1]
  print(chosenK)
  
  templambda <- min(temp1[[chosenK]])#+temp1se[[chosenK]][which(temp1[[chosenK]]==min(temp1[[chosenK]]))][1]
  
  chosenlambda <- lambdas[tail(which(temp1[[chosenK]]<=templambda),1)]
  chosenlambda <- chosenlambda[1]
  print(chosenlambda)
  
  ######### Fit it using CV selected tuning parameters
  Sigmax <- cov(X)
  Sigmafit <- calsigmafit(y,X,nslice=5)
  res <- ssir(Sigmafit,Sigmax,chosenlambda,chosenK,epsilon=1e-04,maxiter=1000,trace=FALSE)
  temp <- eigen(round((res$Pi+t(res$Pi))/2,2))
  est <- temp$vectors[,1:chosenK]%*%t(temp$vectors[,1:chosenK])
  
  cvspec <- sum(abs(diag(est))<1e-5 & beta< 1e-5)/sum(beta==0)
  cvsens <- sum(abs(diag(est))>1e-5 & abs(beta)>1e-5)/sum(beta!=0)
  
  if(chosenK==1){
    abscor <-abs(cor(temp$vectors[,1],beta))
  }
  if(chosenK>1){
    abscor <- apply(abs(cor(temp$vectors[,1:chosenK],beta)),2,max)
  }
  
  
  
  return(list(cvsens=cvsens,cvspec=cvspec,K=chosenK,abscor=abscor,lambda=chosenlambda))
}


############################################
# Simulation for linear regression
############################################
sinfrac <- function(n,p,lambdas,seed){
  
  set.seed(seed)
  beta1 <- rep(0,p)
  beta2 <- rep(0,p)
  beta1[1:3] <- 1
  beta2[4:5] <- 1
  beta <- beta1+beta2	
  # Generate AR-1 type covariance for X
  Sigma <- matrix(0.5,p,p)
  tmpmat <- matrix(0,p,p)
  for(i in 1:(p-1)){
    tmpmat[i,i:p] <- c(0:(length(i:p)-1))
  }
  tmpmat = tmpmat+t(tmpmat)
  Sigma <- Sigma^tmpmat	
  
  X <- mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
  y <- (X[,1]+X[,2]+X[,3])/(0.5+(X[,4]+X[,5]+1.5)^2)+0.1*rnorm(n,0,1)
  
  nfold <- 5
  nslice <- 5
  
  ########## Cv to select lambda #############
  Ks <- 1:3
  nfold <- 5
  nslice <- 5
  #lambdas <- lambdas*sqrt(log(p)/n)
  a <- ssir.cv(X,y,Ks,lambdas,nfold,nslice)
  calmean <- function(mat){ apply(mat,2,mean)}
  calse <- function(mat){ apply(mat,2,sd)/nfold}
  temp1 <- lapply(a,calmean)
  temp1se <- lapply(a,calse)
  cverror <- NULL
  
  for(k in Ks){
    
    tempcv <- min(temp1[[k]])#+temp1se[[k]][which(temp1[[k]]==min(temp1[[k]]))][1]
    
    tempcv2 <- temp1[[k]][tail(which(temp1[[k]]<=tempcv),1)]
    
    cverror <- c(cverror,tempcv2)
  }
  
  #	temp2 <- unlist(lapply(temp1,min))
  chosenK <- which(cverror==min(cverror))
  chosenK <- chosenK[1]
  print(chosenK)
  
  templambda <- min(temp1[[chosenK]])#+temp1se[[chosenK]][which(temp1[[chosenK]]==min(temp1[[chosenK]]))][1]
  
  chosenlambda <- lambdas[tail(which(temp1[[chosenK]]<=templambda),1)]
  chosenlambda <- chosenlambda[1]
  print(chosenlambda)
  
  ######### Fit it using CV selected tuning parameters
  Sigmax <- cov(X)
  Sigmafit <- calsigmafit(y,X,nslice=5)
  res <- ssir(Sigmafit,Sigmax,chosenlambda,chosenK,epsilon=1e-04,maxiter=1000,trace=FALSE)
  temp <- eigen(round((res$Pi+t(res$Pi))/2,2))
  est <- temp$vectors[,1:chosenK]%*%t(temp$vectors[,1:chosenK])
  
  
  if(chosenK==1){
    abscor <-abs(cor((temp$vectors[,1:chosenK]),cbind(beta1,beta2)))
    abscor <- mean(abscor)
  }
  if(chosenK>1){
    abscor <- apply(abs(cor(temp$vectors[,1:chosenK],cbind(beta1,beta2))),2,max)
    abscor <- mean(abscor)
  }
  
  
  
  cvspec <- sum(abs(diag(est))<1e-5 & beta< 1e-5)/sum(beta==0)
  cvsens <- sum(abs(diag(est))>1e-5 & abs(beta)>1e-5)/sum(beta!=0)
  
  
  return(list(cvsens=cvsens,cvspec=cvspec,K=chosenK,abscor=abscor,lambda=chosenlambda))
}



############################################
# Simulation (1) expfro Frobenius Norm
############################################
expfro <- function(ns,p,s,seed){
  
  set.seed(seed)
  beta <- rep(0,p)
  beta[1:s] <- 1/sqrt(s)
  
  # Generate AR-1 type covariance for X
  Sigma <- matrix(0.5,p,p)
  tmpmat <- matrix(0,p,p)
  for(i in 1:(p-1)){
    tmpmat[i,i:p] <- c(0:(length(i:p)-1))
  }
  tmpmat = tmpmat+t(tmpmat)
  Sigma <- Sigma^tmpmat	
  eigenSigma <- eigen(Sigma)
  sqrtSigma <- eigenSigma$vectors%*%diag(sqrt(eigenSigma$values))%*%t(eigenSigma$vectors)	
  
  mse <- number <- NULL
  
  truesubspace <- cbind(beta)%*%t(cbind(beta))
  temp <- svd(truesubspace)
  truesubspace <- temp$u[,1]%*%t(temp$u[,1])
  
  for(n in ns){
    print(n)
    lambdas <- seq(0.2,5,by=0.5)*sqrt(log(p)/n)
    Ks <- 1
    X <- mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
    y <- 1+exp(X%*%beta)+rnorm(n,0,1)
    
    
    chosenlambda <- 2*sqrt(log(p)/n)	
    chosenK <- 1
    
    # Sample covariance
    Sigmafit <- calsigmafit(y,X,nslice=5)
    Sigmax <- cov(X)	
    
    res <- ssir(Sigmafit,Sigmax,chosenlambda,chosenK,epsilon=1e-03,maxiter=1000,trace=FALSE)
    
    Pi <- res$Pi
    temp <- svd(Pi)
    est <- temp$u[,1:chosenK]%*%t(temp$u[,1:chosenK])
    mse <- c(mse,sqrt(sum((truesubspace-est)^2)))
    number <- c(number,n)
    
  }
  
  return(list(mse=mse,number=number))	
}



