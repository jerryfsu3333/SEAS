# AR function
AR <- function(rho, p){
  m <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      m[i,j] <- rho**(abs(i-j))
    }
  }
  return(m)
}

# CS function
CS <- function(rho, p){
  m <- matrix(rho,p,p)
  diag(m) <- 1
  return(m)
}

# Block CS function

CS_blk <- function(rho,p,s){
  block <- CS(rho,s)
  A <- diag(rep(1,p/s),p/s,p/s)
  m <- kronecker(A,block)
  return(m)
}

CS_blk2 <- function(rho,p,s){
  block <- CS(rho,s)
  m <- diag(rep(1,p),p,p)
  m[1:s,1:s] <- block
  return(m)
}

C_IC <- function(mat, all, sig){
  if(is.null(mat)) {return(list(C = NA, IC = NA))}
  tmp <- apply(mat, 1, function(x) any(x!=0))
  C <- sum(which(tmp) %in% sig)/length(sig)
  IC <- sum(which(tmp) %in% setdiff(all, sig))/length(setdiff(all, sig))
  list(C = C, IC = IC)
}

C_IC_cut <- function(mat, all, sig){
  if(is.null(mat)) {return(list(C = NA, IC = NA))}
  mat <- mat %*% t(mat)
  tmp <- diag(mat) > 1e-5
  C <- sum(which(tmp) %in% sig)/length(sig)
  IC <- sum(which(tmp) %in% setdiff(all, sig))/length(setdiff(all, sig))
  list(C = C, IC = IC)
}

# If some whole column is significantly small then we cut the columns to zero
cut_mat <- function(Beta, thrd, rank){
  l <- length(Beta)
  for (i in 1:l){
    if(is.null(Beta[[i]])) next
    mat <- as.matrix(Beta[[i]])
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    r <- rank[[i]]
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
      vec <- as.vector(mat)
      vec[abs(vec) < thrd] <- 0
      Beta[[i]] <- matrix(vec, nobs, nvars)
    }
  }
  return(Beta)
}


eval_val_rmse <- function(Beta, x, y){
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      xnew <- x %*% mat
      fit <- lm(y~xnew)
      rmse <- sqrt(mean((fit$residuals)^2))
      rmse 
    }
  })
  result
}

eval_val_obj <- function(Beta, d, sigma, mu){
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    # for (i in 1:l){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      rank <- d[[i]]
      if(rank != 0){
        tmp <- svd(mat)
        mat <- tmp$u[,1:rank, drop = FALSE] %*% diag(tmp$d[1:rank], rank, rank) %*% t(tmp$v[,1:rank, drop=FALSE])
      }
      # If rank == 0, then it's been cut to zero matrix by cut_mat function
      # result[i] <- sum(diag(t(mat) %*% sigma %*% mat - 2 * t(mu) %*% mat))
      sum(diag(t(mat) %*% sigma %*% mat - 2 * t(mu) %*% mat))
    }
    # }
  })
  return(result)
}

eval_val_dc <- function(Beta, x, y){
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      -dcor(x %*% mat, y)
    }
  })
  return(result)
}

eval_val_ker <- function(Beta, x, y, type=1){
  y <- drop(y)
  nobs <- length(y)
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      if(type == 1){
        
      }else if(type == 2){
        
        rank <- d[[i]]
        if(rank != 0){
          tmp <- svd(mat)
          mat <- tmp$u[,1:rank, drop = FALSE] %*% diag(tmp$d[1:rank], rank, rank) %*% t(tmp$v[,1:rank, drop=FALSE])
        }
        
      }else if(type==3){
        
        rank <- d[[i]]
        if(rank != 0){
          tmp <- svd(mat)
          mat <- tmp$u[,1:rank, drop = FALSE]
        }
        
      }
      
      xnew <- x %*% mat
      bws <- apply(xnew, 2, function(x){IQR(x)/5})
      # fit_bw <- npregbw(xdat = xnew, ydat = y, bws = bws, regtype = 'll', nmulti = 1, bandwidth.compute = FALSE)
      fit_bw <- npregbw(xdat = xnew, ydat = y, bws = bws, regtype = 'lc', bwtype='fixed', ckernel='gaussian', bandwidth.compute = FALSE)
      # fit_bw <- npregbw(xdat = xnew, ydat = y, regtype = 'll', nmulti = 1)
      fit <- npreg(fit_bw)
      rmse <- sqrt(mean((y - fit$mean)^2))
      rmse
    }
  })
  result
}

# Prediction function
lda_pred <- function(xtrain, ytrain, theta, xtest){
  xfit <- xtrain %*% theta
  l <- lda(xfit, ytrain)
  pred <- predict(l, xtest %*% theta)$class
  return(pred)
}

# index of non-zero variables
nz_func <- function(mat){
  if(is.null(mat)){return(NA)}
  if(is.vector(mat)){mat <- as.matrix(mat)}
  s <- which(apply(mat, 1, function(x){any(x!=0)}))
  return(s)
}


orth_mat <- function(Beta, rank){
  l <- length(Beta)
  for (i in 1:l){
    mat <- Beta[[i]]
    r <- rank[i]
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    # If rank is equal to zero, then set Beta to zero
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
      Beta[[i]] <- svd(mat)$u[,1:r, drop=FALSE]
    }
  }
  return(Beta)
}


# input the ssdr returned Beta matrix, returns corresponding predictions
predict_ssdr <- function(x_train, y_train, fit, newx){
  mat <- fit$beta
  rank <- fit$rank
  prior <- sapply(1:K, function(x) sum(y_train==x)/length(y_train))
  n.col <- length(mat)
  n.row <- nrow(newx)
  pred <- matrix(0,n.row,n.col)
  for(i in 1:n.col){
    beta <- mat[[i]]
    nz <- sum(beta[,1] != 0)
    if(is.null(beta) || nz == 0){
      # Sometimes ssdr is not converge, so we leave the corresponding matrix null, in this case we use
      # prior to predict. Or if matrix is a zero matrix then we also use prior.
      pred[,i] <- which.max(prior)
    }else{
      r <- rank[i]
      if(r == 0){
        pred[,i] <- which.max(prior)
      }else{
        subset <- svd(beta)$u[,1:r,drop = FALSE]
        pred[,i] <- lda_pred(x_train,y_train,subset,newx)
      }
    }
  }
  return(pred)
}

prep <- function(x, y, yclass=NULL, H=length(unique(yclass)), type='sir', cut_y=FALSE){
  
  if(is.null(yclass)){
    stop('The class of y is not provided.')
  }
  nclass <- H
  cls <- unique(yclass)
  if(H != length(cls)){stop("H should be equal to the number of classes.")}
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  prior <- sapply(cls, function(i){mean(yclass == i)})
  sigma <- cov(x)
  
  if(type == 'sir'){
    mu <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      mu[, i] <- apply(x[yclass == cls[i], ], 2, mean) - colMeans(x)
    }
  }else if(type == 'pfc'){
    cut_func <- function(x, lb, ub){
      if(x < lb){
        return(lb)
      } else if(x > ub){
        return(ub) 
      } else{
        return(x)
      }
    }
    if(cut_y){
      lb <- quantile(y, 0.2)[[1]]
      ub <- quantile(y, 0.8)[[1]]
      y <- sapply(y, cut_func, lb = lb, ub = ub) 
    }
    Fmat <- cbind(y, y^2, y^3)
    # Fmat <- cbind(y, abs(y))
    # Fmat <- cbind(y, y^2)
    Fmat_c <- scale(Fmat,scale = FALSE)
    x_c <- scale(x, scale = FALSE)
    # invhalf_func <- function(Sigma){
    #   S <- eigen(Sigma)
    #   pos <- 1/sqrt(S$val[S$val > 0.001])
    #   zer <- rep(0,length(S$val < 0.001))
    #   mid <- diag(c(pos,zer), ncol(Sigma), ncol(Sigma))
    #   S$vec%*%mid%*%t(S$vec)
    # }
    # tmp <- svd(t(Fmat_c) %*% Fmat_c)
    # invhalf <- tmp$u %*% diag(1/sqrt(tmp$d)) %*% t(tmp$u)
    # invhalf <- invhalf_func(t(Fmat_c) %*% Fmat_c)
    # mu <- (t(x_c) %*% Fmat_c %*% invhalf)
    mu <- (t(x_c) %*% Fmat_c)
  }else if(type == 'intra'){
    mu <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      y_copy <- y
      y_copy[yclass!=cls[i]] <- 0
      mu[, i] <- cov(y_copy, x)
    }
  }
  
  list(sigma = sigma, mu = mu, nclass = nclass, prior=prior)
}

# rank function

rank_func <- function(B, thrd){
  d <- svd(B)$d
  r <- sum(d >= thrd)
  return(r)
}

rank_func2 <- function(B, thrd){
  d <- svd(B)$d
  if(max(d) == 0){
    r <- 0
  }else{
    r <- sum(d/d[1] >= thrd)
  }
  return(r)
}

# rank_list <- function(Beta, thrd=1e-3){
#   rank <- vector("list", length(Beta))
#   for (i in 1:length(Beta)){
#     if(!is.null(Beta[[i]])){
#       rank[[i]] <- rank_func(Beta[[i]], thrd = thrd)
#     }
#   }
#   return(rank)
# }

# Draw the plot of the ratio of singular values

sv_plot <- function(sv){
  tmp <- rep(0,length(sv)-1)
  for(i in 1:(length(sv)-1)){
    tmp[i] <- (sv[i] + 0.001)/(sv[i+1] + 0.001)
  }
  plot(tmp, main = "Ratio of consecutive singular values", ylab = "ratio")
  plot(sv, main = "Singular values", ylab = "singular values")
}


# subspace function for matrices with the same column dimension
subspace <- function(A,B){
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  d <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*d))
}

# subspace function for matrices with the different column dimension
subspace_2 <- function(Beta, B){
  if(is.vector(Beta)){dim(Beta) <- c(length(Beta), 1)}
  if(is.vector(B)){dim(B) <- c(length(B), 1)}
  if(is.null(B)) {return(NA)}
  if(all(Beta == 0)){
    Pa <- Beta
  }else{
    Pa <- qr.Q(qr(Beta)) 
  }
  Pa <- Pa %*% t(Pa)
  if(all(B == 0)){
    Pb <- B
  }else{
    Pb <- qr.Q(qr(B)) 
  }
  Pb <- Pb %*% t(Pb)
  d <- dim(Beta)[2]
  result <- norm(Pa - Pa %*% Pb %*% Pa, type = 'F')/sqrt(2*d)
  return(result)
}