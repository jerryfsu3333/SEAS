#############  Model I #############
Model1 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 1)
  Beta[1:20,1] <- 1
  Beta <- sqrt(0.5)*Beta/norm(Beta, '2')
  
  nz_vec <- 1:20
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta + 0.5 * rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model II #############
Model2 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 1)
  Beta[1:20,1] <- 1
  Beta <- sqrt(0.5)*Beta/norm(Beta, '2')
  
  nz_vec <- 1:20
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- (x %*% Beta)^3/2 + rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

# #############  (Model III) #############
# Model3 <- function(p=100){
#   Mu <- rep(0,p)
#   Sigma <- AR(0.5, p)
#   
#   # Construct true Beta
#   Beta <- matrix(0, p, 2)
#   Beta[1:6,1] <- 1
#   Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
#   Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
#   
#   nz_vec <- 1:6
#   True_sp <- Beta
#   
#   Data <- function(N){
#     x <- Train(N, Mu, Sigma)
#     nobs <- dim(x)[1]
#     y <- abs((x %*% Beta[,1]) / 4 + 2)^3 * sign(x %*% Beta[,2]) + 0.2 * rnorm(nobs)
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
# 
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  (Model IV) #############
# Model4 <- function(p=100){
#   Mu <- rep(0,p)
#   Sigma <- AR(0.5, p)
#   
#   # Construct true Beta
#   Beta <- matrix(0, p, 2)
#   Beta[1:6,1] <- 1
#   Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
#   Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
#   
#   nz_vec <- 1:6
#   True_sp <- Beta
#   
#   Data <- function(N){
#     x <- Train(N, Mu, Sigma)
#     nobs <- dim(x)[1]
#     y <- abs((x %*% Beta[,1]) / 4 + 2)^3 * sign(x %*% Beta[,2]) + rnorm(nobs)
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
# 
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }

#############  Model III1 (Model V1) #############
Model3_1 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.2 *  rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model (V3) #############
Model3_2 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.6 *  rnorm(nobs))
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model (V5) #############
Model3_3 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model (X1) #############
Model4_1 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- (x %*% Beta[,1]) * exp(x %*% Beta[,2]) + 0.2 * rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model (X3) #############
Model4_2 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- (x %*% Beta[,1]) * exp(x %*% Beta[,2]) + 0.6 * rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model (X5) #############
Model4_3 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- (x %*% Beta[,1]) * exp(x %*% Beta[,2]) + rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

# #############  Model XI(PFC) #############
# # Comment: cant increase too much on lambda.factor and lam_fac_msda, which will lead to a undesirable
# # reduction on the sparsity.
# Model16 <- function(p=100){
# 
# 
#   d <- 2
#   r <- 2
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Beta <- matrix(runif(d*r), d, r)
#   Delta <- 0.1^2*diag(1,p,p)
# 
#   # Construct true Beta
#   nz_vec <- 1:6
#   True_sp <- Gamma
# 
#   Data <- function(N){
#     y <- runif(N,0,4)
#     f <- cbind(y,exp(y))/2
#     eps <- Train(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
# 
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = FALSE)
# 
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############################
# Model16_2 <- function(p=100){
#   
#   
#   d <- 2
#   r <- 2
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   # Beta <- matrix(runif(d*r), d, r)
#   Beta <- diag(c(1,1))
#   Delta <- 2^2*diag(1,p,p)
#   
#   # Construct true Beta
#   nz_vec <- 1:6
#   True_sp <- Gamma
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     f <- cbind(y,y^2)
#     eps <- Train(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  PFC model #############
# Model17 <- function(p=100){
#   
#   d <- 2
#   r <- 3
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Beta <- matrix(rnorm(d*r), d, r)
#   Delta <- diag(1,p,p)
#   
#   # Construct true Beta
#   nz_vec <- 1:6
#   True_sp <- Gamma
#   
#   Data <- function(N){
#     y <- runif(N,0,4)
#     f <- cbind(y,y^2,y^3)
#     eps <- Train(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.9, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = FALSE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model XII (PFC) #############
# Model17_2 <- function(p=100){
#   
#   d <- 2
#   r <- 3
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Beta <- matrix(runif(d*r), d, r)
#   Delta <- diag(1,p,p)
#   
#   # Construct true Beta
#   nz_vec <- 1:6
#   True_sp <- Gamma
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     f <- cbind(y,y^2,y^3)
#     eps <- Train(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + 0.1*eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.9, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model18 <- function(p=100){
#   
#   
#   d <- 2
#   r <- 3
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Beta <- matrix(runif(d*r), d, r)
#   Delta <- 1^2*diag(1,p,p)
#   
#   # Construct true Beta
#   nz_vec <- 1:6
#   True_sp <- Gamma
#   
#   Data <- function(N){
#     y <- runif(N,0,4)
#     f <- cbind(y,y^2,y^3)
#     eps <- Train(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = FALSE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model19 <- function(p=100){
#   
#   
#   d <- 2
#   r <- 3
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Beta <- matrix(runif(d*r), d, r)
#   Delta <- 10^2*diag(1,p,p)
#   
#   # Construct true Beta
#   nz_vec <- 1:6
#   True_sp <- Gamma
#   
#   Data <- function(N){
#     y <- runif(N,0,4)
#     f <- cbind(y,y^2,y^3)
#     eps <- Train(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.6, H = 5)
#   intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.6, H = 5)
#   pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.6, cut_y = FALSE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# ########### PFC-sir model ################
# Model21 <- function(p=100, H = 5){
# 
#   d <- 2
#   r <- H-1
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Beta <- matrix(rnorm(d*r, 0, 1), d, r)
#   Delta <- diag(rep(1,p), p, p)
#   
#   nz_vec <- 1:6
#   True_sp <- Gamma
#   
#   Data <- function(N){
#     y <- runif(N, 0, 4)
#     y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
#     yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
#     f <- c()
#     for (i in 1:(H-1)){
#       tmp <- as.double(yclass == i)
#       f <- cbind(f, tmp)
#     }
#     eps <- mvrnorm(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + 0.2 * eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.1, length.out = 10), cut_y = FALSE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# ########### PFC-sir model ################
# Model21_2 <- function(p=100, H = 5){
#   
#   d <- 2
#   r <- H-1
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Beta <- matrix(runif(d*r), d, r)
#   Delta <- diag(1, p, p)
#   
#   nz_vec <- 1:6
#   True_sp <- Gamma
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
#     yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
#     f <- c()
#     for (i in 1:(H-1)){
#       tmp <- as.double(yclass == i)
#       f <- cbind(f, tmp)
#     }
#     eps <- mvrnorm(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + 0.1*eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.05, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.05, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.7, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.05, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# ########### PFC-sir model ################
# Model21_3 <- function(p=100, H = 5){
# 
#   d <- 2
#   r <- H-1
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Beta <- matrix(runif(d*r), d, r)
#   Delta <- diag(1, p, p)
#   
#   nz_vec <- 1:6
#   True_sp <- Gamma
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     y_breaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
#     yclass <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
#     f <- c()
#     for (i in 1:(H-1)){
#       tmp <- as.double(yclass == i)
#       f <- cbind(f, tmp)
#     }
#     eps <- mvrnorm(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + 0.02 * eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.05, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.05, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.7, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.05, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model22 <- function(p=100){
# 
#   d <- 2
#   r <- 2
#   Gamma <- matrix(0, p, d)
#   Gamma[1:5,1] <- c(1,1,-1,-1,0)
#   Gamma[1:5,2] <- c(1,0,1,0,1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Delta <- diag(1,p,p)
# 
#   nz_vec <- 1:5
#   True_sp <- Gamma
# 
#   Data <- function(N){
#     y <- rnorm(N,0,2)
#     v <- cbind(y, abs(y))
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# 
# #############  Model PFC #############
# Model22_2 <- function(p=100){
#   
#   d <- 2
#   r <- 2
#   tmp <- matrix(0, p, d)
#   tmp[1:6,1] <- 1
#   tmp[1:6,2] <- c(1,-1,1,-1,1,-1)
#   tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
#   tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
#   Delta <- AR(0.1,p)
#   
#   nz_vec <- 1:6
#   Gamma <- Delta %*% tmp
#   True_sp <- tmp
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(y + 1/5*y^2, 2*abs(y))
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model22_3 <- function(p=100){
#   
#   d <- 2
#   r <- 2
#   Gamma <- matrix(0, p, d)
#   Gamma[1:5,1] <- c(1,1,-1,-1,0)
#   Gamma[1:5,2] <- c(1,0,1,0,1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Delta <- AR(0.5,p)
#   
#   nz_vec <- 1:6
#   True_sp <- solve(Delta)%*%Gamma
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(y, abs(y))
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model22_4 <- function(p=100){
#   
#   d <- 2
#   r <- 2
#   Gamma <- matrix(0, p, d)
#   Gamma[1:5,1] <- c(1,1,-1,-1,0)
#   Gamma[1:5,2] <- c(1,0,1,0,1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Delta <- AR(0.5,p)
#   
#   nz_vec <- 1:6
#   True_sp <- solve(Delta) %*% Gamma
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(y + 1/5*y^2, 2*abs(y))
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model22_6 <- function(p=100){
#   
#   d <- 2
#   r <- 2
#   tmp <- matrix(0, p, d)
#   tmp[1:5,1] <- c(1,1,-1,-1,0)
#   tmp[1:5,2] <- c(1,0,1,0,1)
#   tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
#   tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
#   Delta <- AR(0.5,p)
#   
#   nz_vec <- 1:5
#   Gamma <- Delta %*% tmp
#   True_sp <- tmp
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(y, abs(y))
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model22_7 <- function(p=100){
#   
#   d <- 2
#   r <- 2
#   tmp <- matrix(0, p, d)
#   tmp[1:5,1] <- c(1,1,-1,-1,0)
#   tmp[1:5,2] <- c(1,0,1,0,1)
#   tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
#   tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
#   Delta <- AR(0.5,p)
#   
#   nz_vec <- 1:5
#   Gamma <- Delta %*% tmp
#   True_sp <- tmp
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(y, y^2)
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model22_8 <- function(p=100){
#   
#   d <- 2
#   tmp <- matrix(0, p, d)
#   tmp[1:4,1] <- 1
#   tmp[1:4,2] <- c(1,-1,1,-1)
#   tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
#   tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
#   Delta <- AR(0.5,p)
#   
#   nz_vec <- 1:4
#   Gamma <- Delta %*% tmp
#   True_sp <- tmp
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(y, abs(y))
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model intra-slice #############
# Model23 <- function(p=5){
#   
#   Beta1 <- c(1,0,0,0,0)
#   Beta2 <- c(0,1,1,0,0)
#   Beta <- cbind(Beta1, Beta2)
#   
#   # Construct true Beta
#   nz_vec <- 1:3
#   True_sp <- Beta
#   
#   Data <- function(N){
#     v1 <- rt(N,5)
#     v2 <- rt(N,5)
#     v3 <- rt(N,5)
#     w1 <- rgamma(N, shape = 0.2)
#     w2 <- rgamma(N, shape = 0.2)
#     x <- cbind(w1, v1 + w2/2, -v1 + w2/2, v2 + v3, v2 - v3)
#     y <- 1.5 * (5 + x %*% Beta[,1]) * (2 + x %*% Beta[,2]) + 0.5 * rnorm(N)
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.7, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
#   intra_params <- list(lambda.factor = 0.7, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
#   pfc_params <- list(lambda.factor = 0.7, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, cut_y = FALSE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model24 <- function(p=100){
#   
#   d <- 2
#   r <- 2
#   tmp <- matrix(0, p, d)
#   tmp[1:6,1] <- 1
#   tmp[1:6,2] <- c(1,-1,1,-1,1,-1)
#   tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
#   tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
#   Delta <- CS(0.5,p)
# 
#   nz_vec <- 1:6
#   Gamma <- Delta %*% tmp
#   True_sp <- tmp
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(y+1/5*abs(y), y^2/2)
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + 1*eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model25 <- function(p=100){
#   
#   d <- 2
#   r <- 2
#   p0 <- 5
#   Gamma <- matrix(0, p, d)
#   Gamma[1:p0,1] <- 1
#   Gamma[(p0+1):(2*p0),2] <- 1
#   # Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   # Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   
#   # A <- matrix(rnorm(p^2), p, p)
#   # Delta <- t(A) %*% A
#   Delta <- diag(c(rep(2,p0), rep(20,p0), rep(2,p-2*p0)))
#   Beta <- diag(c(2/sqrt(20), 1/sqrt(20)))
#   # Beta <- diag(c(0.2,0.1))
#   
#   # Construct true Beta
#   nz_vec <- 1:(2*p0)
#   True_sp <- solve(Delta) %*% Gamma
#   
#   Data <- function(N){
#     y <- runif(N,0,3)
#     f <- cbind(y, exp(y))
#     eps <- Train(N, rep(0,p), Delta)
#     x <- f %*% t(Beta) %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #############  Model PFC #############
# Model25_2 <- function(p=100){
#   
#   d <- 2
#   r <- 2
#   p0 <- 5
#   Gamma <- matrix(0, p, d)
#   Gamma[1:p0,1] <- 1
#   Gamma[(p0+1):(2*p0),2] <- 1
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   
#   # A <- matrix(rnorm(p^2), p, p)
#   # Delta <- t(A) %*% A
#   Delta <- diag(c(rep(2,p0), rep(20,p0), rep(2,p-2*p0)))
#   # Beta <- diag(c(0.2,0.1))
#   
#   # Construct true Beta
#   nz_vec <- 1:(2*p0)
#   True_sp <- solve(Delta) %*% Gamma
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(2/sqrt(20)*y + 1/sqrt(20)*y^3, 3/sqrt(20)*y^2)
#     eps <- mvrnorm(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + 0.2*eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(0.5,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# 
# #############  Model PFC #############
# # Block CS
# Model26 <- function(p=100){
#   
#   d <- 2
#   r <- 2
#   Gamma <- matrix(0, p, d)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Delta <- diag(1,p,p)
#   Delta[1:6,1:6] <- AR(0.3,6)
#   
#   nz_vec <- 1:6
#   True_sp <- solve(Delta) %*% Gamma
#   
#   Data <- function(N){
#     y <- runif(N,0,3)
#     v <- cbind(y, y^2/2)
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + 0.5*eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.05, length.out = 10), cut_y = FALSE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }

#############  Model (XXVII) #############
# From paper CISE
Model5 <- function(p=100){
  
  d <- 2
  Gamma <- matrix(0, p, d)
  Gamma[1:4,1] <- 1
  Gamma[1:4,2] <- c(1,-1,1,-1)
  Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  Delta <- AR(0.8,p)
  
  nz_vec <- 1:5
  True_sp <- solve(Delta) %*% Gamma
  
  Data <- function(N){
    y <- rnorm(N,0,1)
    v <- cbind(y, y^2)
    eps <- Train(N, rep(0,p), Delta)
    x <- v %*% t(Gamma) + eps
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model (XXVII2) #############
Model6 <- function(p=100){
  
  d <- 2
  Gamma <- matrix(0, p, d)
  Gamma[1:4,1] <- 1
  Gamma[1:4,2] <- c(1,-1,1,-1)
  Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  Delta <- AR(0.8,p)
  
  nz_vec <- 1:5
  True_sp <- solve(Delta)%*%Gamma
  
  Data <- function(N){
    y <- rnorm(N,0,1)
    v <- cbind(y, abs(y))
    eps <- Train(N, rep(0,p), Delta)
    x <- v %*% t(Gamma) + eps
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
  intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.3, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

# #############  Model (XXVII3) #############
# Model7 <- function(p=100){
#   
#   d <- 2
#   Gamma <- matrix(0, p, d)
#   Gamma[1:4,1] <- 1
#   Gamma[1:4,2] <- c(1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   Delta <- AR(0.5,p)
#   
#   nz_vec <- 1:5
#   True_sp <- solve(Delta) %*% Gamma
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(y + 1/5*y^2, 2*abs(y))
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }

# #############  Model (XXVII4) #############
# Model27_4 <- function(p=100){
#   
#   d <- 2
#   tmp <- matrix(0, p, d)
#   tmp[1:4,1] <- 1
#   tmp[1:4,2] <- c(1,-1,1,-1)
#   tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
#   tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
#   Delta <- AR(0.5,p)
#   
#   nz_vec <- 1:4
#   Gamma <- Delta %*% tmp
#   True_sp <- tmp
#   
#   Data <- function(N){
#     y <- rnorm(N,0,1)
#     v <- cbind(y, y^2)
#     eps <- Train(N, rep(0,p), Delta)
#     x <- v %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   intra_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), H = 5)
#   pfc_params <- list(lambda.factor = 0.5, lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.5, length.out = 10), cut_y = TRUE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }



