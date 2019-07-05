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
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
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
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model III #############
Model3 <- function(p=100){
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
    y <- abs((x %*% Beta[,1]) / 4 + 2)^3 * sign(x %*% Beta[,2]) + 0.2 * rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model IV2 #############
Model4 <- function(p=100){
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
    y <- abs((x %*% Beta[,1]) / 4 + 2)^3 * sign(x %*% Beta[,2]) + rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, cut_y = TRUE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VI #############
Model5 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2.5)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.2 *  rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VI2 #############
Model6 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2.5)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.4 *  rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VI3 #############
Model7 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2.5)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.6 *  rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VI4 #############
Model8 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2.5)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.8 *  rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VI5 #############
Model9 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:6,1] <- 1
  Beta[1:6,2] <- c(1,-1,1,-1,1,-1)
  Beta[,1] <- sqrt(0.5)*Beta[,1]/norm(Beta[,1], '2')
  Beta[,2] <- sqrt(2.5)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VIII #############
Model10 <- function(p=100){
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
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VIII2 #############
Model11 <- function(p=100){
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
    y <- (x %*% Beta[,1]) * exp(x %*% Beta[,2]) + 0.4 * rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VIII3 #############
Model12 <- function(p=100){
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
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VIII4 #############
Model13 <- function(p=100){
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
    y <- (x %*% Beta[,1]) * exp(x %*% Beta[,2]) + 0.8 * rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model VIII5 #############
Model14 <- function(p=100){
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
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model X #############
Model15 <- function(p=100){
  Mu <- rep(0,p)
  Sigma <- AR(0.5, p)
  
  # Construct true Beta
  Beta <- matrix(0, p, 2)
  Beta[1:8,1] <- 1
  Beta[9:12,2] <- 1
  
  nz_vec <- 1:12
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- (x %*% Beta[,1]) * (2 + (x %*% Beta[,2]) / 3)^2 + rnorm(nobs)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model PFC #############
Model16 <- function(p=100){

  Delta <- AR(0.5,p)

  d <- 2
  r <- 4
  tmp <- matrix(0, p, d)
  tmp[1:6,1] <- 1
  tmp[1:6,2] <- c(1,-1,1,-1,1,-1)
  tmp[,1] <- tmp[,1]/norm(tmp[,1], '2')
  tmp[,2] <- tmp[,2]/norm(tmp[,2], '2')
  Gamma <- Delta %*% tmp

  Beta <- matrix(runif(d*r), d, r)

  # Construct true Beta
  nz_vec <- 1:6
  True_sp <- tmp

  Data <- function(N){
    y <- runif(N,0,4)
    Fmat <- cbind(y,y^2,y^3,exp(y))
    eps <- Train(N, rep(0,p), Delta)
    # mu <- rnorm(p)
    # mu <- matrix(rep(mu,N), N, p, byrow = TRUE)
    x <- Fmat %*% t(Beta) %*% t(Gamma) + eps
    list(x = x, y = y)
  }

  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.5, lam_fac_ssdr = 0.6, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.5, lam_fac_ssdr = 0.6, H = 5)
  pfc_params <- list(lambda.factor = 0.2, lam_fac_msda = 0.5, lam_fac_ssdr = 0.6, cut_y = FALSE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

# #######################################
# Model17 <- function(p=100){
#   d <- 2
#   r <- 4
#   
#   Delta <- AR(0.5, p)
#   Gamma <- matrix(0, p, 2)
#   Gamma[1:6,1] <- 1
#   Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   
#   Beta <- matrix(runif(d*r), d, r)
#   nz_vec <- 1:6
#   True_sp <- solve(Delta) %*% Gamma
#   
#   Data <- function(N){
#     y <- runif(N,0,4)
#     p <- dim(Gamma)[1]
#     nobs <- length(y)
#     Fmat <- cbind(y, y^2, y^3, exp(y))
#     eps <- mvrnorm(nobs, rep(0,p), Delta)
#     # mu <- rnorm(p)
#     # mu <- matrix(rep(mu,N), N, p, byrow = TRUE)
#     x <- Fmat %*% t(Beta) %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.9, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
#   intra_params <- list(lambda.factor = 0.9, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
#   pfc_params <- list(lambda.factor = 0.9, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, cut_y = FALSE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }
# 
# #######################################
# 
# Model18 <- function(p=100){
#   d <- 2
#   r <- 2
#   
#   # Delta <- AR(0.5, p)
#   Delta <- 10*diag(rep(1,p))
#   Gamma <- matrix(0, p, 2)
#   Gamma[1:4,1] <- c(1,1,-1,-1)
#   Gamma[1:5,2] <- c(1,0,1,0,1,)
#   Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
#   Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
#   
#   Beta <- diag(1,d,d)
#   nz_vec <- 1:5
#   True_sp <- solve(Delta) %*% Gamma
#   
#   Data <- function(N){
#     y <- rnorm(N,0,0.5)
#     nobs <- length(y)
#     Fmat <- cbind(y, abs(y))
#     eps <- mvrnorm(nobs, rep(0,p), Delta)
#     x <- Fmat %*% t(Beta) %*% t(Gamma) + eps
#     list(x = x, y = y)
#   }
#   
#   sir_params <- list(lambda.factor = 0.2, lam_fac_msda = 0.8, lam_fac_ssdr = 0.8, H = 5)
#   intra_params <- list(lambda.factor = 0.2, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
#   pfc_params <- list(lambda.factor = 0.2, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = FALSE)
#   
#   return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
# }