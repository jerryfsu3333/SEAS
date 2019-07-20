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
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.2 *  rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

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
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.4 *  rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.5, cut_y = TRUE)
  
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
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.6 *  rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.6, cut_y = TRUE)
  
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
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + 0.8 *  rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.6, cut_y = TRUE)
  
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
  Beta[,2] <- sqrt(2)*Beta[,2]/norm(Beta[,2], '2')
  
  nz_vec <- 1:6
  True_sp <- Beta
  
  Data <- function(N){
    x <- Train(N, Mu, Sigma)
    nobs <- dim(x)[1]
    y <- x %*% Beta[,1] * exp(x %*% Beta[,2] + rnorm(nobs) )
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.6, cut_y = TRUE)
  
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
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)

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
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
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
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
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
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.6, cut_y = TRUE)
  
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
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.9, lam_fac_ssdr = 0.7, cut_y = TRUE)
  
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

#############  Model XI(PFC) #############
# Comment: cant increase too much on lambda.factor and lam_fac_msda, which will lead to a undesirable
# reduction on the sparsity.
Model16 <- function(p=100){

  # Delta <- AR(0.5,p)
  Delta <- 0.1^2*diag(1,p,p)

  d <- 2
  r <- 2
  Gamma <- matrix(0, p, d)
  Gamma[1:6,1] <- 1
  Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Gamma <- Delta %*% tmp

  Beta <- matrix(runif(d*r), d, r)

  # Construct true Beta
  nz_vec <- 1:6
  True_sp <- Gamma

  Data <- function(N){
    y <- runif(N,0,4)
    Fmat <- cbind(y,exp(y))/2
    eps <- Train(N, rep(0,p), Delta)
    # mu <- rnorm(p)
    # mu <- matrix(rep(mu,N), N, p, byrow = TRUE)
    x <- Fmat %*% t(Beta) %*% t(Gamma) + eps
    list(x = x, y = y)
  }

  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.5, lam_fac_ssdr = 0.5, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.5, lam_fac_ssdr = 0.5, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.5, lam_fac_ssdr = 0.5, cut_y = FALSE)

  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model XII (PFC) #############
Model17 <- function(p=100){
  
  # Delta <- AR(0.5,p)
  Delta <- 0.1^2*diag(1,p,p)
  
  d <- 2
  r <- 3
  Gamma <- matrix(0, p, d)
  Gamma[1:6,1] <- 1
  Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Gamma <- Delta %*% tmp
  # Gamma <- tmp
  
  Beta <- matrix(runif(d*r), d, r)
  
  # Construct true Beta
  nz_vec <- 1:6
  True_sp <- Gamma
  
  Data <- function(N){
    y <- runif(N,0,4)
    Fmat <- cbind(y,y^2,y^3)
    eps <- Train(N, rep(0,p), Delta)
    x <- Fmat %*% t(Beta) %*% t(Gamma) + eps
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.5, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.5, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.5, cut_y = FALSE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model PFC #############
Model18 <- function(p=100){
  
  # Delta <- AR(0.5,p)
  Delta <- 1^2*diag(1,p,p)
  
  d <- 2
  r <- 3
  Gamma <- matrix(0, p, d)
  Gamma[1:6,1] <- 1
  Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Gamma <- Delta %*% tmp
  # Gamma <- tmp
  
  Beta <- matrix(runif(d*r), d, r)
  
  # Construct true Beta
  nz_vec <- 1:6
  True_sp <- Gamma
  
  Data <- function(N){
    y <- runif(N,0,4)
    Fmat <- cbind(y,y^2,y^3)
    eps <- Train(N, rep(0,p), Delta)
    # mu <- rnorm(p)
    # mu <- matrix(rep(mu,N), N, p, byrow = TRUE)
    x <- Fmat %*% t(Beta) %*% t(Gamma) + eps
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.6, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.6, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.6, cut_y = FALSE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model PFC #############
Model19 <- function(p=100){
  
  # Delta <- AR(0.5,p)
  Delta <- 10^2*diag(1,p,p)
  
  d <- 2
  r <- 3
  Gamma <- matrix(0, p, d)
  Gamma[1:6,1] <- 1
  Gamma[1:6,2] <- c(1,-1,1,-1,1,-1)
  Gamma[,1] <- Gamma[,1]/norm(Gamma[,1], '2')
  Gamma[,2] <- Gamma[,2]/norm(Gamma[,2], '2')
  # Gamma <- Delta %*% tmp
  # Gamma <- tmp
  
  Beta <- matrix(runif(d*r), d, r)
  
  # Construct true Beta
  nz_vec <- 1:6
  True_sp <- Gamma
  
  Data <- function(N){
    y <- runif(N,0,4)
    Fmat <- cbind(y,y^2,y^3)
    eps <- Train(N, rep(0,p), Delta)
    # mu <- rnorm(p)
    # mu <- matrix(rep(mu,N), N, p, byrow = TRUE)
    x <- Fmat %*% t(Beta) %*% t(Gamma) + eps
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.6, H = 5)
  intra_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.6, H = 5)
  pfc_params <- list(lambda.factor = 0.6, lam_fac_msda = 0.6, lam_fac_ssdr = 0.6, cut_y = FALSE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model PFC #############
Model20 <- function(p=10){
  
  sigmaY <- 5
  sigma0 <- 1
  sigma <- 1
  Gamma <- matrix(rep(0,p), p, 1)
  Gamma[1] <- 1
  Gamma0 <- qr.Q(qr(Gamma), complete = TRUE)[,2:p]
  
  # Construct true Beta
  nz_vec <- 1
  True_sp <- Gamma
  
  Data <- function(N){
    y <- rnorm(N,0,sigmaY)
    eps <- matrix(rnorm(N*p),N,p)
    x <- y %*% t(Gamma) + sigma0 * eps[,2:p] %*% t(Gamma0) + sigma * eps[,1] %*% t(Gamma)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, cut_y = FALSE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model PFC #############
Model21 <- function(p=10){
  
  sigmaY <- 1
  sigma0 <- 1
  sigma <- 5
  Gamma <- matrix(rep(0,p), p, 1)
  Gamma[1] <- 1
  Gamma0 <- qr.Q(qr(Gamma), complete = TRUE)[,2:p]
  
  # Construct true Beta
  nz_vec <- 1
  True_sp <- Gamma
  
  Data <- function(N){
    y <- rnorm(N,0,sigmaY)
    eps <- matrix(rnorm(N*p),N,p)
    x <- y %*% t(Gamma) + sigma0 * eps[,2:p] %*% t(Gamma0) + sigma * eps[,1] %*% t(Gamma)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, cut_y = FALSE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model PFC #############
Model22 <- function(p=10){
  
  sigmaY <- 1
  sigma0 <- 3
  sigma <- 1
  Gamma <- matrix(rep(0,p), p, 1)
  Gamma[1] <- 1
  Gamma0 <- qr.Q(qr(Gamma), complete = TRUE)[,2:p]
  
  # Construct true Beta
  nz_vec <- 1
  True_sp <- Gamma
  
  Data <- function(N){
    y <- rnorm(N,0,sigmaY)
    eps <- matrix(rnorm(N*p),N,p)
    x <- y %*% t(Gamma) + sigma0 * eps[,2:p] %*% t(Gamma0) + sigma * eps[,1] %*% t(Gamma)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.5, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, cut_y = FALSE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}

#############  Model intra-slice #############
Model23 <- function(p=5){
  
  Beta1 <- c(1,0,0,0,0)
  Beta2 <- c(0,1,1,0,0)
  Beta <- cbind(Beta1, Beta2)
  
  # Construct true Beta
  nz_vec <- 1:3
  True_sp <- Beta
  
  Data <- function(N){
    v1 <- rt(N,5)
    v2 <- rt(N,5)
    v3 <- rt(N,5)
    w1 <- rgamma(N, shape = 0.2)
    w2 <- rgamma(N, shape = 0.2)
    x <- cbind(w1, v1 + w2/2, -v1 + w2/2, v2 + v3, v2 - v3)
    y <- 1.5 * (5 + x %*% Beta[,1]) * (2 + x %*% Beta[,2]) + 0.5 * rnorm(N)
    list(x = x, y = y)
  }
  
  sir_params <- list(lambda.factor = 0.7, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  intra_params <- list(lambda.factor = 0.7, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, H = 5)
  pfc_params <- list(lambda.factor = 0.7, lam_fac_msda = 0.9, lam_fac_ssdr = 0.8, cut_y = FALSE)
  
  return(list(Data = Data, True_sp = True_sp, nz_vec = nz_vec, sir_params = sir_params, intra_params = intra_params, pfc_params = pfc_params))
}
