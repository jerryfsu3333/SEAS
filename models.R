#############  Model I #############
Model1 <- function(p=100, N=500, N_val=500){
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
Model2 <- function(p=100, N=500, N_val=500){
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
Model3 <- function(p=100, N=500, N_val=500){
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
Model4 <- function(p=100, N=500, N_val=500){
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
Model5 <- function(p=100, N=500, N_val=500){
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
Model6 <- function(p=100, N=500, N_val=500){
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
Model7 <- function(p=100, N=500, N_val=500){
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
Model8 <- function(p=100, N=500, N_val=500){
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
Model9 <- function(p=100, N=500, N_val=500){
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
Model10 <- function(p=100, N=500, N_val=500){
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
Model11 <- function(p=100, N=500, N_val=500){
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
Model12 <- function(p=100, N=500, N_val=500){
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
Model13 <- function(p=100, N=500, N_val=500){
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
Model14 <- function(p=100, N=500, N_val=500){
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
Model15 <- function(p=100, N=500, N_val=500){
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

