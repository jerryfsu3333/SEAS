rm(list = ls())
library(msda)
library(MASS)
library(methods)
library(glmnet)
library(energy)
library(caret)
library(LassoSIR)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/my_msda.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")

data <- read.csv('/Users/cengjing/Documents/GitHub/ssdr/wine.data', header = FALSE, col.names = 
                c('class', 'Alc', 'Malic', 'Ash', 'Alka', 'Mag', 'Tot_ph', 'Fla', 'NonFla', 'Proant', 'Col', 'Hue', 'OD', 'Proline'))

y <- data[,1]
x <- as.matrix(data[,-1])
x <- scale(x)
  
fit_lda <- MASS::lda(x,y)
pred <- predict(fit_lda, x)$class
pred <- as.numeric(levels(pred)[pred])
sum(y == pred)/length(y)

fit <- ssdr.cv(x, y, lam1_fac = seq(1,0.01, length.out = 10), categorical=TRUE, type = 'sir')
d <- fit$rank
directions <- svd(fit$mat)$u[,1:d, drop=FALSE]
x_new <- as.matrix(x) %*% directions
plot(x_new[,1], x_new[,2], col=y)

fit_lda <- MASS::lda(x_new,y)
pred <- predict(fit_lda, x_new)$class
pred <- as.numeric(levels(pred)[pred])
sum(y == pred)/length(y)


LassoSIR_fit <- LassoSIR(as.matrix(x), y, categorical = TRUE, no.dim = 2)
d <- LassoSIR_fit$no.dim
x_new <- as.matrix(x) %*% LassoSIR_fit$beta
x_new <- scale(x_new, center = FALSE)
plot(x_new[,1], x_new[,2], col=y)

# fit_lda <- MASS::lda(x_new,y)
# pred <- predict(fit_lda, x_new)$class
# pred <- as.numeric(levels(pred)[pred])
# sum(y == pred)/length(y)
