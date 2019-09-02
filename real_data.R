rm(list = ls())
library(msda)
library(MASS)
library(methods)
library(glmnet)
library(energy)
source("/Users/cengjing/Documents/GitHub/ssdr/msda_prep.R")
source("/Users/cengjing/Documents/GitHub/ssdr/utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/my_msda.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_utility.R")
source("/Users/cengjing/Documents/GitHub/ssdr/ssdr_func.R")

data <- read.csv('/Users/cengjing/Documents/GitHub/ssdr/wine.data', header = FALSE, col.names = 
                c('class', 'Alc', 'Malic', 'Ash', 'Alka', 'Mag', 'Tot_ph', 'Fla', 'NonFla', 'Proant', 'Col', 'Hue', 'OD', 'Proline'))

y <- data[,1]
x <- data[,-1]

ssdr_func <- function(x, y, x_val, y_val, H=5, type = 'sir', lambda.factor=0.5, nlam_msda=10,
                      lam1_fac=seq(1.2,0.01, length.out = 10), lam2_fac=seq(0.001,0.2, length.out = 10),
                      gamma=c(10,30,50), cut_y=TRUE)
