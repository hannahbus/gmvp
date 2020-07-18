# Out of Sample Performance Model KMU ====
# Date: 18.07.20 
# Author: HB
# Contact: hbusshof@smail.uni-koeln.de

library(MASS)
library(stats)

setwd("Desktop/bayesian-gmvp/Computation/")
rm(list = ls())

source("helpers_kmu.R")
source("helpers_kmu.R")
load("returns.Rda")

N <- nrow(returns)
window <- 240 
n_gibbs <- 1000 
burn <- 100
selection <- c("ba", "hd", "ko", "jpm", "ibm")
baseline <- "ba"
stocks <- c(setdiff(selection, baseline))
k <- length(selection)
Q <- 1000 * diag(k)

# Parameters
C_inv <-  diag(k)
b <- rep((1/k), k)
d <- 1 
e <- 1
alpha <- 1/100
lambda <- 0.99 
nu <- 100 

result <- kmu_os_rolling_window(N = N, n_min = 3000,
                                window = window, n_gibbs = n_gibbs, 
                                burn = burn, returns = returns, 
                                baseline = baseline, stocks = stocks, 
                                selection = selection, k = k, Q = Q, 
                                C_inv = C_inv, b = b, d = d, e = e, 
                                alpha = alpha, lambda = lambda, nu = nu) 
result