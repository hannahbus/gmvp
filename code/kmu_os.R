# Out of Sample Performance Model KMU ====
# Date: 18.07.20 
# Author: HB
# Contact: hbusshof@smail.uni-koeln.de

# library(MASS)
# library(stats)

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
N_0 <- 10^(-3) * diag(k)

# Parameters
alpha <- 1/100
lambda <- 0.99 
nu <- 100 

result <- kmu_os(N = N, n_min = 3000,
                                window = window, n_gibbs = n_gibbs, 
                                burn = burn, returns = returns, 
                                baseline = baseline, stocks = stocks, 
                                selection = selection, k = k, Q = Q,  
                                alpha = alpha, lambda = lambda, nu = nu, 
                                N_0 = N_0) 
result