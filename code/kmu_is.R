# Empirical Application ====
# Implementation of Uhlig SV
# In Sample Analysis
# Date: 20.06.20 
# Author: HB
# Contact: hbusshof@smail.uni-koeln.de

library(MASS)
library(stats)
library(coda)
library(mcmcse)

library(tidyverse)
library(reshape2)
library(plotly)

# setwd("~/Desktop/bayesian-gmvp/Computation/")

set.seed(123)
rm(list = ls())
source("helpers_kmu.R")
source("helpers_general.R")

load("returns.Rda")
n_T <- nrow(returns)

# Prepping ====
selection <- c("ba", "hd", "ko", "jpm", "ibm")
# selection <- c(setdiff(names(returns), "date"))
baseline <- "ba"
stocks <- c(setdiff(selection, baseline))
k <- length(selection)
alpha <- 1/100
n_gibbs <- 5000
burn_factor <- 0.5
N_0 <- 10^(-3) * diag(k)
Q <-  1000 * diag(k)

lambda <- 0.99 
nu <- 100 

result_1 <- kmu_is(returns = returns, n_T = n_T, 
                 selection = selection, baseline = baseline, 
                 stocks  = stocks, k = k, Q = Q, N_0 = N_0, alpha = alpha, 
                 lambda = lambda, nu = nu, n_gibbs = n_gibbs, 
                 burn_factor = burn_factor, nu_fix = F)
lambda <- 0.98 
nu <- 50 

result_2 <- kmu_is(returns = returns, n_T = n_T, 
                   selection = selection, baseline = baseline, 
                   stocks  = stocks, k = k, Q = Q, N_0 = N_0, alpha = alpha, 
                   lambda = lambda, nu = nu, n_gibbs = n_gibbs, 
                   burn_factor = burn_factor, nu_fix = F)

lambda <- 0.9 
nu <- 10 

result_3 <- kmu_is(returns = returns, n_T = n_T, 
                   selection = selection, baseline = baseline, 
                   stocks  = stocks, k = k, Q = Q, N_0 = N_0, alpha = alpha, 
                   lambda = lambda, nu = nu, n_gibbs = n_gibbs, 
                   burn_factor = burn_factor, nu_fix = F)

performance_naive <- var(rowMeans(returns[1:n_T, names(returns) %in% 
                                            selection]))

a1 <- analysis_is_performance(result_1$betas, n_gibbs, burn_factor, returns, 
                              choice, plot = F)
a2 <- analysis_is_performance(result_2$betas, n_gibbs, burn_factor, returns, 
                              choice, plot = F)
a3 <- analysis_is_performance(result_3$betas, n_gibbs, burn_factor, returns, 
                              choice, plot = F)

performances <- matrix(0, nrow = 3, ncol = 3)
colnams(performances) <- c("lambda", "nu", "performance")
performances[, 1] <- c(0.99, 0.98, 0.9)
performances[, 2] <- c(100, 50, 10 )
performances[, 3] <- c(a1, a2, a3)
performances