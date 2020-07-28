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

returns <- returns[1:n_T, names(returns) %in% selection]
a1 <- analysis_is_performance(result_1$betas, n_gibbs, burn_factor, returns, 
                              choice, plot = F)
a2 <- analysis_is_performance(result_2$betas, n_gibbs, burn_factor, returns, 
                              choice, plot = F)
a3 <- analysis_is_performance(result_3$betas, n_gibbs, burn_factor, returns, 
                              choice, plot = F)

estimate_mc <- function(object, n_gibbs, burn_factor, scalar = T){
  if (scalar) {
    return(mean(object[(n_gibbs * burn_factor + 1): n_gibbs, 1]))
  } else {
    return(apply(object[, , (n_gibbs * burn_factor + 1): n_gibbs], c(1, 2), mean))
  }
}

effective_size_mc <- function(object, n_gibbs, burn_factor){
  n_gibbs / spectrum0(object[(n_gibbs * burn_factor + 1): n_gibbs, 1])$spec
}

nu1 <- estimate_mc(result_1$nu, n_gibbs = n_gibbs, burn_factor = burn_factor)
nu2 <- estimate_mc(result_2$nu, n_gibbs = n_gibbs, burn_factor = burn_factor)
nu3 <- estimate_mc(result_3$nu, n_gibbs = n_gibbs, burn_factor = burn_factor)

ef1 <- effective_size_mc(result_1$nu,  n_gibbs = n_gibbs, burn_factor = burn_factor)
ef2 <- effective_size_mc(result_2$nu,  n_gibbs = n_gibbs, burn_factor = burn_factor)
ef3 <- effective_size_mc(result_3$nu,  n_gibbs = n_gibbs, burn_factor = burn_factor)

performances <- matrix(0, nrow = 3, ncol = 6)
colnames(performances) <- c("lambda", "nu", "nu_gibbs", "effective_size" ,"performance", 
                            "time")
performances[, 1] <- c(0.99, 0.98, 0.9)
performances[, 2] <- c(100, 50, 10 )
performances[, 3] <- c(nu1, nu2, nu3)
performances[, 4] <- c(ef1, ef2, ef3)
performances[, 5] <- c(a1, a2, a3)
performances[, 6] <- c(result_1$duration, result_2$duration, result_3$duration)


performances

# Plots nu
plot(1:n_gibbs, result_1$nu, "l")
acf(result_1$nu)

plot(1:n_gibbs, result_2$nu, "l")
acf(result_2$nu)

plot(1:n_gibbs, result_3$nu, "l")
acf(result_3$nu)

# Plots h 
h1 <- apply(result_1$h, 1, mean)
h2 <- apply(result_2$h, 1, mean)
h3 <- apply(result_3$h, 1, mean)

ggplot() + geom_line(aes(x = 1:4447, y = 1 / h1), color = "blue") + 
              geom_line(aes(x = 1:4447, y = 1 / h2), color = "red") + 
                geom_line(aes(x = 1:4447, y = 1 / h3), color = "black") + 
  xlab("Day") + 
  ylab("Volatility")

# Plots selected assets 

beta_result <- apply(result[, ,(n_gibbs * burn_factor + 1):(n_gibbs)], 
                     c(1, 2), mean) 
beta_result <- as.data.frame(beta_result)
names(beta_result) <- c("beta_0", stocks)
beta_result$date <- returns[1:n_T, "date"]
df <- beta_result %>% gather(key = "variable", value = "value", -date)
df_plot <- df %>% filter(variable %in% choice)
p1 <- ggplot(df_plot, aes(x = date, y = value)) + 
  geom_line(aes(color = variable))
