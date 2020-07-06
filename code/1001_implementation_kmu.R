# Empirical Application ====
# Implementation of Uhlig SV
# Date: 20.06.20 
# Author: HB
# Contact: hbusshof@smail.uni-koeln.de

library(doParallel)

library(MASS)
library(stats)

library(tidyverse)
library(reshape2)
library(plotly)

setwd("~/Desktop/bayesian-gmvp/Computation/")

rm(list = ls())
source("helpers_kmu.R")
source("helpers_general.R")

load("returns.Rda")
n_max <- 2000
mcoptions <- list(preschedule = FALSE, set.seed= TRUE)

# Prepping ====
selection <- c("ba", "hd", "ko", "jpm", "ibm")
# selection <- c(setdiff(names(returns), "date")) 
baseline <- "ba"
df <- returns[1:n_max, names(returns) %in% selection]
stocks <- c(setdiff(selection, "ba"))
df[1:n_max, stocks] <- - returns[1:n_max, stocks] + returns[1:n_max, "ba"]
model_1 <- lm(data = df, ba ~ .)
summary(model_1)

X_series <- as.matrix(cbind(1, df[, stocks]), ncol = k)
y_series <- as.matrix(df[, "ba"], ncol = k)
n_T <- horizon <- nrow(df)
k <- length(selection)

# Gibbs sampler ====

# Initialize 
b_10 <- matrix(c(model_1$coefficients), nrow = k) 
S_10 <- h_0 <- 1.5  
N_10 <- Q <- diag(k) 
nu <- 100 
lambda <- nu / (nu + 1)
S0_inv <-  k * solve(t(X_series) %*% X_series)
C_inv <-  diag(k)
b <- rep((1/k), k)
d <- 1 
e <- 1
alpha <- 1/200
lambda <- 0.98 
nu <- 50 

n_gibbs <- 1000
density_lambda <- h0_mcmc <- lambda_mcmc <- nu_mcmc <- matrix(0, nrow  = n_gibbs)
beta_mcmc <-  array(0, dim = c(n_max, k, n_gibbs))
Q_mcmc <- array(0, dim = c(k, k, n_gibbs))
h_mcmc <- array(0, dim = c(n_max, ncol = 1, n_gibbs))

start <- Sys.time()

for (i in 1:n_gibbs){
  print(i)
  # FFBS in states 
  result_filter <- filter_forward(X = X_series, y = y_series, nu = nu, Q = Q, 
                                  lambda = lambda, N_10 = N_10, S_10 = S_10, 
                                  b_10 = b_10, n_T = n_T, k = k) 
  result_bsample <- backward_sample(S_m = result_filter$S_m, 
                                    b_m = result_filter$beta_m, 
                                    N_m = result_filter$N_m, nu = nu, Q = Q, 
                                    lambda = lambda, n_T = n_T, k = k)
  beta_mcmc[ , , i] <- result_bsample[ , 2:(k + 1)]
  h_mcmc[ , , i]    <- result_bsample[ , 1]
  print("success FFBS")
  # Q 
  N_10 <- Q <- Q_mcmc[, , i] <- sample_Q(beta_0 = b_10, theta = result_bsample, 
                                         n_T = n_T, k = k, S0_inv = S0_inv)
  print("success Q")
  h_1 <- drop(result_bsample[1, 1])
  h_T <- drop(result_bsample[n_T, 1])
  beta_1 <- matrix(result_bsample[1, 2:(k+1)],nrow = k)
  # nu 
  nu_log <- logp_nu(nu = nu, h_T = h_T, h_0 = h_0, k = k, lambda = lambda, 
                    n_T = n_T, alpha = alpha)
  nu_mcmc[i, 1] <- nu <- nu_sample(nu_old = nu, logdensity_old = nu_log, 
                                   sd = 10, h_T = h_T, 
                                   h_0 = h_0, k = k, lambda = lambda, n_T = n_T,
                                   alpha = alpha)
  print("success nu")
  # lambda
  density_lambda[i, 1] <- logdensity_lambda <- logp_lambda(lambda, nu = nu, 
                                   h = result_bsample[, 1, drop = F], 
                                   h_0 = h_0, k = k, n_T = n_T)
  lambda <- lambda_mcmc[i, 1] <- sample_lambda(lambda, 
                                    logdensity_old = logdensity_lambda,
                                     sd = 0.00000001, nu = nu, 
                                     h = result_bsample[, 1, drop = F], 
                                     h_0 = h_0, k = k, n_T = n_T)
  print("success lambda")
  # beta_0 
  b_10 <- beta_0_sample(C_inv = C_inv, h_1= h_1, Q = Q, beta_1 = beta_1, b = b)
  print("success b_10")
  # h_0 
  log_h0 <- logp_h_0(h_0 = h_0, h_1 = h_1, lambda = lambda, nu = nu, k = k, 
                     d = d, e = e)
  h0_mcmc[i, 1] <- S_10 <- h_0 <- h_0_sample(h0_old = h_0, 
                                             logdensity_old = log_h0, 
                         sd = 0.5, h_1 = h_1, nu = nu, k = k, d = d, e = e)
}
end <- Sys.time()
duration <- end - start 
duration

# Performance ====
beta_result <- apply(beta_mcmc[, ,(n_gibbs * 0.5):(n_gibbs)], c(1, 2), mean)
returns_selected <- returns[1:n_max, names(returns) %in% selection]
performance_1 <- compute_performance_is(beta_smooth = beta_result, k = k, 
                       returns = returns_selected) 
performance_2 <- var(rowMeans(returns[1:n_max, names(returns) %in% selection]))
performance_1
performance_2

# Sanity check 
beta_variance <- apply(beta_mcmc[], c(1, 2), var)
plot_ly(x = paste("beta", 0:(5-1), sep = "_"),
        y = returns$date, z = beta_variance, type = "heatmap", 
        colors = colorRamp(c("white", "black")))
Q_estim <- apply(Q_mcmc, c(1, 2), mean)
plot_ly(x = paste("beta", 0:(k-1), sep = "_"),
        y = paste("beta", 0:(k-1), sep = "_"),  z = solve(Q_estim), 
        type = "heatmap", 
        colors = colorRamp(c("white", "black")))

# Should we just assume homoskedasticity or should we also allow
# for heteroskedasticity! consider HPD interval for the off-diagonal elements. 

# Outcomes I ====
# Plot betas. 
beta_result <- as.data.frame(beta_result)
names(beta_result) <- c("beta_0", stocks)
beta_result$date <- returns[1:n_max, "date"]
df1 <- beta_result %>% gather(key = "variable", value = "value", -date) %>% 
          filter(variable %in% c("hd", "ko", "jpm", "ibm"))
p1 <- ggplot(df1, aes(x = date, y = value)) + geom_line(aes(color = variable))
ggplotly(p1) %>% layout(title = "Betas")

# Plot latent volatility.
h <- as.data.frame(apply(h_mcmc[, , (0.5 * n_gibbs):n_gibbs], c(1,2), mean))
names(h) <- "h"
h$date <- returns[1:n_max, "date"]
p2 <- ggplot(h, aes(x = date, y = h)) + geom_line()
ggplotly(p2) %>% layout(title = "Latent volatility")

# Outcomes II ====
# Trace plots of hyper parameters 
parameters_gibbs <- data.frame(h0_mcmc, lambda_mcmc, nu_mcmc, c(1:n_gibbs))
names(parameters_gibbs) <- c("h0", "lambda", "nu", "round")

p3 <- ggplot(parameters_gibbs) + geom_line(aes(x = round, y = h0)) + 
                labs(x = "round", y = "h0")
p4 <- ggplot(parameters_gibbs) + geom_line(aes(x = round, y = lambda)) + 
                labs(x = "round", y = "lambda") 
p5 <- ggplot(parameters_gibbs) + geom_line(aes(x = round, y = nu)) + 
                labs(x = "round", y = "nu") 
p3
p4
p5

# Watch the evolution of the weights
mean(beta_mcmc[682,44, ])
mean(beta_mcmc[11,3,])
var(beta_mcmc[11,3,])

plot(1:n_gibbs, h_mcmc[2000,1, ], "l")
plot(1:39, beta_mcmc[2,1, 1:39], "l")
plot(1:n_gibbs, beta_mcmc[300,2, ], "l")
plot(1:n_gibbs, beta_mcmc[500, ], "l")
View(density_lambda)
