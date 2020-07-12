# Implementation Taylor 

library(doParallel)

library(MASS)
library(stats)

library(tidyverse)
library(reshape2)
library(plotly)

setwd("~/Desktop/bayesian-gmvp/Computation/")

set.seed(123)
rm(list = ls()) 

source("helpers_kmt.R")
source("helpers_general.R")

load("returns.Rda")
n_max <- 1000 #nrow(returns)
mcoptions <- list(preschedule = FALSE, set.seed= TRUE)

# Prepping ====
selection <- c("ba", "hd", "ko", "jpm", "ibm")
# selection <- c(setdiff(names(returns), "date")) 
baseline <- "ba"
df <- returns[1:n_max, names(returns) %in% selection]
stocks <- c(setdiff(selection, baseline))
df[1:n_max, stocks] <- - returns[1:n_max, stocks] + returns[1:n_max, baseline]
model_1 <- lm(data = df, ba ~ .)
summary(model_1)

X_series <- as.matrix(cbind(1, df[, stocks]), ncol = k)
y_series <- as.matrix(df[, baseline], ncol = k)
n_T <-  nrow(df)
k <- length(selection)

# Global information 
no_mixture <- 10 
n_gibbs <- 10000 
weights <- k_10[, 1, drop = FALSE] 
weights = matrix(rep(weights, n_T), nrow = n_T, byrow = T)

# Initialize 
beta_0 <- matrix(c(model_1$coefficients), nrow = k)
h_series <-  matrix(rep(1, n_T), nrow = n_T, ncol = 1)
s <- matrix(sample(1:no_mixture, n_T, replace = TRUE), nrow = n_T)

a_mcmc <- array(0, dim = c(n_T, 1, n_gibbs)) 
beta_mcmc <- array(0, dim = c(n_T, k, n_gibbs))
mu_h_mcmc <- sigma_vartheta_mcmc <- phi_mcmc <- matrix(0, nrow = n_gibbs, ncol = 1)

# "Hyper parameters" to sample parameters of interest. 

mu_h <- 0

a <- 5/2
b <- 0.05/2 

alpha_phi <- 20
beta_phi  <- 1.5

mu_h_0 <- 0 
sd_h_0 <- 10  

initial <- c(mu_h, 0.1) 
phi <- 0.98 
sigma_vartheta <- 0.1
sigma_eta <-  solve(t(X_series) %*% X_series)
S0_inv <- k * sigma_eta
A_0 <- diag(k)
beta_0_prior <- beta_0

for (i in 1:n_gibbs){
  print(i)
  # Update beta 
  filtered <- filter_taylor(k = k, n_T = n_T, beta_0 = beta_0, 
                            X_series = X_series, y_series = y_series,
                            h_series = h_series, sigma_eta = sigma_eta)
  print("Success Filter!")
  beta_mcmc[, , i] <- beta_series <- smooth_taylor(n_T = n_T, k = k, 
                                                   beta_m = filtered$beta_m, 
                               S_m = filtered$S_m, 
                               S_p = filtered$S_p, beta_p = filtered$beta_p)
  print("Success Smooth!")
  # Update transformed h_series 
  y_star <- log((y_series - matrix(rowSums((X_series * beta_series)), 
                                   nrow = n_T))^2 + 10^(-4))
  result <- ff_a_taylor(y_star = y_star, s = s, phi = phi, mu_h = mu_h, 
                        sigma = sigma_vartheta, initial = initial, n_T = n_T, 
                        k = k_10)
  a_mcmc[, , i] <- a_smooth <- bs_a_taylor(phi, a_measure = result$a_measure, 
                                         a_predict = result$a_predict, 
                                         S_measure = result$S_measure, 
                                         S_predict = result$S_predict, 
                                         n_T = n_T)
  h_series <- exp(a_smooth * 0.5) 
  print("success h_series")
  # Update indicator 
  s <- matrix(c(update_indicator(y_star, a_smooth, k = k_10, no_mixture,
                                 weights = weights, n_T = n_T)), ncol = 1)
  print("success indicator update")
  # Update hyper parameters 
  # sigma_vartheta 
  sigma_vartheta <- sigma_vartheta_mcmc[i, 1] <- sample_sigma_vartheta(a = a, 
                                                              n_T = n_T, b = b, 
                                phi = phi, mu_h = mu_h, logh_series = a_smooth)
  print("Success sigma_vartheta!")
  # phi 
  phi <- phi_mcmc[i, 1] <- sample_phi(logh_series = a_smooth, 
                                      mu_h = mu_h, 
                                      sigma_vartheta = sigma_vartheta, 
                                      alpha_phi = alpha_phi, 
                                      beta_phi = beta_phi, 
                                      phi_old = phi, n_T = n_T)
  print("Success phi!")
  # mu_h 
  mu_h_mcmc[i, 1] <- mu_h <- sample_mu_h(mu_h_0, sd_h_0, logh_series = a_smooth, 
                                         phi = phi, 
                                         sigma_vartheta = sigma_vartheta, 
                                         n_T = n_T)
  sigma_eta <- sample_sigma_eta(n_T = n_T, k = k, beta = beta_series, 
                                beta_0 = beta_0, S0_inv = S0_inv)
  # sample_logh0(mu_h, phi, logh_1, sigma_vartheta)
  beta_0 <-  sample_beta_0(sigma_eta = sigma_eta, A_0 = A_0, 
                           beta_1 = t(beta_series[1, ,drop = F]), 
                           beta_0_prior = beta_0_prior)
}

beta_result <- apply(beta_mcmc[, ,(n_gibbs * 0.5):(n_gibbs)], c(1, 2), mean)
returns_selected <- returns[1:n_max, names(returns) %in% selection]
performance_1 <- compute_performance_is(beta_smooth = beta_result, k = k, 
                                        returns = returns_selected) 
performance_2 <- var(rowMeans(returns[1:n_max, names(returns) %in% selection]))
performance_1
performance_2

beta_variance <- apply(beta_mcmc[, ,(n_gibbs * 0.5):(n_gibbs)], c(1, 2), var)
plot_ly(x = paste("beta", 0:(5-1), sep = "_"),
        y = returns$date, z = beta_variance, type = "heatmap", 
        colors = colorRamp(c("white", "black")))


beta_result <- as.data.frame(beta_result)
names(beta_result) <- c("beta_0", stocks)
beta_result$date <- returns[1:n_max, "date"]
df1 <- beta_result %>% gather(key = "variable", value = "value", -date) %>% 
  filter(variable %in% c("hd", "ko", "jpm", "ibm"))
p1 <- ggplot(df1, aes(x = date, y = value)) + geom_line(aes(color = variable))
ggplotly(p1) %>% layout(title = "Betas")

parameters_gibbs <- data.frame(mu_h_mcmc, sigma_vartheta_mcmc, phi_mcmc, 
                               c(1:n_gibbs))
names(parameters_gibbs) <- c("mu_h", "sigma_vartheta", "phi", "round")

p3 <- ggplot(parameters_gibbs) + geom_line(aes(x = round, y = mu_h)) + 
  labs(x = "round", y = "mu_h")
p4 <- ggplot(parameters_gibbs[, ]) + geom_line(aes(x = round, y = sigma_vartheta)) + 
  labs(x = "round", y = "sigma_vartheta") 
p5 <- ggplot(parameters_gibbs) + geom_line(aes(x = round, y = phi)) + 
  labs(x = "round", y = "phi") 

p3
p4
p5 