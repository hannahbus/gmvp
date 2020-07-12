library(doParallel)
library(MASS)
library(stats) 

library(tidyverse)
library(reshape2)
library(plotly)
library(ggplot2)
library(grid)
library(gridExtra)
setwd("~/Desktop/bayesian-gmvp/Computation/")

set.seed(123)
rm(list = ls()) 

source("helpers_kmt.R")
source("helpers_general.R")
source("simulate_taylor.R")


mcoptions <- list(preschedule = FALSE, set.seed= TRUE)

set.seed(123) 
n_T <- 1000 
k <- 3 

gamma <- 0
phi <- 0.9
logh_0 <- gamma / (1 - phi) # Initialize to unconditional mean

sigma_vartheta <- 0.1 # Initialize to value by Platanioti (2005, p. 17) 

sigma_X <-  diag(k) # should not matter too much as observed
X_series <- mvrnorm(n_T, mu = rep(0, k), Sigma = sigma_X)
# X_series <- matrix(1, nrow = n_T, ncol = k)

beta_0 <- runif(k, min = 0, max = 1)
sigma_eta <- 0.000000000001 * diag(k) #x * diag(k)

simulation <- generate_data(n_T = n_T, sigma_vartheta = sigma_vartheta, 
                            logh_0 = logh_0, gamma = gamma, phi = phi, 
                            sigma_eta = sigma_eta, X_series = X_series, 
                            beta_0 = beta_0)

y_series <- simulation[[1]]
h_series <- simulation[[2]]^2
logh_series <- log((h_series))
beta_series <- simulation[[3]]

# Global information 
no_mixture <- 10 
n_gibbs <- 5000  
weights <- k_10[, 1, drop = FALSE] # seem to be correct
weights = matrix(rep(weights, n_T), nrow = n_T, byrow = T)

# Initialize 
beta_0 <- beta_0
S_0_beta <- diag(k)
h0_0 <-  gamma / (1 - phi)
S_0_h0 <- sigma_vartheta/(1 -  phi^2)
S0_inv <- k * sigma_eta
h_series_true <- h_series
beta_series_true <- beta_series


s <- matrix(sample(1:no_mixture, n_T, replace = T), nrow = n_T)

a_mcmc <- array(0, dim = c(n_T, 1, n_gibbs)) 
beta_mcmc <- array(0, dim = c(n_T, k, n_gibbs))

beta_0_mcmc <- mu_h_mcmc <- sigma_vartheta_mcmc <- phi_mcmc <- matrix(0, nrow = n_gibbs, 
                                                       ncol = 1)

# "Hyper parameters" to sample parameters of interest. 

a <- 5/2
b <- 0.05/2 
mu_h <- gamma / (1 - phi)
alpha_phi <- 20
beta_phi  <- 1.5
mu_h_0 <- mu_h 
sd_h_0 <- 1
A_0 <- diag(k)
beta_0_prior <- beta_0 


for (i in 1:n_gibbs){
  print(i)
  # Update beta 
  filtered <- filter_taylor(k = k, n_T = n_T, beta_0 = beta_0, 
                            S_0_beta = S_0_beta, X_series = X_series, 
                            y_series = y_series, h_series = h_series, 
                            sigma_eta = sigma_eta)
  print("Success Filter!")
  beta_mcmc[, , i] <- beta_series <- smooth_taylor(n_T = n_T, k = k, 
                                                   beta_m = filtered$beta_m, 
                                                   S_m = filtered$S_m, 
                                                   S_p = filtered$S_p, 
                                                   beta_p = filtered$beta_p)
  print("Success Smooth!")
  # Update transformed h_series 
  y_star <- log((y_series - matrix(rowSums((X_series * beta_series)), 
                                   nrow = n_T))^2 + 10^(-4))
  result <- ff_a_taylor(y_star = y_star, s = s, phi = phi, mu_h = mu_h, 
                        sigma = sigma_vartheta, h0_0 = h0_0, S_0_h0 = S_0_h0, 
                        n_T = n_T, k = k_10)
  a_mcmc[, , i] <- a_smooth <- bs_a_taylor(phi = phi, 
                                           a_measure = result$a_measure, 
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
  #Update hyper parameters
  sigma_vartheta <- sigma_vartheta_mcmc[i, 1] <- sample_sigma_vartheta(a = a,
                                                                      n_T = n_T,
                                                                      b = b,
                                                                      phi = phi,
                                                                      mu_h = mu_h,
                                                      logh_series = a_smooth)
  print("Success sigma_vartheta!")
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
  #h0_0 <- sample_logh0(mu_h, phi, logh_1, sigma_vartheta)
  beta_0 <-  sample_beta_0(sigma_eta = sigma_eta, A_0 = A_0, 
                            beta_1 = t(beta_series[1, ,drop = F]), 
                            beta_0_prior = beta_0_prior)
}



parameters_gibbs <- data.frame(mu_h_mcmc, sigma_vartheta_mcmc, phi_mcmc, 
                               c(1:n_gibbs))
names(parameters_gibbs) <- c("mu_h", "sigma_vartheta", "phi", "round")

p3 <- ggplot(parameters_gibbs) + geom_line(aes(x = round, y = mu_h)) + 
  labs(x = "round", y = "mu_h")
p4 <- ggplot(parameters_gibbs) + geom_line(aes(x = round, y = sigma_vartheta)) + 
  labs(x = "round", y = "sigma_vartheta") 
p5 <- ggplot(parameters_gibbs) + geom_line(aes(x = round, y = phi)) + 
  labs(x = "round", y = "phi") 

p3
p4
p5 
plot(1:n_T, logh_series, "l")
var(logh_series)
sigma_vartheta
mean(parameters_gibbs[, "sigma_vartheta"])
mean(parameters_gibbs[(0.5 * n_gibbs +1):n_gibbs, "mu_h"])


# Now, check how a_mcmc performs 

loghseries_estimate <- apply(a_mcmc, c(1,2), mean)
plot(loghseries_estimate, log(h_series_true))

df <- as.data.frame(cbind(loghseries_estimate, log(h_series_true), c(1:n_T)))
names(df) <- c("inferred", "true", "time")

# Actually, does not look so bad! 

ggplot(df) + 
  geom_line(aes(x = time, y = true, color = "true")) + 
  geom_line(aes(x = time, y = inferred , color = "filtered")) +  
  xlab("Day") +
  ylab("h") +
  scale_color_manual(name = "value", breaks = c("filtered", "true"),
                     values = c("filtered" = "blue", "true" = "red"))


# Same consideration for the filtered betas. 

beta_result <- apply(beta_mcmc, c(1,2), mean) 

p1 <- ggplot() + 
  geom_line(aes(x = 1:n_T, y = beta_result[, 1]), color = "blue") +
  geom_line(aes(x = 1:n_T, y = beta_series_true[, 1]), color = "red") + 
  xlab("Day") +
  ylab("beta_1")

p2 <- ggplot() + 
  geom_line(aes(x = 1:n_T, y = beta_result[, 2]), color = "blue") +
  geom_line(aes(x = 1:n_T, y = beta_series_true[, 2]), color = "red") +
  xlab("Day") +
  ylab("beta_2")

p3 <-  ggplot() + 
  geom_line(aes(x = 1:n_T, y = beta_result[, 3], color = "filtered")) +
  geom_line(aes(x = 1:n_T, y = beta_series_true[, 3], color = "true")) + 
  xlab("Day") +
  ylab("beta_3") +
  scale_color_manual(name = "value", breaks = c("filtered", "true"),
                     values = c("filtered" = "blue", "true" = "red"))


grid.arrange(p1, p2, p3, nrow = 3, top = "Check")

y_star_real <- log((y_series - matrix(rowSums((X_series * beta_series_true)), 
                                      nrow = n_T))^2 + 10^(-4))

df <- as.data.frame(cbind(y_star, y_star_real, c(1:n_T)))
names(df) <- c("inferred", "true", "time")

# Actually, does not look so bad! 

ggplot(df) + 
  geom_line(aes(x = time, y = true, color = "true")) + 
  geom_line(aes(x = time, y = inferred , color = "filtered")) +  
  xlab("Day") +
  ylab("y_star") +
  scale_color_manual(name = "value", breaks = c("filtered", "true"),
                     values = c("filtered" = "blue", "true" = "red"))

