# Use Inference for simulated data to detect possible code issues. 


setwd("gmvp/code/")

library(stats)      # Gamma, Wishart
library(MASS)       # MVNORM
library(ggplot2)
library(grid)
library(gridExtra)
# Loading required packages for doParallel library
library(doParallel) # Parallel Computing. Requires foreach, iterators, parallel

rm(list = ls())
source("helpers_kmu.R")
source("helpers_general.R")
source("simulate_kmu.R")
set.seed(345)

# Simulate data. ====
k <- 3 
h_0 <- 1 
Q_true <- Q <- solve(matrix(c(1, 0.4, 0.5, 
                        0.4, 1.5, 0, 
                        0.5, 0, 2), nrow = k, ncol = k, byrow = TRUE))
nu <- 200
lambda <- 0.99
n_T <- n_max <- 200
data <- generate_data(n = n_T, k = k, 
                      lambda = lambda, h_0 = h_0, 
                      Q = Q, nu = nu)

y_series <- data[[1]]
X_series <- data[[2]]
b <- b_10 <- data[[3]]
beta_true <- data[[4]]
h_true <- data[[5]]

remove(data)

# Forward Filtering. ====
N_10 <-  Q
S_10 <- h_0 


# Initialize 
S0_inv <-  k * solve(t(X_series) %*% X_series)
C_inv <-  diag(k)
d <- 1 
e <- 1
alpha <- 1/200

n_gibbs <- 10000
density_lambda <- h0_mcmc <- lambda_mcmc <- nu_mcmc <- matrix(0, nrow = n_gibbs)
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
  # nu_log <- logp_nu(nu = nu, h_T = h_T, h_0 = h_0, k = k, lambda = lambda, 
  #                   n_T = n_T, alpha = alpha)
  # nu_mcmc[i, 1] <- nu <- nu_sample(nu_old = nu, logdensity_old = nu_log, 
  #                                  sd = 10, h_T = h_T, 
  #                                  h_0 = h_0, k = k, lambda = lambda, n_T = n_T,
  #                                  alpha = alpha)
  # print("success nu")
  # lambda
  density_lambda[i, 1] <- logdensity_lambda <- logp_lambda(lambda, nu = nu, 
                                                           h = result_bsample[, 1, drop = F], 
                                                           h_0 = h_0, k = k, n_T = n_T)
  lambda <-  0.99 # lambda_mcmc[i, 1] <- (nu/(nu+1))  # sample_lambda(lambda, 
  #                                   logdensity_old = logdensity_lambda,
  #                                    sd = 0.00000001, nu = nu, 
  #                                    h = result_bsample[, 1, drop = F], 
  #                                    h_0 = h_0, k = k, n_T = n_T)
  print("success lambda")
  # beta_0 
  b_10 <- beta_0_sample(C_inv = C_inv, h_1 = h_1, Q = Q, beta_1 = beta_1, b = b)
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

## Model Parameters ====
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

## States  ====

### TV coefficients 

beta_result <- apply(beta_mcmc[, ,(n_gibbs * 0.5):(n_gibbs), drop = F], c(1, 2), mean) 

p4 <- ggplot() + 
  geom_line(aes(x = 1:n_T, y = beta_result[, 1]), color = "blue") +
  geom_line(aes(x = 1:n_T, y = unlist(beta_true[1, ])), color = "red") + 
  xlab("Day") +
  ylab("beta_1") 

p5 <- ggplot() + 
  geom_line(aes(x = 1:n_T, y = beta_result[, 2]), color = "blue") +
  geom_line(aes(x = 1:n_T, y = unlist(beta_true[2, ])), color = "red") + 
  xlab("Day") +
  ylab("beta_2") 

p6 <- ggplot() + 
  geom_line(aes(x = 1:n_T, y = beta_result[, 3], 
                color = "filtered")) +
  geom_line(aes(x = 1:n_T, y = unlist(beta_true[3, ]), 
                color = "true")) + 
  xlab("Day") +
  ylab("beta_3") +
  scale_color_manual(name = "value", breaks = c("filtered", "true"),
                     values = c("filtered" = "blue", "true" = "red")) 

grid.arrange(p4, p5, p6, nrow = 3, top = "Backward Sampling", 
             bottom = textGrob(paste("Simulated data: nu = ",
                                     nu, "lambda = ", lambda, ", h_0 =", 
                                     h_0),
                               gp = gpar(fontface = 3, fontsize = 9),
                               hjust = 1,
                               x = 1 ))

### Volatility 
h <- as.data.frame(apply(h_mcmc[, , (0.5 * n_gibbs + 1):n_gibbs, drop = F], c(1,2), mean))
dim(h_true)
ggplot() + 
  geom_line(aes(x = 1:n_T, y = unlist(h_true[,1]), color = "true")) + 
  geom_line(aes(x = 1:n_T, y = unlist(h), color = "filtered")) +  
  xlab("Day") +
  ylab("h") +
  scale_color_manual(name = "value", breaks = c("filtered", "true"),
                     values = c("filtered" = "blue", "true" = "red"))


# Matrices 

Q_estim <- apply(Q_mcmc[, , (0.5 * n_gibbs + 1):n_gibbs, drop = F], c(1,2), mean)
Q_cond <- apply(Q_mcmc[, , (0.5 * n_gibbs + 1):n_gibbs, drop = F], c(3), det)

plot(1:(0.5 *  n_gibbs), Q_cond, "l")
Q_estim
solve(Q_estim) - solve(Q_true)