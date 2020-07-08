# Title: Simulation
# Date: 20.01.19 
# Author: Hannah Busshoff
# Email: hbusshof@smail.uni-koeln.de

# Short description: This program tests the forward-filter and backward sampling  
# algorithm for simulated data. 

setwd("~/Desktop/bayesian-gmvp/Computation/") 

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
source("simulation_kmu.R")
set.seed(345)

# Simulate data. ====
k <- 3 
h_0 <- 1 
Q <- solve(100*matrix(c(1, 0.4, 0.5, 
                    0.4, 1.5, 0, 
                    0.5, 0, 2), nrow = k, ncol = k, byrow = TRUE))
nu <- 200
lambda <- 0.98
n_T <- 200
data <- generate_data(n = n_T, k = k, 
                          lambda = lambda, h_0 = h_0, 
                          Q = Q, nu = nu)

y <- data[[1]]
X <- data[[2]]
b_10 <- data[[3]]
beta_true <- data[[4]]
h_true <- data[[5]]

remove(data)

# Forward Filtering. ====
N_10 <-  Q 
S_10 <- h_0 

result_filter <- filter_forward(X = X, y = y, nu = nu, Q = Q, lambda = lambda,
                                N_10 = N_10, S_10 = S_10, b_10 = b_10, 
                                n_T = n_T, k = k) 

beta_filter <- result_filter$beta_m 

# Backward Sampling. ==== 

D <- 1000 
b_m <- beta_filter

mcoptions <- list(preschedule=FALSE, set.seed= TRUE) 

result_bsample <- backwardsample_it(D = D, S_m = result_filter$S_m, b_m = b_m, 
                                  N_m = result_filter$N_m, nu = nu, Q = Q, 
                                  lambda = lambda, n_T = n_T, k = k)

# Plotting Results. ==== 

p1 <- ggplot() + 
        geom_line(aes(x = 1:n_T, y = beta_filter[1, ]), color = "blue") +
        geom_line(aes(x = 1:n_T, y = unlist(beta_true[1, ])), color = "red") + 
        xlab("Day") +
        ylab("beta_1")

p2 <- ggplot() + 
        geom_line(aes(x = 1:n_T, y = beta_filter[2, ]), color = "blue") +
        geom_line(aes(x = 1:n_T, y = unlist(beta_true[2, ])), color = "red") +
        xlab("Day") +
        ylab("beta_2")

p3 <-  ggplot() + 
       geom_line(aes(x = 1:n_T, y = beta_filter[3, ], color = "filtered")) +
       geom_line(aes(x = 1:n_T, y = unlist(beta_true[3, ]), color = "true")) + 
       xlab("Day") +
       ylab("beta_3") +
       scale_color_manual(name = "value", breaks = c("filtered", "true"),
                           values = c("filtered" = "blue", "true" = "red"))
        
p4 <- ggplot() + 
        geom_line(aes(x = 1:n_T, y = result_bsample[, 2]), color = "blue") +
        geom_line(aes(x = 1:n_T, y = unlist(beta_true[1, ])), color = "red") + 
        xlab("Day") +
        ylab("beta_1") 

p5 <- ggplot() + 
        geom_line(aes(x = 1:n_T, y = result_bsample[, 3]), color = "blue") +
        geom_line(aes(x = 1:n_T, y = unlist(beta_true[2, ])), color = "red") + 
        xlab("Day") +
        ylab("beta_2") 

p6 <- ggplot() + 
        geom_line(aes(x = 1:n_T, y = result_bsample[, 4], 
                      color = "filtered")) +
        geom_line(aes(x = 1:n_T, y = unlist(beta_true[3, ]), 
                      color = "true")) + 
        xlab("Day") +
        ylab("beta_3") +
        scale_color_manual(name = "value", breaks = c("filtered", "true"),
                           values = c("filtered" = "blue", "true" = "red")) 

grid.arrange(p1, p2, p3, nrow = 3, top = "Forward Filtering", 
             bottom = textGrob(paste("Simulated data: nu = ",
                                     nu, "lambda = ", lambda, ", h_0 =", 
                                     h_0),
                               gp = gpar(fontface = 3, fontsize = 9),
                               hjust = 1,
                               x = 1 ))

grid.arrange(p4, p5, p6, nrow = 3, top = "Backward Sampling", 
             bottom = textGrob(paste("Simulated data: nu = ",
                                     nu, "lambda = ", lambda, ", h_0 =", 
                                     h_0),
                     gp = gpar(fontface = 3, fontsize = 9),
                     hjust = 1,
                     x = 1 ))

## Plot for h ==== 

ggplot() + 
        geom_line(aes(x = 1:n_T, y = unlist(h_true[,1]), color = "true")) + 
        geom_line(aes(x = 1:n_T, y = result_bsample[,1], color = "filtered")) +  
        xlab("Day") +
        ylab("h") +
        scale_color_manual(name = "value", breaks = c("filtered", "true"),
                           values = c("filtered" = "blue", "true" = "red"))

# This gives some indication that there may be a problem with the formula or 
# with the code on h!