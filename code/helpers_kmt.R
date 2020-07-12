# FFBS  beta ====

measure_taylor <- function(X, y, h, beta_p, S_p){
  R = drop((X %*% S_p %*% t(X) +  h))
  i = drop((y - X %*% beta_p))
  beta_m = beta_p + S_p %*% t(X) * (1/R) * i
  S_m = S_p - (1/R) * S_p  %*% t(X) %*% X %*% S_p
  return(list(beta_m, S_m))
}

predict_taylor <- function(sigma_eta, beta_m, S_m){
  beta_p = beta_m 
  S_p = S_m + sigma_eta
  return(list(beta_p, S_p))
}

filter_taylor <- function(k, n_T, beta_0, S_0_beta, X_series, y_series, 
                          h_series, sigma_eta){
  beta_m <- beta_p <-  matrix(NA, ncol = k, nrow = n_T)
  S_m <- S_p <- array(0, dim = c(k, k, n_T))
  predict <- predict_taylor(sigma_eta = sigma_eta, beta_0, S_0_beta)
  for (i in 1:n_T){
    measure = measure_taylor(X_series[i, , drop = F], y_series[i, 1], 
                             h_series[i, 1], predict[[1]], predict[[2]])
    beta_m[i, ] <- measure[[1]]
    S_m[, , i] <- measure[[2]]
    predict = predict_taylor(sigma_eta = sigma_eta, measure[[1]], measure[[2]])
    beta_p[i, ] <- predict[[1]]
    S_p[, , i] <- predict[[2]]
  }
  result <- list(beta_m, beta_p, S_m, S_p)
  names(result) <- c("beta_m", "beta_p", "S_m", "S_p")
  return(result)
}


smooth_taylor <- function(n_T, k, beta_m, S_m, 
                           S_p, beta_p){
  beta_smooth <- matrix(NA, nrow = n_T, ncol = k)
  S_smooth <- array(0, dim = c(k, k, n_T))
  beta_smooth[1, ] <- mvrnorm(1, mu = beta_m[n_T, ], 
                              Sigma = S_m[, ,n_T])
  for (i in 1:(n_T-1)){
    S_pt <- S_p[, ,n_T - i]
    S_mt <- S_m[, ,n_T - i]
    mu <- beta_m[n_T - i, ] + S_mt %*% solve(S_pt) %*% (beta_smooth[i, ] -
                                                              beta_p[n_T - i, ]) 
    aux <- S_mt - S_mt %*% solve(S_pt) %*% S_mt
    S_smooth[, , i + 1] <-  aux
    beta_smooth[i + 1, ] <- mvrnorm(1, mu = mu, Sigma = aux)
  }
  apply(beta_smooth, 2, rev)
}

# FFBS a := log(h) ====

measure_a_taylor <- function(a_p, S_p, v, m, y_star){
  i <- y_star - a_p - m 
  R <- S_p + v 
  a_m <- a_p + S_p * (1/R) * i 
  S_m <- S_p - S_p * (1/R) * S_p
  return(list(a_m, S_m))
}

predict_a_taylor <- function(gamma, phi, sigma, a_m, S_m){
  a_p <- gamma + phi * a_m 
  S_p <- phi^2*S_m + sigma
  return(list(a_p, S_p))
}

smoothing_a_taylor <- function(phi, a_m, S_m, S_p, a_plus, a_pred){
  a_s <- a_m + phi * S_m * (1/S_p) * (a_plus - a_p)
  S_s <- S_m - phi^2 * S_m * (1 / S_p) * S_m
  return(list(a_s, S_s))
}

ff_a_taylor <- function(y_star, s, phi, mu_h, sigma, h0_0, S_0_h0, n_T, k){
  gamma <- mu_h * (1 - phi)
  a_measure <- S_measure <-  matrix(0, nrow = n_T, ncol = 1)
  a_predict <- S_predict <- matrix(0, nrow = (n_T + 1), ncol = 1)
  # Initialize 
  initial <- predict_a_taylor(gamma = gamma, phi = phi, sigma = sigma, 
                              a_m = h0_0, S_m = S_0_h0)
  a_predict[1, 1] <- initial[[1]]
  S_predict[1, 1] <- initial[[2]]
  # Begin filtering and prediction round 
  for (i in 1:n_T){
    # Filter
    measure <- measure_a_taylor(a_p = a_predict[i, 1], 
                               S_p = S_predict[i, 1], 
                               m = k[s[i, 1], 2], 
                               v = k[s[i, 1], 3], 
                               y_star = y_star[i, 1]) 
    
    a_measure[i, 1] <- measure[[1]]
    S_measure[i, 1] <- measure[[2]]
    # Predict
    predict <- predict_a_taylor(gamma = gamma, phi = phi, sigma = sigma, 
                                a_m = measure[[1]], S_m = measure[[2]])
    a_predict[i + 1, 1] <- predict[[1]]
    S_predict[i + 1, 1] <- predict[[2]]
    
  }
  result <- list(a_measure, a_predict[1:n_T, 1, drop = F], S_measure,
                 S_predict[1:n_T, 1, drop = F])
  names(result) <- c("a_measure", "a_predict", "S_measure", "S_predict")
  return(result)
}


bs_a_taylor <- function(phi, a_measure, a_predict, S_measure, S_predict, n_T){
  a_T <- rnorm(1, mean = a_measure[n_T, 1], 
               sd = sqrt(S_measure[n_T, 1]))
  a_smooth <-  matrix(0, nrow = n_T, ncol = 1)
  a_smooth[n_T, 1] <- a_T
  for (i in 1:(n_T - 1)){
    var <- S_measure[n_T - i, 1] - phi^2 * S_measure[n_T - i, 1]^2 * 
      (1 / S_predict[(n_T + 1 - i), 1])  
    mu <- a_measure[n_T - i, 1] + phi * S_measure[n_T - i, 1] * 
      (1/S_predict[(n_T + 1 - i), 1]) * 
      (a_smooth[(n_T + 1 - i), 1] -  a_predict[(n_T + 1 - i), 1])
    a_smooth[(n_T  - i), 1] <- rnorm(1,  mean = mu, sd = sqrt(var))
  }
  return(a_smooth)
}

# Indicator  ==== 

density_s <- function(y_star, a, m, v){
  apply(((y_star - a - m)), 1, function(x){dnorm(x, mean = 0, sd = sqrt(v))})
}

update_indicator <- function(y_star, a, k, no_mixture, weights, n_T){
  probs <- matrix(0, nrow = n_T, ncol = no_mixture)
  for (s in 1:no_mixture){
    probs[, s] <-  density_s(y_star, a,
                             m = k[s, 2], v = k[s, 3])
  }
  unnormal <- probs * weights 
  factor <- rowSums(unnormal)
  factor <- matrix(rep(factor, no_mixture), nrow = n_T, byrow = F) 
  normal <- unnormal / factor
  apply(normal, 1, function(x){sample(1:no_mixture, 
                                      1, prob = c(x))})
}

# Hierarchical Bayes ====

logprior_phi <- function(phi, alpha_phi, beta_phi){
  log(((phi + 1) / 2))*(alpha_phi - 1) + log(((1 - phi) / 2))*(beta_phi - 1)
}

moments_phi <- function(logh_series, mu_h, sigma_vartheta, n_T){
  nominator <- sum((logh_series[2:n_T, ] - mu_h) * 
                     (logh_series[1:(n_T - 1)] - mu_h))
  denominator <- sum((logh_series[1:(n_T - 1)] - mu_h)^2) 
  mu_phi <- nominator / denominator
  sigma_phi <- sigma_vartheta / denominator
  return(list(mu_phi, sigma_phi))
}

helper_phi <- function(phi, sigma_vartheta, mu_h, logh_series){
  log(1 - phi^2)*(1/2) - (logh_series[1, 1] - mu_h)^2 * (1 - phi^2) / (2 * sigma_vartheta)
}

accept_phi <- function(phi_prop, phi, a, b, mu_h, sigma_vartheta, 
                       logh_series){ 
  term_1 <- logprior_phi(phi_prop, a, b) + helper_phi(phi_prop, sigma_vartheta, mu_h, logh_series)
  term_2 <- logprior_phi(phi, a, b) + helper_phi(phi, sigma_vartheta, mu_h, logh_series)
  if (!is.numeric(term_1)){
    warning("proposal kernel not numeric")
  }
  if (!is.numeric(term_2)){
    warning("kernel with init. value is not numeric")
  }
  term_1 - term_2 
}

sample_phi <- function(phi_old, logh_series, mu_h, sigma_vartheta, alpha_phi, 
                       beta_phi, n_T){
  moments <- moments_phi(logh_series = logh_series, mu_h = mu_h, 
                         sigma_vartheta = sigma_vartheta, n_T = n_T)
  phi_prop <- rnorm(1, mean = moments[[1]], sd = moments[[2]]) 
  if (abs(phi_prop) >= 1){
    return(phi_old)
  } else {
    a_prob <- accept_phi(phi_prop = phi_prop, phi = phi_old, a = alpha_phi, 
                         b = beta_phi, mu_h = mu_h, 
                         sigma_vartheta = sigma_vartheta, 
                         logh_series = logh_series)
    if (!is.numeric(a_prob)){
      return(list(phi_prop, phi_old))
      break
    } else if (a_prob < 0){ 
      if (runif(1, min = 0, max = 1) <= exp(a_prob)){
        phi_old <- phi_prop
      } 
    } else if (a_prob >= 0) {
      phi_old <- phi_prop
    }
    return(phi_old)
  }
} 

sample_sigma_vartheta <- function(a, n_T, b, phi, mu_h, 
                              logh_series){
  a_update <- (2*a + n_T)/2 
  b_1 <- 2*b + (1 - phi^2) * (logh_series[1, 1] - mu_h)^2   
  b_2 <- sum((logh_series[2 : n_T, ] - phi * logh_series[1:n_T - 1, ] - 
                (1 - phi) * mu_h)^2)
  b_update = (b_1 + b_2)/2 
  1 / rgamma(1, shape = a_update, rate = b_update)
}

sample_mu_h <- function(mu_h_0, sd_h_0, logh_series, phi, sigma_vartheta, n_T){
  mu_h <- ((1 + phi)*logh_series[1,1] + 
            sum(logh_series[2 : nrow(logh_series), ] - 
              phi * logh_series[1:(nrow(logh_series) - 1), ])) / 
          ((1 + phi) + (n_T - 1) * (1 - phi) + (sigma_vartheta) / 
            ((1 - phi) * sd_h_0))
  s_h <- (((1 - phi^2) + (n_T - 1) * (1 - phi)^2) / 
            (sigma_vartheta) + 1 / (sd_h_0))^(-1)
  rnorm(1, mean = mu_h, sd = sqrt(s_h))
}

sample_beta_0 <- function(sigma_eta, A_0, beta_1, beta_0_prior){
  sigma_eta_inv <- solve(sigma_eta)
  sigma_udpate <-  solve(sigma_eta_inv + A_0)
  mu <- sigma_udpate %*% (sigma_eta_inv %*% beta_1 + A_0 %*% beta_0_prior)
  mvrnorm(1, mu = mu, Sigma = sigma_udpate)
} 

sample_logh0 <- function(mu_h, phi, logh_1, sigma_vartheta){
  mean <- mu_h + phi * (logh_1 - mu_h)
  rnorm(1, mean = mean, sd = sqrt(sigma_vartheta))
}

sample_sigma_eta <- function(n_T, k, beta, beta_0, S0_inv){
  beta_lag <- rbind(t(beta_0), beta[1:(n_T - 1), ])
  beta_tilde <- beta - beta_lag
  G <- matrix(0, ncol = k, nrow = k)
  for (i in 1:n_T){
    G <- G + t(beta_tilde[i, , drop = F]) %*% beta_tilde[i, , drop = F]
  }
  solve(drop(rWishart(1, n_T + k + 1, solve(G + S0_inv))))
}

# Normal mixtures ====

k_10 <- matrix(c(
  0.00609,
  0.04775,
  0.13057,
  0.20674,
  0.22715,
  0.18842,
  0.12047,
  0.05591,
  0.01575,
  0.00115,
  1.92677,
  1.34744,
  0.73504,
  0.02266,
  -0.85173,
  -1.97278,
  -3.46788,
  -5.55246,
  -8.68384,
  -14.65000,
  0.11265,
  0.17788,
  0.26768,
  0.40611,
  0.62699,
  0.98583,
  1.57469,
  2.54498,
  4.16591,
  7.33342
), byrow = F, ncol = 3)  


k_7 <- matrix(c(0.04395,
                0.24566,
                0.34001,
                0.25750,
                0.10556,
                0.00002,
                0.00730,
                1.50746,
                0.52478,
                -0.65098,
                -2.35859,
                -5.24321,
                -9.83726,
                -11.40039,
                0.16735,
                0.34023,
                0.64009,
                1.26261,
                2.61369,
                5.17950,
                5.79596), byrow = F, ncol = 3)

colnames(k_7) <- colnames(k_10) <- c("p_j", "m_j", "v_j")