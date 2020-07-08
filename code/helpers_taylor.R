# Kalman-filter  ====

measure_taylor <- function(X, y, h, beta_p, S_p){
  R = drop((X %*% S_p %*% t(X) +  h))
  error = drop((y - X %*% beta_p))
  beta_m = beta_p + S_p %*% t(X) * (1/R) * error
  S_m = S_p - (1/R) * S_p  %*% t(X) %*% X %*% S_p
  return(list(beta_m, S_m))
}

predict_taylor <- function(sigma_eta, beta_m, S_m){
  beta_p = beta_m 
  S_p = S_m + sigma_eta
  return(list(beta_p, S_p))
}

filter_taylor <- function(k, horizon, beta_0, X_series, y_series, h_series, 
                          sigma_eta){
  beta_m <- beta_p <-  matrix(NA, ncol = k, nrow = horizon)
  S_m <- S_p <- matrix(NA, ncol = (k * (k+1) * 0.5), nrow = horizon)
  predict <- list(matrix(beta_0, nrow = k), sigma_eta)
  for (i in 1:horizon){
    measure = measure_taylor(X_series[i, , drop = FALSE], y_series[i, 1], 
                             h_series[i, 1], predict[[1]], predict[[2]])
    beta_m[i, ] <- measure[[1]]
    S_m[i, ] <- measure[[2]][lower.tri(measure[[2]], diag = TRUE)]
    predict = predict_taylor(sigma_eta = sigma_eta, measure[[1]], measure[[2]])
    beta_p[i, ] <- predict[[1]]
    S_p[i, ] <- predict[[2]][lower.tri(predict[[2]], diag = TRUE)]
  }
  result <- list(beta_m, beta_p, S_m, S_p)
  names(result) <- c("beta_m", "beta_p", "S_m", "S_p")
  return(result)
}

# Smoothing ==== 

smooth_taylor <- function(horizon, k, beta_m, S_m, 
                           S_p, beta_p){
  beta_smooth <- matrix(NA, nrow = horizon, ncol = k)
  S_smooth <- matrix(NA, nrow = horizon, ncol = (k * (k + 1) * 0.5))
  beta_smooth[1, ] <- mvrnorm(1, mu = beta_m[horizon, ], 
                              Sigma = convert_2_cov(S_m[horizon, ], k = k))
  for (i in 1:(horizon-1)){
    S_pt <- convert_2_cov(S_p[horizon - i, ], k = k)
    S_mt <- convert_2_cov(S_m[horizon - i, ], k = k)
    mu <- beta_m[horizon - i, ] + S_mt %*% solve(S_pt) %*% (beta_smooth[i, ] -
                                                              beta_p[horizon - i, ]) 
    aux <- S_mt + S_mt %*% solve(S_pt) %*% S_mt
    S_smooth[i + 1, ] <-  aux[lower.tri(aux, diag = TRUE)]
    beta_smooth[i + 1, ] <- mvrnorm(1, mu = mu, Sigma = aux)
  }
  apply(beta_smooth, 2, rev)
}

# Hierarchical Bayes for phi ====

prior_phi <- function(phi, alpha_phi, beta_phi){
  ((phi + 1) / 2)^(alpha_phi - 1) * ((1 - phi) / 2)^(beta_phi - 1)
}

moments_phi <- function(logh_series, mu_h, sigma_vartheta){
  nominator <- sum((logh_series[2:nrow(logh_series), ] - mu_h) * 
                     (logh_series[1:(nrow(logh_series) - 1)] - mu_h))
  denominator <- sum((logh_series[1:(nrow(logh_series) - 1)] - mu_h)^2) 
  mu_phi <- nominator / denominator
  sigma_phi <- sigma_vartheta / denominator
  return(list(mu_phi, sigma_phi))
}

helper_phi <- function(phi, sigma_vartheta, mu_h, logh_series){
  (1 - phi) * exp(- (logh_series[1, 1] - mu_h)^2 * (1 - phi^2) / (2 * sigma_vartheta))
}

accept_phi <- function(phi_prop, phi, a, b, mu_h, sigma_vartheta, 
                       logh_series){ 
  nominator <- prior_phi(phi_prop, a, b) * helper_phi(phi_prop, sigma_vartheta, 
                                                      mu_h, logh_series)
  denominator <- prior_phi(phi, a, b) * helper_phi(phi, sigma_vartheta, 
                                                   mu_h, logh_series)
  return <- nominator / denominator
}

phi_hb <- function(mcmc_rounds, logh_series, mu_h, 
                     sigma_vartheta, alpha_phi, beta_phi){
  phi_mcmc <- matrix(rep(0, (mcmc_rounds + 1)), ncol = 1)
  phi_0 <- runif(1, min = -1, max =  1)
  while(abs(phi_0) ==  1 ){
    phi_0 <- runif(1, min = -1, max = 1)
  }
  phi_mcmc[1, 1] <- phi_0
  for (i in 1:mcmc_rounds){
    moments <- moments_phi(logh_series = logh_series, mu_h = mu_h, 
                           sigma_vartheta = sigma_vartheta)
    phi_prop <- rnorm(1, mean = moments[[1]], sd = moments[[2]]) 
    a_prob <- accept_phi(phi_prop = phi_prop, phi = phi_0, a = alpha_phi, 
                         b = beta_phi, mu_h = mu_h, 
                         sigma_vartheta = sigma_vartheta, 
                         logh_series = logh_series)
    if (a_prob < 1){ 
      if (a_prob <= runif(1, min = 0, max = 1)){
        phi_0 <- phi_prop
      } 
    } else {
      phi_0 <- phi_prop
    }
    phi_mcmc[i + 1, 1] <- phi_0
  } 
  return(mean(phi_mcmc))
}

# Hierarchical Bayes for sigma_vartheta ====

sigma_vartheta_hb <- function(a, horizon, b, phi, mu_h, 
                              logh_series, n_mc = T){
  a_update <- a + horizon 
  b_1 <- b + (1 - phi^2) * (logh_series[1, 1] - (mu_h))^2   
  b_2 <- sum((logh_series[2 : nrow(logh_series), ] - 
                phi * logh_series[1:(nrow(logh_series) - 1), ] - 
                (1 - phi) * mu_h)^2)
  b_update = b_1 + b_2 
  if (n_mc == T){
    estimate <- (b_update / 2) / ((a_update / 2) - 1)
  } else { 
    estimate <- mean(1/ rgamma(1, a_update / 2, b_update/ 2))
  }
  return(estimate)
}

# Hierarchical Bayes for mu_h ====

mu_h_hb <- function(mu_phi, sd_phi, logh_series, phi, sigma_vartheta, horizon, 
                   n_mc = T){
mu_h <- ((1 + phi)*logh_series[1,1] + sum(logh_series[2 : nrow(logh_series), ] - 
                                           phi * logh_series[1:(nrow(logh_series) - 1), ])) / 
  ((1 + phi) + (horizon - 1) * (1 - phi) + (sigma_vartheta) / 
     ((1 - phi) * sd_phi))
if (n_mc == T){
 estimate <- mu_h  
} else {
  s_h <- (((1 - phi^2) + (horizon - 1) * (1 - phi)^2) / (sigma_vartheta) + 1 / (sd_phi))
  estimate <- mean(rnorm(n_mc, mean = mu_h,sd = s_h))
  }
}

# Hierarchical Bayes for Sigma_eta ==== 

compute_beta_shock <- function(beta, beta_0, n_T){
  beta_lag <- matrix(cbind(beta_0, beta[, 1:(n_T - 1)]), ncol = n_T, byrow = F)
  as.matrix(c(rowSums((beta - beta_lag)^2)), ncol = 1)
}

sample_sigma_eta <- function(beta_shock, nu_lambda, n_T){
  shape <- (nu_lambda + n_T) * 0.5 
  scale <-  2 / (beta_shock + nu_lambda) 
  precision <- as.matrix(apply(scale, 1, function(x){rgamma(1, shape = shape, scale = x)}), 
            ncol = 1)
  solve(diag(c(precision)))
}

# Hierarchical Bayes for omega ====

sample_omega <- function(beta, beta_0, n_T, lambdas, b, a, n_N, sample){
  beta_lag <- matrix(cbind(beta_0, beta[, 1:(n_T - 1)]), ncol = n_T, byrow = F)
  term_1 <- sum((rowSums((beta - beta_lag)^2) * lambdas))  
  nu <- n_T * n_N
  mu <- (b + nu)/(term_1 + b/a)
  mean(rgamma(sample, shape = nu/2, scale = 2*mu/nu))
}

# Hierarchical Bayes for beta_0 ====

sample_beta_0 <- function(sigma_eta, A_0, beta_1, beta_0_prior, sample, 
                          collapse = T){
  sigma_eta_inv <- solve(sigma_eta)
  sigma_udpate <-  solve(sigma_eta_inv + A_0)
  mu <- sigma_udpate %*% (sigma_eta_inv %*% beta_1 + A_0 %*% beta_0_prior )
  if (collapse == T){
    return(colMeans(mvrnorm(sample, mu = mu, Sigma = sigma_udpate)))
  } else {
    return(mvrnorm(sample, mu = mu, Sigma = sigma_udpate))
  }
} 

# Hierarchical Bayes for h_0 ==== 

sample_logh0 <- function(mu_h, phi, logh_1, sigma_vartheta, sample, 
                      collapse = T){
  mean <- mu_h + phi * (logh_1 - mu_h)
  result <- rnorm(sample, mean = mean, sd = sqrt(sigma_vartheta))
  if (collapse ==  T){
    return(mean(result))
  } else { 
    return(result)
  }
}

# Hierarchical Bayes for nu_lambda ====

logp_nu <- function(nu_0, nu_lambda, n_N, precision){
  term_1  <- 0.5 * n_N * nu_lambda * log(nu_lambda/2)
  term_2  <- (-n_N) * log(gamma((nu_lambda)/2)) 
  upsilon <- 1/(nu_0) + 0.5 * sum(log(1/precision) + precision)
  term_3  <- (-upsilon) * nu_lambda 
  term_1 + term_2 + term_3
}

nu_lambda_hb <- function(mcmc_rounds, nu_0, n_N, precision, nu_lambda_0, c, 
                         collapse = T){
  accept <- 0 
  nu_lambda_mcmc <- matrix(rep(0, (mcmc_rounds)), ncol = 1)
  nu_lambda_old <- nu_lambda_0
  log_old <- logp_nu(nu_0, nu_lambda_old, n_N, precision)
  for (i in 1:mcmc_rounds){
    nu_lambda_prop <- nu_lambda_old + rnorm(1, mean = 0, sd = c)
    while (nu_lambda_prop <= 0){
      nu_lambda_prop <- nu_lambda_old + rnorm(1, mean = 0, sd = c)
    }  
    log_prop <- logp_nu(nu_0, nu_lambda_prop, n_N, precision)
    if ((log_prop - log_old) < 0){
      if (runif(1, min = 0, max = 1) <= exp((log_prop - log_old))){
        nu_lambda_old <- nu_lambda_prop 
        accept <- accept + 1
      }
    } else {
      nu_lambda_old <- nu_lambda_prop
      accept <- accept + 1
    }
    nu_lambda_mcmc[i, 1] <- nu_lambda_old
  } 
  a_rate <- accept/mcmc_rounds
  if (collapse == T){
    return(list((mean(nu_lambda_mcmc[(mcmc_rounds * 0.5) : (mcmc_rounds)])), 
                a_rate))
  }
  else {
    return(list(a_rate, nu_lambda_mcmc))
  }
}

# Mixtures for approximation ====
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

# Mixture simulator Kim ==== 

# FFBS ====

filter_a_taylor <- function(a_p, S_p, v, m, y_star){
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
  return(list(a_p, S_p))
}

ff_a_taylor <- function(y_star, s, phi, gamma, sigma, initial, 
                        n_T, k){
  a_measure <- S_measure <-  matrix(0, nrow = n_T, ncol = 1)
  a_predict <- S_predict <- matrix(0, nrow = (n_T + 1), ncol = 1)
  # Initialize 
  a_predict[1, 1] <- initial[1]
  S_predict[1, 1] <- initial[2]
  # Begin filtering and prediction round 
  for (i in 1:n_T){
    # Filter
    measure <- filter_a_taylor(a_p = a_predict[i, 1], 
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
    var <- S_measure[n_T - i, 1] - phi^2 * S_measure[n_T - i, 1]^2 * (1 / S_predict[(n_T + 1 - i), 1])  
    mu <- a_measure[n_T - i, 1] + phi * S_measure[n_T - i, 1] * (1/S_predict[(n_T + 1 - i), 1]) * (a_smooth[(n_T + 1 - i), 1] - 
                                                                                                     a_predict[(n_T + 1 - i), 1])
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