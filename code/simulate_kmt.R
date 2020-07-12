# Simulation of data in accordance with Taylor ====

generate_data <- function(n_T, sigma_vartheta, logh_0, gamma, 
                          phi, sigma_eta, X_series, 
                          beta_0){
  h_series <- generate_h(n_T, sigma_vartheta, logh_0, gamma, phi)
  beta_series <- generate_beta(n_T, sigma_eta, beta_0)
  epsilon <- matrix(rnorm(n_T, mean = 0, sd = 1), nrow = n_T)
  y_series <- apply(X_series * beta_series, 1, sum) + h_series * epsilon
  return(list(y_series, h_series, beta_series))
}

generate_h <- function(n_T, sigma_vartheta, logh_0, gamma, phi){
  logh_series <- matrix(NA, ncol = 1, nrow = (n_T + 1))
  vartheta <- matrix(rnorm(n_T, mean = 0, sd = sqrt(sigma_vartheta)), 
                     nrow = n_T)
  logh_series[1, 1] <- logh_0
  for (i in 1:n_T){ 
    logh_series[i + 1, 1] <- gamma + phi * logh_series[i, 1] + vartheta[i, 1]
  }
  sqrt(exp(logh_series))[2:nrow(logh_series), ,drop = FALSE] 
}

generate_beta <- function(n_T, sigma_eta, beta_0){ 
  eta_series <- mvrnorm(n_T, mu = rep(0, k), Sigma = sigma_eta)
  eta_series[1, ] <- eta_series[1, ] + beta_0
  apply(eta_series, 2, cumsum)
}