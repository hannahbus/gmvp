# FFBS ==== 

measure <- function(beta_p, N_p, nu, S_p, X, y, k){
  e   = y - X %*% beta_p 
  N_m = N_p + t(X) %*% X
  N_m_inv <- solve(N_m)
  beta_m = N_m_inv %*% (N_p %*% beta_p + y*t(X))  
  S_m = (nu/(nu + 1)) * S_p +
    1/(nu+1) * e * (1 - X %*% N_m_inv %*% t(X)) * e
  return(list(beta_m, N_m, S_m))
} 

predict <- function(beta_m, N_m, nu, S_m, lambda, Q){ 
  beta_p = beta_m 
  N_p = solve((solve(Q) + solve(lambda * N_m)))
  S_p = lambda * (nu + 1) / (nu) * S_m 
  return(list(beta_p, N_p, S_p))
}

filter_forward <- function(X, y, nu, Q, lambda, N_10, S_10, b_10, n_T, k){ 
  # Initialize. 
  S_m <- matrix(0, nrow = n_T, ncol = 1)
  S_p <- matrix(0, nrow = 1, ncol = 1) 
  b_m <- matrix(0, nrow = k, ncol = n_T)
  b_p <- matrix(0, nrow = k, ncol = 1)
  N_m <- array(0, dim = c(k, k, n_T))
  N_p <- N_10 
  S_p  <- S_10
  b_p[, 1] <- b_10
  for (i in 1:n_T){ 
    # Measure.
    result_m <- measure(beta_p = b_p, N_p = N_p, S_p = S_p, 
                        X = X[i, , drop = F], y = y[i, 1], k = k, nu = nu)
    b_m[, i] <- result_m[[1]]
    N_m[, , i] <- result_m[[2]]
    S_m[i, 1] <- result_m[[3]] 
    # Predict.
    result_p <- predict(beta_m = b_m[, i], N_m = N_m[, , i], 
                        nu = nu, S_m = S_m[i, 1], 
                        lambda = lambda, Q = Q)
    b_p[, 1] <- result_p[[1]]
    N_p <- result_p[[2]]
    S_p <- result_p[[3]]
  }
  result <- list(S_m, b_m, N_m)
  names(result) <- c("S_m", "beta_m", "N_m")
  return(result)
}

backward_sample <- function(S_m, b_m, N_m, nu, Q, lambda, n_T, k){ 
  # Initialize result vectors. 
  b_mT = b_m[, n_T, drop = F]
  N_mT = N_m[, , n_T]
  beta  = matrix(nrow = n_T, ncol = k)
  h = matrix(nrow = n_T, ncol = 1)
  # Sample theta_T 
  h[1, 1] = rgamma(1, shape = (nu + 1)/2, scale = 2/(S_m[n_T, 1] * (nu + 1)))
  beta[1, ] = mvrnorm(n = 1, t(b_mT), solve(N_mT * h[1,1]))
  for (i in 1:(n_T-1)){ 
    z = rgamma(1, shape = 0.5, scale = 2 / ((nu + 1) * S_m[n_T - i, 1]))
    h[1 + i, 1] = lambda * h[i, 1] + z
    C = h[i + 1, 1] * N_m[, , n_T - i]
    D = h[i, 1] * Q 
    E = solve(C + D)
    mu = E %*% (C %*% b_m[, n_T - i] + D %*% beta[i, ])
    beta[1 + i, ] = mvrnorm(n = 1, mu, E)
  }
  return(cbind(rev(h), apply(beta, 2, rev)))
}

backwardsample_it <- function(D, S_m, b_m, N_m, nu, Q, lambda, n_T, k,  
                  mcoptions = list(preschedule = F, set.seed= T)){
  gibbs_result <- foreach(icount(D), .options.multicore=mcoptions) %dopar% { 
                      backward_sample(S_m, b_m, N_m, nu, Q, lambda, n_T, k)
                      }
  result <- matrix(0, nrow = nrow(S_m), ncol = (nrow(b_m) + 1))
  colnames(result) <- c("h",  paste("beta", 0:(nrow(b_m) - 1), sep = "_"))
  for (i in 1:(k+1)){ 
    marginal <- sapply(gibbs_result, function(x) x[, i])
    result[, i] <- apply(marginal, 1, mean)
  }
  return(result)
}

# Hierarchical Bayes ==== 

compute_G <- function(beta_0, theta, n_T, k){ 
  beta <- t(theta[, 2:(k + 1)])
  beta_lag <- cbind(beta_0, beta[, 1:(n_T - 1)])
  h <- theta[, 1]
  beta_tilde <- beta - beta_lag 
  G <- matrix(0, ncol = k, nrow = k)
  for (i in 1:n_T){
   G <- G + h[i] * beta_tilde[, i] %*% t(beta_tilde[, i])
  }
  if (!isSymmetric(G)){
    warning("A non-symmetric matrix was generated.")
  }
  return(G)
} 

sample_Q <- function(beta_0, theta, n_T, k, S0_inv){
  G_exp <- compute_G(beta_0, theta, n_T, k)
  scale_W <- 1/(n_T) * (G_exp + S0_inv)
  cond <- kappa(scale_W)
  if (cond > 10^6){
    warning(print("The conditioning number is > 10^6."))
  }
  drop(rWishart(1, df = (n_T + k), Sigma = 1/(n_T) * solve(scale_W)))
}

logp_nu <- function(nu, h_T, h_0, k, lambda, n_T, alpha){
  if (nu > 0){
    logp_nu <- (n_T * log(lambda) + log(h_T / h_0)) * 
      ((nu + k) / 2 - 1) - n_T * log(beta((nu + k) * 0.5, 0.5)) - alpha * nu
    return(logp_nu)
  } else if (nu <= 0){ 
    return(NA)
  }
}

nu_sample <- function(nu_old, logdensity_old, sd, h_T, h_0, k, lambda, n_T, 
                      alpha){
  nu_prop <- nu_old + rnorm(1, mean = 0, sd = sd)
  logdensity_prop <- logp_nu(nu_prop, h_T = h_T, h_0 = h_0, k = k, 
                             lambda = lambda, n_T = n_T, alpha = alpha)
  if (is.na(logdensity_prop)){
    return(nu_old)
  } else if (logdensity_prop - logdensity_old <= 0) {
    acceptance <- exp(logdensity_prop - logdensity_old) 
    if (runif(1, min = 0, max = 1) <= acceptance){ 
      nu_old <- nu_prop 
    } 
  } else if (logdensity_prop - logdensity_old > 0) {
    nu_old <- nu_prop
  }
  return(nu_old)
}

beta_0_sample <- function(C_inv, h_1, Q, beta_1, b){
  sigma = solve(C_inv + h_1 * Q)
  mu_beta0 = sigma %*% (h_1 * Q %*% beta_1 + C_inv %*% b) 
  mvrnorm(1, mu = mu_beta0, Sigma = sigma)
}

# h_0 
logp_h_0 <- function(h_0, h_1, lambda, nu, k, d, e){
  if (h_0 > 0){
  term_1 = 1 - lambda * h_1 / h_0 
  if (term_1 > 0 ){
    return(-0.5 * log(term_1) + (-0.5 * (nu + k) + d - 1) * log(h_0) - e * h_0) 
  } else if (term_1 < 0){
    return(NA)
    } else if (term_1 == 0) {
      return(Inf)
    }
  } else {
    return(NA)
  }
}

h_0_sample <- function(h0_old, logdensity_old, sd, h_1, nu, k, d, e){ 
  h0_prop <- h0_old + rnorm(1, mean = 0, sd = sd)
  logdensity_prop <- logp_h_0(h_0 = h0_prop, h_1 = h_1, 
                                lambda = lambda, nu = nu, k = k, d = d, e = e)
  if (is.na(logdensity_prop)){
    return(h0_old)
  } else {
    if (is.na(logdensity_old)){
      h0_old <- h0_prop
    } else if (logdensity_prop - logdensity_old < 0){
      acceptance <- exp(logdensity_prop - logdensity_old) 
      if (runif(1, min = 0, max = 1) <= acceptance){ 
        h0_old <- h0_prop 
      }  
    } else if (logdensity_prop - logdensity_old >= 0) { 
      h0_old <- h0_prop
    }
    return(h0_old)
  } 
}


logp_lambda <- function(lambda, nu, h, h_0, k, n_T){ 
  h_lag = rbind(h_0, h[1:(n_T - 1), 1, drop = FALSE]) 
  ratio =  h/h_lag
  if (sum((1 - lambda * ratio) < 0) > 0){ 
    return(NA)
  } else if (lambda <= 0) {
     return(NA)
  } else if (sum((1 - lambda * ratio) < 0) == 0) {
    logp_lambda = n_T * (nu + k) * 0.5 * log(lambda) -  
                  0.5 * sum(log(1 - lambda * ratio))
    return(logp_lambda)
  } 
}

sample_lambda <- function(lambda_old, logdensity_old, sd, nu, h, h_0, k, n_T){
  lambda_prop <- lambda_old + rnorm(1, mean = 0, sd = sd)
  logdensity_prop <- logp_lambda(lambda_prop, nu, h, h_0, k, n_T)
  if (is.na(logdensity_prop)){
    return(lambda_old)
  } else if (is.na(logdensity_old)) { 
    if (is.na(logdensity_prop)) {
      return(lambda_old)
    } else {
      return(lambda_prop)
      } 
    } else if (logdensity_prop - logdensity_old < 0){
      acceptance <- exp(logdensity_prop - logdensity_old) 
      if (runif(1, min = 0, max = 1) <= acceptance){ 
        lambda_old <- lambda_prop 
        } 
      } else if (logdensity_prop - logdensity_old >= 0){ 
      lambda_old <- lambda_prop
      }
    return(lambda_old)
}

# OS Experiment ====

kmu_os_rolling_window <- function(N, n_min, window, n_gibbs, burn, 
                                  returns, baseline, stocks, selection,
                                  k, Q, C_inv, b, d, e, alpha, lambda, 
                                  nu, sd = 10){
  rounds <- floor((N - n_min) / window)  
  # Initialize objects 
  h0_mcmc <- nu_mcmc <- S_T_mcmc <- matrix(0, ncol = 1, nrow = (n_gibbs))
  beta_mcmc <-  array(0, dim = c(n_min, k, n_gibbs))
  h_mcmc <- array(0, dim = c(n_min, 1, n_gibbs))
  N_T_mcmc <- array(0, dim = c(k, k, n_gibbs))
  returns_pred <- list()
  starts <- window * c(0:(rounds - 1))
  ends  <-  n_min + starts
  N_10 <- Q
  for (round in 1:rounds){
    print(round)
    if (round == rounds){
      window <- N - ends[rounds]
    }
    df_is <- - returns[starts[round]:ends[round], stocks] + returns[1:n_min, 
                                                                    baseline]
    df_is[, "baseline"] <- returns[starts[round]:ends[round], baseline]
    df_os <- - returns[(ends[round] + 1):(ends[round] + window), 
                       stocks] + returns[(ends[round] + 1):(ends[round] + window), 
                                         baseline]
    df_os <- - returns[(ends[round] + 1):(ends[round] + window), 
                       stocks] + returns[(ends[round] + 1):(ends[round] + window), 
                                         baseline]
    df_os[, "baseline"] <- returns[(ends[round] + 1):(ends[round] + window), 
                                   baseline]
    # Initialize
    model_1 <- lm(data = df_is, baseline ~ .)
    b_10 <- matrix(c(model_1$coefficients), nrow = k) 
    h_0 <- 1 / sum((df_is$baseline - model_1$fitted.values)^2) / (n_min - k)
    S_10  <- lambda * ((nu + 1) / nu) * (1 / h_0) 
    X_series <- as.matrix(cbind(1, df_is[, stocks]), ncol = k)
    y_series <- as.matrix(df_is[, "baseline"], ncol = k)
    remove(model_1)
    remove(df_is)
    for (i in 1:n_gibbs){
      print(i)
      # FFBS in states 
      result_filter <- filter_forward(X = X_series, y = y_series, nu = nu, Q = Q, 
                                      lambda = lambda, N_10 = N_10, S_10 = S_10, 
                                      b_10 = b_10, n_T = n_min, k = k) 
      result_bsample <- backward_sample(S_m = result_filter$S_m, 
                                        b_m = result_filter$beta_m, 
                                        N_m = result_filter$N_m, nu = nu, Q = Q, 
                                        lambda = lambda, n_T = n_min, k = k)
      beta_mcmc[ , , i] <- result_bsample[ , 2:(k + 1)]
      h_mcmc[ , , i]    <- result_bsample[ , 1]
      h_1 <- drop(result_bsample[1, 1])
      h_T <- drop(result_bsample[n_min, 1])
      beta_1 <- matrix(result_bsample[1, 2:(k+1)],nrow = k)
      h_0_help <- lambda * h_1 + rgamma(1, shape = 0.5, 
                                        scale = 2 / ((nu + 1) * S_10))
      # nu 
      nu_log <- logp_nu(nu = nu, h_T = h_T, h_0 = h_0_help, k = k, 
                        lambda = lambda, 
                        n_T = n_min, alpha = alpha)
      nu_mcmc[i, 1] <- nu <- nu_sample(nu_old = nu, logdensity_old = nu_log, 
                                       sd = sd, h_T = h_T, 
                                       h_0 = h_0_help, k = k, lambda = lambda, 
                                       n_T = n_min,
                                       alpha = alpha)
      # beta_0 
      b_10 <- beta_0_sample(C_inv = C_inv, h_1 = h_1, Q = Q, beta_1 = beta_1, 
                            b = b)
      # h_0 
      log_h0 <- logp_h_0(h_0 = h_0_help, h_1 = h_1, lambda = lambda, nu = nu, 
                         k = k, d = d, e = e)
      h0_mcmc[i, 1] <- h_0 <- h_0_sample(h0_old = h_0_help, 
                                         logdensity_old = log_h0, 
                                         sd = 0.5, h_1 = h_1, nu = nu, k = k, 
                                         d = d, e = e)
      S_10 <- (1 / h_0) * lambda * (nu + 1) / nu
      N_T_mcmc[, , i] <- result_filter$N_m[, , n_min, drop = T]
      S_T_mcmc[i, 1] <- result_filter$S_m[n_min, 1]
    }
    beta_result <- apply(beta_mcmc[, ,(n_gibbs - burn + 1):(n_gibbs)], c(1, 2), 
                         mean)
    beta_T1 <- matrix(beta_result[n_min, ], nrow = k)
    nu  <- mean(nu_mcmc[(n_gibbs - burn + 1): n_gibbs, 1])
    S_T <- mean(S_T_mcmc[(n_gibbs - burn + 1): n_gibbs, 1])
    N_T <- apply(N_T_mcmc[, ,(n_gibbs - burn + 1): n_gibbs], c(1, 2), mean)
    X_series <- as.matrix(cbind(1, df_os[, stocks]), ncol = k)
    y_series <- as.matrix(df_os[, "baseline"], ncol = k)
    N_T1 <- solve(solve(Q) + solve(lambda * N_T))
    S_T1 <- lambda * ((nu) / (nu + 1))  * S_T 
    result_filter <- filter_forward(X = X_series, y = y_series, nu = nu, Q = Q, 
                                    lambda = lambda, N_10 = N_T1, S_10 = S_T1, 
                                    b_10 = beta_T1, n_T = window, k = k)
    returns_os <- returns[(ends[round]+ 1):(ends[round] + window), selection]
    weights_pred <- t(cbind(beta_T1, result_filter$beta_m[, 1:(window - 1)]))
    weight_baseline <- 1 - rowSums(weights_pred[, 2:k])
    weights <- matrix(cbind(weight_baseline, weights_pred[, 2:k]), ncol = k)
    returns_pred[[round]] <- matrix(rowSums(returns_os * weights), nrow = 1) 
  }
  return(var(unlist(returns_pred)))
}