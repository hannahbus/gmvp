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

filter_forward <- function(X, y, nu, Q, 
                           lambda, N_10, S_10, b_10, n_T, k){ 
  # Initialize. 
  S_m <- matrix(0, nrow = n_T, ncol = 1)
  S_p <- matrix(0, nrow = 1, ncol = 1) 
  b_m <- matrix(0, nrow = k, ncol = n_T)
  b_p <- matrix(0, nrow = k, ncol = 1)
  N_m <- array(0, dim = c(k, k, n_T))
  N_p  <- N_10
  b_p[, 1]  <- b_10
  S_p[1, 1] <- S_10
  for (i in 1:n_T){ 
    # Measure.
    result_m <- measure(beta_p = b_p, N_p = N_p, S_p = S_p, 
                        X = X[i, , drop = FALSE], 
                        y = y[i, 1], k = k, nu = nu)
    b_m[, i] <- result_m[[1]]
    N_m[, , i] <- result_m[[2]]
    S_m[i, 1] <- result_m[[3]] 
    # Predict.
    result_p <- predict(beta_m = b_m[, i], N_m = N_m[, , i], 
                        nu = nu, S_m = S_m[i, 1], 
                        lambda = lambda, Q = Q)
    b_p[, 1] <- result_p[[1]]
    N_p <- result_p[[2]]
    S_p[1, 1] <- result_p[[3]]
  }
  result <- list(S_m, b_m, N_m)
  names(result) <- c("S_m", "beta_m", "N_m")
  return(result)
}

backward_sample <- function(S_m, b_m, N_m, nu, Q, lambda, n_T, k){ 
  # Initialize result vectors. 
  b_mT = b_m[, n_T, drop = FALSE]
  N_mT = N_m[, , n_T]
  beta  = matrix(nrow = n_T, ncol = k)
  h = matrix(nrow = n_T, ncol = 1)
  # Sample theta_T 
  h[1,1] = rgamma(1, shape = (nu + 1)/2, scale = 2/(S_m[n_T, 1] * (nu + 1)))
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
                  mcoptions = list(preschedule=FALSE, set.seed= TRUE)){
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
  } else if (lambda > 1) { 
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
