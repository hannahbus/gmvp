# Simulate Data ====

h_make <- function(i, vartheta, lambda, h_0){ 
  (prod(vartheta[1:i]) / (lambda^i)) * h_0
} 

eta_make <- function(i, h, k, Q){
  mvrnorm(1, c(rep(0, k)), Sigma = 1/h[i,1] * solve(Q))
}

generate_data <- function(n, k, lambda, h_0, Q, nu){ 
  X <- mvrnorm(n, c(rep(0,k - 1)), diag(k - 1))
  X <- cbind(rep(1, n), X)
  b_10 <- as.matrix(runif(k, min = 0, max = 1), nrow = k, ncol = 1)
  epsilon <- rnorm(n, mean = 0, sd = 1)
  vartheta <- rbeta(n, (nu + k)/2, 1/2)
  h <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:n){ 
    h[i,1] <- h_make(i, vartheta, lambda, h_0)
  }
  eta <- matrix(0, nrow = n, ncol = k)
  for (i in 1:n){ 
    eta[i,] <- eta_make(i, h, k = k, Q = Q)
  }
  shockacc <- apply(eta, 2, cumsum) 
  beta     <- data.frame(apply(shockacc, 1, function(x) x + b_10)) 
  u <- 1 / sqrt(h) * epsilon
  y <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:n){ 
    y[i,1] <- X[i, ] %*% beta[,i]
  }
  result <- list(y, X, b_10, beta, h) 
  names(result) <- c("y", "X", "b_10", "beta", "h")
  return(result)
}