# Matrix manipulations ====

convert_2_cov <- function(lower, k){
  matrix <- diag(k)
  matrix[lower.tri(matrix, diag = TRUE)] <- lower 
  matrix + t(matrix) - diag(diag(matrix))
}

# Performance measures ====

compute_performance_is <- function(beta_smooth, k, returns){
  weight_baseline <- 1 - rowSums(beta_smooth[, 2:k])
  weights <- matrix(cbind(weight_baseline, beta_smooth[, 2:k]), ncol = k)
  var(rowSums(returns * weights))
}


