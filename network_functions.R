# Here are the function to compute EDC and Dn 

# Extreme Downside Correlation (EDC)
EDC <- function(X, tau = 0.10) {
  
  # Strictly enforce matrix data type for linear algebra operations
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  T_obs <- nrow(X)
  d <- ncol(X)
  
  mu_hat <- colMeans(X)
  sigma_hat <- apply(X, 2, sd)
  
  
  # qnorm maps the probability tau to the standard normal quantile (e.g., -1.28 for tau=0.10)
  z_tau <- qnorm(tau) 
  VaR_thresholds <- mu_hat + z_tau * sigma_hat
  X_centered <- scale(X, center = mu_hat, scale = FALSE)
  
  # Construct the Indicator Matrix 
  indicator_matrix <- sweep(X, MARGIN = 2, STATS = VaR_thresholds, FUN = "<")
  
  # Apply the indicator function to the centered returns
  X_trunc <- X_centered * indicator_matrix
  
  # Compute the full d x d downside covariance matrix (Numerator)
  downside_cov <- crossprod(X_trunc) / T_obs
  
  # Extract strictly the isolated downside variances (the main diagonal)
  downside_var <- diag(downside_cov)
  
  normalization_matrix <- sqrt(outer(downside_var, downside_var))
  
  # Calculate the final EDC matrix
  EDC_matrix <- downside_cov / normalization_matrix
  

  diag(EDC_matrix) <- 1
  
  # Preserve column and row names if they exist
  rownames(EDC_matrix) <- colnames(X)
  colnames(EDC_matrix) <- colnames(X)
  
  return(EDC_matrix)
}


#Network Density (D_n)

D_n <- function(EDC){

  d <- ncol(EDC_matrix)
  
  num_edges <- (d * (d - 1)) / 2
  
  unique_edges <- EDC_matrix[upper.tri(EDC_matrix, diag = FALSE)]
  
  sum_absolute_weights <- sum(abs(unique_edges))
  
  density <- sum_absolute_weights / num_edges
  
  return(density)
  
}
