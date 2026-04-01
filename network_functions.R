# Here are the function to compute EDC and Dn 
#Antonio

# Extreme Downside Correlation
EDC <- function(X, tau = 0.10) {
  if (!is.matrix(X)) X <- as.matrix(X)
  T_obs <- nrow(X)
  d     <- ncol(X)
  
  mu_hat    <- colMeans(X)
  sigma_hat <- apply(X, 2, sd)
  
  VaR_thresholds <- mu_hat + qnorm(tau) * sigma_hat
  X_centered     <- scale(X, center = mu_hat, scale = FALSE)
  
  indicator_matrix <- sweep(X, MARGIN = 2, STATS = VaR_thresholds, FUN = "<")
  X_trunc          <- X_centered * indicator_matrix
  
  downside_cov <- crossprod(X_trunc) / T_obs
  downside_var <- diag(downside_cov)
  
  normalization_matrix <- sqrt(outer(downside_var, downside_var))
  
  normalization_matrix[normalization_matrix == 0] <- 1
  
  EDC_matrix           <- downside_cov / normalization_matrix
  diag(EDC_matrix)     <- 1
  rownames(EDC_matrix) <- colnames(X)
  colnames(EDC_matrix) <- colnames(X)
  return(EDC_matrix)
}

# Network Density.
D_n <- function(EDC_mat) {                        
  d         <- ncol(EDC_mat)                      
  num_edges <- d * (d - 1) / 2
  unique_edges         <- EDC_mat[upper.tri(EDC_mat, diag = FALSE)]
  sum_absolute_weights <- sum(abs(unique_edges))
  density              <- sum_absolute_weights / num_edges
  return(density)
}

# Local density wrapper for rollapply
compute_local_density <- function(X_slice) {
  # FIX 1: tau must be in (0,1). Original code had tau = 2 (invalid).
  local_EDC <- EDC(X_slice, tau = 0.10)
  return(D_n(local_EDC))
}
