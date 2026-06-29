library(zoo)
library(MASS)
library(numDeriv)

# ─────────────────────────────────────────────
# Estimating Equations
# ─────────────────────────────────────────────
H_eqn <- function(X_mat, theta, tau, gamma) {
  N_win <- nrow(X_mat)
  d <- ncol(X_mat) 
  pairs <- combn(d, 2)
  num_pairs <- ncol(pairs)
  
  mu_x   <- theta[1:d]
  sig2_x <- theta[(d + 1):(2 * d)]
  mu_y   <- theta[(2 * d + 1):(3 * d)]
  sig2_y <- theta[(3 * d + 1):(4 * d)]
  cov_pairs <- theta[(4 * d + 1):(4 * d + num_pairs)]
  sedcd  <- theta[length(theta)] 
  
  z <- sweep(X_mat, 2, mu_x, FUN="-")
  z <- sweep(z, 2, sqrt(pmax(sig2_x, 1e-6)), FUN="/")
  y <- z * (1 / (1 + exp(-(-tau - z)/ gamma)))
  
  h_mu_x   <- sweep(X_mat, 2, mu_x, FUN="-")
  h_sig2_x <- sweep(h_mu_x^2, 2, sig2_x, FUN="-")
  h_mu_y   <- sweep(y, 2, mu_y, FUN="-")
  h_sig2_y <- sweep(h_mu_y^2, 2, sig2_y, FUN="-")
  
  
  # Estimating equation for Cov(Yi, Yj) which is γ_ij in the paper
  # H_cov = Y_i * Y_j - E[Y_i]E[Y_j] - Cov(Y_i, Y_j)
  # Since E[Y_i] = mu_y[i] and E[Y_j] = mu_y[j], and we are at the true parameter,
  # the equation simplifies. However, for the H_eqn we use the general form.
  y_prod <- y[, pairs[1,], drop=FALSE] * y[, pairs[2,], drop=FALSE]
  mu_y_prod <- mu_y[pairs[1,]] * mu_y[pairs[2,]]
  h_cov <- y_prod - matrix(mu_y_prod, nrow=N_win, ncol=num_pairs, byrow=TRUE) - matrix(cov_pairs, nrow=N_win, ncol=num_pairs, byrow=TRUE)
  
  # Estimating equation for SEDCD
  var_prods <- sig2_y[pairs[1,]] * sig2_y[pairs[2,]]
  cor_ij <- cov_pairs / sqrt(pmax(var_prods, 1e-12))
  sedcd_formula <- mean(abs(cor_ij))
  h_sedcd <- rep(sedcd - sedcd_formula, N_win)
  
  return(cbind(h_mu_x, h_sig2_x, h_mu_y, h_sig2_y, h_cov, h_sedcd))
}

solve_theta <- function(X_window, tau = 3, gamma=0.1) {
  d <- ncol(X_window)
  mu_x <- colMeans(X_window)
  sig2_x <- pmax(apply(X_window, 2, var), 1e-6) 
  
  Z_window <- sweep(X_window, 2, mu_x, FUN="-")
  Z_window <- sweep(Z_window, 2, sqrt(sig2_x), FUN="/")
  Y_window <- Z_window * (1 / (1 + exp(-(-tau -Z_window) / gamma)))
  
  mu_y <- colMeans(Y_window)
  sig2_y <- pmax(apply(Y_window, 2, var), 1e-6)
  
  cov_y <- cov(Y_window)
  pairs_idx <- t(combn(d, 2))
  cov_pairs <- cov_y[pairs_idx]
  
  # Calculate SEDCD from the estimated covariances and variances
  var_prods <- sig2_y[pairs_idx[, 1]] * sig2_y[pairs_idx[, 2]]
  cor_ij <- cov_pairs / sqrt(pmax(var_prods, 1e-12))
  sedcd <- mean(abs(cor_ij))
  
  return(c(mu_x, sig2_x, mu_y, sig2_y, cov_pairs, sedcd))
}

get_theta <- function(x, k, tau=3, gamma=0.1) {
  N <- nrow(x)
  k <- ceiling(k) # bandwidth
  d <- ncol(x)
  P <- 4 * d + ncol(combn(d, 2)) + 1 
  res <- matrix(0, nrow=N, ncol=P)
  
  for(t in k:N) {
    res[t, ] <- solve_theta(x[(t-k+1):t, , drop=FALSE], tau, gamma)
  }
  return(res)
}

# ─────────────────────────────────────────────
# 2. Integration and Variance Estimation
# ─────────────────────────────────────────────

DH_mean_multi <- function(X_window, theta, tau, gamma) {
  mean_H <- function(th) colMeans(H_eqn(X_window, th, tau, gamma))
  return(jacobian(mean_H, theta, method="simple"))
}

int.par.implicit.multi <- function(x, thetahat, lag=1, cutoff=1, k, tau, gamma) {
  N <- nrow(x)
  P <- ncol(thetahat)
  spot_estimators <- matrix(0, nrow=N, ncol=P)
  ridge_penalty <- diag(1e-5, P) 
  
  for(t in (cutoff+1):N) {
    th_hat <- thetahat[t - lag, ]
    window_start <- max(1, t - lag - k + 1)
    X_window <- x[window_start:(t - lag), , drop=FALSE]
    
    DH_mat <- DH_mean_multi(X_window, th_hat, tau, gamma) + ridge_penalty
    DH_inv <- solve(DH_mat)
    
    H_val <- as.vector(H_eqn(x[t, , drop=FALSE], th_hat, tau, gamma))
    spot_estimators[t, ] <- th_hat - as.vector(DH_inv %*% H_val)
  }
  return(apply(spot_estimators, 2, cumsum) / N)
}


estimate_var <- function(x, thetahat, block=1, lag=1, cutoff=1, k, tau, gamma) {
  N <- nrow(x)
  P <- ncol(thetahat)

  
  q_val <- matrix(0, nrow=N, ncol=P)
  ridge_penalty <- diag(1e-5, P)
  
  for(t in (cutoff+1):N) {
    th_hat <- thetahat[t - block - lag, ]
    window_start <- max(1, t - block - lag - k + 1)
    X_window <- x[window_start:(t - block - lag), , drop=FALSE]
    
    DH_mat <- DH_mean_multi(X_window, th_hat, tau, gamma) + ridge_penalty
    DH_inv <- solve(DH_mat)
    
    blocksum <- rep(0, P)
    for(i in (t - block + 1):t) {
      H_val <- as.vector(H_eqn(x[i, , drop=FALSE], th_hat, tau, gamma))
      blocksum <- blocksum + H_val
    }
    influence_sum <- -as.vector(DH_inv %*% blocksum)
    q_val[t, ] <- (influence_sum^2) / block
  }
  return(apply(q_val, 2, cumsum) / N)
}
