# ==============================================================================
# File: sedcd_mscale.R
# Description: Refactored SEDCD Estimating Equations (Option B: Robust M-Scale)
#              for Bivariate (d = 2) Locally Stationary Time Series.
#
# Mathematical Specification:
#   Parameter vector theta in R^11:
#     theta = (mu_1, mu_2, sigma_1, sigma_2, nu_1, nu_2,
#              v_11^Y, v_12^Y, v_21^Y, v_22^Y, SEDCD_tau)^T
#
#   Estimating function H(X, theta) in R^11:
#     1-2:   X_i - mu_i
#     3-4:   chi_c((X_i - mu_i) / sigma_i)
#     5-6:   Y_i - nu_i
#     7-10:  (Y_i - nu_i)(Y_j - nu_j) - v_ij^Y
#     11:    SEDCD_tau - |v_12^Y| / sqrt(v_11^Y * v_22^Y)
#
#   where chi_c(u) = 1 - exp(-(u/c)^2) - delta (Dennis-Welsch M-scale)
#   and delta = 1 - (1 + 2/c^2)^(-1/2) ensures Fisher consistency under N(0, 1).
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("MASS", quietly = TRUE)) library(MASS)
  if (requireNamespace("numDeriv", quietly = TRUE)) library(numDeriv)
})

# ------------------------------------------------------------------------------
# 1. Dennis-Welsch M-Scale Dispersion & Calibration
# ------------------------------------------------------------------------------

#' Dennis-Welsch Fisher Consistency Constant under Standard Normality
#'
#' @param c Numeric tuning parameter (default 2.985 for 95% Gaussian efficiency)
#' @return Numeric scalar delta
welsch_delta_constant <- function(c = 2.985) {
  1.0 - (1.0 + 2.0 / (c^2))^(-0.5)
}

#' Centered Dennis-Welsch M-Scale Score Function chi_c(u)
#'
#' @param u Numeric vector of standardized residuals
#' @param c Numeric tuning constant (default 2.985)
#' @return Numeric vector chi_c(u) bounded in [-delta, 1 - delta]
score_welsch_chi <- function(u, c = 2.985) {
  delta <- welsch_delta_constant(c)
  rho <- 1.0 - exp(-(u / c)^2)
  res <- rho - delta
  attr(res, "delta") <- delta
  attr(res, "c") <- c
  return(res)
}

#' Derivative of chi_c(u) with respect to u
#' chi_c'(u) = (2 * u / c^2) * exp(-(u/c)^2)
score_welsch_chi_prime <- function(u, c = 2.985) {
  (2.0 * u / (c^2)) * exp(-(u / c)^2)
}

#' Smooth Logistic Downside Transformation
#'
#' @param z Standardized variable (X - mu) / sigma
#' @param tau Threshold parameter (typically < 0, e.g. lower 10th percentile)
#' @param gamma Smoothing parameter (default 0.05)
#' @return Numeric vector S_gamma(tau - z) in [0, 1]
sigmoid_smooth <- function(z, tau, gamma = 0.05) {
  w <- (tau - z) / gamma
  # Numerical stability for exp: clamp between -50 and 50
  w_clamped <- pmax(pmin(w, 50), -50)
  1.0 / (1.0 + exp(-w_clamped))
}

#' Downside Variable Y = Z * S_gamma(tau - Z)
#'
#' @param z Standardized variable
#' @param tau Threshold parameter
#' @param gamma Smoothing parameter
#' @return Numeric vector Y
downside_transform_z <- function(z, tau, gamma = 0.05) {
  z * sigmoid_smooth(z, tau = tau, gamma = gamma)
}

# ------------------------------------------------------------------------------
# 2. Pilot 1D Robust M-Scale Root Finder
# ------------------------------------------------------------------------------

#' Solve 1D Robust M-Scale Equation for a Single Marginal Series
#'
#' Finds sigma > 0 such that (1/N) * sum(chi_c((x - mu) / sigma)) = 0.
#' Monotonically decreasing in sigma > 0, ensuring a unique root.
#'
#' @param x Numeric vector of observations
#' @param mu Center location (scalar)
#' @param c Tuning parameter (default 2.985)
#' @param tol Convergence tolerance (default 1e-6)
#' @return Numeric scalar sigma_hat > 0
solve_pilot_mscale <- function(x, mu = mean(x), c = 2.985, tol = 1e-6) {
  x <- as.numeric(x)
  diff_vec <- x - mu
  n_obs <- length(diff_vec)
  
  if (n_obs == 0) return(NA_real_)
  
  # Initial estimate via normalized MAD
  mad_val <- median(abs(diff_vec - median(diff_vec))) / 0.6745
  if (is.na(mad_val) || mad_val < 1e-4) {
    mad_val <- sd(diff_vec)
  }
  if (is.na(mad_val) || mad_val < 1e-4) {
    mad_val <- 1.0
  }
  
  # Objective function: mean(chi_c(diff / sigma))
  obj_fn <- function(sigma) {
    if (sigma <= 1e-8) return(1.0)
    u <- diff_vec / sigma
    mean(score_welsch_chi(u, c = c))
  }
  
  # Determine search bracket [lower, upper]
  lower <- max(1e-5, 0.02 * mad_val)
  upper <- max(1.0, 50.0 * mad_val)
  
  f_lower <- obj_fn(lower)
  f_upper <- obj_fn(upper)
  
  # Expand lower bound if needed
  expand_count <- 0
  while (f_lower <= 0 && expand_count < 10) {
    lower <- lower / 5.0
    f_lower <- obj_fn(lower)
    expand_count <- expand_count + 1
  }
  
  # Expand upper bound if needed
  expand_count <- 0
  while (f_upper >= 0 && expand_count < 10) {
    upper <- upper * 5.0
    f_upper <- obj_fn(upper)
    expand_count <- expand_count + 1
  }
  
  if (f_lower * f_upper >= 0) {
    # Fallback to MAD if root could not be bracketed (degenerate data)
    return(mad_val)
  }
  
  root_res <- tryCatch({
    uniroot(obj_fn, interval = c(lower, upper), tol = tol)$root
  }, error = function(e) {
    mad_val
  })
  
  return(max(1e-5, root_res))
}

# ------------------------------------------------------------------------------
# 3. Local Pilot Estimation Algorithm (Rolling Window)
# ------------------------------------------------------------------------------

#' Estimate 11-Dimensional Parameter Vector theta on a Window
#'
#' Parameters in order:
#' 1: mu_1
#' 2: mu_2
#' 3: sigma_1 (robust M-scale)
#' 4: sigma_2 (robust M-scale)
#' 5: nu_1
#' 6: nu_2
#' 7: v_11^Y
#' 8: v_12^Y
#' 9: v_21^Y
#' 10: v_22^Y
#' 11: SEDCD_tau = |v_12^Y| / sqrt(v_11^Y * v_22^Y)
#'
#' @param X_window Matrix of dimensions M_n x 2
#' @param tau Threshold parameter
#' @param gamma Smoothing parameter (default 0.05)
#' @param c Tuning parameter (default 2.985)
#' @return Named numeric vector of length 11
estimate_local_theta <- function(X_window, tau, gamma = 0.05, c = 2.985) {
  X_window <- as.matrix(X_window)
  if (ncol(X_window) != 2) {
    stop(sprintf("estimate_local_theta expects bivariate data (d = 2), got %d columns", ncol(X_window)))
  }
  M_n <- nrow(X_window)
  
  # Step 1: Locations
  mu1 <- mean(X_window[, 1])
  mu2 <- mean(X_window[, 2])
  
  # Step 2: Robust M-Scales
  sigma1 <- solve_pilot_mscale(X_window[, 1], mu = mu1, c = c)
  sigma2 <- solve_pilot_mscale(X_window[, 2], mu = mu2, c = c)
  
  # Step 3: Downside Variables
  z1 <- (X_window[, 1] - mu1) / sigma1
  z2 <- (X_window[, 2] - mu2) / sigma2
  
  Y1 <- downside_transform_z(z1, tau = tau, gamma = gamma)
  Y2 <- downside_transform_z(z2, tau = tau, gamma = gamma)
  
  # Step 4: Downside Moments (nu)
  nu1 <- mean(Y1)
  nu2 <- mean(Y2)
  
  # Step 5: Downside Covariances (v_ij^Y dividing by M_n)
  dY1 <- Y1 - nu1
  dY2 <- Y2 - nu2
  
  v11 <- mean(dY1^2)
  v12 <- mean(dY1 * dY2)
  v21 <- v12
  v22 <- mean(dY2^2)
  
  # Step 6: Smoothed Extremal Downside Correlation Density
  denom <- sqrt(max(1e-12, v11 * v22))
  sedcd_val <- abs(v12) / denom
  
  theta <- c(
    mu1       = mu1,
    mu2       = mu2,
    sigma1    = sigma1,
    sigma2    = sigma2,
    nu1       = nu1,
    nu2       = nu2,
    v11_Y     = v11,
    v12_Y     = v12,
    v21_Y     = v21,
    v22_Y     = v22,
    SEDCD_tau = sedcd_val
  )
  
  return(theta)
}

# ------------------------------------------------------------------------------
# 4. 11-Dimensional Estimating Function H(X_t, theta)
# ------------------------------------------------------------------------------

#' Compute Estimating Function Vector or Matrix H(X, theta)
#'
#' @param X_obs Matrix (N x 2) or vector of length 2
#' @param theta Vector of length 11
#' @param tau Threshold parameter
#' @param gamma Smoothing parameter (default 0.05)
#' @param c Tuning parameter (default 2.985)
#' @return Matrix of dimensions N x 11 (or numeric vector of length 11)
compute_H_vector <- function(X_obs, theta, tau, gamma = 0.05, c = 2.985) {
  is_vec <- is.null(dim(X_obs))
  if (is_vec) {
    X_mat <- matrix(X_obs, nrow = 1)
  } else {
    X_mat <- as.matrix(X_obs)
  }
  
  N <- nrow(X_mat)
  mu1    <- theta[1]
  mu2    <- theta[2]
  sigma1 <- max(1e-6, theta[3])
  sigma2 <- max(1e-6, theta[4])
  nu1    <- theta[5]
  nu2    <- theta[6]
  v11    <- max(1e-12, theta[7])
  v12    <- theta[8]
  v21    <- theta[9]
  v22    <- max(1e-12, theta[10])
  sedcd  <- theta[11]
  
  # Standardized variables
  z1 <- (X_mat[, 1] - mu1) / sigma1
  z2 <- (X_mat[, 2] - mu2) / sigma2
  
  # Downside variables
  Y1 <- downside_transform_z(z1, tau = tau, gamma = gamma)
  Y2 <- downside_transform_z(z2, tau = tau, gamma = gamma)
  
  # Centered downside variables
  dY1 <- Y1 - nu1
  dY2 <- Y2 - nu2
  
  # 1-2: Location equations
  h1 <- X_mat[, 1] - mu1
  h2 <- X_mat[, 2] - mu2
  
  # 3-4: M-Scale equations (Option B: Bounded Dennis-Welsch)
  h3 <- score_welsch_chi(z1, c = c)
  h4 <- score_welsch_chi(z2, c = c)
  
  # 5-6: Downside mean equations
  h5 <- dY1
  h6 <- dY2
  
  # 7-10: Downside covariance equations
  h7  <- dY1^2 - v11
  h8  <- dY1 * dY2 - v12
  h9  <- dY2 * dY1 - v21
  h10 <- dY2^2 - v22
  
  # 11: Target SEDCD density equation
  denom <- sqrt(v11 * v22)
  sedcd_implied <- abs(v12) / denom
  h11 <- rep(sedcd - sedcd_implied, N)
  
  H_mat <- cbind(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11)
  colnames(H_mat) <- c("H_mu1", "H_mu2", "H_sigma1", "H_sigma2",
                       "H_nu1", "H_nu2", "H_v11", "H_v12", "H_v21", "H_v22", "H_sedcd")
  
  if (is_vec) return(as.vector(H_mat[1, ]))
  return(H_mat)
}

# ------------------------------------------------------------------------------
# 5. Empirical Jacobian DH in R^(11 x 11)
# ------------------------------------------------------------------------------

#' Compute Empirical Jacobian DH = (1/M_n) * sum (dH(X_s, theta) / dtheta^T)
#'
#' Utilizes closed-form analytical derivatives of the block-triangular structure,
#' with optional verification via finite differences (numDeriv).
#'
#' @param X_window Matrix M_n x 2
#' @param theta Vector of length 11
#' @param tau Threshold parameter
#' @param gamma Smoothing parameter (default 0.05)
#' @param c Tuning parameter (default 2.985)
#' @param method "analytic" (default, fast) or "numDeriv"
#' @return Matrix 11 x 11
compute_jacobian_DH <- function(X_window, theta, tau, gamma = 0.05, c = 2.985,
                                method = c("analytic", "numDeriv")) {
  method <- match.arg(method)
  X_window <- as.matrix(X_window)
  M_n <- nrow(X_window)
  
  if (method == "numDeriv") {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("Package 'numDeriv' required for method = 'numDeriv'")
    }
    mean_H <- function(th) {
      colMeans(compute_H_vector(X_window, th, tau = tau, gamma = gamma, c = c))
    }
    return(numDeriv::jacobian(mean_H, theta, method = "Richardson"))
  }
  
  # Analytical closed-form calculation
  mu1    <- theta[1]
  mu2    <- theta[2]
  sigma1 <- max(1e-6, theta[3])
  sigma2 <- max(1e-6, theta[4])
  nu1    <- theta[5]
  nu2    <- theta[6]
  v11    <- max(1e-12, theta[7])
  v12    <- theta[8]
  v21    <- theta[9]
  v22    <- max(1e-12, theta[10])
  sedcd  <- theta[11]
  
  z1 <- (X_window[, 1] - mu1) / sigma1
  z2 <- (X_window[, 2] - mu2) / sigma2
  
  # Logistic weights and derivatives
  S1 <- sigmoid_smooth(z1, tau = tau, gamma = gamma)
  S2 <- sigmoid_smooth(z2, tau = tau, gamma = gamma)
  
  # Derivative dY_i / dZ_i
  # Y_i = Z_i * S_gamma(tau - Z_i)
  # dY_i / dZ_i = S_i - (Z_i / gamma) * S_i * (1 - S_i)
  g1_prime <- S1 - (z1 / gamma) * S1 * (1.0 - S1)
  g2_prime <- S2 - (z2 / gamma) * S2 * (1.0 - S2)
  
  Y1 <- z1 * S1
  Y2 <- z2 * S2
  dY1 <- Y1 - nu1
  dY2 <- Y2 - nu2
  
  # Pre-allocate 11x11 Jacobian matrix
  DH <- matrix(0, nrow = 11, ncol = 11)
  
  # Block 1 (Rows 1-2): H_mu = X - mu
  DH[1, 1] <- -1.0
  DH[2, 2] <- -1.0
  
  # Block 2 (Rows 3-4): H_sigma = chi_c(Z)
  # chi_c'(Z) = (2 * Z / c^2) * exp(-(Z/c)^2)
  chi1_p <- score_welsch_chi_prime(z1, c = c)
  chi2_p <- score_welsch_chi_prime(z2, c = c)
  
  # dH_3/dmu1 = -(1/sigma1) * chi_c'(Z1)
  DH[3, 1] <- mean(-chi1_p / sigma1)
  # dH_3/dsigma1 = -(Z1/sigma1) * chi_c'(Z1)
  DH[3, 3] <- mean(-z1 * chi1_p / sigma1)
  
  # dH_4/dmu2 = -(1/sigma2) * chi_c'(Z2)
  DH[4, 2] <- mean(-chi2_p / sigma2)
  # dH_4/dsigma2 = -(Z2/sigma2) * chi_c'(Z2)
  DH[4, 4] <- mean(-z2 * chi2_p / sigma2)
  
  # Block 3 (Rows 5-6): H_nu = Y - nu
  # dY_i / dmu_i = -g_i' / sigma_i
  # dY_i / dsigma_i = -Z_i * g_i' / sigma_i
  DH[5, 1] <- mean(-g1_prime / sigma1)
  DH[5, 3] <- mean(-z1 * g1_prime / sigma1)
  DH[5, 5] <- -1.0
  
  DH[6, 2] <- mean(-g2_prime / sigma2)
  DH[6, 4] <- mean(-z2 * g2_prime / sigma2)
  DH[6, 6] <- -1.0
  
  # Block 4 (Rows 7-10): Covariances
  # Row 7: H_v11 = (Y1 - nu1)^2 - v11
  DH[7, 1] <- mean(2.0 * dY1 * (-g1_prime / sigma1))
  DH[7, 3] <- mean(2.0 * dY1 * (-z1 * g1_prime / sigma1))
  DH[7, 5] <- mean(-2.0 * dY1)
  DH[7, 7] <- -1.0
  
  # Row 8: H_v12 = (Y1 - nu1)(Y2 - nu2) - v12
  DH[8, 1] <- mean(dY2 * (-g1_prime / sigma1))
  DH[8, 2] <- mean(dY1 * (-g2_prime / sigma2))
  DH[8, 3] <- mean(dY2 * (-z1 * g1_prime / sigma1))
  DH[8, 4] <- mean(dY1 * (-z2 * g2_prime / sigma2))
  DH[8, 5] <- mean(-dY2)
  DH[8, 6] <- mean(-dY1)
  DH[8, 8] <- -1.0
  
  # Row 9: H_v21 = (Y2 - nu2)(Y1 - nu1) - v21
  DH[9, 1] <- DH[8, 1]
  DH[9, 2] <- DH[8, 2]
  DH[9, 3] <- DH[8, 3]
  DH[9, 4] <- DH[8, 4]
  DH[9, 5] <- DH[8, 5]
  DH[9, 6] <- DH[8, 6]
  DH[9, 9] <- -1.0
  
  # Row 10: H_v22 = (Y2 - nu2)^2 - v22
  DH[10, 2] <- mean(2.0 * dY2 * (-g2_prime / sigma2))
  DH[10, 4] <- mean(2.0 * dY2 * (-z2 * g2_prime / sigma2))
  DH[10, 6] <- mean(-2.0 * dY2)
  DH[10, 10] <- -1.0
  
  # Block 5 (Row 11): H_sedcd = sedcd - |v12| / sqrt(v11 * v22)
  denom <- sqrt(v11 * v22)
  R_val <- abs(v12) / denom
  sign_v12 <- if (abs(v12) < 1e-12) 0.0 else sign(v12)
  
  DH[11, 7]  <- 0.5 * R_val / v11
  DH[11, 8]  <- -sign_v12 / denom
  DH[11, 10] <- 0.5 * R_val / v22
  DH[11, 11] <- 1.0
  
  return(DH)
}

# ------------------------------------------------------------------------------
# 6. Recursive Rolling-Window CUSUM Pipeline
# ------------------------------------------------------------------------------

#' Compute Local Pilot Theta Across Rolling Windows
#'
#' @param X Data matrix N x 2
#' @param k Rolling window bandwidth
#' @param tau Threshold parameter
#' @param gamma Smoothing parameter
#' @param c Tuning parameter
#' @return Matrix N x 11
get_local_theta_series <- function(X, k, tau, gamma = 0.05, c = 2.985) {
  X <- as.matrix(X)
  N <- nrow(X)
  res <- matrix(NA_real_, nrow = N, ncol = 11)
  
  for (t in k:N) {
    win_data <- X[(t - k + 1):t, , drop = FALSE]
    res[t, ] <- estimate_local_theta(win_data, tau = tau, gamma = gamma, c = c)
  }
  return(res)
}

#' Linearized Recursive Estimator Path for Parameter 11 (SEDCD)
#'
#' @param X Data matrix N x 2
#' @param teta_mat Matrix N x 11 of pilot estimates
#' @param lag Decoupling lag L_n
#' @param cutoff Starting cutoff
#' @param k Window size M_n
#' @param tau Threshold parameter
#' @param gamma Smoothing parameter
#' @param c Tuning parameter
#' @return Numeric vector of length N representing cumulative linearized sum
int_par_sedcd_mscale <- function(X, teta_mat, lag = 1, cutoff = NULL, k,
                                 tau, gamma = 0.05, c = 2.985) {
  N <- nrow(X)
  if (is.null(cutoff)) cutoff <- k + lag
  if (cutoff < k + lag) cutoff <- k + lag
  
  lin_teta <- matrix(0.0, nrow = N, ncol = 11)
  
  for (t in cutoff:N) {
    param_idx <- t - lag
    teta_lag  <- teta_mat[param_idx, ]
    
    ws <- max(1, param_idx - k + 1)
    X_window <- X[ws:param_idx, , drop = FALSE]
    
    DH <- compute_jacobian_DH(X_window, teta_lag, tau = tau, gamma = gamma, c = c)
    DH_inv <- tryCatch(solve(DH), error = function(e) MASS::ginv(DH))
    
    H_val <- compute_H_vector(X[t, ], teta_lag, tau = tau, gamma = gamma, c = c)
    correction <- as.vector(DH_inv %*% H_val)
    lin_teta[t, ] <- teta_lag - correction
  }
  
  # Target functional is SEDCD (index 11)
  sedcd_lin <- lin_teta[, 11]
  return(cumsum(sedcd_lin) / N)
}

#' Integrated Long-Run Variance Estimation for SEDCD
#'
#' @param X Data matrix N x 2
#' @param teta_mat Matrix N x 11 of pilot estimates
#' @param block Block size b_n
#' @param lag Decoupling lag L_n
#' @param cutoff Starting cutoff
#' @param k Window size M_n
#' @param tau Threshold parameter
#' @param gamma Smoothing parameter
#' @param c Tuning parameter
#' @return Numeric vector Q_n(u) for parameter 11 (SEDCD)
var_est_sedcd_mscale <- function(X, teta_mat, block = 1, lag = 1, cutoff = NULL, k,
                                 tau, gamma = 0.05, c = 2.985) {
  N <- nrow(X)
  min_cutoff <- k + lag + block
  if (is.null(cutoff) || cutoff < min_cutoff) cutoff <- min_cutoff
  
  q_sq <- numeric(N)
  
  for (t in cutoff:N) {
    param_idx <- t - lag - block
    teta_lag  <- teta_mat[param_idx, ]
    
    ws <- max(1, param_idx - k + 1)
    X_window <- X[ws:param_idx, , drop = FALSE]
    
    DH <- compute_jacobian_DH(X_window, teta_lag, tau = tau, gamma = gamma, c = c)
    DH_inv <- tryCatch(solve(DH), error = function(e) MASS::ginv(DH))
    
    # Block estimating vector sum
    block_data <- X[(t - block + 1):t, , drop = FALSE]
    H_block <- compute_H_vector(block_data, teta_lag, tau = tau, gamma = gamma, c = c)
    H_sum <- colSums(H_block)
    
    q_par <- -as.vector(DH_inv %*% H_sum)
    # Target parameter 11
    q_sq[t] <- q_par[11]^2
  }
  
  q_scaled <- q_sq / block
  Q_est <- cumsum(q_scaled) / N
  return(Q_est)
}

#' CUSUM Change-Point Test for SEDCD with Multiplier Bootstrap
#'
#' @param X Data matrix N x 2
#' @param tau Threshold parameter
#' @param k Bandwidth window size M_n (default floor(N^0.5))
#' @param lag Decoupling lag L_n (default max(1, floor(0.1 * (log(N))^2)))
#' @param block Block size b_n (default max(1, floor(0.1 * (log(N))^2)))
#' @param gamma Smoothing parameter (default 0.05)
#' @param c Tuning parameter (default 2.985)
#' @param MC Number of multiplier bootstrap replications (default 500)
#' @param plotting Logical, whether to produce diagnostic CUSUM plot
#' @return List containing test_stat, p_value, q95, max_index, Tu, Qn
cusum_sedcd_test <- function(X, tau, k = NULL, lag = NULL, block = NULL,
                             gamma = 0.05, c = 2.985, MC = 500, plotting = FALSE) {
  X <- as.matrix(X)
  N <- nrow(X)
  d <- ncol(X)
  
  if (d > 2) {
    return(cusum_sedcd_multivariate(X = X, tau = tau, k = k, lag = lag, block = block,
                                     gamma = gamma, c = c, MC = MC, plotting = plotting))
  }
  
  if (is.null(k)) k <- max(10, floor(N^0.5))
  if (is.null(lag)) lag <- max(1, floor(0.1 * (log(N))^2))
  if (is.null(block)) block <- max(1, floor(0.1 * (log(N))^2))
  
  cutoff <- k + lag + block
  
  # 1. Pilot estimation
  teta_mat <- get_local_theta_series(X, k = k, tau = tau, gamma = gamma, c = c)
  
  # 2. Linearized cumulative process
  Mn <- int_par_sedcd_mscale(X, teta_mat, lag = lag, cutoff = cutoff, k = k,
                             tau = tau, gamma = gamma, c = c)
  
  # 3. Variance process
  Qn <- var_est_sedcd_mscale(X, teta_mat, block = block, lag = lag, cutoff = cutoff, k = k,
                             tau = tau, gamma = gamma, c = c)
  
  qn <- pmax(1e-8, diff(c(0, Qn)))
  
  # 4. CUSUM Bridge Statistic
  valid_range <- (cutoff + 1):N
  denom_bridge <- N - cutoff
  bridge <- Mn[valid_range] - ((1:length(valid_range)) / denom_bridge) * Mn[N]
  Tu_valid <- sqrt(N) * bridge
  Tu <- c(rep(0.0, cutoff), Tu_valid)
  
  test_stat <- max(abs(Tu_valid))
  max_idx <- cutoff + which.max(abs(Tu_valid))
  
  # 5. Multiplier Gaussian Bootstrap
  Z_mc <- numeric(MC)
  for (m in 1:MC) {
    xi <- rnorm(N)
    BM <- cumsum(xi * sqrt(qn))
    bridge_mc <- (BM[valid_range] - BM[cutoff]) - ((1:length(valid_range)) / denom_bridge) * (BM[N] - BM[cutoff])
    Z_mc[m] <- max(abs(bridge_mc))
  }
  
  p_val <- mean(Z_mc >= test_stat)
  q90 <- quantile(Z_mc, 0.90)
  q95 <- quantile(Z_mc, 0.95)
  
  if (plotting) {
    u_grid <- (1:N) / N
    plot(u_grid, Tu, type = "l", col = "blue", lwd = 1.5,
         xlab = "Rescaled Time u = t/n", ylab = "T_n(u)",
         main = sprintf("CUSUM SEDCD Test (p-val = %.3f)", p_val))
    abline(h = c(q95, -q95), col = "red", lty = 2)
    abline(v = max_idx / N, col = "darkgreen", lty = 3)
  }
  
  return(list(
    test_stat = test_stat,
    p_value   = p_val,
    q95       = as.numeric(q95),
    q90       = as.numeric(q90),
    max_index = max_idx,
    tau_idx   = max_idx,
    break_u   = max_idx / N,
    tau_u     = max_idx / N,
    Tu        = Tu,
    Qn        = Qn,
    teta_mat  = teta_mat,
    sedcd     = teta_mat[, 11],
    reject_95 = (p_val < 0.05),
    reject_90 = (p_val < 0.10)
  ))
}

# ------------------------------------------------------------------------------
# 7. General Multivariate (d >= 2) Linearized SEDCD CUSUM Pipeline
# ------------------------------------------------------------------------------

#' Fast 1D M-Scale Root Finder for Rolling Windows
solve_mscale_fast <- function(x_diff, mad_val, c = 2.985, tol = 1e-5) {
  delta_val <- 1.0 - (1.0 + 2.0 / (c^2))^(-0.5)
  diff_sq_over_c2 <- (x_diff^2) / (c^2)
  low <- max(1e-4, 0.05 * mad_val)
  high <- max(1.0, 20.0 * mad_val)
  
  f_low <- mean(1.0 - exp(-diff_sq_over_c2 / (low^2))) - delta_val
  f_high <- mean(1.0 - exp(-diff_sq_over_c2 / (high^2))) - delta_val
  
  if (f_low * f_high >= 0) return(mad_val)
  
  uniroot(function(s) mean(1.0 - exp(-diff_sq_over_c2 / (s^2))) - delta_val,
          lower = low, upper = high, tol = tol)$root
}

#' Multivariate Local Pilot Estimation for Arbitrary Dimension d
estimate_local_theta_multi <- function(X_win, tau, gamma = 0.05, c = 2.985) {
  X_win <- as.matrix(X_win)
  d <- ncol(X_win)
  n_win <- nrow(X_win)
  mu_vec <- colMeans(X_win)
  sigma_vec <- numeric(d)
  
  for (j in 1:d) {
    diff_j <- X_win[, j] - mu_vec[j]
    mad_j <- median(abs(diff_j - median(diff_j))) / 0.6745
    if (is.na(mad_j) || mad_j < 1e-4) mad_j <- sd(diff_j)
    sigma_vec[j] <- solve_mscale_fast(diff_j, mad_j, c = c)
  }
  
  Z <- sweep(sweep(X_win, 2, mu_vec, "-"), 2, sigma_vec, "/")
  Y <- Z * sigmoid_smooth(Z, tau = tau, gamma = gamma)
  nu_vec <- colMeans(Y)
  dY <- sweep(Y, 2, nu_vec, "-")
  V_mat <- crossprod(dY) / n_win
  
  sd_vec <- sqrt(pmax(1e-12, diag(V_mat)))
  denom <- sd_vec %*% t(sd_vec)
  cor_mat <- V_mat / denom
  sedcd_val <- (sum(abs(cor_mat)) - d) / (d * (d - 1))
  
  return(c(sedcd_val, mu_vec, sigma_vec, nu_vec, as.vector(V_mat)))
}

#' Multivariate Estimating Function H(X, theta) in R^(d^2 + 3d + 1)
compute_H_multi <- function(X_obs, theta, tau, gamma = 0.05, c = 2.985) {
  if (is.null(dim(X_obs))) X_obs <- matrix(X_obs, nrow = 1)
  d <- ncol(X_obs)
  n_obs <- nrow(X_obs)
  
  sedcd     <- theta[1]
  mu_vec    <- theta[2:(1 + d)]
  sigma_vec <- pmax(1e-6, theta[(2 + d):(1 + 2 * d)])
  nu_vec    <- theta[(2 + 2 * d):(1 + 3 * d)]
  V_mat     <- matrix(theta[(2 + 3 * d):length(theta)], nrow = d, ncol = d)
  
  sd_vec <- sqrt(pmax(1e-12, diag(V_mat)))
  denom <- sd_vec %*% t(sd_vec)
  cor_mat <- V_mat / denom
  smooth_abs_eps <- 1e-8
  cor_abs <- sqrt(cor_mat^2 + smooth_abs_eps)
  sedcd_implied <- (sum(cor_abs) - sum(diag(cor_abs))) / (d * (d - 1))
  h_sedcd <- matrix(sedcd - sedcd_implied, nrow = n_obs, ncol = 1)
  
  h_means <- sweep(X_obs, 2, mu_vec, "-")
  Z <- sweep(sweep(X_obs, 2, mu_vec, "-"), 2, sigma_vec, "/")
  h_sigmas <- score_welsch_chi(Z, c = c)
  
  Y <- Z * sigmoid_smooth(Z, tau = tau, gamma = gamma)
  h_nus <- sweep(Y, 2, nu_vec, "-")
  
  dY <- sweep(Y, 2, nu_vec, "-")
  h_covs <- matrix(0, nrow = n_obs, ncol = d * d)
  col_idx <- 1
  for (j in 1:d) {
    for (i in 1:d) {
      h_covs[, col_idx] <- dY[, i] * dY[, j] - V_mat[i, j]
      col_idx <- col_idx + 1
    }
  }
  
  return(cbind(h_sedcd, h_means, h_sigmas, h_nus, h_covs))
}

#' Multivariate Linearized SEDCD CUSUM Test with Robust Dennis-Welsch M-Scale (Option B)
#'
#' @param X Data matrix N x d
#' @param tau Threshold parameter
#' @param k Rolling window bandwidth M_n
#' @param lag Separation lag L_n
#' @param block Block size b_n
#' @param gamma Smoothing parameter
#' @param c Tuning parameter
#' @param MC Number of bootstrap replications
#' @param dh_step Step caching for Jacobian computation (default 10)
#' @param plotting Logical, whether to produce diagnostic CUSUM plot
#' @return List containing test results
cusum_sedcd_multivariate <- function(X, tau, k = NULL, lag = NULL, block = NULL,
                                    gamma = 0.05, c = 2.985, MC = 500,
                                    dh_step = 10, plotting = FALSE) {
  X <- as.matrix(X)
  N <- nrow(X)
  d <- ncol(X)
  p <- 1 + 3 * d + d^2
  
  if (is.null(k)) k <- max(10, floor(N^0.5))
  if (is.null(lag)) lag <- max(1, floor(0.1 * (log(N))^2))
  if (is.null(block)) block <- max(1, floor(0.1 * (log(N))^2))
  cutoff <- k + lag + block
  
  # 1. Pilot estimation
  teta_mat <- matrix(NA_real_, nrow = N, ncol = p)
  for (t in k:N) {
    win_data <- X[(t - k + 1):t, , drop = FALSE]
    teta_mat[t, ] <- estimate_local_theta_multi(win_data, tau = tau, gamma = gamma, c = c)
  }
  
  # 2. Linearized cumulative process
  Lin_sedcd <- numeric(N)
  DH_inv_row1 <- NULL
  for (t in cutoff:N) {
    t_lag <- t - lag
    teta_lag <- teta_mat[t_lag, ]
    
    if (is.null(DH_inv_row1) || (t - cutoff) %% dh_step == 0) {
      ws <- max(1, t_lag - k + 1)
      X_win <- X[ws:t_lag, , drop = FALSE]
      mean_H_fn <- function(param) colMeans(compute_H_multi(X_win, param, tau = tau, gamma = gamma, c = c))
      DH <- jacobian(mean_H_fn, teta_lag, method = "simple")
      DH_inv_row1 <- ginv(DH)[1, ]
    }
    
    H_t <- compute_H_multi(X[t, , drop = FALSE], teta_lag, tau = tau, gamma = gamma, c = c)
    Lin_sedcd[t] <- teta_lag[1] - sum(DH_inv_row1 * H_t)
  }
  Mn <- cumsum(Lin_sedcd) / N
  
  # 3. Long-run variance process
  q_sq <- numeric(N)
  DH_inv_row1 <- NULL
  for (t in cutoff:N) {
    param_idx <- t - lag - block
    teta_lag <- teta_mat[param_idx, ]
    
    if (is.null(DH_inv_row1) || (t - cutoff) %% dh_step == 0) {
      ws <- max(1, param_idx - k + 1)
      X_win <- X[ws:param_idx, , drop = FALSE]
      mean_H_fn <- function(param) colMeans(compute_H_multi(X_win, param, tau = tau, gamma = gamma, c = c))
      DH <- jacobian(mean_H_fn, teta_lag, method = "simple")
      DH_inv_row1 <- ginv(DH)[1, ]
    }
    
    H_block <- colSums(compute_H_multi(X[(t - block + 1):t, , drop = FALSE], teta_lag, tau = tau, gamma = gamma, c = c))
    q_val <- -sum(DH_inv_row1 * H_block)
    q_sq[t] <- (q_val^2) / block
  }
  Qn <- cumsum(q_sq) / N
  qn <- pmax(1e-8, diff(c(0, Qn)))
  
  # 4. CUSUM Bridge Statistic
  valid_range <- (cutoff + 1):N
  denom_bridge <- N - cutoff
  bridge <- Mn[valid_range] - ((1:length(valid_range)) / denom_bridge) * Mn[N]
  Tu_valid <- sqrt(N) * bridge
  Tu <- c(rep(0.0, cutoff), Tu_valid)
  
  test_stat <- max(abs(Tu_valid))
  max_idx <- cutoff + which.max(abs(Tu_valid))
  
  # 5. Multiplier Gaussian Bootstrap
  Z_mc <- numeric(MC)
  for (m in 1:MC) {
    xi <- rnorm(N)
    BM <- cumsum(xi * sqrt(qn))
    bridge_mc <- (BM[valid_range] - BM[cutoff]) - ((1:length(valid_range)) / denom_bridge) * (BM[N] - BM[cutoff])
    Z_mc[m] <- max(abs(bridge_mc))
  }
  
  p_val <- mean(Z_mc >= test_stat)
  q90 <- as.numeric(quantile(Z_mc, 0.90))
  q95 <- as.numeric(quantile(Z_mc, 0.95))
  
  if (plotting) {
    u_grid <- (1:N) / N
    plot(u_grid, Tu, type = "l", col = "blue", lwd = 1.5,
         xlab = "Rescaled Time u = t/n", ylab = "T_n(u)",
         main = sprintf("Multivariate CUSUM SEDCD Test (d = %d, p-val = %.3f)", d, p_val))
    abline(h = c(q95, -q95), col = "red", lty = 2)
    abline(v = max_idx / N, col = "darkgreen", lty = 3)
  }
  
  return(list(
    test_stat = test_stat,
    p_value   = p_val,
    q95       = q95,
    q90       = q90,
    max_index = max_idx,
    tau_idx   = max_idx,
    break_u   = max_idx / N,
    tau_u     = max_idx / N,
    Tu        = Tu,
    Qn        = Qn,
    teta_mat  = teta_mat,
    sedcd     = teta_mat[, 1],
    reject_95 = (p_val < 0.05),
    reject_90 = (p_val < 0.10)
  ))
}

