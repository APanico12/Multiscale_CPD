# =============================================================================
# Robust Recursive M-Estimation and CUSUM Change-Point Testing
#
# Implements:
#   - Tukey biweight loss, weight, and derivatives
#   - Recursive robust regression estimator (IRLS or Newton-Raphson)
#   - Linearized (bias-corrected) recursive estimator
#   - Robust variance estimation
#   - CUSUM-type test for change in regression coefficients, with Monte Carlo critical values
# =============================================================================
# -----------------------------------------------------------------------------
# 1. Loss, Score, Weight and Derivative Functions
# -----------------------------------------------------------------------------
library(robustbase)
library(MASS)

# Welsh (Welsch) functions
Welsh.rho <- function(u, c) {
  1.0 - exp(-(u / c)^2)
}
Welsh.psi <- function(u, c) {
  u * exp(-(u / c)^2)
}
Welsh.psi.prime <- function(u, c) {
  (1.0 - 2.0 * (u / c)^2) * exp(-(u / c)^2)
}
Welsh.weight <- function(u, c) {
  exp(-(u / c)^2)
}

# Tukey biweight functions
tukey_weight <- function(r, c) {
  ifelse(abs(r) <= c, (1.0 - (r / c)^2)^2, 0.0)
}
tukey_loss_derivative <- function(r, c) {
  ifelse(abs(r) <= c, r * (1.0 - (r / c)^2)^2, 0.0)
}
tukey_loss_2nd_derivative <- function(r, c) {
  ifelse(abs(r) <= c, (1.0 - (r / c)^2) * (1.0 - 5.0 * (r / c)^2), 0.0)
}

# -----------------------------------------------------------------------------
# 2. Single-window M-estimator of regression coefficients (IRWLS)
# -----------------------------------------------------------------------------
#' Solve robust regression M-estimate on a rolling window using IRWLS
#' initialized with Least Median of Squares (LMS)
#' @param Y_win numeric vector (window response, length k)
#' @param X_win numeric matrix (window regressors, shape k x d)
#' @param tol convergence tolerance
#' @param max_iter maximum IRWLS iterations
#' @param loss loss type: "Welsh" or "Tukey"
#' @param c tuning constant (default: 2.985 for Welsh, 4.685 for Tukey)

solve_theta <- function(Y_win, X_win, tol = 1e-6, max_iter = 50, 
                        loss = "Welsh", c = NULL, init_beta = NULL) {
  X_win <- as.matrix(X_win)
  Y_win <- as.vector(Y_win)
  k <- length(Y_win)
  d <- ncol(X_win)
  
  if (loss %in% c("L2", "l2", "OLS", "ols")) {
    beta_ols <- tryCatch({
      as.vector(solve(crossprod(X_win) + 1e-4 * diag(d), crossprod(X_win, Y_win)))
    }, error = function(e) rep(0, d))
    return(beta_ols)
  }
  
  if (loss == "Welsh") {
    if (is.null(c)) c <- 2.985  # 95% Gaussian efficiency for Welsh
    weight_fn <- Welsh.weight
  } else if (loss == "Tukey") {
    if (is.null(c)) c <- 4.685  # 95% Gaussian efficiency for Tukey
    weight_fn <- tukey_weight
  } else {
    stop("Unknown loss function: ", loss)
  }

  # Initial estimate: Use warm-start init_beta if provided, else LMS with fallback to regularized OLS
  if (!is.null(init_beta) && !any(is.na(init_beta))) {
    beta_curr <- init_beta
  } else {
    beta_curr <- tryCatch({
      suppressWarnings(as.vector(coef(MASS::lmsreg(X_win, Y_win, intercept = FALSE))))
    }, error = function(e) NULL)
    
    if (is.null(beta_curr) || any(is.na(beta_curr))) {
      beta_curr <- tryCatch({
        as.vector(solve(crossprod(X_win) + 1e-4 * diag(d), crossprod(X_win, Y_win)))
      }, error = function(e) rep(0, d))
    }
  }

  # IRWLS loop: w(u) = psi(u) / u >= 0
  for (iter in 1:max_iter) {
    r_curr <- as.vector(Y_win - X_win %*% beta_curr)
    sigma_hat <- median(abs(r_curr)) / 0.6745 # MADN
    if (is.na(sigma_hat) || sigma_hat < 1e-5) sigma_hat <- 1.0
    
    # Non-negative weights w(u)
    u_curr <- r_curr / sigma_hat
    w_curr <- weight_fn(u_curr, c = c)
    
    # Weighted least squares update without allocating k x k diagonal matrix
    XtWX <- crossprod(X_win, w_curr * X_win)
    XtWY <- crossprod(X_win, w_curr * Y_win)
    
    beta_new <- tryCatch({
      as.vector(solve(XtWX + 1e-6 * diag(d), XtWY))
    }, error = function(e) beta_curr)
    
    if (any(is.na(beta_new))) break
    
    if (max(abs(beta_new - beta_curr)) < tol * max(sigma_hat, 1e-4)) {
      beta_curr <- beta_new
      break
    }
    beta_curr <- beta_new
  }

  return(as.vector(beta_curr))
}

#' Rolling window estimator for time-varying regression beta(t/n)
#' @param Y numeric vector of responses (length N)
#' @param X numeric matrix of regressors (N x d)
#' @param k rolling window size (rate exponent if < 1, or absolute integer)
#' @param c tuning constant (default 2.985 for Welsh)
#' @param loss loss type: "Welsh" or "Tukey"

get_theta <- function(Y, X, k = 0.65, c = 2.985, loss = "Welsh") {

  N <- length(Y)
  X <- as.matrix(X)
  d <- ncol(X)
  if (k < 1) k <- floor(N^k) else k <- floor(k)
  
  betahat <- matrix(0, nrow = N, ncol = d)
  
  # Rolling window of length k using warm-started IRWLS (LMS cold-start on first window)
  betahat[k, ] <- solve_theta(Y[1:k], X[1:k, , drop = FALSE], c = c, loss = loss, init_beta = NULL)
  
  if (k < N) {
    for (t in (k + 1):N) {
      win_idx <- (t - k + 1):t
      betahat[t, ] <- solve_theta(Y[win_idx], X[win_idx, , drop = FALSE], 
                                  c = c, loss = loss, init_beta = betahat[t - 1, ])
    }
  }
  
  # Pad initial warm-up window with first available estimate
  if (k > 1) {
    for (t in 1:(k - 1)) {
      betahat[t, ] <- betahat[k, ]
    }
  }
  
  return(betahat)
}

# -----------------------------------------------------------------------------
# 3b. Cross-validation for optimal rolling bandwidth k
# -----------------------------------------------------------------------------

#' Cross-Validation for Optimal Rolling Bandwidth k
#' Implements Equation (861) from paper.tex / SEDCD:
#' Lambda(k) = sum_{t=k}^{N - L_n} || H(D_{t+L_n}, \hat{beta}_t(k)) ||^2
#'
#' @param Y       numeric vector of responses (length N)
#' @param X       numeric matrix of regressors (N x d)
#' @param k_grid  candidate bandwidth vector (integers, or rate exponents < 1)
#' @param lag     decoupling lag L_n (default: ceiling(N^0.1))
#' @param loss    loss function ("Welsh" or "Tukey")
#' @param c       loss tuning constant
#' @return list(k_opt, k_opt_rate, lag, cv_losses, k_grid)
cv_optimal_bandwidth_regression <- function(Y, X, k_grid = NULL, lag = NULL, 
                                            loss = "Welsh", c = NULL) {
  Y <- as.numeric(Y)
  N <- length(Y)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  d <- ncol(X)
  
  is_l2 <- loss %in% c("L2", "l2", "OLS", "ols")
  if (is_l2) {
    c <- 0
    psi_fn <- function(r, c = NULL) r
  } else if (is.null(c)) {
    c <- if (loss == "Welsh") 2.985 else 4.685
    psi_fn <- if (loss == "Welsh") Welsh.psi else tukey_loss_derivative
  } else {
    psi_fn <- if (loss == "Welsh") Welsh.psi else tukey_loss_derivative
  }
  
  # Decoupling lag L_n = N^0.1 as specified in the paper
  if (is.null(lag)) {
    lag <- max(1, ceiling(N^0.1))
  } else if (lag < 1) {
    lag <- max(1, ceiling(N^lag))
  }
  
  # Default search grid: [N^0.35, N^0.65] with 5 candidate evaluation points
  if (is.null(k_grid)) {
    k_min  <- max(d + 2, floor(N^0.35))
    k_max  <- min(floor(N - lag - 2), floor(N^0.65))
    k_grid <- unique(pmax(d + 2, floor(seq(k_min, k_max, length.out = 5))))
  } else {
    k_grid <- sapply(k_grid, function(kv) if (kv < 1) floor(N^kv) else floor(kv))
    k_grid <- unique(pmax(d + 2, k_grid))
  }
  
  cv_losses <- numeric(length(k_grid))
  names(cv_losses) <- paste0("k=", k_grid)
  
  # Robust residual scale for score standardization
  pilot_ols <- tryCatch(as.vector(solve(crossprod(X) + 1e-4 * diag(d), crossprod(X, Y))), 
                        error = function(e) rep(0, d))
  global_sig <- median(abs(Y - as.vector(X %*% pilot_ols))) / 0.6745
  if (is.na(global_sig) || global_sig < 1e-5) global_sig <- 1.0
  
  for (i in seq_along(k_grid)) {
    k_val <- k_grid[i]
    betahat_k <- get_theta(Y, X, k = k_val, c = c, loss = loss)
    
    t_start <- k_val
    t_end   <- N - lag
    
    if (t_start <= t_end) {
      t_eval  <- t_start:t_end
      fut_idx <- t_eval + lag
      
      X_fut <- X[fut_idx, , drop = FALSE]
      Y_fut <- Y[fut_idx]
      
      # Vectorized future prediction: r_fut = Y_{t+lag} - X_{t+lag} %*% beta_t
      pred_fut <- rowSums(X_fut * betahat_k[t_eval, , drop = FALSE])
      r_fut    <- Y_fut - pred_fut
      
      # Score weights psi(r_fut / sig)
      psi_vals <- if (is_l2) r_fut else psi_fn(r_fut / global_sig, c)
      
      # Squared norm ||H_{t+lag}||^2 = || X_{t+lag} * psi ||^2
      norm_X_sq <- rowSums(X_fut^2)
      H_sq      <- norm_X_sq * (psi_vals^2)
      
      cv_losses[i] <- mean(H_sq)
    } else {
      cv_losses[i] <- Inf
    }
  }
  
  best_idx  <- which.min(cv_losses)
  best_k    <- k_grid[best_idx]
  best_rate <- log(best_k) / log(N)
  
  return(list(
    k_opt      = best_k,
    k_opt_rate = round(best_rate, 4),
    lag        = lag,
    cv_losses  = round(cv_losses, 4),
    k_grid     = k_grid
  ))
}

# -----------------------------------------------------------------------------
# 4. Linearized (bias-corrected) recursive estimator for regression
# -----------------------------------------------------------------------------

#' Linearized recursive M-estimator for time-varying regression beta(t/n).
#'
#' At each t, starts from the pilot estimate betahat[t - lag, ] (which does NOT
#' contain (X[t, ], Y[t])) and applies a single Newton-Raphson correction using (X[t, ], Y[t]).
#'
#' @param Y       numeric vector of responses (length N)
#' @param X       numeric matrix of regressors (N x d)
#' @param betahat rolling M-estimates from get_theta() (N x d)
#' @param k       rolling window size (proportion or absolute) if< 1 N^k, else integer
#' @param c       Welsh tuning constant (default 2.985)
#' @param lag     lag between pilot estimate and future observation (default 1)
#' @param C_mat   contrast / selection matrix (s x d). Default is diag(d)
#' @return list(Lin.beta = N x d matrix, Mn = N x s cumulative process)

int.par.regression <- function(Y, X, betahat, k = 0.45, lag = NULL, c = NULL, loss = "Tukey", C_mat = NULL) {
  
  N <- length(Y)
  X <- as.matrix(X)
  d <- ncol(X)

  if (is.null(C_mat)) C_mat <- diag(d)
  if (is.vector(C_mat)) C_mat <- matrix(C_mat, nrow = 1, ncol = d)
  C_mat <- as.matrix(C_mat)
  s <- nrow(C_mat)
  
  is_l2 <- loss %in% c("L2", "l2", "OLS", "ols")
  if (is_l2) {
    psi_fn <- function(r, c = NULL) r
    psi_prime_fn <- function(r, c = NULL) rep(1, length(r))
  } else if (loss == "Welsh") {
    if (is.null(c)) c <- 2.985
    psi_fn <- Welsh.psi
    psi_prime_fn <- Welsh.psi.prime
  } else {
    if (is.null(c)) c <- 4.685
    psi_fn <- tukey_loss_derivative
    psi_prime_fn <- tukey_loss_2nd_derivative
  }
  
  if (k < 1) k_win <- max(d + 2, floor(N^k)) else k_win <- floor(k)
  if (is.null(lag)) lag <- max(1, floor(N^0.15))
  if (lag<1) lag <- floor(N^lag)

  Lin.beta <- matrix(0, nrow = N, ncol = d)
  ws <- k_win + lag
  
  
  ridge_penalty <- 1e-4 * diag(d)
  
  for (t in ws:N) {
    pilot_beta <- betahat[t - lag, ]
    start_idx  <- max(1, t - k_win - lag + 1)
    end_idx    <- t - lag
    
    X_hist <- X[start_idx:end_idx, , drop = FALSE]
    Y_hist <- Y[start_idx:end_idx]
    
    if (is_l2) {
      DH_mat   <- crossprod(X_hist) / nrow(X_hist) + ridge_penalty
      r_future <- Y[t] - sum(X[t, ] * pilot_beta)
      H_future <- X[t, ] * r_future
    } else {
      r_hist <- as.vector(Y_hist - X_hist %*% pilot_beta)
      sigma_hat <- median(abs(r_hist)) / 0.6745
      if (is.na(sigma_hat) || sigma_hat < 1e-5) sigma_hat <- 1.0
      
      # Local Jacobian DH_mat = 1/k * X_hist' W_psi' X_hist
      psi_prime <- psi_prime_fn(r_hist / sigma_hat, c) / sigma_hat
      DH_mat    <- crossprod(X_hist * sqrt(pmax(psi_prime, 0))) / length(r_hist) + ridge_penalty
      
      # Future observation score vector at t
      r_future <- Y[t] - sum(X[t, ] * pilot_beta)
      H_future <- X[t, ] * psi_fn(r_future, c) #there is no scale because it is a single observation
    }
    
    # One-step Newton correction: pilot + DH^-1 * H_future
    correction <- tryCatch({
      solve(DH_mat, H_future)
    }, error = function(e) rep(0, d))
    
    Lin.beta[t, ] <- pilot_beta + correction
  }
  
  # Projected cumulative process Mn(u) = 1/N * sum_{i=1}^t C_mat %*% beta_i
  C_Lin_beta <- matrix(t(C_mat %*% t(Lin.beta)), nrow = N, ncol = s)   # N x s
  Mn         <- matrix(apply(C_Lin_beta, 2, cumsum) / N, nrow = N, ncol = s) # N x s
  
  return(list(Lin.beta = Lin.beta, Mn = Mn))
}

# -----------------------------------------------------------------------------
# 5. Robust variance estimation for regression
# -----------------------------------------------------------------------------

#' Blockwise robust variance estimation for the recursive regression estimator.
#'
#' @param Y       numeric vector of responses (length N)
#' @param X       numeric matrix of regressors (N x d)
#' @param betahat rolling M-estimates from get_theta() (N x d)
#' @param k       rolling window size
#' @param block   block size for aggregating variance contributions
#' @param lag     decoupling lag
#' @param c       tuning constant (default 4.685 for Tukey, 2.985 for Welsh)
#' @param loss    loss function type: "Tukey" or "Welsh"
#' @param C_mat   contrast / selection matrix (s x d)
#' @return list(qn = N x s matrix of incremental variances, Q_est = cumulative variance)
var.est.regression <- function(Y, X, betahat, k = 0.45, block = NULL, lag = NULL, 
                               c = NULL, loss = "Tukey", C_mat = NULL) {
  N <- length(Y)
  X <- as.matrix(X)
  d <- ncol(X)
  if (is.null(C_mat)) C_mat <- diag(d)
  if (is.vector(C_mat)) C_mat <- matrix(C_mat, nrow = 1, ncol = d)
  C_mat <- as.matrix(C_mat)
  s <- nrow(C_mat)
  
  is_l2 <- loss %in% c("L2", "l2", "OLS", "ols")
  if (is_l2) {
    psi_fn <- function(r, c = NULL) r
    psi_prime_fn <- function(r, c = NULL) rep(1, length(r))
  } else if (loss == "Welsh") {
    if (is.null(c)) c <- 2.985
    psi_fn <- Welsh.psi
    psi_prime_fn <- Welsh.psi.prime
  } else {
    if (is.null(c)) c <- 4.685
    psi_fn <- tukey_loss_derivative
    psi_prime_fn <- tukey_loss_2nd_derivative
  }
  
  if (k < 1) k_win <- max(d + 2, floor(N^k)) else k_win <- floor(k)
  if (is.null(block)) block <- max(1, floor(N^0.25))
  if (is.null(lag))   lag   <- max(1, floor(N^0.15))
  cutoff <- k_win + block + lag
  
  ridge_penalty <- 1e-4 * diag(d)
  q_mat <- matrix(0, nrow = N, ncol = s)
  
  for (t in cutoff:N) {
    param_idx  <- t - lag - block
    pilot_beta <- betahat[param_idx, ]
    
    ws <- max(1, param_idx - k_win + 1)
    X_hist <- X[ws:param_idx, , drop = FALSE]
    Y_hist <- Y[ws:param_idx]
    
    if (is_l2) {
      DH_mat  <- crossprod(X_hist) / nrow(X_hist) + ridge_penalty
      blk_idx <- (t - block + 1):t
      X_blk   <- X[blk_idx, , drop = FALSE]
      Y_blk   <- Y[blk_idx]
      r_blk   <- as.vector(Y_blk - X_blk %*% pilot_beta)
      H_blk   <- colSums(X_blk * r_blk)
    } else {
      r_hist <- as.vector(Y_hist - X_hist %*% pilot_beta)
      sigma_hat <- median(abs(r_hist)) / 0.6745
      if (is.na(sigma_hat) || sigma_hat < 1e-5) sigma_hat <- 1.0
      
      psi_prime <- psi_prime_fn(r_hist / sigma_hat, c) / sigma_hat
      DH_mat    <- crossprod(X_hist * sqrt(pmax(psi_prime, 0))) / length(r_hist) + ridge_penalty
      
      # Block score sum
      blk_idx <- (t - block + 1):t
      X_blk   <- X[blk_idx, , drop = FALSE]
      Y_blk   <- Y[blk_idx]
      r_blk   <- as.vector(Y_blk - X_blk %*% pilot_beta)
      H_blk   <- colSums(X_blk * psi_fn(r_blk / sigma_hat, c))
    }
    
    DH_inv_H <- tryCatch({
      solve(DH_mat, H_blk)
    }, error = function(e) rep(0, d))
    
    q_mat[t, ] <- as.vector(-C_mat %*% DH_inv_H)
  }
  
  # Incremental variance contributions matching 1/sqrt(N) scaling of CUSUM
  qn    <- (q_mat^2) / (block * N) # N x s
  # Boundary padding: under local stationarity, pad warm-up window with initial variance estimate
  if (cutoff > 1) {
    for (j in 1:s) {
      qn[1:(cutoff - 1), j] <- qn[cutoff, j]
    }
  }
  Q_est <- matrix(apply(qn, 2, cumsum), nrow = N, ncol = s)
  
  return(list(qn = qn, Q_est = Q_est))
}

# -----------------------------------------------------------------------------
# 6. Proposed CUSUM test for change in regression coefficients
# -----------------------------------------------------------------------------

#' Proposed Linearized CUSUM test for time-varying regression coefficients.
#'
#' @param Y          numeric vector of responses (length N)
#' @param X          numeric matrix of regressors (N x d)
#' @param C_mat      contrast / selection matrix (s x d). Defaults to diag(d)
#' @param k          rolling window size (proportion e.g. 0.45 or absolute integer)
#' @param lag        forward decoupling lag (default: floor(N^0.15))
#' @param block      block size for robust variance estimate (default: floor(N^0.25))
#' @param cutoff     trimming threshold for boundary cutoffs
#' @param c          tuning constant (default: 4.685 for Tukey, 2.985 for Welsh)
#' @param loss       loss type: "Tukey" or "Welsh"
#' @param MC         number of Monte Carlo replications for multiplier bootstrap
#' @param linearized use the linearized (bias-corrected) recursive estimator if TRUE
#' @param plotting   produce a diagnostic plot if TRUE
#' @return list(test_stat, p_value, max_index, break_u, crit_value, reject_95)

CUSUM.regression <- function(Y, X, C_mat = NULL, betahat = NULL, k = 0.65, lag = NULL, block = NULL,
                             cutoff = NULL, c = NULL, loss = "Welsh", MC = 1000, B = 1000,
                             linearized = TRUE, use_cv = FALSE, plotting = FALSE) {
  if (!is.null(B)) MC <- B
  N <- length(Y)
  X <- as.matrix(X)
  d <- ncol(X)
  if (is.null(C_mat)) C_mat <- diag(d)
  if (is.vector(C_mat)) C_mat <- matrix(C_mat, nrow = 1, ncol = d)
  C_mat <- as.matrix(C_mat)
  s <- nrow(C_mat)
  
  is_l2 <- loss %in% c("L2", "l2", "OLS", "ols")
  if (is_l2) {
    c <- 0
  } else if (is.null(c)) {
    if (loss == "Welsh") c <- 2.985 else c <- 4.685
  }
  
  # Cross-validation for optimal bandwidth k if requested
  cv_info <- NULL
  if (isTRUE(use_cv) || (is.character(k) && tolower(k) == "cv")) {
    cv_info <- cv_optimal_bandwidth_regression(Y, X, lag = max(1, ceiling(N^0.1)), loss = loss, c = c)
    k <- cv_info$k_opt
  }
  
  if (k < 1) k_win <- max(d + 2, floor(N^k)) else k_win <- floor(k)
  if (is.null(block)) block <- max(1, floor(N^0.25))
  if (is.null(lag))   lag   <- max(1, floor(N^0.15))
  if (lag < 1) lag <- floor(N^lag)
  min_cutoff <- k_win + block + lag
  if (is.null(cutoff) || cutoff < min_cutoff) cutoff <- min_cutoff
  
  # Step 1: Rolling window pilot estimates
  if (is.null(betahat)) {
    betahat <- get_theta(Y, X, k = k, c = c, loss = loss)
  }
  
  # Step 2: Linearized bias-corrected estimator or plug-in
  if (linearized) {
    lin_res <- int.par.regression(Y, X, betahat, k = k, c = c, loss = loss, lag = lag, C_mat = C_mat)
    Mn      <- lin_res$Mn
  } else {
    C_beta  <- matrix(t(C_mat %*% t(betahat)), nrow = N, ncol = s)
    Mn      <- matrix(apply(C_beta, 2, cumsum) / N, nrow = N, ncol = s)
  }
  
  # Step 3: Robust variance estimation
  var_res <- var.est.regression(Y, X, betahat, k = k, block = block, lag = lag, c = c, loss = loss, C_mat = C_mat)
  qn      <- var_res$qn
  Qn      <- var_res$Q_est

  # Step 4: CUSUM test process Tu(u) = sqrt(N) * (Mn(u) - u * Mn(1))
  Tu <- matrix(0, nrow = N, ncol = s)
  for (j in 1:s) {
    Tu[(cutoff + 1):N, j] <- sqrt(N) * (Mn[(cutoff + 1):N, j] - (1:(N - cutoff)) / (N - cutoff) * Mn[N, j])
  }
  
  # Assign descriptive names to coordinates
  param_names <- if (!is.null(rownames(C_mat))) {
    rownames(C_mat)
  } else if (s == d && !is.null(colnames(X))) {
    colnames(X)
  } else if (s == 1 && sum(C_mat != 0) == 1 && any(abs(C_mat) == 1)) {
    j_idx <- which(C_mat[1, ] != 0)
    if (!is.null(colnames(X))) colnames(X)[j_idx] else paste0("beta_", j_idx)
  } else {
    paste0("Contrast_", 1:s)
  }
  colnames(Tu) <- param_names
  colnames(Mn) <- param_names
  colnames(qn) <- param_names
  colnames(Qn) <- param_names

  # Coordinate-wise statistics
  Z.vec           <- apply(abs(Tu), 2, max)
  max_index_coord <- apply(abs(Tu), 2, which.max)
  break_u_coord   <- max_index_coord / N
  
  # Overall joint maximum
  max_coord <- which.max(Z.vec)
  max_index <- max_index_coord[max_coord]
  break_u   <- break_u_coord[max_coord]
  Z         <- max(Z.vec)
  
  # Step 5: Multiplier bootstrap critical values
  Z_mc <- matrix(0, nrow = MC, ncol = s)

  for (j in 1:s) {
    sd_inc <- sqrt(pmax(qn[, j], 0)) # N x 1
    for (m in 1:MC) {
      BM <- cumsum(rnorm(N) * sd_inc)
      BB <- BM[(cutoff + 1):N] - (1:(N - cutoff)) / (N - cutoff) * BM[N]
      Z_mc[m, j] <- max(abs(BB))
    }
  }

  # Coordinate-wise critical values and p-values
  q10 <- apply(Z_mc, 2, quantile, probs = 0.90)
  q5  <- apply(Z_mc, 2, quantile, probs = 0.95)
  p_val_coord <- numeric(s)
  for (j in 1:s) {
    p_val_coord[j] <- mean(Z_mc[, j] >= Z.vec[j])
  }
  names(p_val_coord)     <- param_names
  names(q10)             <- param_names
  names(q5)              <- param_names
  names(Z.vec)           <- param_names
  names(max_index_coord) <- param_names
  names(break_u_coord)   <- param_names

  # Joint test statistic and critical values (supremum norm across coordinates)
  Z_joint_mc    <- apply(Z_mc, 1, max)
  crit_value_90 <- as.numeric(quantile(Z_joint_mc, 0.90))
  crit_value_95 <- as.numeric(quantile(Z_joint_mc, 0.95))
  p_value_joint <- mean(Z_joint_mc >= Z)
  
  # Parameter estimates and long-run variance reporting
  param_series   <- if (linearized) lin_res$Lin.beta else betahat
  C_param_series <- matrix(t(C_mat %*% t(param_series)), nrow = N, ncol = s)
  
  param_mean <- colMeans(C_param_series[(cutoff + 1):N, , drop = FALSE])
  param_pre  <- numeric(s)
  param_post <- numeric(s)
  for (j in 1:s) {
    t_star <- max_index_coord[j]
    param_pre[j]  <- mean(C_param_series[(cutoff + 1):t_star, j])
    param_post[j] <- if (t_star < N) mean(C_param_series[(t_star + 1):N, j]) else param_pre[j]
  }
  
  param_var <- Qn[N, ] # total integrated variance at u = 1
  names(param_mean) <- param_names
  names(param_pre)  <- param_names
  names(param_post) <- param_names
  names(param_var)  <- param_names
  
  param_summary <- data.frame(
    Parameter  = param_names,
    Mean_Value = round(param_mean, 4),
    Pre_Break  = round(param_pre, 4),
    Post_Break = round(param_post, 4),
    Variance   = round(param_var, 5),
    Test_Stat  = round(Z.vec, 4),
    Crit_90    = round(q10, 4),
    Crit_95    = round(q5, 4),
    p_value    = round(p_val_coord, 4),
    Break_t    = max_index_coord,
    Break_u    = round(break_u_coord, 3),
    Reject_95  = (Z.vec > q5),
    stringsAsFactors = FALSE
  )
  
  res <- list(
    stat             = Z.vec,
    test_stat        = Z,
    p_value          = p_value_joint,
    p_val_coord      = p_val_coord,
    crit_value       = crit_value_95,
    crit_value_95    = crit_value_95,
    crit_value_90    = crit_value_90,
    q5               = q5,
    q10              = q10,
    max_index        = max_index,
    break_u          = break_u,
    max_index_coord  = max_index_coord,
    break_u_coord    = break_u_coord,
    reject           = (p_value_joint < 0.05),
    reject_95        = (p_value_joint < 0.05),
    reject_coord     = (Z.vec > q5),
    param_names      = param_names,
    param_mean       = param_mean,
    param_pre        = param_pre,
    param_post       = param_post,
    param_var        = param_var,
    param_summary    = param_summary,
    Tu               = Tu,
    Mn               = Mn,
    qn               = qn,
    Qn               = Qn,
    beta             = betahat,
    Lin.beta         = if (linearized) lin_res$Lin.beta else betahat,
    C_mat            = C_mat,
    cutoff           = cutoff,
    N                = N,
    s                = s,
    k                = k_win,
    k_rate           = round(log(k_win) / log(N), 4),
    cv_info          = cv_info
  )
  class(res) <- "cusum_regression"
  
  if (plotting) {
    plot_cusum_regression(res)
  }
  
  return(res)
}


#' Plotting function for CUSUM regression test using zoo
#'
#' @param res object of class 'cusum_regression' returned by CUSUM.regression
#' @param ... additional arguments passed to plot.zoo
plot_cusum_regression <- function(res, ...) {
  if (!inherits(res, "cusum_regression")) {
    stop("Object must be of class 'cusum_regression'")
  }
  
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("Package 'zoo' is required for plot_cusum_regression")
  }
  
  N <- res$N
  s <- res$s
  u <- (1:N) / N
  Tu_zoo <- zoo::zoo(res$Tu, order.by = u)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  if (s == 1) {
    par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.8, 1.5))
  } else if (s == 2) {
    par(mfrow = c(2, 1), mar = c(4.2, 4.5, 3.2, 1.5))
  } else {
    par(mfrow = c(s, 1), mar = c(3.8, 4.5, 2.8, 1.5))
  }
  
  for (j in 1:s) {
    p_name   <- res$param_names[j]
    zj       <- res$stat[j]
    q5_j     <- res$q5[j]
    q10_j    <- res$q10[j]
    pj       <- res$p_val_coord[j]
    idx_j    <- res$max_index_coord[j]
    u_j      <- res$break_u_coord[j]
    val_j    <- res$param_mean[j]
    pre_j    <- res$param_pre[j]
    post_j   <- res$param_post[j]
    var_j    <- res$param_var[j]
    is_rej   <- res$reject_coord[j]
    
    y_limit  <- max(q5_j, zj) * 1.15
    
    status_str <- if (is_rej) "REJECT H0 (Break detected)" else "Accept H0 (No break)"
    main_title <- paste0("CUSUM Process for ", p_name, 
                         " | Z = ", round(zj, 3), 
                         " (p = ", format.pval(pj, digits = 3, eps = 0.001), ") - ", status_str)
    sub_title  <- paste0("Est. Param Value: Mean = ", round(val_j, 3), 
                         " [Pre: ", round(pre_j, 3), ", Post: ", round(post_j, 3), "]",
                         " | Long-Run Variance Q(1) = ", round(var_j, 5))
    
    # Plot using zoo
    zoo::plot.zoo(Tu_zoo[, j], xlab = "u = t / N", ylab = paste0("T(u) [", p_name, "]"),
                  type = "l", ylim = c(-y_limit, y_limit),
                  main = main_title, col = "black", lwd = 1.6, 
                  cex.main = 0.92, cex.lab = 0.9, cex.axis = 0.85, ...)
    
    mtext(sub_title, side = 3, line = 0.35, cex = 0.76, col = "darkblue")
    
    # Critical boundaries
    abline(h = q10_j, lty = 2, col = "gray40", lwd = 1.2)
    abline(h = -q10_j, lty = 2, col = "gray40", lwd = 1.2)
    abline(h = q5_j, lty = 3, col = "blue", lwd = 1.4)
    abline(h = -q5_j, lty = 3, col = "blue", lwd = 1.4)
    
    # Break indicator
    abline(v = u_j, col = "red", lty = 2, lwd = 1.5)
    text(u_j, 0, labels = paste0("t* = ", idx_j), 
         pos = ifelse(u_j > 0.8, 2, 4), col = "red", font = 2, cex = 0.9)
    
    # Legend with translucent background
    legend("topleft", 
           legend = c("T(u)", 
                      paste0("95% Crit (q5 = ", round(q5_j, 2), ")"), 
                      paste0("90% Crit (q10 = ", round(q10_j, 2), ")"),
                      paste0("Break t* = ", idx_j)),
           col = c("black", "blue", "gray40", "red"),
           lty = c(1, 3, 2, 2),
           lwd = c(1.6, 1.4, 1.2, 1.5),
           bg = rgb(1, 1, 1, 0.85), box.col = "gray80", cex = 0.75)
  }
  
  invisible(Tu_zoo)
}

plot.cusum_regression <- function(x, ...) {
  plot_cusum_regression(x, ...)
}

print.cusum_regression <- function(x, ...) {
  cat("\n================================================================================\n")
  cat("   Robust CUSUM Test for Structural Changes in Regression Coefficients\n")
  cat("================================================================================\n")
  cat(sprintf("Sample Size (N)     : %d\n", as.integer(x$N)))
  cat(sprintf("Warm-up Cutoff      : %d (u = %.3f)\n", as.integer(x$cutoff), x$cutoff / x$N))
  cat(sprintf("Monitored Dimensions: %d\n", as.integer(x$s)))
  cat(sprintf("Overall Test Stat Z : %.4f\n", x$test_stat))
  cat(sprintf("Joint 95%% Crit Value: %.4f\n", x$crit_value_95))
  cat(sprintf("Joint 90%% Crit Value: %.4f\n", x$crit_value_90))
  cat(sprintf("Joint p-value       : %s\n", format.pval(x$p_value, digits = 4, eps = 0.0001)))
  cat(sprintf("Decision (alpha=.05): %s\n", ifelse(x$reject, "REJECT H0 (Structural change detected)", "FAIL TO REJECT H0 (Stable coefficients)")))
  cat(sprintf("Estimated Break Loc : t* = %d (u* = %.3f)\n", as.integer(x$max_index), x$break_u))
  cat("--------------------------------------------------------------------------------\n")
  cat("Parameter Breakdown (Value, Variance, Test Statistic, and Decision):\n\n")
  print(x$param_summary, row.names = FALSE)
  cat("================================================================================\n\n")
  invisible(x)
}




# ==============================================================================
# 9. Example Usage: Structural Change Point in Beta (Single Parameter Monitoring)
# ==============================================================================

if (!exists("generate_tv_regression_dgp")) {
  if (file.exists("DGP.R")) {
    source("DGP.R")
  } else if (file.exists("Multiscale_CPD/Robust_regression/DGP.R")) {
    source("Multiscale_CPD/Robust_regression/DGP.R")
  }
}

# 1. Generate LMHC data with an abrupt change point in beta_1 at u = 0.5 (delta = 1.0)
#    and constant beta_2 = 2.0 under 5% Additive Outliers (AO)
d_cp <- generate_tv_regression_dgp(n = 1000, model_type = "LMHC", hp_scenario = "H1",
                                   delta = 1, contamination = "AO",
                                   K = 10, epsilon = 0.05, seed = 123)

# 2. Pilot robust M-estimates
beta <- get_theta(d_cp$Y, d_cp$X, k = 0.65, c = 2.985, loss = "Welsh")
#plot true beta and estimated beta
plot(1:nrow(d_cp$X), d_cp$beta_true[, 1], type = "l", col = "blue", lwd = 2, ylim = range(c(d_cp$beta_true[, 1], beta[, 1])),
     xlab = "Time Index", ylab = "Beta Coefficient", main = "True vs Estimated Beta_1")
lines(1:nrow(d_cp$X), beta[, 1], type = "l", col = "red", lwd = 2)
legend("topright", legend = c("True", "Estimated"), col = c("blue", "red"), lty = 1, lwd = 2)

# 3. Test specifically for the parameter with the change point: beta_1 (1 parameter, 1 image)
c_beta1 <- matrix(c(1, 0), nrow = 1)
rownames(c_beta1) <- "beta_1"

cusum_res <- CUSUM.regression(d_cp$Y, d_cp$X, C_mat = c_beta1, betahat = beta, 
                              k = 0.65, c = 2.985, loss = "Welsh", lag = 0.1, block = 0.25, 
                              MC = 500, linearized = TRUE, plotting = TRUE)

print(cusum_res)
