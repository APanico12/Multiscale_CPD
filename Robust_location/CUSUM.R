# =============================================================================
# Robust Recursive M-Estimation and CUSUM Change-Point Testing (Location)
#
# Implements:
#   - Welsh (Welsch) and Tukey biweight loss, score, weight, and derivative functions
#   - Robust scale standardization via normalized MAD
#   - Single-window M-estimator of location (IRLS or Newton-Raphson)
#   - Recursive robust location estimator (rolling window with warm-starting)
#   - Out-of-sample forward Cross-Validation for optimal bandwidth k selection
#   - Linearized (bias-corrected) recursive estimator
#   - Robust variance estimation
#   - CUSUM-type test for change in location, with Monte Carlo critical values
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Loss, score, weight, and derivative functions
# -----------------------------------------------------------------------------

# Welsh (Welsch) functions: rho(u) = 1 - exp(-(u/c)^2)
Welsh.rho <- function(u, c = 2.985) {
  1.0 - exp(-(u / c)^2)
}

#' Welsh psi function (first derivative of rho): psi(u) = u * exp(-(u/c)^2)
Welsh.psi <- function(u, c = 2.985) {
  u * exp(-(u / c)^2)
}

#' Welsh psi prime (second derivative of rho / derivative of psi)
Welsh.psi.prime <- function(u, c = 2.985) {
  (1.0 - 2.0 * (u / c)^2) * exp(-(u / c)^2)
}

#' Welsh weight function: w(u) = psi(u) / u = exp(-(u/c)^2)
Welsh.weight <- function(u, c = 2.985) {
  exp(-(u / c)^2)
}

# Tukey biweight functions
#' Tukey biweight weight function w(r) = (1 - (r/c)^2)^2 for |r| <= c, else 0
tukey_weight <- function(r, c = 4.685) {
  ifelse(abs(r) <= c,
         (1.0 - (r / c)^2)^2,
         0.0)
}

#' Tukey biweight psi function (first derivative of rho): psi(r) = r * w(r)
tukey_loss_derivative <- function(r, c = 4.685) {
  ifelse(abs(r) <= c,
         r * (1.0 - (r / c)^2)^2,
         0.0)
}

#' Derivative of psi (used as the Jacobian / weight in Newton steps)
tukey_loss_2nd_derivative <- function(r, c = 4.685) {
  ifelse(abs(r) <= c,
         (1.0 - (r / c)^2) * (1.0 - 5.0 * (r / c)^2),
         0.0)
}

# -----------------------------------------------------------------------------
# 2. Single-window M-estimator of location (IRLS or Newton-Raphson)
# -----------------------------------------------------------------------------

#' Solve for the robust location M-estimate on a fixed window of data
#'
#' @param x          numeric vector (one window of data)
#' @param c          tuning constant (default: 2.985 for Welsh, 4.685 for Tukey)
#' @param loss       loss function: "Welsh", "Tukey", or "L2"
#' @param method     "IRLS" or "Newton"
#' @param tol        convergence tolerance
#' @param max_iter   maximum number of iterations
#' @param init_theta optional initial value for warm-starting
#' @return scalar robust location estimate
solve_teta <- function(x, c = NULL, loss = c("Welsh", "Tukey", "L2"),
                       method = c("IRLS", "Newton"),
                       tol = 1e-6, max_iter = 100, init_theta = NULL) {

  loss <- match.arg(loss)
  method <- match.arg(method)

  if (loss == "L2") return(mean(x))

  if (is.null(c)) {
    c <- if (loss == "Welsh") 2.985 else 4.685
  }

  weight_fn <- if (loss == "Welsh") Welsh.weight else tukey_weight
  psi_fn    <- if (loss == "Welsh") Welsh.psi else tukey_loss_derivative
  psi_p_fn  <- if (loss == "Welsh") Welsh.psi.prime else tukey_loss_2nd_derivative

  # Robust starting value (warm-start if provided)
  theta_current <- if (!is.null(init_theta) && is.finite(init_theta)) init_theta else median(x)

  for (iter in 1:max_iter) {

    r <- x - theta_current
    sigma_hat <- median(abs(r)) / 0.6745
    if (is.na(sigma_hat) || sigma_hat < 1e-5) sigma_hat <- 1.0

    u <- r / sigma_hat

    if (method == "IRLS") {

      w <- weight_fn(u, c)
      sum_w <- sum(w)

      if (sum_w == 0) {
        warning("IRLS: All weights are zero. Algorithm terminated early.")
        break
      }

      theta_new <- sum(w * x) / sum_w

    } else if (method == "Newton") {

      H_sum  <- sum(psi_fn(u, c))
      DH_sum <- sum(psi_p_fn(u, c)) / sigma_hat

      if (DH_sum == 0) {
        warning("Newton: Jacobian sum is zero. Algorithm terminated early.")
        break
      }

      theta_new <- theta_current + (H_sum / DH_sum)
    }

    if (abs(theta_new - theta_current) < tol * max(sigma_hat, 1e-4)) {
      return(theta_new)
    }

    theta_current <- theta_new
  }

  return(theta_current)
}

# -----------------------------------------------------------------------------
# 3. Recursive (rolling-window) robust location estimator
# -----------------------------------------------------------------------------

#' Compute the recursive robust location estimate teta[t] over a rolling
#' window of size k ending at each time t.
#'
#' @param X        numeric data vector
#' @param k        window size; if k < 1, treated as a proportion rate of length(X) (e.g. 0.45 or 0.65)
#' @param c        tuning constant
#' @param loss     "Welsh", "Tukey", or "L2"
#' @param method   "IRLS" or "Newton"
#' @param tol      convergence tolerance
#' @param max_iter maximum number of iterations
#' @return numeric vector of length(X), with teta[t] = 0 for t < k
get_teta <- function(X, k = 0.45, c = NULL, loss = c("Welsh", "Tukey", "L2"),
                     method = c("IRLS", "Newton"),
                     tol = 1e-6, max_iter = 100) {

  loss <- match.arg(loss)
  method <- match.arg(method)

  n_obs <- length(X)
  if (!is.numeric(X)) stop("X must be a numeric vector.")
  if (k < 1) k <- floor(n_obs^k) else k <- floor(k)
  k <- max(3, k)

  teta <- numeric(n_obs)

  # Cold start at window k
  teta[k] <- solve_teta(X[1:k], c = c, loss = loss, method = method,
                        tol = tol, max_iter = max_iter, init_theta = NULL)

  # Rolling warm-started updates
  if (k < n_obs) {
    for (t in (k + 1):n_obs) {
      teta[t] <- solve_teta(X[(t - k + 1):t], c = c, loss = loss, method = method,
                            tol = tol, max_iter = max_iter, init_theta = teta[t - 1])
    }
  }

  return(teta)
}

# -----------------------------------------------------------------------------
# 3b. Cross-validation for optimal rolling bandwidth k
# -----------------------------------------------------------------------------

#' Cross-Validation for Optimal Rolling Bandwidth k (Location)
#'
#' Evaluates candidate rolling bandwidths using out-of-sample forward evaluation
#' with decoupling lag L_n = max(1, floor(0.1 * (log(N))^2)):
#'   Lambda(k) = 1 / (N - k - L_n + 1) * sum_{t=k}^{N - L_n} psi^2((X_{t+L_n} - hat{theta}_t(k)) / sigma_hat)
#'
#' @param x       numeric vector (length N)
#' @param k_grid  candidate bandwidth vector (integers, or rate exponents < 1)
#' @param lag     decoupling lag L_n (default: floor(0.1 * (log(N))^2))
#' @param loss    loss function ("Welsh", "Tukey", or "L2")
#' @param c       loss tuning constant
#' @return list(k_opt, k_opt_rate, lag, cv_losses, k_grid)
cv_optimal_bandwidth_location <- function(x, k_grid = NULL, lag = NULL,
                                          loss = c("Welsh", "Tukey", "L2"), c = NULL) {
  loss <- match.arg(loss)
  x <- as.numeric(x)
  N <- length(x)

  is_l2 <- (loss == "L2")
  if (is_l2) {
    c <- 0
    psi_fn <- function(r, c = NULL) r
  } else if (is.null(c)) {
    c <- if (loss == "Welsh") 2.985 else 4.685
    psi_fn <- if (loss == "Welsh") Welsh.psi else tukey_loss_derivative
  } else {
    psi_fn <- if (loss == "Welsh") Welsh.psi else tukey_loss_derivative
  }

  # Decoupling lag L_n = max(1, floor(0.1 * (log(N))^2))
  if (is.null(lag)) {
    lag <- max(1, floor(0.1 * (log(N))^2))
  } else if (lag < 1) {
    lag <- max(1, ceiling(N^lag))
  }

  # Default search grid: [N^0.35, N^0.65] with 5 candidate evaluation points
  if (is.null(k_grid)) {
    k_min  <- max(5, floor(N^0.35))
    k_max  <- min(floor(N - lag - 2), floor(N^0.65))
    k_grid <- unique(pmax(5, floor(seq(k_min, k_max, length.out = 5))))
  } else {
    k_grid <- sapply(k_grid, function(kv) if (kv < 1) floor(N^kv) else floor(kv))
    k_grid <- unique(pmax(5, k_grid))
  }

  cv_losses <- numeric(length(k_grid))
  names(cv_losses) <- paste0("k=", k_grid)

  # Robust residual scale for score standardization
  global_sig <- median(abs(x - median(x))) / 0.6745
  if (is.na(global_sig) || global_sig < 1e-5) global_sig <- 1.0

  for (i in seq_along(k_grid)) {
    k_val <- k_grid[i]
    teta_k <- get_teta(x, k = k_val, c = c, loss = loss)

    t_start <- k_val
    t_end   <- N - lag

    if (t_start <= t_end) {
      t_eval  <- t_start:t_end
      fut_idx <- t_eval + lag

      r_fut    <- x[fut_idx] - teta_k[t_eval]
      psi_vals <- if (is_l2) r_fut else psi_fn(r_fut / global_sig, c)

      cv_losses[i] <- mean(psi_vals^2)
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
# 4. Linearized (bias-corrected) recursive estimator
# -----------------------------------------------------------------------------

#' Linearized one-step-ahead update of the recursive location M-estimator.
#'
#' At each t, starts from the PILOT estimate teta[t - lag] (which does NOT
#' contain x[t]) and applies a single Newton correction using x[t]. This
#' avoids double-counting x[t]'s influence and eliminates the first-order
#' smoothing/transition bias of the rolling estimator.
#'
#' @param x    numeric data vector
#' @param teta recursive M-estimate from get_teta()
#' @param k    window size (rate exponent if < 1, or absolute integer)
#' @param c    tuning constant
#' @param lag  forward lag for future observation (default: max(1, floor(0.1 * (log(N))^2)))
#' @param loss "Welsh", "Tukey", or "L2"
#' @return numeric vector: cumulative mean of the linearized estimates
int.par.mean <- function(x, teta, k, c = NULL, lag = NULL, loss = c("Welsh", "Tukey", "L2")) {

  loss <- match.arg(loss)
  n_obs <- length(x)
  if (!is.numeric(x)) stop("x must be a numeric vector.")

  if (k < 1) k <- floor(n_obs^k) else k <- floor(k)
  if (k < 1) stop("Window size 'k' must be at least 1.")

  if (is.null(lag)) {
    lag <- max(1, floor(0.1 * (log(n_obs))^2))
  } else if (lag < 1) {
    lag <- max(1, ceiling(n_obs^lag))
  }

  is_l2 <- (loss == "L2")
  if (is_l2) {
    c <- 0
    psi_fn <- function(r, c = NULL) r
    psi_p_fn <- function(r, c = NULL) rep(1, length(r))
  } else if (is.null(c)) {
    c <- if (loss == "Welsh") 2.985 else 4.685
    psi_fn <- if (loss == "Welsh") Welsh.psi else tukey_loss_derivative
    psi_p_fn <- if (loss == "Welsh") Welsh.psi.prime else tukey_loss_2nd_derivative
  } else {
    psi_fn <- if (loss == "Welsh") Welsh.psi else tukey_loss_derivative
    psi_p_fn <- if (loss == "Welsh") Welsh.psi.prime else tukey_loss_2nd_derivative
  }

  Lin.teta <- numeric(n_obs)

  # Loop starts at k + lag (full historical window) and runs to n_obs
  ws <- k + lag
  for (t in ws:n_obs) {
    pilot <- teta[t - lag]

    # 1. Historical window for the Jacobian, using the PILOT estimate teta[t - lag]
    start_idx <- max(1, t - k - lag + 1)
    end_idx   <- t - lag
    r_hist <- x[start_idx:end_idx] - pilot

    if (is_l2) {
      DH_rob_loc <- 1.0
      r_future <- x[t] - pilot
      H_future <- r_future
    } else {
      sigma_hat <- median(abs(r_hist)) / 0.6745
      if (is.na(sigma_hat) || sigma_hat < 1e-5) sigma_hat <- 1.0

      DH_rob_loc <- mean(psi_p_fn(r_hist / sigma_hat, c)) / sigma_hat
      if (is.na(DH_rob_loc) || DH_rob_loc < 1e-4) DH_rob_loc <- 1e-4

      # 2. Newton correction using the FUTURE observation x[t], relative to pilot
      r_future <- x[t] - pilot
      H_future <- psi_fn(r_future / sigma_hat, c)
    }

    # 3. Linearized estimator: pilot + one-step Newton correction
    Lin.teta[t] <- pilot + (H_future / DH_rob_loc)
  }

  return(cumsum(Lin.teta) / n_obs)
}

# -----------------------------------------------------------------------------
# 5. Robust variance estimation
# -----------------------------------------------------------------------------

#' Blockwise robust variance estimate for the recursive location estimator.
#'
#' @param x     numeric data vector
#' @param teta  recursive M-estimate from get_teta()
#' @param k     window size (rate exponent if < 1, or absolute integer)
#' @param block block size for aggregating the variance contribution
#' @param lag   lag between parameter estimate and evaluation block
#' @param c     tuning constant
#' @param loss  "Welsh", "Tukey", or "L2"
#' @return numeric vector: cumulative integrated variance process Q_n(u)
var.est.mean <- function(x, teta, k = 0.45, block = NULL, lag = NULL, c = NULL,
                         loss = c("Welsh", "Tukey", "L2")) {

  loss <- match.arg(loss)
  N <- length(x)
  if (k < 1) k <- floor(N^k) else k <- floor(k)

  if (is.null(lag)) lag <- max(1, floor(0.1 * (log(N))^2))
  if (lag < 1) lag <- max(1, ceiling(N^lag))

  if (is.null(block)) block <- max(1, floor(0.1 * (log(N))^2))
  if (block < 1) block <- max(1, ceiling(N^block))

  is_l2 <- (loss == "L2")
  if (is_l2) {
    c <- 0
    psi_fn <- function(r, c = NULL) r
    psi_p_fn <- function(r, c = NULL) rep(1, length(r))
  } else if (is.null(c)) {
    c <- if (loss == "Welsh") 2.985 else 4.685
    psi_fn <- if (loss == "Welsh") Welsh.psi else tukey_loss_derivative
    psi_p_fn <- if (loss == "Welsh") Welsh.psi.prime else tukey_loss_2nd_derivative
  } else {
    psi_fn <- if (loss == "Welsh") Welsh.psi else tukey_loss_derivative
    psi_p_fn <- if (loss == "Welsh") Welsh.psi.prime else tukey_loss_2nd_derivative
  }

  cutoff <- k + block + lag
  q_sq <- numeric(N)

  for (t in cutoff:N) {

    # Parameter estimate from the past, ensuring near-independence
    param_idx <- t - lag - block
    teta_t_lag <- teta[param_idx]

    ws <- max(1, param_idx - k + 1)
    X_window <- x[ws:param_idx]
    r_win <- X_window - teta_t_lag

    if (is_l2) {
      DH_inv <- 1.0
      H_val_block <- x[(t - block + 1):t] - teta_t_lag
    } else {
      sigma_hat <- median(abs(r_win)) / 0.6745
      if (is.na(sigma_hat) || sigma_hat < 1e-5) sigma_hat <- 1.0

      DH <- mean(psi_p_fn(r_win / sigma_hat, c)) / sigma_hat
      if (is.na(DH) || DH < 1e-4) DH <- 1e-4
      DH_inv <- 1 / DH

      r_block <- x[(t - block + 1):t] - teta_t_lag
      H_val_block <- psi_fn(r_block / sigma_hat, c)
    }

    q_par <- -DH_inv * sum(H_val_block)
    q_sq[t] <- q_par^2
  }

  q_scaled <- q_sq / block
  Q.est <- cumsum(q_scaled) / N

  return(Q.est)
}

# -----------------------------------------------------------------------------
# 6. CUSUM test for change in location
# -----------------------------------------------------------------------------

#' CUSUM-type test for a change in the robust location parameter.
#'
#' @param x           numeric data vector
#' @param teta        optional precomputed recursive M-estimate; computed if NULL
#' @param lag         forward lag used in the linearized estimator
#' @param block       block size for the variance estimate
#' @param cutoff      minimum starting index for the test statistic
#' @param k           window size (proportion rate or absolute integer)
#' @param c           tuning constant
#' @param loss        "Welsh", "Tukey", or "L2"
#' @param MC          number of Monte Carlo replications for critical values
#' @param linearized  use the linearized (bias-corrected) estimator if TRUE
#' @param plotting    produce a diagnostic plot if TRUE
#' @return list(p_value, test_stat, max_index)
CUSUM.mean <- function(x, teta = NULL, lag = NULL, block = NULL, cutoff = 1,
                       k = 0.45, c = NULL, loss = c("Welsh", "Tukey", "L2"),
                       MC = 1000, linearized = TRUE, plotting = FALSE) {

  loss <- match.arg(loss)
  N <- length(x)

  if (k < 1) k_win <- floor(N^k) else k_win <- floor(k)
  if (is.null(lag)) lag <- max(1, floor(0.1 * (log(N))^2))
  if (is.null(block)) block <- max(1, floor(0.1 * (log(N))^2))

  Check.cutoff <- k_win + lag + block
  if (cutoff < Check.cutoff) cutoff <- Check.cutoff

  # If teta is not provided, compute it (main parameter estimation step)
  if (is.null(teta)) teta <- get_teta(x, k = k, c = c, loss = loss)

  if (linearized) {
    Mn <- int.par.mean(x = x, teta = teta, k = k, c = c, lag = lag, loss = loss)
  } else {
    Mn <- cumsum(teta) / N
  }

  Qn <- var.est.mean(x = x, teta = teta, k = k, block = block, lag = lag, c = c, loss = loss)
  qn <- diff(c(0, Qn))

  Tu <- sqrt(N) * (Mn[(cutoff + 1):N] -
                     (1:(N - cutoff)) / (N - cutoff) * Mn[length(Mn)])
  Tu <- c(rep(0, cutoff), Tu)

  Z <- max(abs(Tu))
  max_index <- which.max(abs(Tu))

  # Monte Carlo critical values via a Brownian-bridge-type construction
  Z.mc <- numeric(MC)
  for (i in 1:MC) {
    BM <- cumsum(rnorm(N) * sqrt(qn))
    Z.mc[i] <- max(abs(BM[(cutoff + 1):N] -
                          (1:(N - cutoff)) / (N - cutoff) * BM[N]))
  }

  q10 <- quantile(Z.mc, 0.90)
  q5  <- quantile(Z.mc, 0.95)

  if (plotting) {
    plot((1:N) / N, Tu, xlab = "u", ylab = "T(u)", type = "l",
         ylim = c(-max(q5, Z) * 1.1, max(q5, Z) * 1.1))
    abline(h = q10, lty = 2)
    abline(h = -q10, lty = 2)
    abline(h = q5, lty = 3)
    abline(h = -q5, lty = 3)
    abline(v = max_index / N, col = "red", lty = 2)
    text(max_index / N, 0, labels = paste(max_index), pos = 3, col = "red")
  }

  res <- list(p_value = mean(Z.mc > Z), test_stat = Z, max_index = max_index)
  return(res)
}

# ######################
# #Example usage:
# ########################

# source("DGP.R")
# n=1000
# ar_params <- c(0.2, -0.1)
# ma_params <- c(0.2)
# hp_scenario <- "H1"
# var_scenario <- "i"

# seed <- 123


# d_ao <- ARMA_mu(n = n, ar_coeffs = ar_params, ma_coeffs = ma_params,
#                   mu_scenario = hp_scenario, var_scenario = var_scenario, contamination_scenario = "AO",
#                   gamma = 10, epsilon = 0.06, seed = seed)


# res <- CUSUM.mean(x = d_ao$Xt, k = 0.45, loss = "Welsh", c = 2.985, MC = 1000, linearized = TRUE, plotting = TRUE)