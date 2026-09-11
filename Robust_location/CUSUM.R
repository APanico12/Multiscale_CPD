# =============================================================================
# Robust Recursive M-Estimation and CUSUM Change-Point Testing
#
# Implements:
#   - Tukey biweight loss, weight, and derivatives
#   - Recursive robust location estimator (IRLS or Newton-Raphson)
#   - Linearized (bias-corrected) recursive estimator
#   - Robust variance estimation
#   - CUSUM-type test for change in location, with Monte Carlo critical values
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Tukey biweight functions
# -----------------------------------------------------------------------------

#' Tukey biweight weight function w(r) = (1 - (r/c)^2)^2 for |r| <= c, else 0
tukey_weight <- function(r, c) {
  ifelse(abs(r) <= c,
         (1 - (r / c)^2)^2,
         0)
}

#' Tukey biweight psi function (first derivative of rho): psi(r) = r * w(r)
tukey_loss_derivative <- function(r, c) {
  ifelse(abs(r) <= c,
         r * (1 - (r / c)^2)^2,
         0)
}

#' Derivative of psi (used as the Jacobian / weight in Newton steps)
tukey_loss_2nd_derivative <- function(r, c) {
  ifelse(abs(r) <= c,
         (1 - (r / c)^2) * (1 - 5 * (r / c)^2),
         0)
}

# -----------------------------------------------------------------------------
# 2. Single-window M-estimator of location (IRLS or Newton-Raphson)
# -----------------------------------------------------------------------------

#' Solve for the robust location M-estimate on a fixed window of data
#'
#' @param x        numeric vector (one window of data)
#' @param c        Tukey tuning constant
#' @param method   "IRLS" or "Newton"
#' @param tol      convergence tolerance
#' @param max_iter maximum number of iterations
#' @return scalar robust location estimate
solve_teta <- function(x, c = 4.685, method = c("IRLS", "Newton"),
                        tol = 1e-6, max_iter = 1000) {

  method <- match.arg(method)

  # Robust starting value
  theta_current <- median(x)

  for (iter in 1:max_iter) {

    r <- x - theta_current

    if (method == "IRLS") {

      w <- tukey_weight(r, c)
      sum_w <- sum(w)

      if (sum_w == 0) {
        warning("IRLS: All weights are zero. Algorithm terminated early.")
        break
      }

      theta_new <- sum(w * x) / sum_w

    } else if (method == "Newton") {

      H_sum  <- sum(tukey_loss_derivative(r, c))
      DH_sum <- sum(tukey_loss_2nd_derivative(r, c))

      if (DH_sum == 0) {
        warning("Newton: Jacobian sum is zero. Algorithm terminated early.")
        break
      }

      theta_new <- theta_current + (H_sum / DH_sum)
    }

    if (abs(theta_new - theta_current) < tol) {
      return(theta_new)
    }

    theta_current <- theta_new
  }

  warning(paste("Maximum iterations reached without convergence using", method))
  return(theta_current)
}

# -----------------------------------------------------------------------------
# 3. Recursive (rolling-window) robust location estimator
# -----------------------------------------------------------------------------

#' Compute the recursive robust location estimate teta[t] over a rolling
#' window of size k ending at each time t.
#'
#' @param X        numeric data vector
#' @param k        window size; if k < 1, treated as a proportion of length(X)
#' @param c        Tukey tuning constant
#' @param method   "IRLS" or "Newton"
#' @param tol      convergence tolerance
#' @param max_iter maximum number of iterations
#' @return numeric vector of length(X), with teta[t] = 0 for t < k
get_teta <- function(X, k = 0.65, c = 4.685, method = c("IRLS", "Newton"),
                      tol = 1e-6, max_iter = 1000) {

  n_obs <- length(X)
  if (!is.numeric(X)) stop("X must be a numeric vector.")
  if (k < 1) k <- floor(n_obs^k)

  teta <- numeric(n_obs)

  for (t in k:n_obs) {
    teta[t] <- solve_teta(X[(t - k + 1):t], c = c, method = method,
                           tol = tol, max_iter = max_iter)
  }

  return(teta)
}

# -----------------------------------------------------------------------------
# 4. Linearized (bias-corrected) recursive estimator
# -----------------------------------------------------------------------------

#' Linearized one-step-ahead update of the recursive M-estimator.
#'
#' At each t, starts from the PILOT estimate teta[t - lag] (which does NOT
#' contain x[t]) and applies a single Newton correction using x[t]. This
#' avoids double-counting x[t]'s influence (it must not already be baked
#' into the estimate being corrected).
#'
#' @param x    numeric data vector
#' @param teta recursive M-estimate from get_teta()
#' @param k    window size (proportion or absolute)
#' @param c    Tukey tuning constant
#' @param lag  forward lag for the "future" observation used in the correction
#' @return numeric vector: cumulative mean of the linearized estimates

int.par.mean <- function(x, teta, k, c = 4.685, lag = 1) {

  n_obs <- length(x)
  if (!is.numeric(x)) stop("x must be a numeric vector.")

  if (k < 1) k <- floor(n_obs^k)
  if (k < 1) stop("Window size 'k' must be at least 1.")

  Lin.teta <- numeric(n_obs)

  # Loop starts at k + lag (full historical window) and runs to n_obs
  ws <- k + lag
  for (t in ws:n_obs) {

    # 1. Historical window for the Jacobian, using the PILOT estimate teta[t - lag]
    start_idx <- max(1, t - k - lag + 1)
    end_idx   <- t - lag
    r_hist <- x[start_idx:end_idx] - teta[t - lag]
    DH_rob_loc <- mean(tukey_loss_2nd_derivative(r_hist, c))

    if (DH_rob_loc == 0) {
      warning(paste("Jacobian is zero at t =", t, "- Skipping to prevent division by zero."))
      next
    }

    # 2. Newton correction using the FUTURE observation x[t], relative to the pilot
    future_obs <- x[t]
    r_future <- future_obs - teta[t - lag]
    H_future <- tukey_loss_derivative(r_future, c)

    # 3. Linearized estimator: pilot + one-step Newton correction
    Lin.teta[t] <- teta[t - lag] + (DH_rob_loc^-1 * H_future)
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
#' @param k     window size (proportion or absolute)
#' @param block block size for aggregating the variance contribution
#' @param lag   lag between parameter estimate and evaluation block
#' @return numeric vector: cumulative integrated variance process Q_n(u)
var.est.mean <- function(x, teta, k = 0.65, block = 1, lag = 1) {

  N <- length(x)
  if (k < 1) k <- floor(N^k)
  cutoff <- k + block + lag
  q_sq <- numeric(N)

  for (t in cutoff:N) {

    # Parameter estimate from the past, to ensure (near) independence
    param_idx <- t - lag - block
    teta_t_lag <- teta[param_idx]

    ws <- max(1, param_idx - k + 1)
    X_window <- x[ws:param_idx]
    DH <- mean(tukey_loss_2nd_derivative(X_window - teta_t_lag, c = 4.685))
    DH_inv <- 1 / DH

    H_val_block <- tukey_loss_derivative(x[(t - block + 1):t] - teta_t_lag, c = 4.685)
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
#' @param k           window size (proportion or absolute)
#' @param c           Tukey tuning constant
#' @param MC          number of Monte Carlo replications for critical values
#' @param linearized  use the linearized (bias-corrected) estimator if TRUE
#' @param plotting    produce a diagnostic plot if TRUE
#' @return list(p_value, test_stat, max_index)
CUSUM.mean <- function(x, teta = NULL, lag = 1, block = 1, cutoff = 1,
                        k = 0.65, c = 4.685, MC = 1000,
                        linearized = TRUE, plotting = FALSE) {

  N <- length(x)

  Check.cutoff <- N^k + lag + block
  if (cutoff < Check.cutoff) cutoff <- Check.cutoff

  # If teta is not provided, compute it (main parameter estimation step)
  if (is.null(teta)) teta <- get_teta(x, k = k, c = c)

  N <- length(teta)

  if (linearized) {
    Mn <- int.par.mean(x = x, teta = teta, k = k, c = c, lag = lag)
  } else {
    Mn <- cumsum(teta) / N
  }

  Qn <- var.est.mean(x = x, teta = teta, k = k, block = block, lag = lag)
  qn <- diff(c(0, Qn))

  Tu <- sqrt(N) * (Mn[(cutoff + 1):N] -
                     (1:(N - cutoff)) / (N - cutoff) * Mn[length(Mn)])
  Tu <- c(rep(0, cutoff), Tu)   # length now matches N (was cutoff + 1 before)

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

# =============================================================================
# Example usage
# =============================================================================
# set.seed(42)
# X <- rnorm(1000, mean = 2)
# X <- c(X, rnorm(1000, mean = 2))
#
# # Add random (heavy) outliers
# outlier_indices <- rbinom(length(X), 1, 0.05) == 1
# X[outlier_indices] <- X[outlier_indices] * 10
#
# plot(X, type = "l", main = "Simulated Data with Outliers", ylab = "X", xlab = "Time")
#
# teta <- get_teta(X, k = 0.65, c = 4.685, max_iter = 100)
# lin_teta_example <- int.par.mean(x = X, teta = teta, k = 0.65, lag = 1)
#
# plot(cumsum(teta) / length(teta), type = "l",
#      main = "Standard Estimator", ylab = "teta", xlab = "Time")
# plot(lin_teta_example, type = "l",
#      main = "Linearized Estimator", ylab = "Lin.teta", xlab = "Time")
#
# loc.NW <- ksmooth(1:length(X), X, kernel = "box", bandwidth = 100)$y
# plot(loc.NW, type = "l", main = "Nadaraya-Watson Estimator", ylab = "NW Estimate", xlab = "Time")
#
# b <- 10
# var_est_example <- var.est.mean(x = X, teta = teta, k = 0.65, block = b, lag = 1)
# plot(var_est_example, type = "l", main = "Variance Estimate", ylab = "Q_n(u)", xlab = "Time")
#
# result <- CUSUM.mean(x = X, teta = teta, lag = 1, block = b, cutoff = 1,
#                       k = 0.65, c = 4.685, MC = 1000, plotting = TRUE)
# p_value   <- result$p_value
# test_stat <- result$test_stat
# max_index <- result$max_index