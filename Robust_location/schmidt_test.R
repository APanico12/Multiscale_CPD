# ==============================================================================
# File: schmidt_test.R
# Description: Implementation of the Gini-based test for changes in the trend
#              function of heteroscedastic time series from:
#              Sara Kristin Schmidt (2021), "Detecting Changes in the Trend
#              Function of Heteroscedastic Time Series", arXiv:2108.09206v1.
#
# Implements:
#   - Local block sample means (mu_hat) and sample variances (sigma_hat^2)
#   - Gini's mean difference U-statistic U(n)
#   - Centring term C_n
#   - Subsampling long-run variance estimator kappa_hat(n)
#   - Limit variance estimator psi_hat_n^2 via exact folded-normal conditional expectation
#   - Standard normal asymptotic test statistic T_n and p-value
# ==============================================================================

#' Schmidt (2021) Asymptotic Test for a Constant Mean under Heteroskedasticity
#'
#' @param x        Numeric vector, the observed time series.
#' @param s        Block length exponent in (0.5, 1). Default s = 0.7 (page 16).
#' @param q        Subsampling block exponent in (0, s). Default q = 0.4 (page 16).
#' @param c0       Number of surrounding blocks for subsampling centring. Default c0 = 10 (page 16).
#' @param M_psi    Number of Monte Carlo draws for estimating psi^2. Default 400.
#' @param seed     Optional random seed for the Monte Carlo draws of psi.
#'
#' @return A list containing:
#'         - statistic: Test statistic T_n
#'         - p_value:   Asymptotic p-value (1 - Phi(T_n))
#'         - U_n:       Gini mean difference statistic U(n)
#'         - C_n:       Estimated centring term
#'         - kappa_hat: Estimated long-run standard deviation of innovations
#'         - psi_hat:   Estimated standard deviation of limit distribution
#'         - l_n:       Block length ell_n
#'         - b_n:       Number of blocks
#'         - mu_hat:    Block sample means
#'         - sigma_hat: Block sample standard deviations
schmidt_test <- function(x, s = 0.7, q = 0.4, c0 = 10, M_psi = 400, seed = NULL) {
  n <- length(x)
  if (n < 20) stop("Sample size too small for Schmidt test.")

  # ---------------------------------------------------------------------------
  # 1. Blocking: ell_n = n^s, b_n = n^(1-s)
  # ---------------------------------------------------------------------------
  l_n <- max(2, floor(n^s))
  b_n <- floor(n / l_n)
  
  if (b_n < 2) {
    b_n <- 2
    l_n <- floor(n / 2)
  }

  # Effective sample length used for complete blocks
  n_eff <- b_n * l_n
  x_eff <- x[1:n_eff]

  # Reshape into matrix: l_n rows, b_n columns
  x_mat <- matrix(x_eff, nrow = l_n, ncol = b_n)

  # Block sample means mu_hat_j and sample variances sigma_hat_j^2
  mu_hat    <- colMeans(x_mat)
  diff_mat  <- sweep(x_mat, 2, mu_hat, "-")
  sigma_sq  <- colSums(diff_mat^2) / l_n
  sigma_sq  <- pmax(sigma_sq, 1e-6)
  sigma_hat <- sqrt(sigma_sq)

  # ---------------------------------------------------------------------------
  # 2. Gini's Mean Difference Statistic U(n) (eq. 3)
  #    U(n) = 2 / (b_n * (b_n - 1)) * sum_{j < k} |mu_hat_j - mu_hat_k|
  # ---------------------------------------------------------------------------
  diff_mu_pairs <- abs(outer(mu_hat, mu_hat, "-"))
  sum_mu_diff <- sum(diff_mu_pairs[upper.tri(diff_mu_pairs)])
  U_n <- (2 / (b_n * (b_n - 1))) * sum_mu_diff

  # ---------------------------------------------------------------------------
  # 3. Estimated Centring Term C_n (Proposition 2.9 & eq. on page 9)
  #    C_n = 2 / (b_n * (b_n - 1)) * sqrt(2/pi) * sum_{j < k} sqrt(sig_j^2 + sig_k^2)
  # ---------------------------------------------------------------------------
  sig_sq_sum_pairs <- outer(sigma_sq, sigma_sq, "+")
  sig_root_pairs   <- sqrt(sig_sq_sum_pairs)
  sum_sig_pairs    <- sum(sig_root_pairs[upper.tri(sig_root_pairs)])
  C_n <- (2 / (b_n * (b_n - 1))) * sqrt(2 / pi) * sum_sig_pairs

  # ---------------------------------------------------------------------------
  # 4. Long Run Variance Estimator kappa_hat(n) (Section 2.3.1, eq. 5 & 6)
  # ---------------------------------------------------------------------------
  l_sub <- max(1, floor(n^q))
  b_sub <- floor(n / l_sub)
  if (b_sub < 2) {
    b_sub <- 2
    l_sub <- floor(n / 2)
  }

  x_sub_mat <- matrix(x[1:(b_sub * l_sub)], nrow = l_sub, ncol = b_sub)
  sub_sums  <- colSums(x_sub_mat)

  sub_terms <- numeric(b_sub)
  for (j in 1:b_sub) {
    left_idx  <- if (j > 1) max(1, j - c0):(j - 1) else integer(0)
    right_idx <- if (j < b_sub) (j + 1):min(b_sub, j + c0) else integer(0)
    surr_idx  <- c(left_idx, right_idx)
    
    n_surr <- length(surr_idx)
    if (n_surr > 0) {
      surr_mean <- mean(sub_sums[surr_idx])
      prefactor <- sqrt(n_surr / (1 + n_surr))
      diff_val  <- sub_sums[j] - surr_mean
      sub_terms[j] <- prefactor * abs(diff_val)
    } else {
      sub_terms[j] <- abs(sub_sums[j])
    }
  }

  kappa_X_tilde <- (1 / b_sub) * sqrt(pi / 2) * (1 / sqrt(l_sub)) * sum(sub_terms)
  mean_sigma    <- mean(sigma_hat)
  kappa_hat     <- kappa_X_tilde / max(mean_sigma, 1e-4)
  if (is.na(kappa_hat) || kappa_hat < 1e-4) kappa_hat <- 1.0

  # ---------------------------------------------------------------------------
  # 5. Limit Variance Estimator psi_hat_n^2 (Section 2.3.3, eq. 7)
  #    Fast vectorized evaluation using folded normal conditional expectation
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  z_draws <- rnorm(M_psi)
  E_jk_uncond <- sig_root_pairs * sqrt(2 / pi)

  A_j_mat <- matrix(0, nrow = b_n, ncol = M_psi)
  ones_b  <- rep(1, b_n)
  ones_M  <- rep(1, M_psi)

  for (j in 1:b_n) {
    mu_j <- abs(sigma_hat[j] * z_draws)
    u_mat <- outer(1 / sigma_hat, mu_j, "*")
    
    cond_mat <- outer(sigma_hat, ones_M) * sqrt(2 / pi) * exp(-0.5 * u_mat^2) + 
                outer(ones_b, mu_j) * (1 - 2 * pnorm(-u_mat))
    
    diff_mat_k <- cond_mat - outer(E_jk_uncond[j, ], ones_M)
    diff_mat_k[j, ] <- 0
    A_j_mat[j, ] <- colSums(diff_mat_k) / (b_n - 1)
  }

  f_n_vals <- colMeans(A_j_mat^2)
  psi_sq   <- 4 * mean(f_n_vals)
  psi_hat  <- sqrt(max(psi_sq, 1e-4))

  # ---------------------------------------------------------------------------
  # 6. Test Statistic T_n & Asymptotic p-value (Theorem 2.12)
  # ---------------------------------------------------------------------------
  scaled_U  <- (sqrt(l_n) / kappa_hat) * U_n
  test_stat <- (sqrt(b_n) / psi_hat) * (scaled_U - C_n)
  p_val     <- pnorm(test_stat, lower.tail = FALSE)

  return(list(
    statistic = as.numeric(test_stat),
    p_value   = as.numeric(p_val),
    U_n       = as.numeric(U_n),
    C_n       = as.numeric(C_n),
    kappa_hat = as.numeric(kappa_hat),
    psi_hat   = as.numeric(psi_hat),
    l_n       = l_n,
    b_n       = b_n,
    mu_hat    = mu_hat,
    sigma_hat = sigma_hat
  ))
}

