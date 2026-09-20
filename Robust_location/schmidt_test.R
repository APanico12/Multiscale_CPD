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
# ==============================================================================
# 7. Example Usage with DGP (ARMA_mu from DGP.R)
# ==============================================================================

#' Demonstrate Schmidt (2021) Gini Change-Point Test on Simulated DGP Trajectories
#'
#' Evaluates the test across Null (H0), Abrupt Break (H1), Gradual Shift (H2),
#' and Outlier Contamination (AO, RO) scenarios generated by ARMA_mu from DGP.R.
#'
#' @param n                 Sample size. Default 500.
#' @param seed              Random seed for reproducibility. Default 123.
#' @param plot_diagnostics  Logical; if TRUE, generates a 2x2 diagnostic plot
#'                          showing series, block partitions, and local block means.
#'
#' @return Invisible list containing simulated DGP objects and test outputs.
#' @export
example_schmidt_dgp <- function(n = 500, seed = 123, plot_diagnostics = TRUE) {
  # Smart locator for DGP.R (supports working directory at project root or Robust_location)
  dgp_path <- if (file.exists("DGP.R")) {
    "DGP.R"
  } else if (file.exists(file.path("Robust_location", "DGP.R"))) {
    file.path("Robust_location", "DGP.R")
  } else {
    stop("Cannot find DGP.R. Please set your working directory to 'Robust_location' or the project root.")
  }
  source(dgp_path)

  cat("================================================================================\n")
  cat(" Demonstration: Schmidt (2021) Gini Test for Changes in Heteroscedastic Trend   \n")
  cat(sprintf(" Sample Size n = %d, Exponents s = 0.70, q = 0.40, M_psi = 400                  \n", n))
  cat("================================================================================\n\n")

  ar_p <- c(0.2, -0.1)
  ma_p <- c(0.2)

  scenarios <- list(
    list(id = "H0_Clean",  name = "H0: Constant Mean (Clean)", 
         mu = "H0", delta = 0.0, contam = "clean", eps = 0.0,  var = "i"),
    list(id = "H1_Abrupt", name = "H1: Abrupt Break (delta = 0.5, u = 0.5)", 
         mu = "H1", delta = 0.5, contam = "clean", eps = 0.0,  var = "i"),
    list(id = "H2_Gradual", name = "H2: Gradual Linear Trend (delta = 1.0)", 
         mu = "H2", delta = 0.75, contam = "clean", eps = 0.0,  var = "i"),
    list(id = "H1_AO",     name = "H1: Additive Outliers (AO, eps = 0.05)", 
         mu = "H1", delta = 0.75, contam = "AO",    eps = 0.05, var = "i"),
    list(id = "H2_RO",     name = "H2: Replacement Outliers (RO, eps = 0.05)", 
         mu = "H2", delta = 0.75, contam = "RO",    eps = 0.05, var = "i")
  )

  results    <- list()
  table_rows <- list()

  for (sc in scenarios) {
    # 1. Generate series from DGP
    dgp_data <- ARMA_mu(
      n                      = n,
      ar_coeffs              = ar_p,
      ma_coeffs              = ma_p,
      mu_scenario            = sc$mu,
      delta                  = sc$delta,
      var_scenario           = sc$var,
      contamination_scenario = sc$contam,
      epsilon                = sc$eps,
      gamma                  = 10,
      seed                   = seed
    )

    # 2. Run Schmidt (2021) test on observed series Xt
    res <- schmidt_test(dgp_data$Xt, s = 0.7, q = 0.4, c0 = 10, M_psi = 400, seed = seed)

    decision <- ifelse(res$p_value < 0.05, "REJECT H0", "Fail to reject")

    table_rows[[length(table_rows) + 1]] <- data.frame(
      Scenario    = sc$name,
      Statistic   = sprintf("%.3f", res$statistic),
      p_value     = sprintf("%.4f", res$p_value),
      Decision_05 = decision,
      U_n         = sprintf("%.4f", res$U_n),
      C_n         = sprintf("%.4f", res$C_n),
      kappa_hat   = sprintf("%.3f", res$kappa_hat),
      psi_hat     = sprintf("%.3f", res$psi_hat),
      b_n         = res$b_n,
      l_n         = res$l_n,
      stringsAsFactors = FALSE
    )

    results[[sc$id]] <- list(dgp = dgp_data, test = res, name = sc$name)
  }

  df_summary <- do.call(rbind, table_rows)
  print(df_summary, row.names = FALSE)
  cat("\nNote: Rejection region is one-sided upper tail (p < 0.05 implies significant break in trend).\n\n")

  # 3. Diagnostic Visualization (2x2 panels)
  if (plot_diagnostics) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    par(mfrow = c(2, 2), mar = c(4.5, 4.8, 2.8, 1.0), mgp = c(3.0, 1.0, 0), tcl = -0.5)

    # Dynamically pick up to 4 scenarios executed in results
    plot_ids <- names(results)[1:min(4, length(results))]

    for (i in seq_along(plot_ids)) {
      pid  <- plot_ids[i]
      item <- results[[pid]]
      d    <- item$dgp
      tst  <- item$test
      u    <- d$u
      x    <- d$Xt
      l_n  <- tst$l_n
      b_n  <- tst$b_n

      # Block partition boundaries in rescaled time u in [0, 1]
      block_ends   <- (1:b_n) * l_n / n
      block_starts <- c(0, block_ends[-b_n])

      plot_title <- if (!is.null(item$name)) item$name else pid

      plot(u, x, type = "l", col = "gray60", lwd = 1.0,
           xlim = c(0, 1.0), ylim = extendrange(x, f = 0.12),
           xlab = expression(bold(paste("Rescaled Time ", italic(u == t/n)))),
           ylab = expression(bold(paste("Observed ", italic(X[t])))),
           main = sprintf("%s\n(T_n = %.2f, p = %.3f)", plot_title, tst$statistic, tst$p_value),
           cex.lab = 1.2, font.lab = 2, cex.axis = 1.2, cex.main = 1.3, font.main = 2,
           las = 1, xaxt = "n", yaxt = "n", bty = "n")

      # Light dashed grid matching tick levels
      grid(col = "gray88", lty = 2, lwd = 0.9)

      # Custom prominent axes and bounding frame
      axis(1, at = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"),
           cex.axis = 1.2, font.axis = 1, lwd = 0, lwd.ticks = 1.8)
      axis(2, cex.axis = 1.2, font.axis = 1, lwd = 0, lwd.ticks = 1.8, las = 1)
      box(which = "plot", lty = "solid", lwd = 1.8, col = "black")

      # Vertical dashed block partition dividers
      abline(v = block_ends, col = "gray80", lty = 3, lwd = 1.1)

      # Horizontal bar for each Schmidt block sample mean
      for (j in 1:b_n) {
        lines(c(block_starts[j], block_ends[j]), rep(tst$mu_hat[j], 2), 
              col = "#1565c0", lwd = 2.6)
      }

      # True expectation line under hypothesis m(u)
      if (!is.null(d$m_u)) {
        lines(u, d$m_u, col = "#d32f2f", lwd = 1.8, lty = 2)
      }

      # Highlight outlier observations if present
      out_idx <- which(d$outlier_flags == 1)
      if (length(out_idx) > 0) {
        points(u[out_idx], x[out_idx], col = "#d32f2f", pch = 1, lwd = 1.8, cex = 1.2)
      }

      # Panel 1 master legend
      if (i == 1) {
        legend("topleft",
               legend = c("Observed Series", "Schmidt Block Means", "True Mean m(u)", "Outlier"),
               col    = c("gray60", "#1565c0", "#d32f2f", "#d32f2f"),
               lty    = c(1, 1, 2, NA),
               pch    = c(NA, NA, NA, 1),
               pt.lwd = c(NA, NA, NA, 1.8),
               pt.cex = c(NA, NA, NA, 1.2),
               lwd    = c(1.0, 2.6, 1.8, NA),
               bty    = "o", box.col = "gray80", bg = "white",
               cex    = 0.95, inset = c(0.02, 0.03))
      }
    }
  }

  invisible(results)
}

# Execute automatically if run as a script
if (sys.nframe() == 0) {
  example_schmidt_dgp()
}