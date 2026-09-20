# Generate data for change in shift in location
# according to    X_{t,n} = \mu(u) + \sum_{j=1}^p a_j X_{t-j,n} + \sigma(u) e_t + \sum_{k=1}^q b_k \sigma(u_k) e_{t-k}
# ==============================================================================
# File: DGP.R
# Description: Data Generating Process (DGP) for Time-Varying Robust Location:
#              X_{t,n} = mu(u) + sum_{j=1}^p a_j X_{t-j,n} + sigma(u) e_t 
#                               + sum_{k=1}^q b_k sigma(u_k) e_{t-k}
#              Supports Contamination: Clean, AO, IO, and RO (Replacement Outliers).
#              Supports Heteroskedasticity:
#                (i)   Smooth trending: sigma(u) = exp(u)
#                (ii)  Cyclic:          sigma(u) = 1 + sin(2*pi*u)
#                (iii) Abrupt break:    sigma(u) = 1.0 * I(u <= 0.5) + 2.0 * I(u > 0.5)
# ==============================================================================

#' Generate time series from an ARMA(p,q) model with a time-varying intercept,
#' time-varying scale (heteroskedasticity), and optional contamination.
#'
#' @param n Integer, the length of the time series.
#' @param ar_coeffs Numeric vector, the AR(p) coefficients (a_j).
#' @param ma_coeffs Numeric vector, the MA(q) coefficients (b_k).
#' @param mu_scenario String, the scenario for the time-varying intercept μ(u).
#'        One of "H0", "H1", "H2".
#' @param k Numeric, the magnitude of the shift for scenarios H1 and H2.
#' @param var_scenario String, the heteroskedasticity scenario for scale σ(u).
#'        One of:
#'        - "i" (or "trending"): Smooth trending variance σ(u) = exp(u / 2)
#'        - "i" (or "trending"): Smooth trending variance σ(u) = exp(u)
#'        - "ii" (or "cyclic"): Cyclic variance σ(u) = 1 + sin(2 * pi * u)
#'        - "iii" (or "break"): Abrupt variance break σ(u) = 0.5 * I(u <= 0.5) + 1.0 * I(u > 0.5)
#'        - "iii" (or "break"): Abrupt variance break σ(u) = 1.0 * I(u <= 0.5) + 2.0 * I(u > 0.5)
#'        - "constant" (or "none"): Homoskedastic σ(u) = 1.0
#' @param innov_dist String, the distribution of the innovations (ϵ_t).
#'        One of "gaussian" or "t3".
#' @param contamination_scenario String, the contamination scenario.
#'        One of "clean" / "none", "AO" (Additive Outliers), or "IO" (Innovation Outliers).
#' @param epsilon Numeric, the probability of outlier occurrence (Bernoulli parameter).
#' @param gamma Numeric, for "IO" and "AO", the fixed outlier magnitude multiplier.
#'        One of "clean" / "none", "AO" (Additive Outliers), "IO" (Innovation Outliers),
#'        or "RO" (Replacement Outliers via Markov Chain).
#' @param epsilon Numeric, the probability of outlier occurrence.
#' @param gamma Numeric, for "IO", "AO", and "RO", the fixed outlier magnitude multiplier.
#' @param kappa Numeric, Markov chain recovery probability for "RO" (default: 0.50).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list containing:
#'         - Xt: The generated time series (numeric vector).
#'         - Xt_clean: The uncontaminated latent series (numeric vector).
#'         - m_u: The true local expected value m(u) (numeric vector).
#'         - mu_u: The time-varying intercept μ(u) (numeric vector).
#'         - sigma_u: The time-varying scale σ(u) (numeric vector).
#'         - u: Rescaled time vector (1:n)/n.
#'         - outlier_flags: Binary 0/1 vector indicating contaminated time points.
ARMA_mu <- function(n, ar_coeffs = NULL, ma_coeffs = NULL, mu_scenario = "H0", k = 0.5,
                    var_scenario = "i", innov_dist = "gaussian", contamination_scenario = "clean",
                    epsilon = 0.05, gamma = 10) {
                    epsilon = 0.05, gamma = 10, kappa = 0.50, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # 1. Generate the time-varying intercept mu(u) and scale sigma(u)
  u <- (1:n) / n
  mu_u <- numeric(n)
  if (mu_scenario == "H1") {
    mu_u[u > 0.5] <- k
  } else if (mu_scenario == "H2") {
    mu_u <- k * u
  }

  # Scale function sigma(u) under the 3 heteroskedasticity scenarios
  # Scale function sigma(u) under the heteroskedasticity scenarios
  var_scen <- tolower(as.character(var_scenario))
  if (var_scen %in% c("i", "trending", "smooth", "exp")) {
    sigma_u <- exp(u / 2)
    sigma_u <- exp(u)                      # Updated to exp(u)
  } else if (var_scen %in% c("ii", "cyclic", "sin")) {
    sigma_u <- 1 + sin(2 * pi * u)
  } else if (var_scen %in% c("iii", "break", "abrupt")) {
    sigma_u <- ifelse(u <= 0.5, 0.5, 1.0)
    sigma_u <- ifelse(u <= 0.5, 1.0, 2.0)  # Jump of magnitude 1.0 (from 1.0 to 2.0)
  } else if (var_scen %in% c("none", "constant", "homoskedastic")) {
    sigma_u <- rep(1.0, n)
  } else {
    warning(sprintf("Unknown variance scenario '%s', defaulting to scenario (i).", var_scenario))
    sigma_u <- exp(u / 2)
    sigma_u <- exp(u)
  }

  # 2. Calculate the true local expected value m(u)
  sum_ar <- if (!is.null(ar_coeffs)) sum(ar_coeffs) else 0
  m_u <- mu_u / (1 - sum_ar)

  # 3. Generate white noise innovations
  if (innov_dist == "gaussian") {
  innov_norm <- tolower(as.character(innov_dist))
  if (innov_norm %in% c("gaussian", "normal")) {
    innovations <- rnorm(n)
  } else if (innov_dist == "t3") {
  } else if (innov_norm %in% c("t3", "t_3")) {
    innovations <- rt(n, df = 3)
  } else {
    stop("Unknown innovation distribution.")
  }

  # 4. Apply Innovation Outliers (IO) 
  if (contamination_scenario == "IO") {
  # 4. Apply Innovation Outliers (IO)
  contam_norm <- toupper(as.character(contamination_scenario))
  outlier_flags <- integer(n)

  if (contam_norm == "IO") {
    It <- rbinom(n, 1, epsilon)
    Xi <- rbinom(n, 1, 0.5) * 2 - 1  # Generates -1 or 1 with equal probability
    innovations <- innovations + It * gamma * Xi 
    Xi <- rbinom(n, 1, 0.5) * 2 - 1  # -1 or +1 with equal probability
    innovations <- innovations + It * gamma * Xi
    outlier_flags <- It
  }

  # Scale innovations by time-varying scale sigma(u)
  scaled_innovations <- sigma_u * innovations

  # 5. Simulate the ARMA(p,q) process
  p <- if (!is.null(ar_coeffs)) length(ar_coeffs) else 0
  q <- if (!is.null(ma_coeffs)) length(ma_coeffs) else 0
  Xt <- numeric(n)
  
  # Pad innovations and Xt for easier indexing
  padded_innovations <- c(rep(0, q), scaled_innovations)
  padded_Xt <- c(rep(0, p), Xt)

  for (t in 1:n) {
    ar_term <- if (p > 0) sum(ar_coeffs * padded_Xt[(t+p-1):(t)]) else 0
    ma_term <- if (q > 0) sum(ma_coeffs * padded_innovations[(t+q-1):(t)]) else 0
    ar_term <- if (p > 0) sum(ar_coeffs * padded_Xt[(t + p - 1):t]) else 0
    ma_term <- if (q > 0) sum(ma_coeffs * padded_innovations[(t + q - 1):t]) else 0
    
    # Equation: X_{t,n} = mu(u) + sum(a_j * X_{t-j}) + sigma(u)*e_t + sum(b_k * sigma(u_k)*e_{t-k})
    padded_Xt[t + p] <- mu_u[t] + ar_term + padded_innovations[t + q] + ma_term
  }
  Xt <- padded_Xt[(p+1):(n+p)]
  Xt_clean <- padded_Xt[(p + 1):(n + p)]
  Xt <- Xt_clean

  # 6. Apply Additive Outliers (AO) using gamma as the fixed magnitude
  if (contamination_scenario == "AO") {
  # 6. Apply Additive Outliers (AO)
  if (contam_norm == "AO") {
    It <- rbinom(n, 1, epsilon)
    Xi <- rbinom(n, 1, 0.5) * 2 - 1  # Generates -1 or 1 with equal probability
    Xt <- Xt + It * gamma * Xi 
    Xi <- rbinom(n, 1, 0.5) * 2 - 1  # -1 or +1 with equal probability
    Xt <- Xt_clean + It * gamma * Xi
    outlier_flags <- It
  }

  return(list(Xt = Xt, m_u = m_u, mu_u = mu_u, sigma_u = sigma_u))
  # 7. Apply Replacement Outliers (RO) via 2-state Markov Chain
  if (contam_norm == "RO") {
    # Stationary transition probabilities:
    # P(0 -> 1) = p_01, P(1 -> 0) = p_10 = kappa
    # Unconditional probability: P(S_t = 1) = epsilon
    p_10 <- kappa
    p_01 <- min(max((kappa * epsilon) / (1.0 - epsilon), 0.0), 1.0)
    
    S <- integer(n)
    S[1] <- if (runif(1) < epsilon) 1L else 0L
    
    for (t in 2:n) {
      rand_val <- runif(1)
      if (S[t - 1] == 0L) {
        S[t] <- if (rand_val < p_01) 1L else 0L
      } else {
        S[t] <- if (rand_val < p_10) 0L else 1L
      }
    }
    
    Xi <- rbinom(n, 1, 0.5) * 2 - 1
    Xt_star <- Xt_clean + gamma * Xi
    Xt <- (1 - S) * Xt_clean + S * Xt_star
    outlier_flags <- S
  }

  return(list(
    Xt            = Xt,
    Xt_clean      = Xt_clean,
    m_u           = m_u,
    mu_u          = mu_u,
    sigma_u       = sigma_u,
    u             = u,
    outlier_flags = outlier_flags
  ))
}


# ==============================================================================
# Visualization Functions (Publication 2x2 Matrix Style)
# ==============================================================================

#' Plot a single contaminated time series panel
plot_contaminated_location_series <- function(dgp_res, title_text, show_legend = FALSE,
                                              cex_axis = 1.4, cex_lab = 1.5, cex_main = 1.8) {
  y_lims <- extendrange(dgp_res$Xt, f = 0.10)
  
  plot(NA, NA,
       xlim = c(0, 1.0),
       ylim = y_lims,
       xlab = expression(bold(paste("Rescaled Time ", italic(u == t/n)))),
       ylab = expression(bold(paste("Observed ", italic(X[t])))),
       cex.lab = cex_lab,
       font.lab = 2,
       col.lab = "black",
       xaxt = "n", yaxt = "n",
       bty = "n")
  
  # Interior light dashed gridlines matching tick levels
  grid(col = "gray88", lty = 2, lwd = 0.9)
  
  # Observed series path
  lines(dgp_res$u, dgp_res$Xt, col = "gray25", lwd = 1.3)
  
  # Highlight outlier points
  out_idx <- which(dgp_res$outlier_flags == 1)
  if (length(out_idx) > 0) {
    points(dgp_res$u[out_idx], dgp_res$Xt[out_idx],
           col = "#c62828", pch = 1, lwd = 2.0, cex = 1.3)
  }
  
  # Prominent axis ticks
  x_ticks <- seq(0, 1.0, by = 0.25)
  axis(1, at = x_ticks, labels = c("0", "0.25", "0.5", "0.75", "1"),
       cex.axis = cex_axis, font.axis = 1, lwd = 0, lwd.ticks = 2.0)
  axis(2, cex.axis = cex_axis, font.axis = 1, lwd = 0, lwd.ticks = 2.0, las = 1)
  
  # Sharp outer bounding frame
  box(which = "plot", lty = "solid", lwd = 2.0, col = "black")
  
  # Prominent bold title matching paper style
  title(main = title_text, col.main = "black", font.main = 2, cex.main = cex_main, line = 1.1)
  
  # Legend if requested
  if (show_legend) {
    legend("topleft",
           legend  = c("Series", "Outlier"),
           col     = c("gray25", "#c62828"),
           lty     = c(1, NA),
           lwd     = c(1.3, NA),
           pch     = c(NA, 1),
           pt.lwd  = c(NA, 2.0),
           pt.cex  = c(NA, 1.3),
           bty     = "o",
           box.col = "gray80",
           box.lwd = 1.0,
           bg      = "white",
           cex     = 1.2,
           inset   = c(0.02, 0.03))
  }
}

#' Publication figure demonstration: 2x2 matrix of contaminated location series
demo_location_dgp <- function(save_files = TRUE,
                              output_prefix = "X_location",
                              n = 300,
                              seed = 123,
                              cex_axis = 1.4,
                              cex_lab = 1.5,
                              cex_main = 1.8) {
  
  ar_params <- c(0.2, -0.1)
  ma_params <- c(0.2)
  
  # 1. Clean
  d_clean <- ARMA_mu(n = n, ar_coeffs = ar_params, ma_coeffs = ma_params,
                     mu_scenario = "H0", var_scenario = "i", contamination_scenario = "clean",
                     seed = seed)
  
  # 2. Additive Outliers (AO)
  d_ao <- ARMA_mu(n = n, ar_coeffs = ar_params, ma_coeffs = ma_params,
                  mu_scenario = "H0", var_scenario = "i", contamination_scenario = "AO",
                  gamma = 10, epsilon = 0.06, seed = seed)
  
  # 3. Innovation Outliers (IO)
  d_io <- ARMA_mu(n = n, ar_coeffs = ar_params, ma_coeffs = ma_params,
                  mu_scenario = "H0", var_scenario = "i", contamination_scenario = "IO",
                  gamma = 10, epsilon = 0.06, seed = seed)
  
  # 4. Replacement Outliers (RO)
  d_ro <- ARMA_mu(n = n, ar_coeffs = ar_params, ma_coeffs = ma_params,
                  mu_scenario = "H0", var_scenario = "i", contamination_scenario = "RO",
                  gamma = 10, epsilon = 0.06, kappa = 0.40, seed = seed)
  
  draw_panels <- function() {
    par(mfrow = c(2, 2),
        mar   = c(4.5, 4.8, 2.8, 1.0),
        mgp   = c(3.0, 1.0, 0),
        tcl   = -0.5)
    
    plot_contaminated_location_series(d_clean, "Clean (Reference)", show_legend = FALSE,
                                     cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
    plot_contaminated_location_series(d_ao, "Additive Outliers (AO)", show_legend = TRUE,
                                     cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
    plot_contaminated_location_series(d_io, "Innovation Outliers (IO)", show_legend = FALSE,
                                     cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
    plot_contaminated_location_series(d_ro, "Replacement Outliers (RO)", show_legend = FALSE,
                                     cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
  }
  
  if (save_files) {
    # 1. High-resolution PNG (wide rectangular 2x2 layout)
    png_file <- paste0(output_prefix, ".png")
    png(png_file, width = 16, height = 7.5, units = "in", res = 300)
    draw_panels()
    dev.off()
    cat(sprintf("[SUCCESS] Saved PNG: %s\n", png_file))
    
    # 2. Vector PDF for LaTeX inclusion
    pdf_file <- paste0(output_prefix, ".pdf")
    pdf(pdf_file, width = 16, height = 7.5)
    draw_panels()
    dev.off()
    cat(sprintf("[SUCCESS] Saved PDF: %s\n", pdf_file))
    
    # Copy to paper img directory if it exists
    paper_img_dir <- file.path("..", "paper", "img")
    if (dir.exists(paper_img_dir)) {
      file.copy(png_file, file.path(paper_img_dir, basename(png_file)), overwrite = TRUE)
      file.copy(pdf_file, file.path(paper_img_dir, basename(pdf_file)), overwrite = TRUE)
      cat(sprintf("[SUCCESS] Copied %s to %s\n", basename(pdf_file), paper_img_dir))
    }
  } else {
    draw_panels()
  }
  
  invisible(list(clean = d_clean, ao = d_ao, io = d_io, ro = d_ro))
}
