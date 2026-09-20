# ==============================================================================
# File: DGP.R
# Description: Data Generating Processes (DGPs) for Time-Varying Robust M-Regression:
#              - Model I:  Linear Model with Conditional Heteroscedasticity (LMCH)
#                          Y_{t,n} = X_{t,n}' beta(t/n) + Gamma(X_{t,n}) e_{t,n}
#                          Gamma(X_{t,n}) = sqrt(1 + (x1_{t,n})^2 + (x2_{t,n})^2)
#              - Model II: Linear Model with Unconditional Heteroscedasticity (LMUH)
#                          Y_{t,n} = X_{t,n}' beta(t/n) + sigma(t/n) e_{t,n}
#                          sigma(u) = 0.5 + |sin(2*pi*u)|
#              Supports Contamination: Clean, AO, IO, and RO (Replacement Outliers).
# ==============================================================================

# ==============================================================================
# 1. Parameter Specifications under Null (H0) and Alternative (H1/H2) Hypotheses
# ==============================================================================

# Standard bivariate regression parameters for Model I (LMCH) & Model II (LMUH)
# Target parameter: beta_0 = (1.0, 2.0)'
# Perturbation: v = (1.0, 0.0)'
beta_null <- function(u, ...) {
  c(1.0, 2.0)
}

beta_h1_abrupt <- function(u, delta = 0.8, ...) {
  c(1.0 + delta * (u > 0.5), 2.0)
}

beta_h1_gradual <- function(u, delta = 0.8, ...) {
  c(1.0 + delta * u, 2.0)
}

# Legacy trivariate parameters for backwards compatibility with LMAT
beta_null_lmat <- function(u, ...) {
  c(0.4, 1.0, 2.0)
}

beta_h1_abrupt_lmat <- function(u, delta = 0.8, ...) {
  c(0.4, 1.0 + delta * (u > 0.5), 2.0)
}

beta_h1_gradual_lmat <- function(u, delta = 0.8, ...) {
  c(0.4, 1.0 + delta * u, 2.0)
}


# ==============================================================================
# 2. Main Simulation Data Generator
# ==============================================================================

generate_tv_regression_dgp <- function(n = 500,
                                       model_type = "LMCH",
                                       error_dist = "normal",
                                       innov_dist = NULL,
                                       contamination = "Clean",
                                       epsilon = 0.05,
                                       K = 20,
                                       ro_multiplier = 4.0,    # Multiplier for beta effect during RO
                                       kappa = 0.50,           # Markov chain recovery probability
                                       hp_scenario = "H0",
                                       beta_fn = NULL,
                                       delta = 0.8,            # Shift magnitude for H1/H2
                                       add_baseline = FALSE,   # Adds 2*sin(2*pi*u) baseline
                                       burn_in = 150,
                                       seed = NULL,
                                       ...) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(innov_dist)) error_dist <- innov_dist
  
  # Normalize model_type and support aliases
  model_type <- toupper(model_type)
  if (model_type == "LMHC") model_type <- "LMCH"
  if (model_type %in% c("LMUC", "LUHC")) model_type <- "LMUH"
  
  if (!model_type %in% c("LMCH", "LMUH", "LMAT")) {
    stop(paste("Unknown model_type:", model_type, "- Expected 'LMCH' or 'LMUH'"))
  }
  
  # Normalize error distribution
  error_dist <- tolower(error_dist)
  if (error_dist == "gaussian") error_dist <- "normal"
  if (error_dist %in% c("t_3", "t-3")) error_dist <- "t3"
  if (!error_dist %in% c("normal", "t3")) {
    stop(paste("Unknown error_dist:", error_dist, "- Expected 'normal'/'gaussian' or 't3'/'t_3'"))
  }
  
  # Normalize contamination regime
  contamination <- toupper(contamination)
  if (!contamination %in% c("CLEAN", "NONE", "AO", "IO", "RO")) {
    stop(paste("Unknown contamination:", contamination, "- Expected 'Clean', 'AO', 'IO', or 'RO'"))
  }
  
  total_n <- n + burn_in
  u_eval  <- (1:n) / n
  lags    <- 0:burn_in
  
  # ---------------------------------------------------------------------------
  # A. Select Default Beta Functions based on hp_scenario
  # ---------------------------------------------------------------------------
  hp_scenario <- toupper(hp_scenario)
  if (is.null(beta_fn)) {
    if (model_type %in% c("LMCH", "LMUH")) {
      beta_fn <- switch(hp_scenario,
                        "H0" = beta_null,
                        "H1" = beta_h1_abrupt,
                        "H2" = beta_h1_gradual,
                        beta_null)
    } else if (model_type == "LMAT") {
      beta_fn <- switch(hp_scenario,
                        "H0" = beta_null_lmat,
                        "H1" = beta_h1_abrupt_lmat,
                        "H2" = beta_h1_gradual_lmat,
                        beta_null_lmat)
    }
  }
  
  d <- if (model_type %in% c("LMCH", "LMUH")) 2 else 3
  beta_mat <- matrix(0, nrow = n, ncol = d)
  for (t_idx in 1:n) {
    if (!is.null(delta)) {
      beta_mat[t_idx, ] <- beta_fn(u_eval[t_idx], delta = delta, ...)
    } else {
      beta_mat[t_idx, ] <- beta_fn(u_eval[t_idx], ...)
    }
  }
  
  # ---------------------------------------------------------------------------
  # B. Baseline Driving Shocks
  # ---------------------------------------------------------------------------
  eps1 <- rnorm(total_n, mean = 0, sd = 1)
  eps2 <- rnorm(total_n, mean = 0, sd = 1)
  
  if (error_dist == "normal") {
    eta <- rnorm(total_n, mean = 0, sd = 1)
  } else if (error_dist == "t3") {
    eta <- rt(total_n, df = 3) / sqrt(3.0)  # Standardized to unit variance
  }
  
  # ---------------------------------------------------------------------------
  # C. Scenario 2: Innovation Outliers (Pre-filter Contamination)
  # ---------------------------------------------------------------------------
  outlier_flags <- integer(n)
  if (contamination == "IO") {
    mag_io    <- if (is.null(K)) 20.0 else K
    I_eta_obs <- rbinom(n, size = 1, prob = epsilon)
    I_eta     <- c(integer(burn_in), I_eta_obs)
    xi_eta    <- rbinom(total_n, size = 1, prob = 0.5)
    sign_eta  <- ifelse(xi_eta == 1, -1.0, 1.0)
    eta       <- eta + I_eta * sign_eta * mag_io
    outlier_flags <- I_eta_obs
  }
  
  # ---------------------------------------------------------------------------
  # D. MA(infinity) Time-Varying Linear Filters
  # ---------------------------------------------------------------------------
  phi1_u <- 0.5 - 0.5 * u_eval
  phi2_u <- 0.25 + 0.5 * u_eval
  xi_u   <- 0.5 - (u_eval - 0.5)^2
  
  X <- matrix(0, nrow = n, ncol = 2)
  e <- numeric(n)
  
  for (i in 1:n) {
    idx <- i + burn_in
    w1  <- phi1_u[i]^lags
    w2  <- phi2_u[i]^lags
    we  <- xi_u[i]^lags
    
    window_idx <- idx - lags
    X[i, 1] <- sum(w1 * eps1[window_idx])
    X[i, 2] <- sum(w2 * eps2[window_idx])
    e[i]    <- sum(we * eta[window_idx])
  }
  
  # ---------------------------------------------------------------------------
  # E. Structural Signal Synthesis (Clean Path)
  # ---------------------------------------------------------------------------
  baseline_trend <- if (add_baseline) 2.0 * sin(2.0 * pi * u_eval) else numeric(n)
  
  if (model_type == "LMCH") {
    # Model I: Conditional Heteroscedasticity Gamma(X) = sqrt(1 + x1^2 + x2^2)
    gamma_X      <- sqrt(1.0 + X[, 1]^2 + X[, 2]^2)
    signal_clean <- rowSums(X * beta_mat)
    noise_clean  <- gamma_X * e
    
    Y_clean    <- baseline_trend + signal_clean + noise_clean
    regressors <- X
    
  } else if (model_type == "LMUH") {
    # Model II: Unconditional Heteroscedasticity sigma(u) = 0.5 + |sin(2*pi*u)|
    sigma_u      <- 0.5 + abs(sin(2.0 * pi * u_eval))
    signal_clean <- rowSums(X * beta_mat)
    noise_clean  <- sigma_u * e
    
    Y_clean    <- baseline_trend + signal_clean + noise_clean
    regressors <- X
    
  } else if (model_type == "LMAT") {
    # Legacy Model: Autoregressive with Unconditional Heteroscedasticity
    sigma_u <- 0.5 + abs(sin(2.0 * pi * u_eval))
    Y_clean <- numeric(n)
    W       <- matrix(0, nrow = n, ncol = 3)
    
    y_lag <- 0.0
    for (i in 1:n) {
      W[i, ]     <- c(y_lag, X[i, 1], X[i, 2])
      signal     <- baseline_trend[i] + sum(W[i, ] * beta_mat[i, ])
      noise      <- sigma_u[i] * e[i]
      Y_clean[i] <- signal + noise
      y_lag      <- Y_clean[i]
    }
    regressors  <- W
    noise_clean <- sigma_u * e
  }
  
  # ---------------------------------------------------------------------------
  # F. Post-Processing Contamination (AO & Amplified-Beta RO)
  # ---------------------------------------------------------------------------
  Y_obs <- Y_clean
  mag_K <- if (is.null(K)) 4.0 * sd(Y_clean) else K
  
  if (contamination == "AO") {
    I_ao    <- rbinom(n, size = 1, prob = epsilon)
    xi_ao   <- rbinom(n, size = 1, prob = 0.5)
    sign_ao <- ifelse(xi_ao == 1, -1.0, 1.0)
    Y_obs   <- Y_clean + I_ao * sign_ao * mag_K
    outlier_flags <- I_ao
    
  } else if (contamination == "RO") {
    # Two-state stationary Markov Chain S_t in {0, 1}
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
    
    # Replacement Rule: Amplify the beta effect during regime bursts S_t = 1
    if (model_type == "LMCH") {
      signal_amplified <- rowSums(X * (ro_multiplier * beta_mat))
      Y_star <- baseline_trend + signal_amplified + noise_clean
    } else if (model_type == "LMUH") {
      signal_amplified <- rowSums(X * (ro_multiplier * beta_mat))
      Y_star <- baseline_trend + signal_amplified + noise_clean
    } else if (model_type == "LMAT") {
      signal_amplified <- rowSums(W * (ro_multiplier * beta_mat))
      Y_star <- baseline_trend + signal_amplified + noise_clean
    }
    
    Y_obs <- (1 - S) * Y_clean + S * Y_star
    outlier_flags <- S
  }
  
  # Align lagged feedback in legacy LMAT if contaminated
  if (model_type == "LMAT" && !contamination %in% c("CLEAN", "NONE")) {
    regressors[, 1] <- c(0.0, Y_obs[1:(n - 1)])
  }
  
  if (model_type %in% c("LMCH", "LMUH")) {
    colnames(regressors) <- c("beta_1 (X1)", "beta_2 (X2)")
  } else {
    colnames(regressors) <- c("beta_AR (Y_lag)", "beta_1 (X1)", "beta_2 (X2)")
  }
  
  return(list(
    u             = u_eval,
    Y             = Y_obs,
    Y_clean       = Y_clean,
    regressors    = regressors,
    X             = regressors,
    beta_true     = beta_mat,
    outlier_flags = outlier_flags
  ))
}


# ==============================================================================
# 3. Execution & Visualization Example (Paper Style)
# ==============================================================================

plot_contaminated_series <- function(dgp_res, title_text, show_legend = FALSE,
                                     cex_axis = 2, cex_lab = 2.0, cex_main = 2.4) {
  y_lims <- extendrange(dgp_res$Y, f = 0.10)
  
  # Base empty plot matching power curve style
  plot(NA, NA,
       xlim = c(0, 1.0),
       ylim = y_lims,
       xlab = expression(bold(paste("Rescaled Time ", italic(u == t/n)))),
       ylab = expression(bold(paste("Observed ", italic(Y[t])))),
       cex.lab = cex_lab,
       font.lab = 2,
       col.lab = "black",
       xaxt = "n", yaxt = "n",
       bty = "n")
  
  # Interior light dashed gridlines matching tick levels (paper style)
  grid(col = "gray88", lty = 2, lwd = 0.9)
  
  # Latent / observed time series path
  lines(dgp_res$u, dgp_res$Y, col = "gray25", lwd = 1.4)
  
  # Highlight outlier points
  out_idx <- which(dgp_res$outlier_flags == 1)
  if (length(out_idx) > 0) {
    points(dgp_res$u[out_idx], dgp_res$Y[out_idx],
           col = "#c62828", pch = 1, lwd = 2.0, cex = 1.3)
  }
  
  # Prominent axis ticks matching plot_power_matrix (cex.axis = 2.0, lwd.ticks = 2.0)
  x_ticks <- seq(0, 1.0, by = 0.25)
  axis(1, at = x_ticks, labels = c("0", "0.25", "0.5", "0.75", "1"),
       cex.axis = cex_axis, font.axis = 1, lwd = 0, lwd.ticks = 2.0)
  axis(2, cex.axis = cex_axis, font.axis = 1, lwd = 0, lwd.ticks = 2.0, las = 1)
  
  # Sharp outer bounding frame (lwd = 2.0 matching power matrix figures)
  box(which = "plot", lty = "solid", lwd = 2.0, col = "black")
  
  # Prominent bold title matching paper style
  title(main = title_text, col.main = "black", font.main = 2, cex.main = cex_main, line = 1.1)
  
  # Legend if requested
  if (show_legend) {
    legend("topleft",
           legend  = c("Series", "Outlier"),
           col     = c("gray25", "#c62828"),
           lty     = c(1, NA),
           lwd     = c(1.4, NA),
           pch     = c(NA, 1),
           pt.lwd  = c(NA, 2.0),
           pt.cex  = c(NA, 1.3),
           bty     = "o",
           box.col = "gray80",
           box.lwd = 1.0,
           bg      = "white",
           cex     = 1.5,
           inset   = c(0.02, 0.03))
  }
}

# Visual verification demo: publication figure for Model I (LMCH)
# Default: 2x2 matrix of horizontal rectangular panels occupying ~1/3 of a LaTeX page
demo_dgp <- function(save_files = TRUE,
                     output_prefix = "Y_regression",
                     layout = "2x2",
                     n = 300,
                     seed = 123,
                     cex_axis = 2.0,
                     cex_lab = 2.0,
                     cex_main = 2.2) {
  
  # Generate 4 contamination scenarios exclusively for Model I (LMCH)
  # 1. Clean
  d_clean <- generate_tv_regression_dgp(n = n, model_type = "LMCH", contamination = "Clean", seed = seed)
  
  # 2. Additive Outliers (AO)
  d_ao <- generate_tv_regression_dgp(n = n, model_type = "LMCH", contamination = "AO",
                                     K = 20, epsilon = 0.06, seed = seed)
  
  # 3. Innovation Outliers (IO)
  d_io <- generate_tv_regression_dgp(n = n, model_type = "LMCH", contamination = "IO",
                                     K = 20, epsilon = 0.06, seed = seed)
  
  # 4. Replacement Outliers (RO)
  d_ro <- generate_tv_regression_dgp(n = n, model_type = "LMCH", contamination = "RO",
                                     ro_multiplier = 4.0, epsilon = 0.06, kappa = 0.40, seed = seed)
  
  draw_panels <- function() {
    if (layout == "2x2") {
      # 2x2 layout: Wide horizontal rectangular panels (X-axis longer), occupying ~1/3 of page height
      par(mfrow = c(2, 2),
          mar   = c(4.8, 5.4, 2.8, 1.0),
          mgp   = c(3.3, 1.1, 0),
          tcl   = -0.5)
      
      plot_contaminated_series(d_clean, "Clean (Reference)", show_legend = FALSE,
                               cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
      plot_contaminated_series(d_ao, "Additive Outliers (AO)", show_legend = TRUE,
                               cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
      plot_contaminated_series(d_io, "Innovation Outliers (IO)", show_legend = FALSE,
                               cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
      plot_contaminated_series(d_ro, "Replacement Outliers (RO)", show_legend = FALSE,
                               cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
    } else {
      # 1x4 layout: Single row across page
      par(mfrow = c(1, 4),
          mar   = c(4.8, 5.2, 3.5, 1.0),
          oma   = c(0.5, 0.5, 0.5, 0.5),
          mgp   = c(3.3, 1.2, 0),
          tcl   = -0.6)
      
      plot_contaminated_series(d_clean, "Clean (Reference)", show_legend = FALSE,
                               cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
      plot_contaminated_series(d_ao, "Additive Outliers (AO)", show_legend = TRUE,
                               cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
      plot_contaminated_series(d_io, "Innovation Outliers (IO)", show_legend = FALSE,
                               cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
      plot_contaminated_series(d_ro, "Replacement Outliers (RO)", show_legend = FALSE,
                               cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main)
    }
  }
  
  # Save publication figures if requested
  if (save_files) {
    pdf_file <- paste0(output_prefix, ".pdf")
    png_file <- paste0(output_prefix, ".png")
    
    if (layout == "2x2") {
      # Vector PDF & 300 DPI PNG sized for ~1/3 of a LaTeX page height across full linewidth
      pdf(pdf_file, width = 16.0, height = 7.6)
      draw_panels()
      dev.off()
      cat(sprintf("Saved vector PDF: %s\n", pdf_file))
      
      png(png_file, width = 4800, height = 2280, res = 300)
      draw_panels()
      dev.off()
      cat(sprintf("Saved high-res PNG: %s\n", png_file))
    } else {
      pdf(pdf_file, width = 16.0, height = 5.2)
      draw_panels()
      dev.off()
      cat(sprintf("Saved vector PDF: %s\n", pdf_file))
      
      png(png_file, width = 4800, height = 1560, res = 300)
      draw_panels()
      dev.off()
      cat(sprintf("Saved high-res PNG: %s\n", png_file))
    }
    
    # Copy to paper/img/ directory if present
    img_dirs <- c(
      file.path("..", "paper", "img"),
      file.path("Multiscale_CPD", "paper", "img"),
      file.path("paper", "img")
    )
    for (idir in img_dirs) {
      if (dir.exists(idir)) {
        target_png <- file.path(idir, "Y_regression.png")
        target_pdf <- file.path(idir, "Y_regression.pdf")
        file.copy(png_file, target_png, overwrite = TRUE)
        file.copy(pdf_file, target_pdf, overwrite = TRUE)
        cat(sprintf("Synced figure to paper directory: %s\n", target_png))
        break
      }
    }
  }
  
  # Also render in active graphic device if interactive
  if (interactive()) {
    draw_panels()
  }
}

# Run visual verification when executed directly from terminal
if (sys.nframe() == 0) {
  demo_dgp()
}
