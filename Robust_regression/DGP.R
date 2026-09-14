# ==============================================================================
# File: dgp_tv_regression.R
# Description: Comprehensive Data Generating Processes (DGPs) for Time-Varying
#              Robust M-Regression under Non-Stationarity, Heteroscedasticity,
#              Autoregressive Dynamics, and Various Contamination Scenarios.
# ==============================================================================

# ==============================================================================
# 1. Parameter Specifications under Null (H0) and Alternative (H1) Hypotheses
# ==============================================================================

# Model I: LMHC
beta_null_lmhc <- function(u, ...) {
  c(1.0, 2.0)
}

beta_h1_abrupt_lmhc <- function(u, delta = 0.8, ...) {
  c(1.0 + delta * (u > 0.5), 2.0)
}

beta_h1_gradual_lmhc <- function(u, delta = 0.8, ...) {
  c(1.0 + delta * u, 2.0)
}

# Model II: LMAT
beta_null_lmat <- function(u, ...) {
  c(0.4, 1.0, 0.5)
}

beta_h1_abrupt_lmat <- function(u, delta = 0.8, ...) {
  c(0.4, 1.0 + delta * (u > 0.5), 0.5)
}

beta_h1_gradual_lmat <- function(u, delta = 0.8, ...) {
  c(0.4, 1.0 + delta * u, 0.5)
}


# ==============================================================================
# 2. Main Simulation Data Generator
# ==============================================================================

generate_tv_regression_dgp <- function(n = 500,
                                       model_type = "LMHC",
                                       error_dist = "normal",
                                       innov_dist = NULL,
                                       contamination = "Clean",
                                       epsilon = 0.10,
                                       K = 20,
                                       ro_multiplier = 4.0,    # Multiplier for beta effect during RO
                                       kappa = 0.50,
                                       hp_scenario = "H0",
                                       beta_fn = NULL,
                                       delta = 0.8,            # Shift magnitude for H1/H2
                                       add_baseline = FALSE,   # Adds 2*sin(2*pi*u) baseline
                                       burn_in = 150,
                                       seed = NULL,
                                       ...) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(innov_dist)) error_dist <- innov_dist
  
  # Normalize inputs
  model_type    <- toupper(model_type)
  if (!model_type %in% c("LMHC", "LMAT")) {
    stop(paste("Unknown model_type:", model_type, "- Expected 'LMHC' or 'LMAT'"))
  }
  
  error_dist    <- tolower(error_dist)
  if (error_dist == "gaussian") error_dist <- "normal"
  if (error_dist %in% c("t_3", "t-3")) error_dist <- "t3"
  if (!error_dist %in% c("normal", "t3")) {
    stop(paste("Unknown error_dist:", error_dist, "- Expected 'normal'/'gaussian' or 't3'/'t_3'"))
  }
  
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
    if (model_type == "LMHC") {
      beta_fn <- switch(hp_scenario,
                        "H0" = beta_null_lmhc,
                        "H1" = beta_h1_abrupt_lmhc,
                        "H2" = beta_h1_gradual_lmhc,
                        beta_null_lmhc)
    } else {
      beta_fn <- switch(hp_scenario,
                        "H0" = beta_null_lmat,
                        "H1" = beta_h1_abrupt_lmat,
                        "H2" = beta_h1_gradual_lmat,
                        beta_null_lmat)
    }
  }
  
  d <- if (model_type == "LMHC") 2 else 3
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
  
  if (model_type == "LMHC") {
    gamma_X      <- sqrt(1.0 + X[, 1]^2 + X[, 2]^2)
    signal_clean <- rowSums(X * beta_mat)
    noise_clean  <- gamma_X * e
    
    Y_clean    <- baseline_trend + signal_clean + noise_clean
    regressors <- X
    
  } else if (model_type == "LMAT") {
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
    regressors <- W
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
    if (model_type == "LMHC") {
      signal_amplified <- rowSums(X * (ro_multiplier * beta_mat))
      Y_star <- baseline_trend + signal_amplified + noise_clean
    } else {
      signal_amplified <- rowSums(W * (ro_multiplier * beta_mat))
      Y_star <- baseline_trend + signal_amplified + (sigma_u * e)
    }
    
    Y_obs <- (1 - S) * Y_clean + S * Y_star
    outlier_flags <- S
  }
  
  # Align lagged feedback in LMAT if contaminated
  if (model_type == "LMAT" && !contamination %in% c("CLEAN", "NONE")) {
    regressors[, 1] <- c(0.0, Y_obs[1:(n - 1)])
  }
  
  if (model_type == "LMHC") {
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
# 3. Execution & Visualization Example
# ==============================================================================

plot_contaminated_series <- function(dgp_res, title_text) {
  y_lims <- extendrange(dgp_res$Y, f = 0.08)
  plot(dgp_res$u, dgp_res$Y, type = "l", col = "gray25", lwd = 1.2,
       main = title_text, xlab = "Rescaled Time (u = t/n)", ylab = "Observed Y",
       ylim = y_lims)
  
  out_idx <- which(dgp_res$outlier_flags == 1)
  if (length(out_idx) > 0) {
    points(dgp_res$u[out_idx], dgp_res$Y[out_idx],
           col = "firebrick", pch = 1, lwd = 1.8, cex = 1.3)
  }
  
  grid(col = "gray85")
}

# Run visual verification only if executed as the main script
if (sys.nframe() == 0) {
# Run visual verification demo
demo_dgp <- function() {
  par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.5, 1))
  
  # 1. Clean
  d_clean <- generate_tv_regression_dgp(n = 300, model_type = "LMHC", contamination = "Clean", seed = 123)
  plot_contaminated_series(d_clean, "1. Clean (Reference)")
  
  # 2. Additive Outliers (AO)
  d_ao <- generate_tv_regression_dgp(n = 300, model_type = "LMHC", contamination = "AO",
                                     K = 20, epsilon = 0.08, seed = 123)
  plot_contaminated_series(d_ao, "2. Additive Outliers (AO)")
  
  # 3. Innovation Outliers (IO)
  d_io <- generate_tv_regression_dgp(n = 300, model_type = "LMHC", contamination = "IO",
                                     K = 20, epsilon = 0.08, seed = 123)
  plot_contaminated_series(d_io, "3. Innovation Outliers (IO)")
  
  # 4. Replacement Outliers (RO with 4x Beta Amplification)
  d_ro <- generate_tv_regression_dgp(n = 300, model_type = "LMHC", contamination = "RO",
                                     ro_multiplier = 4.0, epsilon = 0.08, kappa = 0.40, seed = 123)
  plot_contaminated_series(d_ro, "4. Replacement Outliers (RO)")
  
  par(mfrow = c(1, 1))
}
}

