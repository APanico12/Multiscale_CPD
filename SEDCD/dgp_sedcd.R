# ==============================================================================
# File: dgp_sedcd.R
# Description: Multivariate Simulation Data Generating Process (DGP) for
#              Smoothed Extremal Downside Correlation Density (SEDCD)
#              Change-Point Detection.
#
# Specifications (Section 5.3):
#   1. Dimension d in {2, 5} (supports arbitrary d >= 2) and sample sizes
#      n in {1000, 2000, 5000}.
#   2. Dynamic univariate conditional margins (Patton, 2006):
#        X_{t,j} = mu_{t,j}(beta_j) + sigma_{t,j}(beta_j) * eps_{t,j}
#      where:
#        - mu_{t,j} follows ARMA(1, 1):
#            mu_{t,j} = phi_{0,j} + phi_{1,j} X_{t-1,j} + theta_{1,j} a_{t-1,j}
#        - sigma_{t,j}^2 follows GJR-GARCH(1, 1) (Glosten et al., 1993):
#            sigma_{t,j}^2 = omega_j + (alpha_j + gamma_j * I(a_{t-1,j} < 0)) a_{t-1,j}^2 + beta_j sigma_{t-1,j}^2
#        - Marginal innovations eps_{t,j} i.i.d. ~ t^*_{nu_j}(0, 1) standardized Student-t:
#            F_j(z; nu_j) = T_{nu_j}(z * sqrt(nu_j / (nu_j - 2)))
#        - Conditionally uniform random variables via PIT:
#            U_{t,j} = F_j(eps_{t,j}; nu_j) ~ U(0, 1).
#
#   3. Conditional copula C(. | F_{t-1}; theta_t) driven by GAS(1, 1)
#      (De Lira Salvatierra & Patton, 2015):
#      (i) Student's t-copula (C^t_{rho_t, nu}):
#            - Exchangeable dynamic linear correlation rho_t in (-1, 1)
#            - Constant degrees of freedom nu > 2
#            - Link: rho_t = Lambda(psi_t) = tanh(psi_t) in (-1, 1)
#      (ii) Clayton copula (C^C_{theta_t}):
#            - Asymmetric lower tail dependence: lambda_l = 2^(-1/theta_t), lambda_u = 0
#            - Link: theta_t = Lambda(psi_t) = exp(psi_t) in (0, Inf)
#
#   4. GAS(1, 1) autoregressive score recursion:
#        psi_t = omega_t + beta * psi_{t-1} + alpha * s_{t-1}
#      where s_{t-1} = S_{t-1} * nabla_{t-1} is the scaled score of the conditional copula log-likelihood:
#        nabla_{t-1} = (d log c(u_{t-1}; theta_{t-1}) / d theta_{t-1}) * Lambda'(psi_{t-1})
#      and directional scaling S_{t-1} = I_{t-1}^{-1/2} (inverse Fisher information square root).
#
#   5. Structural break hypotheses on baseline intercept omega_t:
#        H0: omega_t = omega_0 (null of structural stability)
#        H1: omega_t = omega_pre for t <= floor(u^* n),
#            omega_post = omega_pre + delta_omega for t > floor(u^* n) (abrupt change)
#        H2: omega_t = omega_0 + delta_omega * (t / n) (gradual drift)
# ==============================================================================

#' Simulate Multivariate Time Series with Dynamic GAS Copula & ARMA-GJR-GARCH Margins
#'
#' @param n Sample size (integer, e.g. 1000, 2000, 5000)
#' @param d Multivariate dimension (integer >= 2, default 2; also supports d = 5)
#' @param scenario "H0" (constant), "H1" (abrupt jump at u_star), or "H2" (gradual drift)
#' @param copula_type "clayton" (asymmetric lower tail dependence) or "t" (Student's t copula)
#' @param omega_0 Baseline GAS intercept (scalar, or NULL for automatic stationary calibration)
#' @param delta_omega Shift magnitude for intercept under H1 and H2 (scalar, or NULL for default)
#' @param u_star Break point fraction under H1 (default 0.50)
#' @param alpha_gas GAS score coefficient (default 0.08)
#' @param beta_gas GAS autoregressive persistence coefficient (default 0.90)
#' @param df_copula Degrees of freedom for Student's t copula (default 4)
#' @param df_margins Degrees of freedom for marginal standardized Student-t (scalar or vector of length d, default 5)
#' @param arma_params List containing ARMA(1, 1) parameters: phi0, phi1, theta1 (default list(phi0 = 0.0, phi1 = 0.05, theta1 = -0.02))
#' @param garch_params List containing GJR-GARCH(1, 1) parameters: omega, alpha, gamma, beta (default list(omega = 0.02, alpha = 0.04, gamma = 0.08, beta = 0.85))
#' @param alpha_tail Target quantile probability for adaptive thresholding (default 0.10)
#' @param burn_in Integer burn-in length to eliminate initialization transients (default 500)
#' @param scaling Scaling matrix for score: "fisher" (inverse Fisher information) or "unit" (S_t = 1)
#' @param seed Optional integer seed for reproducibility
#' @param ... Optional arguments for backward compatibility
#' @return List containing:
#'   - X: Observed multivariate series matrix N x d (including ARMA-GJR-GARCH nuisance dynamics)
#'   - X_tilde: Demeaned and volatility-standardized marginal series N x d
#'   - U: Copula uniform realizations N x d
#'   - eps: Standardized Student-t marginal innovations N x d
#'   - sigma_matrix: Conditional volatility matrix N x d
#'   - mu_matrix: Conditional mean matrix N x d
#'   - theta_path: Time-varying copula parameter vector of length N
#'   - lambda_L_path: Theoretical lower tail dependence path of length N
#'   - psi_path: Latent unconstrained parameter path of length N
#'   - score_path: Scaled score innovations s_t of length N
#'   - omega_path: Time-varying intercept omega_t path of length N
#'   - tau: Adaptive threshold scalar (empirical alpha_tail quantile of X)
#'   - n: Sample size
#'   - d: Dimension
#'   - scenario: Hypothesis scenario ("H0", "H1", "H2")
#'   - copula_type: Copula family ("clayton", "t")
simulate_sedcd_dgp <- function(n = 2000,
                               d = 2,
                               scenario = c("H0", "H1", "H2"),
                               copula_type = c("clayton", "t"),
                               omega_0 = NULL,
                               delta_omega = NULL,
                               u_star = 0.50,
                               alpha_gas = 0.08,
                               beta_gas = 0.90,
                               df_copula = 4,
                               df_margins = 5,
                               arma_params = list(phi0 = 0.0, phi1 = 0.05, theta1 = -0.02),
                               garch_params = list(omega = 0.02, alpha = 0.04, gamma = 0.08, beta = 0.85),
                               alpha_tail = 0.10,
                               burn_in = 500,
                               scaling = c("fisher", "unit"),
                               seed = NULL,
                               ...) {
  if (!is.null(seed)) set.seed(seed)
  scenario    <- match.arg(scenario)
  copula_type <- match.arg(copula_type)
  scaling     <- match.arg(scaling)
  
  # Backward compatibility: extract potential legacy parameters
  dots <- list(...)
  if ("marginal" %in% names(dots)) {
    if (dots$marginal == "t5") df_margins <- 5
  }
  
  if (d < 2) stop("Dimension d must be >= 2")
  if (length(df_margins) == 1) {
    df_margins <- rep(df_margins, d)
  }
  if (length(df_margins) != d) {
    stop("Argument 'df_margins' must have length 1 or d")
  }
  if (any(df_margins <= 2)) stop("Marginal degrees of freedom must be strictly > 2")
  
  total_n <- n + burn_in
  
  # ----------------------------------------------------------------------------
  # 1. Baseline Calibration of omega_0 and delta_omega
  # ----------------------------------------------------------------------------
  if (is.null(omega_0)) {
    if (copula_type == "clayton") {
      # Target baseline theta_0 = 1.5 -> lambda_L = 2^(-1/1.5) = 0.630
      target_theta0 <- 1.5
      psi_target    <- log(target_theta0)
      omega_0       <- (1.0 - beta_gas) * psi_target
    } else {
      # Target baseline rho_0 = 0.40 -> moderate positive lower tail dependence
      target_rho0 <- 0.40
      psi_target  <- atanh(target_rho0)
      omega_0     <- (1.0 - beta_gas) * psi_target
    }
  }
  
  if (is.null(delta_omega)) {
    if (copula_type == "clayton") {
      # Target shift: theta jumps from 1.5 to 3.0 -> lambda_L from 0.630 to 0.794
      delta_omega <- (1.0 - beta_gas) * log(2.0)
    } else {
      # Target shift: rho jumps from 0.40 to 0.75 -> substantial tail-dependence shift
      delta_omega <- (1.0 - beta_gas) * (atanh(0.75) - atanh(0.40))
    }
  }
  
  # ----------------------------------------------------------------------------
  # 2. Construct Time-Varying Intercept omega_t Path
  # ----------------------------------------------------------------------------
  omega_path <- numeric(total_n)
  omega_path[1:burn_in] <- omega_0
  
  u_grid <- (1:n) / n
  idx_analysis <- (burn_in + 1):total_n
  
  if (scenario == "H0") {
    omega_path[idx_analysis] <- omega_0
  } else if (scenario == "H1") {
    omega_path[idx_analysis] <- ifelse(u_grid <= u_star, omega_0, omega_0 + delta_omega)
  } else if (scenario == "H2") {
    omega_path[idx_analysis] <- omega_0 + delta_omega * u_grid
  }
  
  # ----------------------------------------------------------------------------
  # 3. Pre-allocate Storage & Initialize State
  # ----------------------------------------------------------------------------
  psi_path      <- numeric(total_n)
  theta_path    <- numeric(total_n)
  lambda_L_path <- numeric(total_n)
  score_path    <- numeric(total_n)
  
  U_mat         <- matrix(0.0, nrow = total_n, ncol = d)
  eps_mat       <- matrix(0.0, nrow = total_n, ncol = d)
  sigma_mat     <- matrix(0.0, nrow = total_n, ncol = d)
  mu_mat        <- matrix(0.0, nrow = total_n, ncol = d)
  a_mat         <- matrix(0.0, nrow = total_n, ncol = d)
  X_mat         <- matrix(0.0, nrow = total_n, ncol = d)
  
  # Stationary GJR-GARCH unconditional variance initialization
  denom_garch <- 1.0 - garch_params$alpha - 0.5 * garch_params$gamma - garch_params$beta
  if (denom_garch <= 0) denom_garch <- 0.05
  sigma2_init <- garch_params$omega / denom_garch
  
  sigma_mat[1, ] <- sqrt(sigma2_init)
  mu_mat[1, ]    <- arma_params$phi0 / (1.0 - arma_params$phi1)
  
  # Stationary unconditional expectation for latent state psi
  psi_path[1] <- omega_0 / (1.0 - beta_gas)
  
  eps_tol <- 1e-9
  
  # ----------------------------------------------------------------------------
  # 4. Sequential Simulation Loop: GAS Dynamics + Copula + ARMA-GJR-GARCH
  # ----------------------------------------------------------------------------
  for (t in 1:total_n) {
    psi_cur <- psi_path[t]
    
    # 4.1 Copula Calibration and Variate Generation
    if (copula_type == "clayton") {
      # Clayton calibration: theta_t = exp(psi_t) in (0, Inf)
      psi_cur_clamped <- max(-4.0, min(4.0, psi_cur))
      theta_cur <- exp(psi_cur_clamped)
      theta_path[t] <- theta_cur
      lambda_L_path[t] <- 2.0^(-1.0 / theta_cur)
      
      # Marshall-Olkin (1988) algorithm for d-dimensional Clayton copula:
      # V_t ~ Gamma(shape = 1 / theta_t, rate = 1)
      # E_{t,j} ~ Exponential(1) i.i.d. for j = 1, ..., d
      # U_{t,j} = (1 + E_{t,j} / V_t)^(-1 / theta_t)
      shape_gamma <- 1.0 / max(1e-4, theta_cur)
      V_t <- rgamma(1, shape = shape_gamma, rate = 1.0)
      V_t <- max(1e-12, V_t)
      E_t <- rexp(d, rate = 1.0)
      
      u_t <- (1.0 + E_t / V_t)^(-1.0 / theta_cur)
      u_t <- pmin(pmax(u_t, eps_tol), 1.0 - eps_tol)
      U_mat[t, ] <- u_t
      
      # Exact Score of Clayton log-copula density wrt theta
      u_pow <- u_t^(-theta_cur)
      A_t <- sum(u_pow) - d + 1.0
      A_t <- max(eps_tol, A_t)
      
      sum_j <- sum((1:(d - 1)) / (1.0 + (1:(d - 1)) * theta_cur))
      term_log_u <- sum(log(u_t))
      term_pow_log <- sum(u_pow * log(u_t))
      
      dlogc_dtheta <- sum_j - term_log_u + (1.0 / (theta_cur^2)) * log(A_t) +
                      (d + 1.0 / theta_cur) * (term_pow_log / A_t)
      
      # Score wrt unconstrained parameter psi: nabla_t = dlogc/dtheta * Lambda'(psi)
      nabla_t <- dlogc_dtheta * theta_cur
      
      # Fisher Information Directional Scaling: S_t = I_t^(-1/2)
      if (scaling == "fisher") {
        I_theta <- (d - 1.0) / ((1.0 + theta_cur) * (1.0 + 2.0 * theta_cur))
        I_psi <- max(1e-4, I_theta * (theta_cur^2))
        S_t <- 1.0 / sqrt(I_psi)
      } else {
        S_t <- 1.0
      }
      
      s_t <- S_t * nabla_t
      # Stabilizing winsorization
      s_t <- max(-10.0, min(10.0, s_t))
      score_path[t] <- s_t
      
    } else {
      # Student's t-copula calibration: rho_t = tanh(psi_t) in (-1, 1)
      psi_cur_clamped <- max(-4.0, min(4.0, psi_cur))
      rho_cur <- tanh(psi_cur_clamped)
      
      # Bound rho within admissible positive definite range for dimension d
      rho_min <- -1.0 / (d - 1.0) + 0.02
      rho_max <- 0.98
      rho_cur <- max(rho_min, min(rho_max, rho_cur))
      
      theta_path[t] <- rho_cur
      # Theoretical symmetric tail dependence: lambda_l = lambda_u
      arg_td <- -sqrt((df_copula + 1.0) * (1.0 - rho_cur) / (1.0 + rho_cur))
      lambda_L_path[t] <- 2.0 * pt(arg_td, df = df_copula + 1.0)
      
      # Generate d-variate Student's t copula variates with exchangeable correlation R_t
      R_mat <- matrix(rho_cur, nrow = d, ncol = d)
      diag(R_mat) <- 1.0
      
      chol_R <- chol(R_mat)
      z_std  <- rnorm(d)
      z_vec  <- as.numeric(z_std %*% chol_R)
      s_chi  <- rchisq(1, df = df_copula)
      xi_vec <- z_vec / sqrt(s_chi / df_copula)
      
      u_t <- pt(xi_vec, df = df_copula)
      u_t <- pmin(pmax(u_t, eps_tol), 1.0 - eps_tol)
      U_mat[t, ] <- u_t
      
      # Exact Score of Student's t log-copula density wrt rho
      xi_t <- qt(u_t, df = df_copula)
      S1_t <- sum(xi_t^2)
      S2_t <- (sum(xi_t))^2
      
      denom_det1 <- 1.0 - rho_cur
      denom_det2 <- 1.0 + (d - 1.0) * rho_cur
      
      Q_t <- (S1_t / denom_det1) - (rho_cur / (denom_det1 * denom_det2)) * S2_t
      Q_t <- max(1e-4, Q_t)
      
      dQ_drho <- (S1_t / (denom_det1^2)) -
                 ((1.0 + (d - 1.0) * (rho_cur^2)) / ((denom_det1^2) * (denom_det2^2))) * S2_t
      
      dlogc_drho <- (d * (d - 1.0) * rho_cur) / (2.0 * denom_det1 * denom_det2) -
                    ((df_copula + d) / (2.0 * (df_copula + Q_t))) * dQ_drho
      
      # Score wrt unconstrained parameter psi: nabla_t = dlogc/drho * (1 - rho^2)
      nabla_t <- dlogc_drho * (1.0 - rho_cur^2)
      
      # Fisher Information Directional Scaling: S_t = I_t^(-1/2)
      if (scaling == "fisher") {
        I_psi <- (d * (d - 1.0) * (1.0 + (d - 1.0) * rho_cur^2) * ((1.0 + rho_cur)^2)) /
                 (2.0 * (denom_det2^2)) * ((df_copula + d) / (df_copula + d + 2.0))
        I_psi <- max(1e-4, I_psi)
        S_t <- 1.0 / sqrt(I_psi)
      } else {
        S_t <- 1.0
      }
      
      s_t <- S_t * nabla_t
      s_t <- max(-10.0, min(10.0, s_t))
      score_path[t] <- s_t
    }
    
    # 4.2 Standardized Student-t Marginal Innovations
    for (j in 1:d) {
      nu_j <- df_margins[j]
      # F_j(z; nu_j) = T_{nu_j}(z * sqrt(nu_j / (nu_j - 2)))
      # eps_{t,j} = sqrt((nu_j - 2) / nu_j) * T_{nu_j}^(-1)(U_{t,j})
      eps_mat[t, j] <- sqrt((nu_j - 2.0) / nu_j) * qt(U_mat[t, j], df = nu_j)
    }
    
    # 4.3 ARMA(1, 1) Mean & GJR-GARCH(1, 1) Volatility Updates
    if (t > 1) {
      for (j in 1:d) {
        prev_a    <- a_mat[t - 1, j]
        prev_sig2 <- sigma_mat[t - 1, j]^2
        
        # GJR-GARCH(1, 1) asymmetric leverage volatility:
        # sigma_{t,j}^2 = omega + (alpha + gamma * I(a_{t-1} < 0)) a_{t-1}^2 + beta * sigma_{t-1}^2
        leverage_ind <- if (prev_a < 0.0) 1.0 else 0.0
        arch_coeff   <- garch_params$alpha + garch_params$gamma * leverage_ind
        sig2_t       <- garch_params$omega + arch_coeff * (prev_a^2) + garch_params$beta * prev_sig2
        sigma_mat[t, j] <- sqrt(max(1e-6, sig2_t))
        
        # ARMA(1, 1) conditional mean:
        # mu_{t,j} = phi0 + phi1 * X_{t-1,j} + theta1 * a_{t-1,j}
        mu_mat[t, j] <- arma_params$phi0 + arma_params$phi1 * X_mat[t - 1, j] +
                        arma_params$theta1 * prev_a
      }
    }
    
    # 4.4 Construct Observed Process: X_{t,j} = mu_{t,j} + sigma_{t,j} * eps_{t,j}
    a_mat[t, ] <- sigma_mat[t, ] * eps_mat[t, ]
    X_mat[t, ] <- mu_mat[t, ] + a_mat[t, ]
    
    # 4.5 Advance GAS State: psi_{t+1} = omega_{t+1} + beta * psi_t + alpha * s_t
    if (t < total_n) {
      psi_path[t + 1] <- omega_path[t + 1] + beta_gas * psi_cur + alpha_gas * s_t
    }
  }
  
  # ----------------------------------------------------------------------------
  # 5. Extraction of Post-Burn-In Sample
  # ----------------------------------------------------------------------------
  keep_idx <- (burn_in + 1):total_n
  
  X_out       <- X_mat[keep_idx, , drop = FALSE]
  U_out       <- U_mat[keep_idx, , drop = FALSE]
  eps_out     <- eps_mat[keep_idx, , drop = FALSE]
  sigma_out   <- sigma_mat[keep_idx, , drop = FALSE]
  mu_out      <- mu_mat[keep_idx, , drop = FALSE]
  theta_out   <- theta_path[keep_idx]
  lambda_out  <- lambda_L_path[keep_idx]
  psi_out     <- psi_path[keep_idx]
  score_out   <- score_path[keep_idx]
  omega_out   <- omega_path[keep_idx]
  
  # Adaptive threshold: empirical lower alpha_tail quantile of X
  tau_empirical <- as.numeric(quantile(X_out, probs = alpha_tail))
  
  colnames(X_out)       <- paste0("X", 1:d)
  colnames(U_out)       <- paste0("U", 1:d)
  colnames(eps_out)     <- paste0("eps", 1:d)
  colnames(sigma_out)   <- paste0("sigma", 1:d)
  colnames(mu_out)      <- paste0("mu", 1:d)
  
  return(list(
    X             = X_out,
    X_tilde       = eps_out,
    U             = U_out,
    eps           = eps_out,
    sigma_matrix  = sigma_out,
    mu_matrix     = mu_out,
    theta_path    = theta_out,
    lambda_L_path = lambda_out,
    psi_path      = psi_out,
    score_path    = score_out,
    omega_path    = omega_out,
    tau           = tau_empirical,
    n             = n,
    d             = d,
    scenario      = scenario,
    copula_type   = copula_type,
    # Backward compatibility aliases
    W             = eps_out[, 1:min(2, d), drop = FALSE],
    V             = U_out[, 1:min(2, d), drop = FALSE],
    scale_matrix  = sigma_out[, 1:min(2, d), drop = FALSE]
  ))
}
