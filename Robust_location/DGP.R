# Generate data for change in shift in location
# according to    X_{t,n} = \mu(u) + \sum_{j=1}^p a_j X_{t-j,n} + \sigma(u) e_t + \sum_{k=1}^q b_k \sigma(u_k) e_{t-k}

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
#'        - "ii" (or "cyclic"): Cyclic variance σ(u) = 1 + sin(2 * pi * u)
#'        - "iii" (or "break"): Abrupt variance break σ(u) = 0.5 * I(u <= 0.5) + 1.0 * I(u > 0.5)
#'        - "constant" (or "none"): Homoskedastic σ(u) = 1.0
#' @param innov_dist String, the distribution of the innovations (ϵ_t).
#'        One of "gaussian" or "t3".
#' @param contamination_scenario String, the contamination scenario.
#'        One of "clean" / "none", "AO" (Additive Outliers), or "IO" (Innovation Outliers).
#' @param epsilon Numeric, the probability of outlier occurrence (Bernoulli parameter).
#' @param gamma Numeric, for "IO" and "AO", the fixed outlier magnitude multiplier.
#'
#' @return A list containing:
#'         - Xt: The generated time series (numeric vector).
#'         - m_u: The true local expected value m(u) (numeric vector).
#'         - mu_u: The time-varying intercept μ(u) (numeric vector).
#'         - sigma_u: The time-varying scale σ(u) (numeric vector).
ARMA_mu <- function(n, ar_coeffs = NULL, ma_coeffs = NULL, mu_scenario = "H0", k = 0.5,
                    var_scenario = "i", innov_dist = "gaussian", contamination_scenario = "clean",
                    epsilon = 0.05, gamma = 10) {

  # 1. Generate the time-varying intercept mu(u) and scale sigma(u)
  u <- (1:n) / n
  mu_u <- numeric(n)
  if (mu_scenario == "H1") {
    mu_u[u > 0.5] <- k
  } else if (mu_scenario == "H2") {
    mu_u <- k * u
  }

  # Scale function sigma(u) under the 3 heteroskedasticity scenarios
  var_scen <- tolower(as.character(var_scenario))
  if (var_scen %in% c("i", "trending", "smooth", "exp")) {
    sigma_u <- exp(u / 2)
  } else if (var_scen %in% c("ii", "cyclic", "sin")) {
    sigma_u <- 1 + sin(2 * pi * u)
  } else if (var_scen %in% c("iii", "break", "abrupt")) {
    sigma_u <- ifelse(u <= 0.5, 0.5, 1.0)
  } else if (var_scen %in% c("none", "constant", "homoskedastic")) {
    sigma_u <- rep(1.0, n)
  } else {
    warning(sprintf("Unknown variance scenario '%s', defaulting to scenario (i).", var_scenario))
    sigma_u <- exp(u / 2)
  }

  # 2. Calculate the true local expected value m(u)
  sum_ar <- if (!is.null(ar_coeffs)) sum(ar_coeffs) else 0
  m_u <- mu_u / (1 - sum_ar)

  # 3. Generate white noise innovations
  if (innov_dist == "gaussian") {
    innovations <- rnorm(n)
  } else if (innov_dist == "t3") {
    innovations <- rt(n, df = 3)
  } else {
    stop("Unknown innovation distribution.")
  }

  # 4. Apply Innovation Outliers (IO) 
  if (contamination_scenario == "IO") {
    It <- rbinom(n, 1, epsilon)
    Xi <- rbinom(n, 1, 0.5) * 2 - 1  # Generates -1 or 1 with equal probability
    innovations <- innovations + It * gamma * Xi 
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
    
    # Equation: X_{t,n} = mu(u) + sum(a_j * X_{t-j}) + sigma(u)*e_t + sum(b_k * sigma(u_k)*e_{t-k})
    padded_Xt[t + p] <- mu_u[t] + ar_term + padded_innovations[t + q] + ma_term
  }
  Xt <- padded_Xt[(p+1):(n+p)]

  # 6. Apply Additive Outliers (AO) using gamma as the fixed magnitude
  if (contamination_scenario == "AO") {
    It <- rbinom(n, 1, epsilon)
    Xi <- rbinom(n, 1, 0.5) * 2 - 1  # Generates -1 or 1 with equal probability
    Xt <- Xt + It * gamma * Xi 
  }

  return(list(Xt = Xt, m_u = m_u, mu_u = mu_u, sigma_u = sigma_u))
}

