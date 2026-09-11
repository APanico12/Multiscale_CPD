# Generate data for change in shift in location
# according to    X_{t,n} = \mu(u) + \sum_{j=1}^p a_j X_{t-j,n} + \epsilon_t + \sum_{k=1}^q b_k \epsilon_{t-k}

#' Generate time series from an ARMA(p,q) model with a time-varying intercept and optional contamination.
#'
#' @param n Integer, the length of the time series.
#' @param ar_coeffs Numeric vector, the AR(p) coefficients (a_j).
#' @param ma_coeffs Numeric vector, the MA(q) coefficients (b_k).
#' @param mu_scenario String, the scenario for the time-varying intercept μ(u).
#'        One of "H0", "H1", "H2".
#' @param k Numeric, the magnitude of the shift for scenarios H1 and H2.
#' @param innov_dist String, the distribution of the innovations (ϵ_t).
#'        One of "gaussian" or "t3".
#' @param contamination_scenario String, the contamination scenario.
#'        One of "none", "AO" (Additive Outliers), or "IO" (Innovation Outliers).
#' @param epsilon Numeric, the probability of outlier occurrence (Bernoulli parameter).
#' @param gamma Numeric, for "IO" (Innovation Outliers) it's the scale parameter for the Cauchy distribution of outliers.
#'        For "AO" (Additive Outliers) it's the fixed value to be added when an outlier occurs.
#'
#' @return A list containing:
#'         - Xt: The generated time series (numeric vector).
#'         - m_u: The true local expected value m(u) (numeric vector).
#'         - mu_u: The time-varying intercept μ(u) (numeric vector).
ARMA_mu <- function(n, ar_coeffs = NULL, ma_coeffs = NULL, mu_scenario = "H0", k = 1,
                    innov_dist = "gaussian", contamination_scenario = "clean",
                    epsilon = 0.05, gamma = 10) {

  # 1. Generate the time-varying intercept mu(u)
  u <- (1:n) / n
  mu_u <- numeric(n)
  if (mu_scenario == "H1") {
    mu_u[u > 0.5] <- k
  } else if (mu_scenario == "H2") {
    mu_u <- k * u
  }

  # 2. Calculate the true local expected value m(u)
  sum_ar <- if (!is.null(ar_coeffs)) sum(ar_coeffs) else 0
  m_u <- mu_u / (1 - sum_ar)

  # 3. Generate innovations
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

  # 5. Simulate the ARMA(p,q) process
  p <- if (!is.null(ar_coeffs)) length(ar_coeffs) else 0
  q <- if (!is.null(ma_coeffs)) length(ma_coeffs) else 0
  Xt <- numeric(n)
  
  # Pad innovations and Xt for easier indexing
  padded_innovations <- c(rep(0, q), innovations)
  padded_Xt <- c(rep(0, p), Xt)

  for (t in 1:n) {
    ar_term <- if (p > 0) sum(ar_coeffs * padded_Xt[(t+p-1):(t)]) else 0
    ma_term <- if (q > 0) sum(ma_coeffs * padded_innovations[(t+q-1):(t)]) else 0
    
    # The model is defined with +epsilon_t and +b_k*epsilon_{t-k}
    padded_Xt[t + p] <- mu_u[t] + ar_term + padded_innovations[t + q] + ma_term
  }
  Xt <- padded_Xt[(p+1):(n+p)]

  # 6. Apply Additive Outliers (AO) using gamma as the fixed magnitude
  if (contamination_scenario == "AO") {
    It <- rbinom(n, 1, epsilon)
    Xi <- rbinom(n, 1, 0.5) * 2 - 1  # Generates -1 or 1 with equal probability
    Xt <- Xt + It * gamma * Xi 
  }

  return(list(Xt = Xt, m_u = m_u, mu_u = mu_u))
}

# # Example Usage:
# #
# set.seed(173)
# # Set parameters
# n_obs <- 2000
# ar_params <- c(0.5, -0.1) # AR(2)
# ma_params <- c(0.4)       # MA(1)
# shift_k <- 2

# # H0 with Gaussian innovations
# h0_data <- ARMA_mu(n = n_obs, ar_coeffs = ar_params, ma_coeffs = ma_params)

# # H1 with t3 innovations
# h1_t_data <- ARMA_mu(n = n_obs, ar_coeffs = ar_params, ma_coeffs = ma_params,
#                      contamination_scenario = "AO",  mu_scenario = "H1", k = shift_k, innov_dist = "t3", epsilon = 0.05, gamma = 10)

# # H2 with Additive Outliers
# h2_ao_data <- ARMA_mu(n = n_obs, ar_coeffs = ar_params, ma_coeffs = ma_params,
#                       mu_scenario = "H2", k = shift_k, innov_dist = "t3",
#                       contamination_scenario = "IO", epsilon = 0.05, gamma = 10)

# # Plot the results
# par(mfrow=c(3,1), mar=c(4,4,2,1))
# plot(h0_data$Xt, type='l', main="H0: No Shift, Gaussian Innovations", ylab="Xt")
# lines(h0_data$m_u, col='red', lwd=2)

# plot(h1_t_data$Xt, type='l', main="H1: Abrupt Shift, t3 Innovations", ylab="Xt")
# lines(h1_t_data$m_u, col='red', lwd=2)

# plot(h2_ao_data$Xt, type='l', main="H2: Gradual Shift, Additive Outliers", ylab="Xt")
# lines(h2_ao_data$m_u, col='red', lwd=2)

# #plot true integrated parameters
# plot(cumsum(h2_ao_data$m_u)/n_obs,main="H2: Gradual Shift, Additive Outliers integrated mu")
# plot(cumsum(h1_t_data$m_u)/n_obs,main="H1: Abrupt Shift, t3 Innovations integrated mu")  

# #compare estimators at 1 
# normal = int_teta[n_obs]
# lin = int_lin_teta[n_obs]
# oracle_vec = cumsum(h1_t_data$m_u)/n_obs
# oracle = oracle_vec[n_obs] # Safely extracting the scalar value
# cat(normal,lin,oracle)