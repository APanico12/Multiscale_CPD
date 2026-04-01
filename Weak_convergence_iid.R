# Author: Antonio Panico
# Simulates the null distribution of the CUSUM test for Network Density
# and checks whether sqrt(n) * T*_n converges to the KS (Brownian Bridge) distribution.

library(MASS)
#install.packages("zoo")
library(zoo)
#install.packages("mvpd")
library(mvpd) #contains Kolmogorov distribution
#install.packages("cumulcalib")
library(cumulcalib) # contains qKolmogorov to generate quantile of Kolmogorov distribution
source("network_functions.R")

# ─────────────────────────────────────────────
# Simulation settings
# ─────────────────────────────────────────────
M   <- 1000
n   <- c(100, 200, 500, 1000, 2000, 5000)
d   <- 5

mu    <- rep(0, d)
sigma <- diag(d)

# Results matrix: one column per sample size
results <- matrix(NA, nrow = M, ncol = length(n))

# ─────────────────────────────────────────────
# Simulation loop
# ─────────────────────────────────────────────
for (i in seq_along(n)) {
  
  cat(sprintf("Running n = %d ...\n", n[i]))
  
  W <- floor(n[i]^(2/3))   # local window width
  
  for (j in 1:M) {
    
    X_sim <- mvrnorm(n = n[i], mu = mu, Sigma = sigma)
    
    D_array_raw <- rollapply(
      data       = X_sim,
      width      = W,
      FUN        = compute_local_density,
      by.column  = FALSE,
      align      = "right",
      fill       = NA
    )
    
    #rollapply returns NAs for the first (W-1) positions.
    D_array <- D_array_raw[!is.na(D_array_raw)]
    
    steps   <- length(D_array)
             
    # CUMSUM test statistic 
    
    m_n_full  <- cumsum(D_array) / steps
    total_m_n <- tail(m_n_full, 1)
    u_seq     <- seq_len(steps) / steps
  
    T_n_u <- m_n_full - u_seq * total_m_n
    
    # FIX Look here with Fabian
    test_stat <- max(abs(T_n_u))
    
    results[j, i] <- sqrt(n[i]) * test_stat
  }
}

# ─────────────────────────────────────────────
# QQ plots: compare empirical quantiles against KS distribution
# ─────────────────────────────────────────────
#generate Kolmogorov distribution quantiles using rkolm  to approximate them 

p <- (1:M - 0.5) / M
theoretical_quantiles <- sapply(p, qKolmogorov)

# ─────────────────────────────────────────────
# Plotting the Q-Q plot for each sample size
# ─────────────────────────────────────────────
# sort results for each column so that we can build the emipirical quantiles
sorted_results <- apply(results, 2, sort)

for (i in seq_along(n)) {
  # Linear regression between theoretical and empirical quantiles
  fit <- lm(sorted_results[, i] ~ theoretical_quantiles)
  
  qqplot(theoretical_quantiles, sorted_results[, i],
         main = sprintf("Q-Q Plot: n = %d", n[i]),
         xlab = "Sample Data",
         ylab = "Empirical Quantiles",
         pch  = 19, col = "darkblue")
 
  abline(fit, col = "green", lwd = 2, lty = 2)  # regression line
}

 