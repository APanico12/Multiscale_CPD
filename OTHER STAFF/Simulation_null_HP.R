# Author: Antonio Panico
# This script runs the simulation under the null hypothesis of no change point,
#for different DGPs and sample sizes, and produces QQ plots 
#comparing the empirical distribution of the test statistic to the theoretical Kolmogorov distribution. 
#The DGPs include iid, shifted mean, drift, heteroskedasticity, and mixed cases.
#The test statistic is computed using a rolling window approach to estimate local density changes. 

library(MASS)
library(zoo)
library(cumulcalib)
require(doParallel)
require(foreach)
require(SymTS)

source("network_functions.R")
source("DGP.R")
source("Plot_functions.R")
# ─────────────────────────────────────────────
# Simulation settings
# ─────────────────────────────────────────────
MC         <- 1e5
MC.Bootstrap <- 1000
n_values  <- c(100,500, 1000, 2000, 5000)
do = 0.5   
#A <- diag(5)
A <- matrix(c(1, 1, 1, 1, 2,
              1, 1, 0, 1, 0,
              1, 0, 1, 0, 0,
              1, 1, 0, 1, 1,
              2, 0, 0, 1, 1),
            nrow = 5, ncol = 5, byrow = TRUE)

sigma <- Vectorize(function(u) abs(sin(u * 2*pi)^2) + 0.5)
alpha <- Vectorize(function(u) 1 + 10 * (u > 0.5))
df    <- Vectorize(function(u) 2 + 10 * (u > 0.5))
a = function(u) 0.2 + 0.1*u

#list of different DGPs to simulate source: DGP.R
dgp_functions <- list(
  iid     = gen_iid,
  iid_change = gen_iid_with_change,
  shifted = gen_shifted,
  drift   = gen_drift,
  hetero  = gen_hetero,
  mixed   = gen_mixed
 )

# ─────────────────────────────────────────────
# Single replication function
# ─────────────────────────────────────────────
run_one <- function(n, dgp_name, dgp_functions, A, sigma, alpha,a) {
  
  X_sim <- dgp_functions[[dgp_name]](A, n)
  W     <- floor(n^(1/2))
  
  D_array_raw <- rollapply(
    data      = X_sim,
    width     = W,
    FUN       = compute_local_density,
    by.column = FALSE,
    align     = "right",
    fill      = NA
  )
  
  D_array <- D_array_raw[!is.na(D_array_raw)]
  steps   <- length(D_array)
  
  m_n_full  <- cumsum(D_array) / steps
  total_m_n <- tail(m_n_full, 1)
  u_seq     <- seq_len(steps) / steps
  T_n_u     <- m_n_full - u_seq * total_m_n
  test_stat <- max(abs(T_n_u))
  return(sqrt(n) * test_stat)
}


# ─────────────────────────────────────────────
# Simulation grid
# ─────────────────────────────────────────────
sim_grid <- expand.grid(
  iteration = 1:M,
  n         = n_values,
  dgp_name  = names(dgp_functions),
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────
# Parallel backend
# ─────────────────────────────────────────────
set.seed(123, kind = "L'Ecuyer-CMRG")

cores <- detectCores() 
cl    <- makeCluster(cores)
registerDoParallel(cl)

clusterExport(cl, c("dgp_functions", "A", "sigma", "alpha",
                    "gen_iid", "gen_shifted", "gen_drift",
                    "gen_hetero", "gen_mixed",
                    "compute_local_density", "run_one",
                    "do","a","rAR1"))      

clusterEvalQ(cl, {
  source("network_functions.R")   
  library(zoo)
  library(SymTS)
})

# ─────────────────────────────────────────────
# Parallel simulation
# ─────────────────────────────────────────────
cat("Running", nrow(sim_grid), "simulations across", cores, "cores...\n")

sim_results <- foreach(
  row      = iter(sim_grid, by = "row"),
  .combine = rbind
) %dopar% {
  stat <- run_one(
    n         = row$n,
    dgp_name  = row$dgp_name,
    dgp_functions = dgp_functions,
    A         = A,
    sigma     = sigma,
    alpha     = alpha
  )
  data.frame(
    iteration = row$iteration,
    n         = row$n,
    dgp_name  = row$dgp_name,
    stat      = stat
  )
}

stopCluster(cl)

# ─────────────────────────────────────────────
# QQ plots
# ─────────────────────────────────────────────

p                     <- (1:M - 0.5) / M
theoretical_quantiles <- sapply(p, qKolmogorov)
dgp_names <- names(dgp_functions)
n_rows <- length(dgp_names)
n_cols <- length(n_values)

png("qqplots_all.png",
    width  = 300 * n_cols,
    height = 350 * n_rows,
    res    = 96)

par(mfrow = c(n_rows, n_cols),
    mar   = c(3, 3, 3, 1),
    oma   = c(0, 3, 3, 0))   # outer margins for row/column labels

for (dgp in dgp_names) {
  for (n_val in n_values) {
    
    stats <- sim_results$stat[
      sim_results$dgp_name == dgp & sim_results$n == n_val
    ]
    empirical_quantiles <- sort(stats)
    fit <- lm(empirical_quantiles ~ theoretical_quantiles)
    
    qqplot(theoretical_quantiles, empirical_quantiles,
           main = sprintf("n = %d", n_val),
           xlab = "Kolmogorov quantiles",
           ylab = "",#if (n_val == n_values[1]) dgp else "",
           pch  = 19, col = "darkblue", cex = 0.5,
           axes = TRUE)
    
    abline(fit, col = "red", lwd = 2, lty = 2)
    
    # Row label on the left for the first column only
    if (n_val == n_values[1]) {
      mtext(dgp, side = 2, line = 3, cex = 0.85, font = 2)
    }
  }
}

# Overall title
mtext("QQ plots: empirical vs Kolmogorov quantiles",
      side = 3, outer = TRUE, line = 1, cex = 1.1, font = 2)

dev.off()
cat("Saved: qqplots_all.png\n")