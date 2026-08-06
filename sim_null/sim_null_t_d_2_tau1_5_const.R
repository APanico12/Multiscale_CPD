# Author: Antonio Panico
library(MASS)
# Force underlying math libraries to use only 1 thread per R worker
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(NUMEXPR_NUM_THREADS = 1)

require(doParallel)
require(foreach)


# ─────────────────────────────────────────────
# Source Dependences
# ─────────────────────────────────────────────
source("sedcd_try.R")
source("DGP.R")
#source("Plot_functions.R")


# ─────────────────────────────────────────────
# Simulation settings
# ─────────────────────────────────────────────
MC.Bootstrap <- 1000
MC.levels <- 5000
#Set up DGP
n_values = c(2000)
d_values = c(2)
rho = 0 # pairwisecorr sigmat_matrix
nu = 5 #dgf t
theta.Clayton = 0
nu_param = 5  #df t marginals

#tau and L_n types
tau_values = c(-1.5)
L_n_types = c("1(log(n))^2","0.75(log(n))^2","0.5(log(n))^2","0.25(log(n))^2","0.1(log(n))^2")

#Implement the DGP
dgp_functions <- list(
  
  # Base Copulas (Constant Marginals)
  "T_Copula" = function(n_val, d_val) {
    Gen.from.t(n = n_val, d = d_val, rho = rho, nu = nu,
               margins_df = nu_param, time_varying = FALSE)
  }
)

# ─────────────────────────────────────────────
# Build the Grid
# ─────────────────────────────────────────────

grid_list <- lapply(MC.levels, function(m) {
  expand.grid(
    iteration = 1:m,
    MC_level  = m,
    n         = n_values,
    d         = d_values,
    dgp_name  = names(dgp_functions),
    tau       = tau_values,
    L_n_type  = L_n_types,
    stringsAsFactors = FALSE
  )
})

# Combine the list of grids into one 
sim_grid <- do.call(rbind, grid_list)

# ─────────────────────────────────────────────
# Optimal Bandwidth Selection function
# ─────────────────────────────────────────────
find_k_opt <- function(X, L_n, tau, gamma) {
  N <- nrow(X)
  k_min <- max(1, floor(N^0.35))
  k_max <- floor(N^0.75)
  # evaluate on a small grid of 5 points to save computation time
  k_grid <- unique(floor(seq(k_min, k_max, length.out = 5)))
  
  best_k <- k_max
  best_val <- Inf
  
  for (k in k_grid) {
    theta_hat <- get.teta.sedcd(X, k, tau, gamma)
    
    val <- 0
    for (t in k:(N - L_n)) {
      h_val <- Heq.sedcd(X[t + L_n, , drop = FALSE], theta_hat[t, ], tau, gamma)
      val <- val + sum(h_val^2)
    }
    
    if (val < best_val) {
      best_val <- val
      best_k <- k
    }
  }
  return(best_k)
}

# ─────────────────────────────────────────────
# Single replication function
# ─────────────────────────────────────────────
run_one <- function(n, d, dgp_name, dgp_functions, tau, L_n_type) {
  
  # Generate data from the specified DGP
  X_sim <- dgp_functions[[dgp_name]](n, d)

  # Robustly parse the lag multiplier from the L_n_type string
  lag_multiplier <- as.numeric(regmatches(L_n_type, regexpr("^[0-9.]+", L_n_type)))
  lag <- max(1, floor(lag_multiplier * log(n)^2))
  
  gamma <- 0.1
  
  k_opt <- find_k_opt(X_sim, lag, tau, gamma)
  
  b      <- max(1,lag)
  cutoff <- k_opt + b + lag
  
  res <- CUSUM.sedcd(
    x = X_sim, 
    k = k_opt, 
    cutoff = cutoff, 
    lag = lag, 
    block = b, 
    tau = tau, 
    gamma = gamma, 
    MC = MC.Bootstrap)
  return(res$p_value)
}

# ─────────────────────────────────────────────
# Parallel backend
# ─────────────────────────────────────────────

set.seed(123, kind = "L'Ecuyer-CMRG")
# HPC-aware core detection.
cores_str <- Sys.getenv("SLURM_NTASKS")
if (nchar(cores_str) > 0) {
  cores <- as.integer(cores_str)
  cat("SLURM environment detected. Using allocated cores:", cores, "\n")
} else {
  cores <- detectCores()
  cat("SLURM_NTASKS not found. Using all available local cores:", cores, "\n")
}
cl <- makeCluster(cores)

registerDoParallel(cl)
clusterSetRNGStream(cl, 123)
clusterExport(cl, c("dgp_functions", "run_one", "find_k_opt", "rho", "nu",
                    "theta.Clayton", "nu_param", "MC.Bootstrap"))


clusterEvalQ(cl, {
  # Prevent thread explosion on the worker nodes
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
  source("DGP.R")
  source("sedcd_try.R")
  library(MASS)
  library(numDeriv)
  library(copula)
  library(evd)
})

# ─────────────────────────────────────────────
# Parallel simulation
# ─────────────────────────────────────────────

cat("Running", nrow(sim_grid), "simulations across", cores, "cores...\n")

sim_results <- foreach(
  row      = iter(sim_grid, by = "row"),
  .combine = rbind,
  .packages = c("MASS", "numDeriv", "copula", "evd")
) %dopar% {
  stat <- run_one(
    n         = row$n,
    d         = row$d,
    dgp_name  = row$dgp_name,
    dgp_functions = dgp_functions,
    tau       = row$tau,
    L_n_type  = row$L_n_type
  )
  data.frame(
    iteration = row$iteration,
    level     = row$MC_level,
    d         = row$d,
    dgp_name  = row$dgp_name,
    n_values  = row$n,
    tau       = row$tau,
    L_n_type  = row$L_n_type,
    stat      = stat
  )
}

stopCluster(cl)

write.csv(sim_results, "sim_null_t_d_2_tau1_5_const.csv", row.names = FALSE)
cat("Simulation complete. Results saved to sim_null_t_d_2_tau1_5_const.csv\n")
