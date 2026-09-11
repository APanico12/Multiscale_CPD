# This script runs a Monte Carlo simulation to evaluate the bias and MSE of
# robust location estimators (standard and linearized) under various scenarios.
# It is designed to run in parallel and saves results to a CSV file.

library(doParallel)
library(foreach)

# Source the required function files
source("DGP.R")
source("CUSUM.R")

# --- Simulation Parameters ---
# NOTE: Reduced numbers for quick local testing.

MC.simulations <- 1000 
n_values <- c(100,500,1000,5000) #Sample sizes to test
shift_k_opts <- c(2)   # Shift magnitudes for H1 and H2 scenarios

# --- Model Parameters ---
ar_params <- c(0.2, -0.1) # AR(2)
ma_params <- c(0.2)       # MA(1)
epsilon <- c(0.05)      # Contamination proportion for AO and IO
# --- Contamination and Model Scenarios ---
scenarios.contamination <-c("clean", "AO", "IO") # Contamination scenarios
innov_dist_opts <- c("gaussian", "t3") # Innovation distributions to test
mu_scenarios <- c("H0", "H1", "H2") # Hypotheses scenarios

# --- 2. Simulation Grid Generation ---

cat("Generating simulation grid...\n")

# Base grid with parameters that apply to all scenarios
grid_list <- lapply(MC.simulations, function(m) {
  expand.grid(
    iteration = 1:m,
    MC_level  = m,
    n         = n_values,
    ar_params = list(ar_params),
    ma_params = list(ma_params),
    hp_scenario = mu_scenarios,
    percentage_contamination = epsilon,
    contamination_scenario = scenarios.contamination,
    innov_dist = innov_dist_opts,
    shift_k = shift_k_opts,
    stringsAsFactors = FALSE
  )
})

# Combine the list of grids into one 
sim_grid <- do.call(rbind, grid_list)

# --- 3. Parallel Simulation Execution ---

# This function runs one iteration of the simulation for a given parameter set
run_one <- function(n_val, ar_params, ma_params, innov_dist, shift_k_val, mu_scenario, contamination_scenario, percentage_contamination) {
  
  # Generate the time series data based on the specified DGP and parameters
  ts_data <- ARMA_mu(n = n_val, ar_coeffs = ar_params, ma_coeffs = ma_params,
                     mu_scenario = mu_scenario, k = shift_k_val,
                     innov_dist = innov_dist, contamination_scenario = contamination_scenario, 
                     epsilon = percentage_contamination, gamma = 10)

                     
  # FIXED: The true mean at the end of the series (at u=1) reduced to a scalar
  true_mean_at_end <- (cumsum(ts_data$m_u) / n_val)[n_val]
  
  # Estimate the location parameter
  teta <- get_teta(ts_data$Xt, k = 0.65, c = 4.685)
  
  # Calculate the integrated estimators
  int_teta <- cumsum(teta) / n_val
  int_lin_teta <- int.par.mean(x = ts_data$Xt, teta = teta, k = 0.65, c = 4.685, lag = ceiling(n_val^0.1))
  
  # Get the final value of the estimators (at u=1)
  est_teta_final <- int_teta[n_val]
  est_lin_teta_final <- int_lin_teta[n_val]
  
  # Handle potential non-finite results from estimators
  if (!is.finite(est_teta_final)) est_teta_final <- NA
  if (!is.finite(est_lin_teta_final)) est_lin_teta_final <- NA

  # Calculate error and squared error for bias and MSE calculation (Now safely scalar operations)
  error_teta <- est_teta_final - true_mean_at_end
  error_lin_teta <- est_lin_teta_final - true_mean_at_end
  
  res <- list()
  res$error_teta = error_teta
  res$sq_error_teta = error_teta^2
  res$error_lin_teta = error_lin_teta
  res$sq_error_lin_teta = error_lin_teta^2
  
  # Return a named list of results
  return(res)
}

# --- Setup and run the parallel computation ---
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
clusterExport(cl, c("run_one", "ARMA_mu", "get_teta", "int.par.mean", "tukey_loss_2nd_derivative", "tukey_loss_derivative", "tukey_weight"))

clusterEvalQ(cl, {
  # Prevent thread explosion on the worker nodes
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
  source("DGP.R")
  source("CUSUM.R")
})

# ─────────────────────────────────────────────
# Parallel simulation
# ─────────────────────────────────────────────

cat("Running", nrow(sim_grid), "simulations across", cores, "cores...\n")

sim_results <- foreach(
    row      = iter(sim_grid, by = "row"),
    .combine = rbind
) %dopar% {
  res <- run_one(
    n_val = row$n,
    ar_params = unlist(row$ar_params),
    ma_params = unlist(row$ma_params),
    innov_dist = row$innov_dist,
    shift_k_val = row$shift_k,
    mu_scenario = row$hp_scenario,
    contamination_scenario = row$contamination_scenario,
    percentage_contamination = row$percentage_contamination
  )
  data.frame(
    iteration = row$iteration,
    replication = row$MC_level,                 
    scenario = row$hp_scenario,
    contamination = row$contamination_scenario,
    percentage_contamination = row$percentage_contamination,
    shift_k = row$shift_k,
    innov_dist = row$innov_dist,
    n = row$n,
    ar_params = paste(unlist(row$ar_params), collapse = ","),  
    ma_params = paste(unlist(row$ma_params), collapse = ","),
    bias_teta = res[["error_teta"]],             
    bias_lin_teta = res[["error_lin_teta"]],
    sq_error_teta = res[["sq_error_teta"]],
    sq_error_lin_teta = res[["sq_error_lin_teta"]],
    stringsAsFactors = FALSE
  )
}

stopCluster(cl)

write.csv(sim_results, "sim_results.csv", row.names = FALSE)
cat("Simulation complete. Results saved to sim_results.csv.\n")