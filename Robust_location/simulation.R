# ==============================================================================
# Monte Carlo Simulation Study: Robust Location Estimation
# Comparing Standard Recursive M-Estimator vs Linearized Recursive Estimator
# Supports Welsh and Tukey loss functions, and Cross-Validation for bandwidth selection
# ==============================================================================

suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(parallel)
})

# Source required function files
source("DGP.R")
source("CUSUM.R")

# ------------------------------------------------------------------------------
# 1. Configuration & Parameters
# ------------------------------------------------------------------------------
# You can set parameters here directly OR override them via command line arguments:
#   Rscript simulation.R --use_cv=TRUE --loss=Welsh --reps=1000 --epsilon=0.10
#   Rscript simulation.R --quick --use_cv=TRUE
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
is_quick    <- any(c("--quick", "--test", "-q") %in% args)
has_cv_flag <- any(c("--cv", "--use-cv") %in% args)

parse_arg <- function(arg_name, default_val) {
  match_idx <- grep(paste0("^--", arg_name, "="), args)
  if (length(match_idx) > 0) {
    val_str <- sub(paste0("^--", arg_name, "="), "", args[match_idx[1]])
    return(val_str)
  }
  return(default_val)
}

# --- PRIMARY CONTROLS ---
use_cv_default  <- FALSE  # Set to TRUE to enable Cross-Validation, or use --use_cv=TRUE or --cv
loss_default    <- "Welsh" # "Welsh" or "Tukey"
k_fixed_default <- 0.65   # Bandwidth rate exponent when use_cv is FALSE (e.g. 0.65 -> k = floor(n^0.65))

loss_arg         <- parse_arg("loss", loss_default)
use_cv_arg       <- as.logical(parse_arg("use_cv", as.character(use_cv_default)))
if (has_cv_flag) use_cv_arg <- TRUE
k_fixed_arg      <- as.numeric(parse_arg("k", as.character(k_fixed_default)))
var_scenario_arg <- parse_arg("var_scenario", "i")
shift_arg        <- as.numeric(parse_arg("shift", "0.5"))

if (is_quick) {
  cat("========================================================\n")
  cat("RUNNING IN QUICK / TEST MODE (Reduced Grid & Reps)\n")
  cat("========================================================\n")
  MC.simulations          <- 5
  n_values                <- c(100, 500)
  shift_k_opts            <- c(shift_arg)
  epsilon                 <- c(0.10)
  scenarios.contamination <- c("clean", "AO", "IO")
  innov_dist_opts         <- c("gaussian", "t3")
  mu_scenarios            <- c("H0", "H1", "H2")
} else {
  reps_str        <- parse_arg("reps", Sys.getenv("MC_REPS", "1000"))
  MC.simulations  <- as.integer(reps_str)
  
  n_str           <- parse_arg("n", "100,500,1000,5000")
  n_values        <- as.integer(strsplit(n_str, ",")[[1]])
  
  eps_str         <- parse_arg("epsilon", "0.10")
  epsilon         <- as.numeric(strsplit(eps_str, ",")[[1]])
  
  shift_k_opts    <- c(shift_arg)
  scenarios.contamination <- c("clean", "AO", "IO")
  innov_dist_opts <- c("gaussian", "t3")
  mu_scenarios    <- c("H0", "H1", "H2")
}

# Model ARMA parameters
ar_params <- c(0.2, -0.1) # AR(2)
ma_params <- c(0.2)       # MA(1)

cat(sprintf("Configuration: Loss = %s | Bandwidth Mode = %s (fixed k = %.2f) | Var Scenario = %s | Shift = %.2f | Reps = %d | Epsilon = %s\n",
            loss_arg, if (use_cv_arg) "Cross-Validation (CV)" else "Fixed Rate",
            k_fixed_arg, var_scenario_arg, shift_arg, MC.simulations, paste(epsilon, collapse = ",")))

# ------------------------------------------------------------------------------
# 2. Simulation Grid Generation
# ------------------------------------------------------------------------------

cat("Generating simulation grid...\n")

base_design <- expand.grid(
  n                        = n_values,
  ar_params                = list(ar_params),
  ma_params                = list(ma_params),
  hp_scenario              = mu_scenarios,
  var_scenario             = var_scenario_arg,
  percentage_contamination = epsilon,
  contamination_scenario   = scenarios.contamination,
  innov_dist               = innov_dist_opts,
  shift_k                  = shift_k_opts,
  stringsAsFactors         = FALSE
)

grid_list <- lapply(1:MC.simulations, function(m) {
  df <- base_design
  df$iteration <- m
  df$MC_level  <- MC.simulations
  df
})

sim_grid <- do.call(rbind, grid_list)

# ------------------------------------------------------------------------------
# 3. Single Iteration Worker Function
# ------------------------------------------------------------------------------

run_one <- function(n_val, ar_params, ma_params, innov_dist, shift_k_val,
                    mu_scenario, contamination_scenario, percentage_contamination,
                    var_scenario = "i", loss_type = "Welsh", use_cv = FALSE, k_fixed = 0.65) {

  # 1. Generate the time series data based on specified DGP
  ts_data <- ARMA_mu(
    n                      = n_val,
    ar_coeffs              = ar_params,
    ma_coeffs              = ma_params,
    mu_scenario            = mu_scenario,
    k                      = shift_k_val,
    var_scenario           = var_scenario,
    innov_dist             = innov_dist,
    contamination_scenario = contamination_scenario,
    epsilon                = percentage_contamination,
    gamma                  = 10
  )

  # True integrated mean at u = 1
  true_mean_at_end <- (cumsum(ts_data$m_u) / n_val)[n_val]

  # Bandwidth selection: CV or fixed rate
  if (use_cv) {
    cv_res   <- cv_optimal_bandwidth_location(ts_data$Xt, loss = loss_type)
    k_choice <- cv_res$k_opt
    k_rate   <- cv_res$k_opt_rate
  } else {
    k_choice <- k_fixed
    k_rate   <- if (k_fixed < 1) k_fixed else log(k_fixed) / log(n_val)
  }

  # Decoupling lag L_n = max(1, floor(0.1 * (log(n))^2))
  lag_val <- max(1, floor(0.1 * (log(n_val))^2))

  # 2. Standard recursive M-estimator
  teta <- get_teta(ts_data$Xt, k = k_choice, loss = loss_type)
  int_teta <- cumsum(teta) / n_val

  # 3. Linearized recursive estimator
  int_lin_teta <- int.par.mean(
    x     = ts_data$Xt,
    teta  = teta,
    k     = k_choice,
    lag   = lag_val,
    loss  = loss_type
  )

  # Final values at u = 1
  est_teta_final     <- int_teta[n_val]
  est_lin_teta_final <- int_lin_teta[n_val]

  if (!is.finite(est_teta_final)) est_teta_final <- NA
  if (!is.finite(est_lin_teta_final)) est_lin_teta_final <- NA

  # Calculate error and squared error
  error_teta     <- est_teta_final - true_mean_at_end
  error_lin_teta <- est_lin_teta_final - true_mean_at_end

  list(
    error_teta        = error_teta,
    sq_error_teta     = error_teta^2,
    error_lin_teta    = error_lin_teta,
    sq_error_lin_teta = error_lin_teta^2,
    k_choice          = k_choice,
    k_rate            = k_rate
  )
}

# ------------------------------------------------------------------------------
# 4. Setup and Run Parallel Computation
# ------------------------------------------------------------------------------

set.seed(123, kind = "L'Ecuyer-CMRG")

# HPC / Multi-core detection
cores_str <- Sys.getenv("SLURM_NTASKS")
if (nchar(cores_str) > 0) {
  cores <- as.integer(cores_str)
  cat("SLURM environment detected. Using allocated cores:", cores, "\n")
} else {
  cores <- detectCores()
  cat("Using all available local cores:", cores, "\n")
}

cl <- makeCluster(cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 123)

# Export functions and configuration to workers
clusterExport(cl, c(
  "run_one", "ARMA_mu", "get_teta", "solve_teta", "int.par.mean",
  "cv_optimal_bandwidth_location",
  "Welsh.rho", "Welsh.psi", "Welsh.psi.prime", "Welsh.weight",
  "tukey_weight", "tukey_loss_derivative", "tukey_loss_2nd_derivative",
  "loss_arg", "use_cv_arg", "k_fixed_arg"
))

invisible(clusterEvalQ(cl, {
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
  source("DGP.R")
  source("CUSUM.R")
}))

cat("Running", nrow(sim_grid), "simulation tasks across", cores, "cores...\n")
t_start <- proc.time()

sim_results <- foreach(
  row      = iter(sim_grid, by = "row"),
  .combine = rbind
) %dopar% {
  res <- run_one(
    n_val                    = row$n,
    ar_params                = unlist(row$ar_params),
    ma_params                = unlist(row$ma_params),
    innov_dist               = row$innov_dist,
    shift_k_val              = row$shift_k,
    mu_scenario              = row$hp_scenario,
    contamination_scenario   = row$contamination_scenario,
    percentage_contamination = row$percentage_contamination,
    var_scenario             = row$var_scenario,
    loss_type                = loss_arg,
    use_cv                   = use_cv_arg,
    k_fixed                  = k_fixed_arg
  )

  data.frame(
    iteration                = row$iteration,
    replication              = row$MC_level,
    scenario                 = row$hp_scenario,
    contamination            = row$contamination_scenario,
    percentage_contamination = row$percentage_contamination,
    shift_k                  = row$shift_k,
    var_scenario             = row$var_scenario,
    innov_dist               = row$innov_dist,
    n                        = row$n,
    loss                     = loss_arg,
    bias_teta                = res[["error_teta"]],
    bias_lin_teta            = res[["error_lin_teta"]],
    sq_error_teta            = res[["sq_error_teta"]],
    sq_error_lin_teta        = res[["sq_error_lin_teta"]],
    k_rate                   = res[["k_rate"]],
    stringsAsFactors         = FALSE
  )
}

stopCluster(cl)

elapsed <- proc.time() - t_start
cat(sprintf("Computation finished in %.2f seconds.\n", elapsed[3]))

write.csv(sim_results, "sim_results.csv", row.names = FALSE)
cat("Simulation complete. Results saved to sim_results.csv.\n")