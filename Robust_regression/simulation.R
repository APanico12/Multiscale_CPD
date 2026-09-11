#!/usr/bin/env Rscript
# ==============================================================================
# Monte Carlo Simulation Study: Robust TV Regression Change-Point Detection
# Comparing Proposed Linearized CUSUM, Plug-in CUSUM, and Liu & Zhou (2024) SCB
# ==============================================================================

suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(parallel)
})

# Source required DGP and testing modules
source("DGP.R")
source("CUSUM.R")

# ------------------------------------------------------------------------------
# 1. Configuration & Argument Parsing
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
is_quick <- any(c("--quick", "--test", "-q") %in% args)

parse_arg <- function(arg_name, default_val) {
  match_idx <- grep(paste0("^--", arg_name, "="), args)
  if (length(match_idx) > 0) {
    val_str <- sub(paste0("^--", arg_name, "="), "", args[match_idx[1]])
    return(val_str)
  }
  return(default_val)
}

if (is_quick) {
  cat("========================================================\n")
  cat("RUNNING IN QUICK / TEST MODE (Reduced Grid & Reps)\n")
  cat("========================================================\n")
  MC_reps     <- 2
  B_boot      <- 50
  n_values    <- c(200)
  model_types <- c("LMHC")
  innov_dists <- c("normal")
  hp_list     <- c("H0", "H1")
  delta_shift <- 0.8
  contam_grid <- data.frame(
    contamination = c("Clean", "AO"),
    epsilon       = c(0.00, 0.05),
    stringsAsFactors = FALSE
  )
} else {
  mc_str   <- parse_arg("reps", Sys.getenv("MC_REPS", "1000"))
  MC_reps  <- as.integer(mc_str)
  
  b_str    <- parse_arg("boot", Sys.getenv("B_BOOT", "200"))
  B_boot   <- as.integer(b_str)
  
  n_str    <- parse_arg("n", "500,1000")
  n_str    <- parse_arg("n", "500,1000,5000")
  n_values <- as.integer(strsplit(n_str, ",")[[1]])
  
  model_types <- c("LMHC", "LMAT")
  innov_dists <- c("normal", "t3")
  hp_list     <- c("H0", "H1", "H2")
  delta_shift <- 0.8
  
  # Contamination grid aligning with paper.tex Section 5.3
  contam_grid <- data.frame(
    contamination = c("Clean", "AO", "AO", "IO", "IO", "RO", "RO"),
    epsilon       = c(0.00,    0.05, 0.10, 0.05, 0.10, 0.05, 0.10),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. Build Simulation Design Grid
# ------------------------------------------------------------------------------

cat("Constructing simulation parameter grid...\n")

design_list <- list()
for (m in model_types) {
  for (nv in n_values) {
    for (idist in innov_dists) {
      for (cg in seq_len(nrow(contam_grid))) {
        contam <- contam_grid$contamination[cg]
        eps    <- contam_grid$epsilon[cg]
        for (hp in hp_list) {
          cur_delta <- if (hp == "H0") 0.0 else delta_shift
          design_list[[length(design_list) + 1]] <- data.frame(
            model_type    = m,
            n             = nv,
            innov_dist    = idist,
            contamination = contam,
            epsilon       = eps,
            hp_scenario   = hp,
            delta         = cur_delta,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

base_design <- do.call(rbind, design_list)

# Expand across Monte Carlo replications
full_grid <- do.call(rbind, lapply(1:MC_reps, function(rep_id) {
  df <- base_design
  df$iteration <- rep_id
  df
}))

# Randomize row order to balance core loads
set.seed(42)
full_grid <- full_grid[sample(nrow(full_grid)), ]
rownames(full_grid) <- NULL

cat(sprintf("Total simulation tasks: %d (%d design settings x %d replications)\n",
            nrow(full_grid), nrow(base_design), MC_reps))

# ------------------------------------------------------------------------------
# 3. Single Iteration Worker Function
# ------------------------------------------------------------------------------

run_one_simulation <- function(model_type, n, innov_dist, contamination, epsilon,
                               hp_scenario, delta, B_boot) {
  # 1. Generate data according to design
  dgp_out <- generate_tv_regression_dgp(
    n             = n,
    model_type    = model_type,
    hp_scenario   = hp_scenario,
    delta         = delta,
    innov_dist    = innov_dist,
    contamination = contamination,
    epsilon       = epsilon,
    K             = 10
  )
  
  Y <- dgp_out$Y
  X <- dgp_out$regressors
  d <- ncol(X)
  C_mat <- diag(d)
  
  # 2. Proposed Linearized CUSUM Test
  t0 <- proc.time()
  res_lin <- tryCatch({
    CUSUM.regression(
      Y          = Y,
      X          = X,
      C_mat      = C_mat,
      k          = 0.45,
      lag        = NULL,
      block      = NULL,
      MC         = B_boot,
      linearized = TRUE
    )
  }, error = function(e) {
    list(stat = NA, crit_value = NA, p_value = NA, break_u = NA, reject = NA)
  })
  time_lin <- (proc.time() - t0)[3]
  
  # 3. Standard Plug-in CUSUM Test (Ablation baseline)
  t0 <- proc.time()
  res_plug <- tryCatch({
    CUSUM.regression(
      Y          = Y,
      X          = X,
      C_mat      = C_mat,
      k          = 0.45,
      lag        = NULL,
      block      = NULL,
      MC         = B_boot,
      linearized = FALSE
    )
  }, error = function(e) {
    list(stat = NA, crit_value = NA, p_value = NA, break_u = NA, reject = NA)
  })
  time_plug <- (proc.time() - t0)[3]
  
  # 4. Liu & Zhou (2024) Self-Convolved Bootstrap Test
  t0 <- proc.time()
  res_scb <- tryCatch({
    SCB.regression(
      Y     = Y,
      X     = X,
      C_mat = C_mat,
      b     = 0.20,
      B     = B_boot,
      loss  = "tukey"
    )
  }, error = function(e) {
    list(stat = NA, crit_value = NA, p_value = NA, reject = NA)
  })
  time_scb <- (proc.time() - t0)[3]
  
  return(list(
    # Linearized
    rej_lin    = as.integer(isTRUE(res_lin$reject)),
    stat_lin   = res_lin$stat,
    pval_lin   = res_lin$p_value,
    break_lin  = res_lin$break_u,
    time_lin   = time_lin,
    # Plug-in
    rej_plug   = as.integer(isTRUE(res_plug$reject)),
    stat_plug  = res_plug$stat,
    pval_plug  = res_plug$p_value,
    break_plug = res_plug$break_u,
    time_plug  = time_plug,
    # Liu & Zhou (2024) SCB
    rej_scb    = as.integer(isTRUE(res_scb$reject)),
    stat_scb   = res_scb$stat,
    pval_scb   = res_scb$p_value,
    time_scb   = time_scb
  ))
}

# ------------------------------------------------------------------------------
# 4. Cluster Initialization & Parallel Execution
# ------------------------------------------------------------------------------

cores_str <- Sys.getenv("SLURM_NTASKS")
if (nchar(cores_str) == 0) cores_str <- Sys.getenv("SLURM_CPUS_PER_TASK")

if (nchar(cores_str) > 0) {
  cores <- as.integer(cores_str)
  cat(sprintf("SLURM cluster detected: Allocating %d worker processes.\n", cores))
} else {
  cores <- max(1, parallel::detectCores() - 1)
  cat(sprintf("Local environment detected: Utilizing %d cores.\n", cores))
}

if (is_quick && cores > 4) cores <- 4

cl <- makeCluster(cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = 2026)

clusterExport(cl, "run_one_simulation")

invisible(clusterEvalQ(cl, {
  source("DGP.R")
  source("CUSUM.R")
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
  Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
  Sys.setenv(NUMEXPR_NUM_THREADS = "1")
}))

cat(sprintf("Executing simulation across %d cores...\n", cores))
sim_start_time <- proc.time()

results_list <- foreach(
  row_idx  = seq_len(nrow(full_grid)),
  .combine = rbind,
  .errorhandling = "pass",
  .packages = c()
) %dopar% {
  row_cfg <- full_grid[row_idx, ]
  
  res <- run_one_simulation(
    model_type    = row_cfg$model_type,
    n             = row_cfg$n,
    innov_dist    = row_cfg$innov_dist,
    contamination = row_cfg$contamination,
    epsilon       = row_cfg$epsilon,
    hp_scenario   = row_cfg$hp_scenario,
    delta         = row_cfg$delta,
    B_boot        = B_boot
  )
  
  data.frame(
    iteration     = row_cfg$iteration,
    model_type    = row_cfg$model_type,
    n             = row_cfg$n,
    innov_dist    = row_cfg$innov_dist,
    contamination = row_cfg$contamination,
    epsilon       = row_cfg$epsilon,
    hp_scenario   = row_cfg$hp_scenario,
    delta         = row_cfg$delta,
    # Linearized proposed
    rej_lin       = res$rej_lin,
    stat_lin      = res$stat_lin,
    pval_lin      = res$pval_lin,
    break_lin     = res$break_lin,
    time_lin      = res$time_lin,
    # Plug-in baseline
    rej_plug      = res$rej_plug,
    stat_plug     = res$stat_plug,
    pval_plug     = res$pval_plug,
    break_plug    = res$break_plug,
    time_plug     = res$time_plug,
    # Liu & Zhou (2024) SCB
    rej_scb       = res$rej_scb,
    stat_scb      = res$stat_scb,
    pval_scb      = res$pval_scb,
    time_scb      = res$time_scb,
    stringsAsFactors = FALSE
  )
}

stopCluster(cl)
total_elapsed <- (proc.time() - sim_start_time)[3]
cat(sprintf("\nSimulation completed in %.2f minutes (%.1f seconds).\n",
            total_elapsed / 60, total_elapsed))

# ------------------------------------------------------------------------------
# 5. Persist Raw & Aggregated Results
# ------------------------------------------------------------------------------

raw_csv <- "sim_results_regression.csv"
write.csv(results_list, raw_csv, row.names = FALSE)
cat(sprintf("Raw replication results saved to: %s (%d rows)\n", raw_csv, nrow(results_list)))

# Compute aggregated rejection frequencies (Size under H0, Power under H1/H2)
agg_summary <- aggregate(
  cbind(rej_lin, rej_plug, rej_scb, time_lin, time_plug, time_scb) ~
    model_type + n + innov_dist + contamination + epsilon + hp_scenario,
  data = results_list,
  FUN = function(x) round(mean(x, na.rm = TRUE), 4)
)

names(agg_summary)[names(agg_summary) == "rej_lin"]  <- "rate_proposed_lin"
names(agg_summary)[names(agg_summary) == "rej_plug"] <- "rate_plugin"
names(agg_summary)[names(agg_summary) == "rej_scb"]  <- "rate_liu_scb"

summary_csv <- "sim_summary_regression.csv"
write.csv(agg_summary, summary_csv, row.names = FALSE)
cat(sprintf("Aggregated summary metrics saved to: %s\n", summary_csv))

cat("\n=========================================================================\n")
cat("SIMULATION SUMMARY: Empirical Rejection Frequencies (alpha = 0.05)\n")
cat("=========================================================================\n")
print(agg_summary)
cat("=========================================================================\n")