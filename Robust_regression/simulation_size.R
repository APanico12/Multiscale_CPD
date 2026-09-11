#!/usr/bin/env Rscript
# ==============================================================================
# File: simulation_size.R
# Description: Monte Carlo Simulation Study for Empirical Size (under H0).
#              Generates results for Table 2 across Model I (LMHC) and Model II (LMAT),
#              evaluating sample sizes n in {500, 1000, 5000} under various contamination
#              scenarios (Clean, AO, IO, RO) and innovation distributions (normal, t3).
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(foreach)
})

# Source required modules
source("DGP.R")
source("CUSUM.R")
source("read_sim_res_size.R")

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
  cat("======================================================================\n")
  cat(" RUNNING SIZE SIMULATION IN QUICK / TEST MODE\n")
  cat("======================================================================\n")
  MC_reps     <- 2
  B_boot      <- 30
  n_values    <- c(150)
  model_types <- c("LMHC")
  innov_dists <- c("normal")
  loss_choice <- "Welsh"
  contam_grid <- data.frame(
    contamination = c("Clean", "AO"),
    epsilon       = c(0.00, 0.05),
    stringsAsFactors = FALSE
  )
} else {
  mc_str      <- parse_arg("reps", Sys.getenv("MC_REPS", "1000"))
  MC_reps     <- as.integer(mc_str)
  
  b_str       <- parse_arg("boot", Sys.getenv("B_BOOT", "100"))
  B_boot      <- as.integer(b_str)
  
  n_str       <- parse_arg("n", Sys.getenv("N_VALS", "500,1000,5000"))
  n_values    <- as.integer(strsplit(n_str, ",")[[1]])
  
  mod_str     <- parse_arg("models", "LMHC,LMAT")
  model_types <- strsplit(mod_str, ",")[[1]]
  
  innov_str   <- parse_arg("innov", "normal,t3")
  innov_dists <- strsplit(innov_str, ",")[[1]]
  
  loss_choice <- parse_arg("loss", "Welsh")
  
  # Contamination grid aligning with Table 2 in paper.tex
  contam_grid <- data.frame(
    contamination = c("Clean", "AO", "AO", "IO", "IO", "RO", "RO"),
    epsilon       = c(0.00,    0.05, 0.10, 0.05, 0.10, 0.05, 0.10),
    stringsAsFactors = FALSE
  )
}

use_cv      <- any(c("--cv", "--use_cv") %in% args) || (tolower(parse_arg("bandwidth", "0.45")) == "cv")
k_param     <- if (use_cv) "CV" else as.numeric(parse_arg("bandwidth", "0.45"))

# ------------------------------------------------------------------------------
# 2. Build Simulation Design Grid
# ------------------------------------------------------------------------------

cat("Constructing size simulation parameter grid (H0: no structural break)...\n")

design_list <- list()
for (m in model_types) {
  for (nv in n_values) {
    for (idist in innov_dists) {
      for (cg in seq_len(nrow(contam_grid))) {
        contam <- contam_grid$contamination[cg]
        eps    <- contam_grid$epsilon[cg]
        design_list[[length(design_list) + 1]] <- data.frame(
          model_type    = m,
          n             = nv,
          innov_dist    = idist,
          contamination = contam,
          epsilon       = eps,
          hp_scenario   = "H0",
          delta         = 0.0,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

base_design <- do.call(rbind, design_list)

full_grid <- do.call(rbind, lapply(1:MC_reps, function(rep_id) {
  df <- base_design
  df$iteration <- rep_id
  df
}))

set.seed(2026)
full_grid <- full_grid[sample(nrow(full_grid)), ]
rownames(full_grid) <- NULL

cat(sprintf("Total size simulation tasks: %d (%d design points x %d replications)\n",
            nrow(full_grid), nrow(base_design), MC_reps))
cat(sprintf("Bootstrap iterations B: %d | Loss: %s\n", B_boot, loss_choice))
cat(sprintf("Bootstrap iterations B: %d | Loss: %s | Bandwidth k: %s\n", 
            B_boot, loss_choice, as.character(k_param)))

# ------------------------------------------------------------------------------
# 3. Worker Function
# ------------------------------------------------------------------------------

run_one_size_sim <- function(model_type, n, innov_dist, contamination, epsilon,
                             B_boot, loss, rep_seed) {
                             B_boot, loss, rep_seed, k_val = 0.45) {
  dgp_out <- generate_tv_regression_dgp(
    n             = n,
    model_type    = model_type,
    hp_scenario   = "H0",
    delta         = 0.0,
    innov_dist    = innov_dist,
    contamination = contamination,
    epsilon       = epsilon,
    K             = 20,
    ro_multiplier = 4.0,
    kappa         = 0.50,
    seed          = rep_seed
  )
  
  Y <- dgp_out$Y
  X <- dgp_out$regressors
  d <- ncol(X)
  C_mat <- diag(d)
  
  # 1. Proposed Linearized CUSUM Test
  res_lin <- tryCatch({
    CUSUM.regression(
      Y          = Y,
      X          = X,
      C_mat      = C_mat,
      k          = 0.45,
      k          = k_val,
      lag        = NULL,
      block      = NULL,
      B          = B_boot,
      loss       = loss,
      linearized = TRUE,
      plotting   = FALSE
    )
  }, error = function(e) list(stat = NA, crit_value = NA, p_value = NA, reject = NA))
  
  # 2. Classical L2 / OLS CUSUM Test
  res_l2 <- tryCatch({
    CUSUM.regression(
      Y          = Y,
      X          = X,
      C_mat      = C_mat,
      k          = k_val,
      lag        = NULL,
      block      = NULL,
      B          = B_boot,
      loss       = "L2",
      linearized = TRUE,
      plotting   = FALSE
    )
  }, error = function(e) list(stat = NA, crit_value = NA, p_value = NA, reject = NA))
  
  return(list(
    rej_lin  = as.integer(isTRUE(res_lin$reject_95)),
    stat_lin = res_lin$test_stat,
    pval_lin = res_lin$p_value,
    rej_l2   = as.integer(isTRUE(res_l2$reject_95)),
    stat_l2  = res_l2$test_stat,
    pval_l2  = res_l2$p_value
  ))
}

# ------------------------------------------------------------------------------
# 4. Cluster Initialization & Parallel Execution
# ------------------------------------------------------------------------------

cores_str <- Sys.getenv("SLURM_NTASKS")
if (nchar(cores_str) == 0) cores_str <- Sys.getenv("SLURM_CPUS_PER_TASK")

if (nchar(cores_str) > 0) {
  cores <- as.integer(cores_str)
  cat(sprintf("SLURM cluster detected: Allocating %d worker cores.\n", cores))
} else {
  cores <- max(1, parallel::detectCores() - 2)
  cat(sprintf("Local environment detected: Utilizing %d cores.\n", cores))
}

if (is_quick && cores > 4) cores <- 4

cl <- makeCluster(cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = 202609)
clusterExport(cl, c("run_one_size_sim", "loss_choice", "B_boot"))

invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(zoo)
    library(robustbase)
  })
  source("DGP.R")
  source("CUSUM.R")
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
}))

cat(sprintf("Executing SIZE simulation across %d cores...\n", cores))
sim_start_time <- proc.time()

results_df <- foreach(
  row_idx  = seq_len(nrow(full_grid)),
  .combine = rbind,
  .errorhandling = "pass"
) %dopar% {
  row_cfg <- full_grid[row_idx, ]
  r_seed  <- 202600 + row_idx
  
  res <- run_one_size_sim(
    model_type    = row_cfg$model_type,
    n             = row_cfg$n,
    innov_dist    = row_cfg$innov_dist,
    contamination = row_cfg$contamination,
    epsilon       = row_cfg$epsilon,
    B_boot        = B_boot,
    loss          = loss_choice,
    rep_seed      = r_seed
    rep_seed      = r_seed,
    k_val         = k_param
  )
  
  data.frame(
    iteration     = row_cfg$iteration,
    model_type    = row_cfg$model_type,
    n             = row_cfg$n,
    innov_dist    = row_cfg$innov_dist,
    contamination = row_cfg$contamination,
    epsilon       = row_cfg$epsilon,
    hp_scenario   = "H0",
    rej_lin       = res$rej_lin,
    stat_lin      = res$stat_lin,
    pval_lin      = res$pval_lin,
    rej_l2        = res$rej_l2,
    stat_l2       = res$stat_l2,
    pval_l2       = res$pval_l2,
    stringsAsFactors = FALSE
  )
}

stopCluster(cl)
total_elapsed <- (proc.time() - sim_start_time)[3]
cat(sprintf("\nSIZE Simulation completed in %.2f minutes (%.1f seconds).\n",
            total_elapsed / 60, total_elapsed))

# ------------------------------------------------------------------------------
# 5. Persist Raw Results & Generate LaTeX Table 2
# ------------------------------------------------------------------------------

raw_csv <- "sim_results_size.csv"
write.csv(results_df, raw_csv, row.names = FALSE)
cat(sprintf("Raw replication results saved to: %s (%d rows)\n", raw_csv, nrow(results_df)))

# Compute aggregated rejection frequencies (Empirical Size)
summary_size <- aggregate(
  cbind(rej_lin, rej_l2) ~ model_type + n + innov_dist + contamination + epsilon,
  data = results_df,
  FUN = function(x) round(mean(x, na.rm = TRUE), 4)
)

names(summary_size)[names(summary_size) == "rej_lin"] <- "rate_linearized"
names(summary_size)[names(summary_size) == "rej_l2"]  <- "rate_l2"

summary_csv <- "sim_summary_size.csv"
write.csv(summary_size, summary_csv, row.names = FALSE)
cat(sprintf("Aggregated empirical size summary saved to: %s\n", summary_csv))

# Generate LaTeX Table 2
tex_file <- "table_size_output.tex"
generate_latex_size_table(summary_size, tex_file)

cat("\n=========================================================================\n")
cat("SIZE SIMULATION SUMMARY: Empirical Size Rejection Frequencies (alpha = 0.05)\n")
cat("=========================================================================\n")
print(summary_size)
cat("=========================================================================\n")

