#!/usr/bin/env Rscript
# ==============================================================================
# File: simulation_power.R
# Description: Monte Carlo Simulation Study for Empirical Power under H1 (Abrupt)
#              and H2 (Gradual) across a grid of shift magnitudes (delta).
#              Evaluates Model I (LMHC) and Model II (LMAT) for n in {1000, 5000, 10000}.
#              Evaluates Model I (LMHC) for n in {500, 5000, 10000}.
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(foreach)
})

# Source required modules
source("DGP.R")
source("CUSUM.R")

log_msg <- function(...) {
  cat(...)
  flush.console()
}

# ------------------------------------------------------------------------------
# 1. Configuration & Argument Parsing
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
is_quick <- any(c("--quick", "--test", "-q") %in% args)
is_fresh <- any(c("--fresh", "--reset") %in% args)

parse_arg <- function(arg_name, default_val) {
  match_idx <- grep(paste0("^--", arg_name, "="), args)
  if (length(match_idx) > 0) {
    val_str <- sub(paste0("^--", arg_name, "="), "", args[match_idx[1]])
    return(val_str)
  }
  return(default_val)
}

if (is_quick) {
  log_msg("======================================================================\n")
  log_msg(" RUNNING POWER SIMULATION IN QUICK / TEST MODE\n")
  log_msg("======================================================================\n")
  MC_reps     <- 1
  B_boot      <- 20
  n_values    <- c(150)
  model_types <- c("LMHC")
  innov_dists <- c("normal", "t3")
  hp_list     <- c("H1", "H2")
  delta_grid  <- c(0.0, 0.50, 1.00)
  loss_choice <- "Welsh"
  contam_grid <- data.frame(
    contamination = c("Clean", "AO", "IO", "RO"),
    epsilon       = c(0.00,    0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )
} else {
  mc_str      <- parse_arg("reps", Sys.getenv("MC_REPS", "1000"))
  MC_reps     <- as.integer(mc_str)
  
  b_str       <- parse_arg("boot", Sys.getenv("B_BOOT", "100"))
  B_boot      <- as.integer(b_str)
  
  n_str       <- parse_arg("n", Sys.getenv("N_VALS", "500,1000,5000"))
  n_values    <- as.integer(strsplit(n_str, ",")[[1]])
  
  delta_str   <- parse_arg("delta", "0.0,0.25,0.50,0.75,1.00,1.25")
  delta_grid  <- as.numeric(strsplit(delta_str, ",")[[1]])
  
  hyp_str     <- parse_arg("hyp", "H1,H2")
  hp_list     <- strsplit(hyp_str, ",")[[1]]
  
  mod_str     <- parse_arg("models", "LMHC")
  model_types <- strsplit(mod_str, ",")[[1]]
  
  innov_str   <- parse_arg("innov", Sys.getenv("INNOV_DISTS", "normal,t3"))
  innov_str   <- gsub("t_3", "t3", innov_str, ignore.case = TRUE)
  innov_dists <- tolower(trimws(strsplit(innov_str, ",")[[1]]))
  innov_dists <- ifelse(innov_dists %in% c("t_3", "t-3"), "t3", innov_dists)
  
  loss_choice <- parse_arg("loss", "Welsh")
  
  # Default to all contamination regimes (Clean + 5% AO/IO/RO)
  contam_mode <- tolower(parse_arg("contam", "all"))
  if (contam_mode == "clean") {
    contam_grid <- data.frame(
      contamination = c("Clean"),
      epsilon       = c(0.00),
      stringsAsFactors = FALSE
    )
  } else if (contam_mode == "ao") {
    contam_grid <- data.frame(
      contamination = c("Clean", "AO"),
      epsilon       = c(0.00, 0.05),
      stringsAsFactors = FALSE
    )
  } else {
    contam_grid <- data.frame(
      contamination = c("Clean", "AO", "IO", "RO"),
      epsilon       = c(0.00,    0.05, 0.05, 0.05),
      stringsAsFactors = FALSE
    )
  }
}

use_cv      <- any(c("--cv", "--use_cv") %in% args) || (tolower(parse_arg("bandwidth", "0.45")) == "cv")
k_param     <- if (use_cv) "CV" else as.numeric(parse_arg("bandwidth", "0.45"))

# ------------------------------------------------------------------------------
# 2. Build Simulation Design Grid
# ------------------------------------------------------------------------------

cat("Constructing power simulation parameter grid...\n")
cat(sprintf("Delta Grid    : %s\n", paste(delta_grid, collapse = ", ")))
cat(sprintf("Hypotheses    : %s\n", paste(hp_list, collapse = ", ")))
cat(sprintf("Models        : %s\n", paste(model_types, collapse = ", ")))
cat(sprintf("Sample Sizes n: %s\n", paste(n_values, collapse = ", ")))

design_list <- list()
for (m in model_types) {
  for (nv in n_values) {
    for (hp in hp_list) {
      for (d_val in delta_grid) {
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
              hp_scenario   = hp,
              delta         = d_val,
              stringsAsFactors = FALSE
            )
          }
        }
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

set.seed(202620)
full_grid <- full_grid[sample(nrow(full_grid)), ]
rownames(full_grid) <- NULL

cat(sprintf("Total power simulation tasks: %d (%d design points x %d replications)\n",
            nrow(full_grid), nrow(base_design), MC_reps))
cat(sprintf("Bootstrap iterations B: %d | Loss: %s | Bandwidth k: %s\n", 
            B_boot, loss_choice, as.character(k_param)))

# ------------------------------------------------------------------------------
# 3. Worker Function
# ------------------------------------------------------------------------------

run_one_power_sim <- function(model_type, n, innov_dist, contamination, epsilon,
                              hp_scenario, delta, B_boot, loss, rep_seed, k_val = 0.45) {
  # When delta = 0, DGP behaves as H0
  cur_hp <- if (abs(delta) < 1e-6) "H0" else hp_scenario
  
  dgp_out <- generate_tv_regression_dgp(
    n             = n,
    model_type    = model_type,
    hp_scenario   = cur_hp,
    delta         = delta,
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
      k          = k_val,
      lag        = NULL,
      block      = NULL,
      B          = B_boot,
      loss       = loss,
      linearized = TRUE,
      plotting   = FALSE
    )
  }, error = function(e) list(test_stat = NA_real_, p_value = NA_real_, break_u = NA_real_, reject_95 = FALSE))
  
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
  }, error = function(e) list(test_stat = NA_real_, p_value = NA_real_, break_u = NA_real_, reject_95 = FALSE))
  
  return(list(
    rej_lin   = as.integer(isTRUE(res_lin$reject_95)),
    stat_lin  = if (is.null(res_lin$test_stat) || length(res_lin$test_stat) == 0) NA_real_ else res_lin$test_stat,
    pval_lin  = if (is.null(res_lin$p_value) || length(res_lin$p_value) == 0) NA_real_ else res_lin$p_value,
    break_lin = if (is.null(res_lin$break_u) || length(res_lin$break_u) == 0) NA_real_ else as.numeric(res_lin$break_u),
    rej_l2    = as.integer(isTRUE(res_l2$reject_95)),
    stat_l2   = if (is.null(res_l2$test_stat) || length(res_l2$test_stat) == 0) NA_real_ else res_l2$test_stat,
    pval_l2   = if (is.null(res_l2$p_value) || length(res_l2$p_value) == 0) NA_real_ else res_l2$p_value,
    break_l2  = if (is.null(res_l2$break_u) || length(res_l2$break_u) == 0) NA_real_ else as.numeric(res_l2$break_u)
  ))
}

# ------------------------------------------------------------------------------
# 4. Cluster Initialization & Parallel Execution
# ------------------------------------------------------------------------------

cores_str <- Sys.getenv("SLURM_NTASKS")
if (nchar(cores_str) == 0) cores_str <- Sys.getenv("SLURM_CPUS_PER_TASK")

if (nchar(cores_str) > 0) {
  cores <- as.integer(cores_str)
  log_msg(sprintf("SLURM cluster detected: Allocating %d worker cores.\n", cores))
} else {
  cores <- max(1, parallel::detectCores() - 2)
  log_msg(sprintf("Local environment detected: Utilizing %d cores.\n", cores))
}

if (is_quick && cores > 4) cores <- 4

cl <- makeCluster(cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = 202630)
clusterExport(cl, c("run_one_power_sim", "loss_choice", "B_boot", "k_param"))

invisible(clusterEvalQ(cl, {
  source("DGP.R")
  source("CUSUM.R")
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
}))

log_msg(sprintf("Executing POWER simulation across %d cores...\n", cores))

raw_csv <- "sim_results_power.csv"

# Handle fresh start flag (--fresh or --reset)
if (is_fresh && file.exists(raw_csv)) {
  file.remove(raw_csv)
  log_msg("Fresh start flag detected: Removed existing sim_results_power.csv.\n")
}

grid_keys <- paste(full_grid$iteration, full_grid$model_type, full_grid$n,
                   full_grid$innov_dist, full_grid$contamination, full_grid$epsilon,
                   full_grid$hp_scenario, full_grid$delta, sep = "_")

if (file.exists(raw_csv) && file.info(raw_csv)$size > 50) {
  existing_res <- tryCatch(read.csv(raw_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  req_cols <- c("iteration", "model_type", "n", "innov_dist", "contamination", "epsilon", "hp_scenario", "delta")
  if (!is.null(existing_res) && nrow(existing_res) > 0 && all(req_cols %in% names(existing_res))) {
    done_keys <- paste(existing_res$iteration, existing_res$model_type, existing_res$n,
                       existing_res$innov_dist, existing_res$contamination, existing_res$epsilon,
                       existing_res$hp_scenario, existing_res$delta, sep = "_")
    remaining_idx <- which(!grid_keys %in% done_keys)
    log_msg(sprintf(">>> Checkpoint found: %d tasks already finished in %s.\n",
                    nrow(existing_res), raw_csv))
    log_msg(sprintf(">>> Resuming with %d remaining tasks (out of %d total).\n",
                    length(remaining_idx), nrow(full_grid)))
    tasks_to_run <- full_grid[remaining_idx, ]
  } else {
    tasks_to_run <- full_grid
    if (file.exists(raw_csv)) file.remove(raw_csv)
  }
} else {
  tasks_to_run <- full_grid
}

# If CSV does not exist, initialize with header
if (!file.exists(raw_csv)) {
  header_df <- data.frame(
    iteration     = integer(),
    model_type    = character(),
    n             = integer(),
    innov_dist    = character(),
    contamination = character(),
    epsilon       = numeric(),
    hp_scenario   = character(),
    delta         = numeric(),
    rej_lin       = integer(),
    stat_lin      = numeric(),
    pval_lin      = numeric(),
    break_lin     = numeric(),
    rej_l2        = integer(),
    stat_l2       = numeric(),
    pval_l2       = numeric(),
    break_l2      = numeric(),
    stringsAsFactors = FALSE
  )
  write.csv(header_df, raw_csv, row.names = FALSE)
}

total_to_run <- nrow(tasks_to_run)
sim_start_time <- proc.time()

if (total_to_run > 0) {
  chunk_size <- max(min(cores * 5, 250), min(total_to_run, 50))
  num_chunks <- ceiling(total_to_run / chunk_size)
  
  log_msg(sprintf("Batching %d tasks into %d chunks of up to %d tasks...\n",
                  total_to_run, num_chunks, chunk_size))
  
  completed_in_run <- 0
  
  for (chk in seq_len(num_chunks)) {
    idx_start <- (chk - 1) * chunk_size + 1
    idx_end   <- min(chk * chunk_size, total_to_run)
    chunk_indices <- idx_start:idx_end
    
    batch_res <- foreach(
      row_idx  = chunk_indices,
      .combine = rbind,
      .errorhandling = "pass"
    ) %dopar% {
      row_cfg <- tasks_to_run[row_idx, ]
      r_seed  <- 202640 + row_idx + chk * 100000L
      
      res <- run_one_power_sim(
        model_type    = row_cfg$model_type,
        n             = row_cfg$n,
        innov_dist    = row_cfg$innov_dist,
        contamination = row_cfg$contamination,
        epsilon       = row_cfg$epsilon,
        hp_scenario   = row_cfg$hp_scenario,
        delta         = row_cfg$delta,
        B_boot        = B_boot,
        loss          = loss_choice,
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
        hp_scenario   = row_cfg$hp_scenario,
        delta         = row_cfg$delta,
        rej_lin       = res$rej_lin,
        stat_lin      = res$stat_lin,
        pval_lin      = res$pval_lin,
        break_lin     = res$break_lin,
        rej_l2        = res$rej_l2,
        stat_l2       = res$stat_l2,
        pval_l2       = res$pval_l2,
        break_l2      = res$break_l2,
        stringsAsFactors = FALSE
      )
    }
    
    if (is.data.frame(batch_res) && nrow(batch_res) > 0) {
      write.table(batch_res, raw_csv, append = TRUE, sep = ",",
                  col.names = FALSE, row.names = FALSE, qmethod = "double")
    }
    
    completed_in_run <- completed_in_run + length(chunk_indices)
    elapsed_sec <- (proc.time() - sim_start_time)[3]
    rate <- completed_in_run / max(elapsed_sec, 0.1)
    eta_sec <- (total_to_run - completed_in_run) / max(rate, 1e-4)
    overall_done <- nrow(full_grid) - total_to_run + completed_in_run
    overall_pct <- (overall_done / nrow(full_grid)) * 100
    
    log_msg(sprintf("[%s] Chunk %d/%d done | Completed: %d/%d (%.1f%% overall) | Rate: %.2f tasks/s | Elapsed: %.1f min | ETA: %.1f min\n",
                    format(Sys.time(), "%H:%M:%S"), chk, num_chunks,
                    overall_done, nrow(full_grid), overall_pct, rate,
                    elapsed_sec / 60, eta_sec / 60))
  }
} else {
  log_msg("All tasks already completed! Proceeding to summarize and plot.\n")
}

stopCluster(cl)
total_elapsed <- (proc.time() - sim_start_time)[3]
log_msg(sprintf("\nPOWER Simulation session completed in %.2f minutes (%.1f seconds).\n",
                total_elapsed / 60, total_elapsed))

# ------------------------------------------------------------------------------
# 5. Persist Raw Results, LaTeX Table, and Power Curves Plots
# ------------------------------------------------------------------------------

# Load all accumulated records from disk
results_df <- read.csv(raw_csv, stringsAsFactors = FALSE)
log_msg(sprintf("Raw replication results loaded from: %s (%d rows)\n", raw_csv, nrow(results_df)))

# Compute aggregated empirical power and mean break localization
summary_power <- aggregate(
  cbind(rej_lin, rej_l2, break_lin, break_l2) ~
    model_type + n + innov_dist + contamination + epsilon + hp_scenario + delta,
  data = results_df,
  FUN = function(x) round(mean(x, na.rm = TRUE), 4)
)

names(summary_power)[names(summary_power) == "rej_lin"]   <- "rate_linearized"
names(summary_power)[names(summary_power) == "rej_l2"]    <- "rate_l2"
names(summary_power)[names(summary_power) == "break_lin"] <- "mean_loc_linearized"
names(summary_power)[names(summary_power) == "break_l2"]  <- "mean_loc_l2"

summary_csv <- "sim_summary_power.csv"
write.csv(summary_power, summary_csv, row.names = FALSE)
log_msg(sprintf("Aggregated empirical power summary saved to: %s\n", summary_csv))

log_msg("\n=========================================================================\n")
log_msg("POWER SIMULATION COMPLETE: Raw replications and summary CSV saved.\n")
log_msg("=========================================================================\n")
