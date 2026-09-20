# ==============================================================================
# Simulation Study: Robust Change-Point Detection Comparison in Location
# Comparing:
#   1. Our Proposed Linearized CUSUM test with Welsh score function (Our_Welsh)
#   2. Two-sample Hodges-Lehmann test (Hodges_Lehmann via robcp::hl_test)
#   3. Huberized CUSUM test (Huber_CUSUM via robcp::huber_cusum)
#   4. Wilcoxon-Mann-Whitney test (Wilcoxon via robcp::wmw_test)
#
# Design:
#   - Hypotheses: H0 (Size), H1 (Abrupt shift), H2 (Gradual shift)
#   - Sample sizes: n in {200, 500, 1000}
#   - Scenarios: Clean, AO (Additive Outliers), IO (Innovation Outliers)
#   - Innovations: Gaussian and Student-t3
# ==============================================================================

required_pkgs <- c("doParallel", "foreach", "parallel", "robcp")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(sprintf("\n[ERROR] Missing required R package(s): %s\nPlease install them from the login node via:\n  Rscript -e 'install.packages(c(%s), repos=\"https://cloud.r-project.org\")'\n",
               paste(missing_pkgs, collapse = ", "),
               paste(sprintf('"%s"', missing_pkgs), collapse = ", ")))
}

suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(parallel)
  library(robcp)
})

# Source required local functions
source("DGP.R")
source("CUSUM.R")
source("schmidt_test.R")

# ------------------------------------------------------------------------------
# 1. Configuration & Command Line Argument Parsing
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
  MC_reps                 <- 5
  n_values                <- c(200, 500)
  shift_str               <- parse_arg("shift", "0.5")
  shift_mag               <- as.numeric(shift_str)
  epsilon                 <- 0.05
  scenarios_contamination <- c("clean", "AO", "IO", "RO")
  innov_dists             <- c("gaussian", "t3")
  hp_scenarios            <- c("H0", "H1", "H2")
  var_scenarios_str       <- parse_arg("var_scenarios", "i,ii,iii")
  var_scenarios           <- strsplit(var_scenarios_str, ",")[[1]]
  mc_cusum_reps           <- 50
} else {
  reps_str                <- parse_arg("reps", Sys.getenv("MC_REPS", "500"))
  MC_reps                 <- as.integer(reps_str)
  
  n_str                   <- parse_arg("n", "200,500,1000")
  n_values                <- as.integer(strsplit(n_str, ",")[[1]])
  
  eps_str                 <- parse_arg("epsilon", "0.05")
  epsilon                 <- as.numeric(eps_str)
  
  shift_str               <- parse_arg("shift", "0.5")
  shift_mag               <- as.numeric(shift_str)
  
  contam_str              <- parse_arg("contamination", "clean,AO,IO,RO")
  scenarios_contamination <- strsplit(contam_str, ",")[[1]]
  innov_dists             <- c("gaussian", "t3")
  hp_scenarios            <- c("H0", "H1", "H2")
  var_scenarios_str       <- parse_arg("var_scenarios", "i,ii,iii")
  var_scenarios           <- strsplit(var_scenarios_str, ",")[[1]]
  b_str                   <- parse_arg("B", "200")
  mc_cusum_reps           <- as.integer(b_str)
}

# ARMA parameters
ar_params <- c(0.2, -0.1) # AR(2)
ma_params <- c(0.2)       # MA(1)

cat(sprintf("Configuration: Reps = %d | Sample sizes = %s | Shift = %.2f | Epsilon = %.2f | Var Scenarios = %s | B (CUSUM draws) = %d\n",
            MC_reps, paste(n_values, collapse = ","), shift_mag, epsilon, paste(var_scenarios, collapse = ","), mc_cusum_reps))

# ------------------------------------------------------------------------------
# 2. Build Simulation Grid
# ------------------------------------------------------------------------------

cat("Building simulation grid...\n")

base_design <- expand.grid(
  n             = n_values,
  hp_scenario   = hp_scenarios,
  contamination = scenarios_contamination,
  innov_dist    = innov_dists,
  var_scenario  = var_scenarios,
  epsilon       = epsilon,
  shift_k       = shift_mag,
  stringsAsFactors = FALSE
)

grid_list <- lapply(1:MC_reps, function(m) {
  df <- base_design
  df$iteration <- m
  df
})

sim_grid <- do.call(rbind, grid_list)

# Randomize row order to balance core loads across workers
set.seed(42)
sim_grid <- sim_grid[sample(nrow(sim_grid)), ]
rownames(sim_grid) <- NULL

cat(sprintf("Total simulation tasks: %d (%d design settings x %d reps)\n",
            nrow(sim_grid), nrow(base_design), MC_reps))

# ------------------------------------------------------------------------------
# 3. Worker Function: Run All 4 Tests on the Same Data Stream
# ------------------------------------------------------------------------------

run_one_comparison <- function(n_val, hp_scenario, contamination, innov_dist,
                               var_scenario, epsilon_val, shift_val, ar_p, ma_p, b_cusum) {

  # 1. Generate data according to design
  ts_data <- ARMA_mu(
    n                      = n_val,
    ar_coeffs              = ar_p,
    ma_coeffs              = ma_p,
    mu_scenario            = hp_scenario,
    k                      = shift_val,
    var_scenario           = var_scenario,
    innov_dist             = innov_dist,
    contamination_scenario = contamination,
    epsilon                = epsilon_val,
    gamma                  = 10
  )
  
  x <- ts_data$Xt

  # 2. Test 1: Our Proposed Linearized CUSUM test with Welsh score function
  rej_our <- tryCatch({
    res <- CUSUM.mean(x = x, loss = "Welsh", k = 0.45, MC = b_cusum, linearized = TRUE)
    as.integer(res$p_value < 0.05)
  }, error = function(e) NA_integer_)

  # 3. Test 2: Two-sample Hodges-Lehmann test (robcp::hl_test)
  rej_hl <- tryCatch({
    res <- robcp::hl_test(x)
    as.integer(res$p.value < 0.05)
  }, error = function(e) NA_integer_)

  # 4. Test 3: Huberized CUSUM test (robcp::huber_cusum)
  rej_huber <- tryCatch({
    res <- robcp::huber_cusum(x, fun = "HLm")
    as.integer(res$p.value < 0.05)
  }, error = function(e) NA_integer_)

  # 5. Test 4: Wilcoxon-Mann-Whitney test (robcp::wmw_test)
  rej_wmw <- tryCatch({
    res <- robcp::wmw_test(x, h = 1L)
    as.integer(res$p.value < 0.05)
  }, error = function(e) NA_integer_)

  # 6. Test 5: Schmidt (2021) Gini test for heteroscedastic time series
  rej_schmidt <- tryCatch({
    res <- schmidt_test(x, s = 0.7, q = 0.4, c0 = 10, M_psi = 300)
    as.integer(res$p_value < 0.05)
  }, error = function(e) NA_integer_)

  list(
    rej_our     = rej_our,
    rej_hl      = rej_hl,
    rej_huber   = rej_huber,
    rej_wmw     = rej_wmw,
    rej_schmidt = rej_schmidt
  )
}

# ------------------------------------------------------------------------------
# 4. Parallel Setup and Execution
# ------------------------------------------------------------------------------

set.seed(123, kind = "L'Ecuyer-CMRG")

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

clusterExport(cl, c(
  "run_one_comparison", "ARMA_mu", "CUSUM.mean", "get_teta", "solve_teta",
  "int.par.mean", "var.est.mean", "schmidt_test",
  "Welsh.rho", "Welsh.psi", "Welsh.psi.prime", "Welsh.weight",
  "tukey_weight", "tukey_loss_derivative", "tukey_loss_2nd_derivative",
  "ar_params", "ma_params", "mc_cusum_reps"
))

invisible(clusterEvalQ(cl, {
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
  suppressPackageStartupMessages(library(robcp))
  source("DGP.R")
  source("CUSUM.R")
  source("schmidt_test.R")
}))

cat("Running parallel simulation across", cores, "cores...\n")
t_start <- proc.time()

results_list <- foreach(
  row      = iter(sim_grid, by = "row"),
  .combine = rbind
) %dopar% {
  res <- run_one_comparison(
    n_val         = row$n,
    hp_scenario   = row$hp_scenario,
    contamination = row$contamination,
    innov_dist    = row$innov_dist,
    var_scenario  = row$var_scenario,
    epsilon_val   = row$epsilon,
    shift_val     = row$shift_k,
    ar_p          = ar_params,
    ma_p          = ma_params,
    b_cusum       = mc_cusum_reps
  )

  data.frame(
    iteration     = row$iteration,
    hp_scenario   = row$hp_scenario,
    contamination = row$contamination,
    innov_dist    = row$innov_dist,
    var_scenario  = row$var_scenario,
    n             = row$n,
    rej_our       = res$rej_our,
    rej_hl        = res$rej_hl,
    rej_huber     = res$rej_huber,
    rej_wmw       = res$rej_wmw,
    rej_schmidt   = res$rej_schmidt,
    stringsAsFactors = FALSE
  )
}

stopCluster(cl)

elapsed <- (proc.time() - t_start)[3]
cat(sprintf("Simulation finished in %.2f seconds (%.2f minutes).\n", elapsed, elapsed / 60))

# Save raw replication data
write.csv(results_list, "sim_results_cpd_comparison.csv", row.names = FALSE)
cat("Saved raw simulation results to sim_results_cpd_comparison.csv\n")

# Compute and save summary rates
suppressPackageStartupMessages(library(dplyr))
summary_df <- results_list %>%
  group_by(hp_scenario, contamination, innov_dist, var_scenario, n) %>%
  summarise(
    reps         = n(),
    rate_our     = mean(rej_our, na.rm = TRUE),
    rate_hl      = mean(rej_hl, na.rm = TRUE),
    rate_huber   = mean(rej_huber, na.rm = TRUE),
    rate_wmw     = mean(rej_wmw, na.rm = TRUE),
    rate_schmidt = mean(rej_schmidt, na.rm = TRUE),
    .groups      = "drop"
  )
write.csv(summary_df, "sim_summary_cpd_comparison.csv", row.names = FALSE)
cat("Saved summary rates to sim_summary_cpd_comparison.csv\n")


