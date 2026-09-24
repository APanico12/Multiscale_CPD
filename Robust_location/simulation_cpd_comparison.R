# ==============================================================================
# Simulation Study: Robust Change-Point Detection Comparison in Location
# Comparing:
#   1. Our Proposed Linearized CUSUM test with Welsh score function (Our_Welsh)
#   2. Two-sample Hodges-Lehmann test (Hodges_Lehmann via robcp::hl_test)
#   3. Huberized CUSUM test (Huber_CUSUM via robcp::huber_cusum)
#   4. Wilcoxon-Mann-Whitney test (Wilcoxon via robcp::wmw_test)
#   5. Schmidt (2021) Gini test for heteroscedastic time series (Schmidt_Gini)

# Design:
#   - Hypotheses: H0 (Size), H1 (Abrupt shift), H2 (Gradual shift)
#   - Sample sizes: n in {200, 500, 1000}
#   - Scenarios: Clean, AO (Additive Outliers), IO (Innovation Outliers), RO (Replacement Outliers)
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
  MC_reps                 <- as.integer(parse_arg("reps", "5"))
  n_str                   <- parse_arg("n", "200,500")
  n_values                <- as.integer(strsplit(n_str, ",")[[1]])
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

# Bandwidth rate exponent k (default: 0.45 -> k_n = floor(n^0.45))
k_str       <- parse_arg("k", "0.45")
k_bandwidth <- as.numeric(k_str)

# Cross-validation for bandwidth selection (Welsh test)
# Accepts --use_cv=TRUE/FALSE, --cv=TRUE/FALSE, or logical flag
cv_str      <- parse_arg("use_cv", parse_arg("cv", "FALSE"))
use_cv      <- as.logical(toupper(cv_str) %in% c("TRUE", "T", "1", "YES"))

# CV candidate grid (rates or integer window sizes, e.g. "0.45,0.65" or "default")
cv_grid_str <- parse_arg("cv_grid", "0.45,0.65")
cv_grid_vec <- if (tolower(cv_grid_str) %in% c("default", "null", "none")) {
  NULL
} else {
  as.numeric(strsplit(cv_grid_str, ",")[[1]])
}

# Parse active models/tests to run
# Available options: W (Welsh), L or HL (Hodges-Lehmann), H (Huber), WMW (Wilcoxon), S (Schmidt)
# Accepts formats like: --models=L,W,S or --models=c(L,W,S) or --models="W,HL,H,WMW,S"
models_str_raw <- parse_arg("models", "W,L,H,WMW,S")
models_clean   <- gsub("[c()\"' ]", "", models_str_raw)
models_tokens  <- toupper(strsplit(models_clean, ",")[[1]])

# Normalize aliases
active_models <- unique(sapply(models_tokens, function(m) {
  if (m %in% c("L", "HL", "HODGES_LEHMANN", "HODGESLEHMANN")) "L"
  else if (m %in% c("W", "WELSH", "OUR_WELSH")) "W"
  else if (m %in% c("H", "HUBER", "HUBER_CUSUM")) "H"
  else if (m %in% c("WMW", "WILCOXON")) "WMW"
  else if (m %in% c("S", "SCHMIDT", "SCHMIDT_GINI")) "S"
  else m
}))

# ARMA parameters
ar_params <- c(0.2, -0.1) # AR(2)
ma_params <- c(0.2)       # MA(1)

# Output file configuration
raw_out_file     <- parse_arg("output_file", parse_arg("out", parse_arg("output", "sim_results_cpd_comparison.csv")))
summary_out_file <- parse_arg("output_summary", parse_arg("summary_out", "sim_summary_cpd_comparison.csv"))

cv_info_str <- if (use_cv) {
  sprintf("TRUE (Grid: %s)", if (is.null(cv_grid_vec)) "default [N^0.35, N^0.65]" else paste(cv_grid_vec, collapse = ","))
} else {
  "FALSE (Fixed k)"
}

cat(sprintf("Configuration: Reps = %d | Sample sizes = %s | Shift = %.2f | Epsilon = %.2f | Var Scenarios = %s | B (CUSUM draws) = %d | Bandwidth k = %.2f | Use CV = %s | Active Models = %s\n",
            MC_reps, paste(n_values, collapse = ","), shift_mag, epsilon, paste(var_scenarios, collapse = ","), mc_cusum_reps, k_bandwidth, cv_info_str, paste(active_models, collapse = ", ")))
cat(sprintf("Output files: Raw results = %s | Summary table = %s\n", raw_out_file, summary_out_file))

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
# 3. Worker Function: Run Selected Tests on the Same Data Stream
# ------------------------------------------------------------------------------

run_one_comparison <- function(n_val, hp_scenario, contamination, innov_dist,
                               var_scenario, epsilon_val, shift_val, ar_p, ma_p, b_cusum, k_val = 0.45,
                               active_models = c("W", "L", "H", "WMW", "S"),
                               use_cv = FALSE, cv_grid = c(0.45, 0.65)) {

  # 1. Generate data according to design
  ts_data <- ARMA_mu(
    n                      = n_val,
    ar_coeffs              = ar_p,
    ma_coeffs              = ma_p,
    mu_scenario            = hp_scenario,
    delta                  = shift_val,
    var_scenario           = var_scenario,
    innov_dist             = innov_dist,
    contamination_scenario = contamination,
    epsilon                = epsilon_val,
    gamma                  = 10
  )
  
  x <- ts_data$Xt

  # 2. Test 1: Our Proposed Linearized CUSUM test with Welsh score function (W)
  rej_our <- if ("W" %in% active_models) {
    tryCatch({
      if (isTRUE(use_cv)) {
        cv_res    <- cv_optimal_bandwidth_location(x = x, loss = "Welsh", k_grid = cv_grid)
        k_eval    <- cv_res$k_opt
        teta_eval <- if (!is.null(cv_res$teta_opt)) cv_res$teta_opt else NULL
      } else {
        k_eval    <- k_val
        teta_eval <- NULL
      }
      res <- CUSUM.mean(x = x, teta = teta_eval, loss = "Welsh", k = k_eval, MC = b_cusum, linearized = TRUE)
      as.integer(res$p_value < 0.05)
    }, error = function(e) NA_integer_)
  } else NA_integer_

  # 3. Test 2: Two-sample Hodges-Lehmann test (robcp::hl_test) (L / HL)
  rej_hl <- if (any(c("L", "HL") %in% active_models)) {
    tryCatch({
      res <- robcp::hl_test(x)
      as.integer(res$p.value < 0.05)
    }, error = function(e) NA_integer_)
  } else NA_integer_

  # 4. Test 3: Huberized CUSUM test (robcp::huber_cusum) (H)
  rej_huber <- if ("H" %in% active_models) {
    tryCatch({
      res <- robcp::huber_cusum(x, fun = "HLm")
      as.integer(res$p.value < 0.05)
    }, error = function(e) NA_integer_)
  } else NA_integer_

  # 5. Test 4: Wilcoxon-Mann-Whitney test (robcp::wmw_test) (WMW)
  rej_wmw <- if ("WMW" %in% active_models) {
    tryCatch({
      res <- robcp::wmw_test(x, h = 1L)
      as.integer(res$p.value < 0.05)
    }, error = function(e) NA_integer_)
  } else NA_integer_

  # 6. Test 5: Schmidt (2021) Gini test for heteroscedastic time series (S)
  rej_schmidt <- if ("S" %in% active_models) {
    tryCatch({
      res <- schmidt_test(x, s = 0.7, q = 0.4, c0 = 10, M_psi = 300)
      as.integer(res$p_value < 0.05)
    }, error = function(e) NA_integer_)
  } else NA_integer_

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
  "int.par.mean", "var.est.mean", "schmidt_test", "cv_optimal_bandwidth_location",
  "Welsh.rho", "Welsh.psi", "Welsh.psi.prime", "Welsh.weight",
  "tukey_weight", "tukey_loss_derivative", "tukey_loss_2nd_derivative",
  "ar_params", "ma_params", "mc_cusum_reps", "k_bandwidth", "active_models",
  "use_cv", "cv_grid_vec"
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
    b_cusum       = mc_cusum_reps,
    k_val         = k_bandwidth,
    active_models = active_models,
    use_cv        = use_cv,
    cv_grid       = cv_grid_vec
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
write.csv(results_list, raw_out_file, row.names = FALSE)
cat(sprintf("Saved raw simulation results to %s\n", raw_out_file))

# Compute and save summary rates
suppressPackageStartupMessages(library(dplyr))
safe_mean <- function(v) if (all(is.na(v))) NA_real_ else round(mean(v, na.rm = TRUE), 4)

summary_df <- results_list %>%
  group_by(hp_scenario, contamination, innov_dist, var_scenario, n) %>%
  summarise(
    reps         = n(),
    rate_our     = safe_mean(rej_our),
    rate_hl      = safe_mean(rej_hl),
    rate_huber   = safe_mean(rej_huber),
    rate_wmw     = safe_mean(rej_wmw),
    rate_schmidt = safe_mean(rej_schmidt),
    .groups      = "drop"
  )
write.csv(summary_df, summary_out_file, row.names = FALSE)
cat(sprintf("Saved summary rates to %s\n", summary_out_file))


