# Demo power run
source("power_analysis.R")

cat("Running demo power study...\n")
res_df <- run_power_analysis(
  delta_grid = c(0.0, 0.25, 0.50, 0.75, 1.00, 1.25),
  n_values   = c(500 , 1000, 5000),
  MC_reps    = 1000,
  B_boot     = 500,
  C_mat      = c(1, 0), # monitoring beta_1
  n_cores    = max(1, parallel::detectCores() - 1),
  seed       = 42
)

#save results to CSV
write.csv(res_df, file = "power_analysis_results.csv", row.names = FALSE)

