# ==============================================================================
# Empirical Application: 3:2:1 Crack Spread Time Series
# Structural Break Detection via Robust Linearized Wild Binary Segmentation (WBS)
# ==============================================================================

# 1. Load implementations
source("CUSUM.R")
source("RobWBSV2.R")

# 2. Read dataset and extract crack spread & dates
df <- read.csv("crack_spread_dataset.csv")
x <- df$crack_spread
dates <- as.Date(df$date)

# 3. Run Linearized Robust WBS with Cross-Validation and sSIC selection
res <- WBS.mean(
  x              = log(log(x)^2),
  dates          = dates,
  M              = 5000,
  Kmax           = 10,
  selection      = "sSIC",
  loss           = "Welsh",
  c              = 2.985,
  linearized     = TRUE,
  use_cv         = TRUE,
  cv_grid        = c(0.45, 0.50, 0.65),
  segment_level  = "robust",
  flag_outliers  = TRUE,
  outlier_thresh = 4.0,
  plot           = FALSE,
  main           = "Robust Linearized WBS (3:2:1 Crack Spread)",
  ylab           = expression(bold("Crack Spread ($/bbl)"))
)

# 4. Save high-resolution publication figures (Vector PDF & 300 DPI PNG)
plot_wbs(
  res,
  dates         = dates,
  save_files    = TRUE,
  output_prefix = "crack_spread_wbs_results",
  main          = "",
  show_legend    = FALSE,
  ylab          = expression(bold("Crack Spread ($/bbl)"))
)

# 5. Display bandwidth selection & CV summary
cat("\n========================================================\n")
cat("       BANDWIDTH SELECTION & CROSS-VALIDATION           \n")
cat("========================================================\n")
if (!is.null(res$cv_info)) {
  cat(sprintf("Selected optimal bandwidth k* = %d (rate = %.4f)\n", 
              res$cv_info$k_opt, res$cv_info$k_opt_rate))
  cat(sprintf("Decoupling lag L_n            = %d\n", res$cv_info$lag))
  cat("Predictive out-of-sample CV losses:\n")
  print(res$cv_info$cv_losses)
} else {
  cat(sprintf("Fixed bandwidth k = %s\n", as.character(res$bandwidth_k)))
}

# 6. Display detected break dates and sSIC table
cat("\n========================================================\n")
cat("          DETECTED STRUCTURAL BREAKS SUMMARY            \n")
cat("========================================================\n")
if (length(res$cpts) > 0) {
  break_summary <- data.frame(
    Break_ID   = seq_along(res$cpts),
    Index      = res$cpts,
    Date       = format(res$cpts_dates, "%Y-%m-%d"),
    Pre_Level  = round(res$segment_levels[1:length(res$cpts)], 3),
    Post_Level = round(res$segment_levels[2:(length(res$cpts) + 1)], 3),
    Shift      = round(res$segment_levels[2:(length(res$cpts) + 1)] - res$segment_levels[1:length(res$cpts)], 3)
  )
  print(break_summary, row.names = FALSE)
} else {
  cat("No structural breaks detected under the selected criterion.\n")
}

cat("\n========================================================\n")
cat("             sSIC MODEL SELECTION TABLE                 \n")
cat("========================================================\n")
print(res$ssic_table, row.names = FALSE)

