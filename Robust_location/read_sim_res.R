library(dplyr)
library(readr)
library(tidyr)

path <- "sim_results.csv"

if (!file.exists(path)) {
  stop("File sim_results.csv not found in current directory.")
}

# 1. Read and summarize data
df <- read.csv(path)

cat(sprintf("Loaded %d simulation records from %s\n", nrow(df), path))

average_summaryb <- df %>% 
  group_by(scenario, contamination, innov_dist, n) %>%
  summarise(
    mean_bias_teta     = mean(bias_teta, na.rm = TRUE),
    sd_teta            = sd(bias_teta, na.rm = TRUE),      
    mean_bias_lin_teta = mean(bias_lin_teta, na.rm = TRUE),
    sd_lin_teta        = sd(bias_lin_teta, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename(
    bias_teta     = mean_bias_teta,
    bias_lin_teta = mean_bias_lin_teta
  )

# Save the original CSV summary
write.csv(average_summaryb, "final_results.csv", row.names = FALSE)
cat("Saved summary results to final_results.csv\n")

# 2. Reshape data for the LaTeX table layout
df_long <- average_summaryb %>%
  pivot_longer(
    cols = c(bias_teta, sd_teta, bias_lin_teta, sd_lin_teta),
    names_to = "metric_name",
    values_to = "value"
  ) %>%
  mutate(
    # Map raw names to table labels (using SD for Standard Deviation and Bias for mean error)
    metric = ifelse(grepl("bias", metric_name), "Bias", "SD"),
    estimator = ifelse(grepl("lin", metric_name), "lin", "teta")
  ) %>% 
  select(-metric_name)

# Pivot so that (n + estimator) become columns
df_wide <- df_long %>%
  pivot_wider(
    names_from = c(n, estimator),
    values_from = value,
    names_sep = "_"
  )

# 3. Custom formatting function for 3 decimal places
fmt_val <- function(col_name, df_row) {
  if (col_name %in% names(df_row) && !is.null(df_row[[col_name]]) && length(df_row[[col_name]]) > 0 && !is.na(df_row[[col_name]])) {
    sprintf("%.3f", df_row[[col_name]])
  } else {
    "-"
  }
}

# Determine present n values
all_n <- sort(unique(df$n))
target_n <- c(100, 500, 1000, 5000)

# 4. Generate the complete LaTeX table
sink("table_output.tex")
cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{Standard Deviation and Bias of $\\widehat{\\Theta}_n^{\\text{lin}}(1)$ and $\\widetilde{\\Theta}_n(1)$ as estimators of $\\int_0^1 m(u)du$ under structural breaks and data contamination. Values based on Monte Carlo replications.}\n")
cat("\\label{tab:robust_location_bias}\n")
cat("\\resizebox{\\textwidth}{!}{\n")
cat("\\begin{tabular}{clllcc|cc|cc|cc}\n")
cat("\\toprule\n")
cat("& & & & \\multicolumn{2}{c}{$n=100$} & \\multicolumn{2}{c}{$n=500$} & \\multicolumn{2}{c}{$n=1000$} & \\multicolumn{2}{c}{$n=5000$} \\\\\n")
cat("\\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10} \\cmidrule(lr){11-12}\n")
cat("\\textbf{Hypothesis} & \\textbf{Scenario} & \\textbf{Innovations} & \\textbf{Metric} & $\\widetilde{\\Theta}_n$ & $\\widehat{\\Theta}_n^{\\text{lin}}$ & $\\widetilde{\\Theta}_n$ & $\\widehat{\\Theta}_n^{\\text{lin}}$ & $\\widetilde{\\Theta}_n$ & $\\widehat{\\Theta}_n^{\\text{lin}}$ & $\\widetilde{\\Theta}_n$ & $\\widehat{\\Theta}_n^{\\text{lin}}$ \\\\\n")
cat("\\midrule\n")

# Dictionaries for matching layout order
scenarios <- c("H0" = "\\mathcal{H}_0 \\text{ (No Shift)}", 
               "H1" = "\\mathcal{H}_1 \\text{ (Abrupt)}", 
               "H2" = "\\mathcal{H}_2 \\text{ (Gradual)}")
conts <- c("clean" = "Clean", "AO" = "AO", "IO" = "IO")
inns  <- c("gaussian" = "Gaussian", "t3" = "$t_3$")

# Iterate logically through the structure
for (sc in names(scenarios)) {
  cat(sprintf("%% --- %s Section ---\n", sc))
  cat(sprintf("\\multirow{12}{*}{$%s$}\n", scenarios[[sc]]))
  
  for (co in names(conts)) {
    for (inn in names(inns)) {
      
      # Filter for exact row combinations
      row_sd   <- df_wide %>% filter(scenario == sc, contamination == co, innov_dist == inn, metric == "SD")
      row_bias <- df_wide %>% filter(scenario == sc, contamination == co, innov_dist == inn, metric == "Bias")
      
      # Skip iteration if no data exists for this specific combination
      if (nrow(row_sd) == 0 || nrow(row_bias) == 0) next
      
      # Handle \multirow and \cmidrule logic based on column positions
      if (co == "clean" && inn == "gaussian") {
        prefix_sd <- sprintf("& \\multirow{4}{*}{%s} & \\multirow{2}{*}{%s} & SD", conts[[co]], inns[[inn]])
      } else if (inn == "gaussian") {
        cat("\\cmidrule{2-12}\n")
        prefix_sd <- sprintf("& \\multirow{4}{*}{%s} & \\multirow{2}{*}{%s} & SD", conts[[co]], inns[[inn]])
      } else {
        cat("\\cmidrule{3-12}\n")
        prefix_sd <- sprintf("& & \\multirow{2}{*}{%s} & SD", inns[[inn]])
      }
      
      prefix_bias <- "& & & Bias"
      
      # Extract values based on n sizes
      sd_vals <- c(
        fmt_val("100_teta", row_sd),   fmt_val("100_lin", row_sd),
        fmt_val("500_teta", row_sd),   fmt_val("500_lin", row_sd),
        fmt_val("1000_teta", row_sd),  fmt_val("1000_lin", row_sd),
        fmt_val("5000_teta", row_sd),  fmt_val("5000_lin", row_sd)
      )
      
      bias_vals <- c(
        fmt_val("100_teta", row_bias),  fmt_val("100_lin", row_bias),
        fmt_val("500_teta", row_bias),  fmt_val("500_lin", row_bias),
        fmt_val("1000_teta", row_bias), fmt_val("1000_lin", row_bias),
        fmt_val("5000_teta", row_bias), fmt_val("5000_lin", row_bias)
      )
      
      # Print formatted rows
      cat(paste(prefix_sd, "&", paste(sd_vals, collapse = " & "), "\\\\\n"))
      cat(paste(prefix_bias, "&", paste(bias_vals, collapse = " & "), "\\\\\n"))
    }
  }
  
  # Terminate sections cleanly
  if (sc != "H2") {
    cat("\\midrule\n")
  } else {
    cat("\\bottomrule\n")
  }
}

cat("\\end{tabular}\n")
cat("}\n")
cat("\\end{table}\n")
sink()

cat("LaTeX table generated and saved as table_output.tex\n")

