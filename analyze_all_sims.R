# Title: Combined Simulation Results Analysis
# Author: Gemini
# Date: 2024-05-24
# Description: This script loads all simulation results from the 'sim_null_csv' directory,
#              combines them, and computes the empirical size of the test.

# ─────────────────────────────────────────────
# 1. LOAD LIBRARIES
# ─────────────────────────────────────────────
library(dplyr)
library(readr)
library(tidyr)
library(knitr)

# ─────────────────────────────────────────────
# 2. LOAD AND PROCESS DATA
# ─────────────────────────────────────────────

cat("Loading and combining all simulation results from 'sim_null_csv/'...
")

# Define the path to the directory containing the CSV files
results_dir <- "sim_null_csv"

# Get a list of all CSV files in the directory
csv_files <- list.files(results_dir, pattern = "[.]csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No simulation results files found in ", results_dir)
}

# Read all CSV files, add a 'd' column from the filename, and combine them
sim_results_df <- lapply(csv_files, function(file_path) {
  # Extract 'd' value from filename, e.g., '..._d_2_...' -> 2
  d_value <- as.numeric(sub(".*_d_(\\d+)_.*", "\\1", basename(file_path)))

  read_csv(file_path, show_col_types = FALSE) %>%
    mutate(d = d_value)
}) %>%
  bind_rows()

cat("All simulation files have been combined.\n")

# see unique combinations of parameters
unique_params <- sim_results_df %>%
  dplyr::select(level, d, dgp_name, n_values, tau, L_n_type) %>%
  distinct()


# Calculate empirical size (rejection rates) for each combination of parameters
# The 'stat' column contains the p-values from the simulation.
# Empirical size is the proportion of p_values < alpha_level.
# Calculate empirical size (rejection rate) for the 10% significance level
# for each combination of parameters.
empirical_size_summary <- sim_results_df %>%
  group_by(level, d, dgp_name, n_values, tau, L_n_type) %>%
  summarise(
    # Correct calculation: mean of a logical vector gives the proportion of TRUEs.
    `empirical_size_10%` = mean(stat < 0.05, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Rename columns for clarity in the output table
  rename(
    MC_levels = level,
    DGP = dgp_name,
    `L_n` = L_n_type,
    n = n_values
  )

# Print the summary table
print(empirical_size_summary)

# ─────────────────────────────────────────────
# 4. GENERATE LATEX TABLE
# ─────────────────────────────────────────────

cat("\n\n--- Generating LaTeX table for the combined results ---\n")

# Reshape the data to a wide format for the table
# Rows: DGP, d, and tau
# Columns: Sample size (n). We filter for a specific L_n to avoid duplicates.
latex_table_data <- empirical_size_summary %>%
  group_by(DGP, d, tau, n, L_n) %>%
  summarise(
    `empirical_size_10%` = mean(`empirical_size_10%`, na.rm = TRUE), # Average if there are duplicates
    .groups = 'drop'
  ) %>%
  select(DGP, d, tau, n, L_n, `empirical_size_10%`) %>% # Select necessary columns
  pivot_wider(names_from = L_n, values_from = `empirical_size_10%`)

# Create and print the LaTeX table code
latex_table <- kable(latex_table_data, 
                     format = "latex", 
                     booktabs = TRUE, 
                     caption = "Combined Empirical Size at 10% Significance Level.",
                     label = "tab:combined_empirical_size_10",
                     digits = 3,
                     # Align columns: DGP (l), tau (l), n=2000 (c)
                     align = 'llc')

cat(latex_table)

# Save the LaTeX table to a file
writeLines(latex_table, "paper/combined_empirical_size.tex")

cat("

LaTeX table saved to 'paper/combined_empirical_size.tex'.
")
