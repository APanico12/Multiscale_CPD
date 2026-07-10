# Title: Simulation Results Analysis
# Author: Antonio Panico
# Date: 2024-05-24
# Description: This script loads the simulation results from 'simulation_results.csv'
#              and computes the empirical size of the test for various significance levels.

# ─────────────────────────────────────────────
# 1. LOAD LIBRARIES
# ─────────────────────────────────────────────
library(dplyr)
library(readr)
library(tidyr)   # For reshaping data (pivot_wider)
library(knitr) # For creating LaTeX tables

# ─────────────────────────────────────────────
# 2. LOAD AND PROCESS DATA
# ─────────────────────────────────────────────

cat("Loading simulation results from 'simulation_results.csv'...\n")

# Define the file path (assumes it's in the same directory)
results_file <- "simulation_results.csv"

if (!file.exists(results_file)) {
  stop("Simulation results file not found: ", results_file, 
       "\nPlease run the 'New_Simulation_sedcd_null.R' script first.")
}

# Read the CSV data
sim_results_df <- read_csv(results_file, show_col_types = FALSE)

# see unique combinations of parameters
unique_params <- sim_results_df %>%
  select(level, dgp_name, n_values, tau, L_n_type) %>%
  distinct()


# Calculate empirical size (rejection rates) for each combination of parameters
# The 'stat' column contains the p-values from the simulation.
# Empirical size is the proportion of p_values < alpha_level.
# Calculate empirical size (rejection rate) for the 10% significance level
# for each combination of parameters.
empirical_size_summary <- sim_results_df %>%
  group_by(level, dgp_name, n_values, tau, L_n_type) %>%
  summarise(
    # Correct calculation: mean of a logical vector gives the proportion of TRUEs.
    `empirical_size_10%` = mean(stat < 0.10),
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

cat("\n\n--- Generating LaTeX table for the paper ---\n")

# Reshape the data to a wide format for the table
# Rows: DGP and tau
# Columns: Sample size (n)
latex_table_data <- empirical_size_summary %>%
  select(DGP, tau, n, `empirical_size_10%`) %>% # Select only necessary columns
  pivot_wider(names_from = n, values_from = `empirical_size_10%`)

# Create and print the LaTeX table code
latex_table <- kable(latex_table_data, 
                     format = "latex", 
                     booktabs = TRUE, 
                     caption = "Empirical Size at 10% Significance Level.",
                     label = "tab:empirical_size_10",
                     digits = 3,
                     align = 'llccc') # Align columns: left, left, center, center, center

cat(latex_table)

cat("\n\nLaTeX table code generated. You can copy and paste it into your paper.tex file.\n")

source("Plot_functions.R")
plot_pvalue_histograms_fabian(sim_results_df)
