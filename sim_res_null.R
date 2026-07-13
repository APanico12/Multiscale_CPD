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
    `empirical_size_10%` = mean(stat < 0.1 ),
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
# Columns: Sample size (n). We filter for a specific L_n to avoid duplicates.
latex_table_data <- empirical_size_summary %>%
  # The warning occurs because there are multiple L_n types for each (DGP, tau, n) combo.
  # We select one L_n type to display in the table, as is common for paper summaries.
  filter(`L_n` == "0.25(log(n))^2") %>%
  select(DGP, tau, n, `empirical_size_10%`) %>% # Select only necessary columns
  pivot_wider(names_from = n, values_from = `empirical_size_10%`)

# Create and print the LaTeX table code
latex_table <- kable(latex_table_data, 
                     format = "latex", 
                     booktabs = TRUE, 
                     caption = "Empirical Size at 10% Significance Level.",
                     label = "tab:empirical_size_10",
                     digits = 3,
                     # Align columns: DGP (l), tau (l), n=2000 (c)
                     align = 'llc')

cat(latex_table)

cat("\n\nLaTeX table code generated. You can copy and paste it into your paper.tex file.\n")

# ─────────────────────────────────────────────
# 5. PLOT P-VALUE DISTRIBUTIONS
# ─────────────────────────────────────────────

source("Plot_functions.R")

# Ensure the output directory exists
img_dir <- "paper/img"
if (!dir.exists(img_dir)) {
  dir.create(img_dir, recursive = TRUE)
}

# --- Generate the main multi-panel histogram plot ---
cat("\n\n--- Generating main p-value distribution plot ---\n")
pdf(file = file.path(img_dir, "pvalue_hist_main.pdf"), width = 10, height = 12)
plot_pvalue_histograms_fabian(sim_results_df)
dev.off()
cat(sprintf("Saved main p-value histogram grid to '%s'\n", file.path(img_dir, "pvalue_hist_main.pdf")))


# --- Generate a plot for a single, specific combination ---
cat("\n\n--- Generating p-value distribution plot for a single combination ---\n")

# Define the specific combination you want to plot
single_n <- 2000
single_dgp <- "Clayton"
single_tau <- -2
single_Ln <- "0.25(log(n))^2"

# Save the single plot to a PDF
pdf(file = file.path(img_dir, "pvalue_hist_single.pdf"), width = 6, height = 6)
plot_pvalue_histograms_fabian(
  sim_results_df,
  n_to_plot = single_n,
  dgp_to_plot = single_dgp,
  tau_to_plot = single_tau,
  Ln_to_plot = single_Ln
)
dev.off()

cat(sprintf("Saved single p-value histogram to '%s'\n", file.path(img_dir, "pvalue_hist_single.pdf")))
