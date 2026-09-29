# ==============================================================================
# Script: plot_electricity_empirical.R
# Purpose: Generate publication-ready figures for the European Electricity Market
#          empirical application in SEDCD:
#          1. electricity_ts_returns (6x2 multi-panel time series & daily returns)
#          2. electricity_residuals_bivariate (3x5 pairwise residual scatter with time hue)
#
# Outputs: Saved in PDF (vector) and PNG (300 DPI) format in:
#          - Multiscale_CPD/SEDCD/
#          - Multiscale_CPD/paper/img/
# ==============================================================================

cat("\n===================================================================\n")
cat(" Generating Empirical Application Figures for Electricity Market\n")
cat("===================================================================\n\n")

# Set working directory to SEDCD if running interactively
if (!file.exists("Plot_functions.R") && file.exists("Multiscale_CPD/SEDCD/Plot_functions.R")) {
  setwd("Multiscale_CPD/SEDCD")
}

# 1. Source required modules
source("sedcd_mscale.R")
source("Plot_functions.R")

# 2. Load dataset
data_file <- "prices.csv"
if (!file.exists(data_file)) {
  stop("Error: prices.csv not found in current directory!")
}

df <- read.csv(data_file, stringsAsFactors = FALSE)
cat(sprintf(" Loaded dataset: %s (%d rows, %d columns)\n", data_file, nrow(df), ncol(df)))

core_countries <- c("DE_LU", "IT_NORD", "FR", "BE", "NL", "PL")
cp_date <- as.Date("2022-07-07")

# Target directories
sedcd_dir <- "."
paper_img_dir <- "../paper/img"
art_dir <- "C:/Users/Antonio/.gemini/antigravity/brain/57f79b8c-2c63-47f6-9959-52420cf4e82f"

if (!dir.exists(paper_img_dir)) {
  dir.create(paper_img_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# FIGURE 1: Electricity Time Series & Daily Returns (6x2 Grid)
# ------------------------------------------------------------------------------
cat("\n[1/2] Generating Electricity Time Series and Returns figure (6x2)...\n")

file_ts_pdf_sedcd <- file.path(sedcd_dir, "electricity_ts_returns.pdf")
file_ts_png_sedcd <- file.path(sedcd_dir, "electricity_ts_returns.png")
file_ts_pdf_paper <- file.path(paper_img_dir, "electricity_ts_returns.pdf")
file_ts_png_paper <- file.path(paper_img_dir, "electricity_ts_returns.png")

plot_electricity_ts_returns(df, core_countries = core_countries,
                            filename = file_ts_pdf_sedcd, is_pdf = TRUE)
plot_electricity_ts_returns(df, core_countries = core_countries,
                            filename = file_ts_png_sedcd, is_pdf = FALSE)

# Sync to paper/img/
file.copy(file_ts_pdf_sedcd, file_ts_pdf_paper, overwrite = TRUE)
file.copy(file_ts_png_sedcd, file_ts_png_paper, overwrite = TRUE)

# Sync to brain artifact if exists
if (dir.exists(art_dir)) {
  file.copy(file_ts_png_sedcd, file.path(art_dir, "electricity_ts_returns.png"), overwrite = TRUE)
}

cat("  -> Saved:", file_ts_pdf_sedcd, "\n")
cat("  -> Saved:", file_ts_png_sedcd, "\n")
cat("  -> Synced to:", file_ts_pdf_paper, "\n")
cat("  -> Synced to:", file_ts_png_paper, "\n")

# ------------------------------------------------------------------------------
# FIGURE 2: Pairwise Bivariate Residual Scatter with Time Progression Hue (3x5 Grid)
# ------------------------------------------------------------------------------
cat("\n[2/2] Generating Pairwise Bivariate Residual Scatter figure (3x5)...\n")

file_res_pdf_sedcd <- file.path(sedcd_dir, "electricity_residuals_bivariate.pdf")
file_res_png_sedcd <- file.path(sedcd_dir, "electricity_residuals_bivariate.png")
file_res_pdf_paper <- file.path(paper_img_dir, "electricity_residuals_bivariate.pdf")
file_res_png_paper <- file.path(paper_img_dir, "electricity_residuals_bivariate.png")

plot_residual_bivariate_scatter(df, core_countries = core_countries,
                                residuals_type = "mscale", point_alpha = 0.45,
                                filename = file_res_pdf_sedcd, is_pdf = TRUE)
plot_residual_bivariate_scatter(df, core_countries = core_countries,
                                residuals_type = "mscale", point_alpha = 0.45,
                                filename = file_res_png_sedcd, is_pdf = FALSE)

# Sync to paper/img/
file.copy(file_res_pdf_sedcd, file_res_pdf_paper, overwrite = TRUE)
file.copy(file_res_png_sedcd, file_res_png_paper, overwrite = TRUE)

# Sync to brain artifact if exists
if (dir.exists(art_dir)) {
  file.copy(file_res_png_sedcd, file.path(art_dir, "electricity_residuals_bivariate.png"), overwrite = TRUE)
}

cat("  -> Saved:", file_res_pdf_sedcd, "\n")
cat("  -> Saved:", file_res_png_sedcd, "\n")
cat("  -> Synced to:", file_res_pdf_paper, "\n")
cat("  -> Synced to:", file_res_png_paper, "\n")

cat("\n===================================================================\n")
cat(" All figures generated and synchronized successfully!\n")
cat("===================================================================\n\n")

