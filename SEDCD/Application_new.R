# ==============================================================================
# Application_new.R
# Purpose: Master Empirical Application Script
#          Multivariate Robust M-Scale SEDCD CUSUM Test for Core 6 European Markets
#          (Germany/Luxembourg, Italy, France, Belgium, Netherlands, Poland)
#          Square Side-by-Side Diagnostics: SEDCD (Left) & CUSUM (Right)
# ==============================================================================

cat("\n================================================================================\n")
cat("   MULTIVARIATE ROBUST M-SCALE SEDCD TEST: CORE 6 EUROPEAN MARKETS\n")
cat("   (Germany/Luxembourg, Italy, France, Belgium, Netherlands, Poland)\n")
cat("================================================================================\n\n")

suppressPackageStartupMessages({
  library(zoo)
  library(MASS)
  library(numDeriv)
})

source("sedcd_mscale.R")
source("Plot_functions.R")

# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------
file_path <- "C:/Users/Antonio/Desktop/prices.csv"
file_path <- if (file.exists("prices.csv")) "prices.csv" else "C:/Users/Antonio/Desktop/prices.csv"
if (!file.exists(file_path)) stop("Data file not found at: ", file_path)

df <- read.csv(file_path)
calendar_dates <- as.Date(substr(df$Date, 1, 10))

# Core 6 Continental European Markets
core_countries <- c("DE_LU", "IT_NORD", "FR", "BE", "NL", "PL")
all_cols <- c("Date", c("DE_LU", "IT_NORD", "FR", "AT", "PL", "ES", "BE", "NL"))
colnames(df) <- all_cols

country_names <- c(
  "DE_LU"   = "Germany / Luxembourg",
  "IT_NORD" = "Italy (North)",
  "FR"      = "France",
  "BE"      = "Belgium",
  "NL"      = "Netherlands",
  "PL"      = "Poland"
)

price_mat <- as.matrix(df[, core_countries])
ret_mat <- na.omit(diff(price_mat))
ret_dates <- calendar_dates[-1]
N <- nrow(ret_mat)
d <- ncol(ret_mat)
p <- 1 + 3 * d + d^2

cat(sprintf("1. Data Loaded: %d observations (%s to %s)\n",
            N, format(ret_dates[1], "%Y-%m-%d"), format(tail(ret_dates, 1), "%Y-%m-%d")))
cat(sprintf("   Core 6 Markets (d = %d, p = %d): %s\n\n",
            d, p, paste(core_countries, collapse = ", ")))

# ------------------------------------------------------------------------------
# 2. Parameters Calibration
# ------------------------------------------------------------------------------
alpha_tail  <- 0.10                       # 10% lower downside quantile
tau_val     <- qnorm(alpha_tail)          # tau = Phi^-1(alpha_tail) = -1.28155
gamma_val   <- 0.05                       # Logistic smoother bandwidth
c_scale     <- 2.985                      # Dennis-Welsch M-scale tuning (95% efficiency)
k_val       <- 50                        # Backward window bandwidth M_n
L_n         <- max(1, floor(0.1 * (log(N))^2)) # Separation lag (6)
block_size  <- L_n                        # Variance block size (6)
cutoff      <- k_val + L_n + block_size   # Boundary cutoff (152)
MC_reps     <- 500                        # Multiplier bootstrap replications
dh_step_val <- 10                         # Step caching for Jacobian updates

cat("2. Theoretical Calibration:\n")
cat(sprintf("   - Downside Tail:        alpha_tail = %.2f -> tau = Phi^-1(alpha) = %.4f\n", alpha_tail, tau_val))
cat(sprintf("   - Smoothing Bandwidth:  gamma = %.2f\n", gamma_val))
cat(sprintf("   - Dennis-Welsch M-Scale: c = %.3f (delta = %.5f)\n", c_scale, welsch_delta_constant(c_scale)))
cat(sprintf("   - Window & Cutoff:      k = %d, L_n = %d, b = %d, cutoff = %d\n", k_val, L_n, block_size, cutoff))
cat(sprintf("   - Bootstrap:            MC = %d replications\n\n", MC_reps))

# ------------------------------------------------------------------------------
# 3. Execute Multivariate Test on Core 6 Countries
# ------------------------------------------------------------------------------
cat("3. Executing Multivariate Robust M-Scale SEDCD CUSUM Test (d = 6, p = 55)...\n")
t0 <- Sys.time()
res_core6 <- cusum_sedcd_multivariate(
  X        = ret_mat,
  tau      = tau_val,
  gamma    = gamma_val,
  c        = c_scale,
  k        = k_val,
  lag      = L_n,
  block    = block_size,
  MC       = MC_reps,
  dh_step  = dh_step_val,
  plotting = FALSE
)
t_elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cp_idx <- res_core6$tau_idx
cp_date <- ret_dates[cp_idx]
cp_u <- res_core6$tau_u

sedcd_s <- res_core6$sedcd
pre_mean  <- mean(sedcd_s[(cutoff + 1):cp_idx], na.rm = TRUE)
post_mean <- mean(sedcd_s[(cp_idx + 1):N], na.rm = TRUE)
overall_mean <- mean(sedcd_s[(cutoff + 1):N], na.rm = TRUE)
shift_val <- post_mean - pre_mean
pct_shift <- (shift_val / pre_mean) * 100

cat("\n================================================================================\n")
cat("      MULTIVARIATE ROBUST M-SCALE SEDCD TEST RESULTS: CORE 6 MARKETS            \n")
cat("================================================================================\n")
cat(sprintf("Markets Analyzed (6):        %s\n", paste(core_countries, collapse = ", ")))
cat(sprintf("Observations (N):            %d days (%s to %s)\n", N, format(ret_dates[1], "%Y-%m-%d"), format(tail(ret_dates, 1), "%Y-%m-%d")))
cat(sprintf("Test Statistic (Z):          %.4f\n", res_core6$test_stat))
cat(sprintf("90%% Critical Value (q90):    %.4f\n", res_core6$q90))
cat(sprintf("95%% Critical Value (q95):    %.4f\n", res_core6$q95))
cat(sprintf("Bootstrap P-Value:           %.4f\n", res_core6$p_value))

dec_str <- if (res_core6$p_value < 0.05) {
  "REJECT H0 at 5% level (Statistically Significant Break)"
} else if (res_core6$p_value < 0.10) {
  "REJECT H0 at 10% level"
} else {
  "FAIL TO REJECT H0"
}
cat(sprintf("Decision:                    %s\n", dec_str))
cat(sprintf("Detected Break Date:         %s (Index: %d, Fractional time u = %.4f)\n",
            format(cp_date, "%Y-%m-%d"), cp_idx, cp_u))
cat(sprintf("Overall Mean SEDCD:          %.4f\n", overall_mean))
cat(sprintf("Pre-Crisis Mean SEDCD:       %.4f\n", pre_mean))
cat(sprintf("Post-Crisis Mean SEDCD:      %.4f\n", post_mean))
cat(sprintf("Structural Shift:            %+.4f (%+.1f%% jump in downside contagion)\n", shift_val, pct_shift))
cat(sprintf("Execution Time:              %.2f seconds\n", t_elapsed))
cat("================================================================================\n\n")

# ------------------------------------------------------------------------------
# 4. Publication Plot: Square Side-by-Side (SEDCD Left, CUSUM Right)
# ------------------------------------------------------------------------------
cat("4. Generating Square Side-by-Side Publication Diagnostic Figure...\n")

plot_core6_sedcd <- function(filename = NULL, is_pdf = FALSE) {
  if (!is.null(filename)) {
    if (is_pdf) {
      pdf(filename, width = 11.0, height = 5.0)
      on.exit(dev.off(), add = TRUE)
    } else {
      png(filename, width = 3300, height = 1500, res = 300)
      on.exit(dev.off(), add = TRUE)
    }
  }

  # Square panels side by side (1 row, 2 columns, pty = "s")
  par(mfrow = c(1, 2), pty = "s", oma = c(0.8, 1.0, 0.8, 1.0), mar = c(3.0, 4.4, 1.0, 1.2),
      mgp = c(2.4, 0.6, 0), tcl = -0.35)

  dates_v <- ret_dates[(cutoff + 1):N]

  # Common yearly calendar ticks
  year_dates <- as.Date(c("2020-01-01", "2021-01-01", "2022-01-01", "2023-01-01", 
                          "2024-01-01", "2025-01-01", "2026-01-01"))
  year_labels <- c("2020", "2021", "2022", "2023", "2024", "2025", "2026")

  # ----------------------------------------------------------------------------
  # LEFT PANEL: Panel A - Aggregate Downside Correlation Density SEDCD
  # ----------------------------------------------------------------------------
  sedcd_plot <- sedcd_s
  sedcd_plot[1:cutoff] <- NA
  ymax_s <- max(c(sedcd_plot, 0.65), na.rm = TRUE) * 1.15

  plot(ret_dates, sedcd_plot, type = "n",
       xlim = range(ret_dates), ylim = c(0, ymax_s),
       xlab = "", ylab = expression(widehat(SEDCD)[tau](t)),
       main = "", # Title removed as requested
       las = 1, xaxt = "n")

  grid(col = "gray90", lty = 2, lwd = 0.6)
  axis.Date(1, at = year_dates, labels = year_labels, cex.axis = 0.85)

  # Trajectory line
  lines(ret_dates, sedcd_plot, col = "#2e7d32", lwd = 1.1)

  # Break line and pre/post regime means
  abline(v = cp_date, col = "#b71c1c", lty = 2, lwd = 1.0)
  segments(ret_dates[cutoff + 1], pre_mean, cp_date, pre_mean, col = "#1565c0", lwd = 1.6)
  segments(cp_date, post_mean, tail(ret_dates, 1), post_mean, col = "#c62828", lwd = 1.6)

  # Legend removed from Left Panel as requested

  # ----------------------------------------------------------------------------
  # RIGHT PANEL: Panel B - Linearized CUSUM Process
  # ----------------------------------------------------------------------------
  Tu_v <- res_core6$Tu[(cutoff + 1):N]
  q95 <- res_core6$q95
  q90 <- res_core6$q90
  ymax_t <- max(abs(c(Tu_v, q95))) * 1.25
  fill_col <- adjustcolor("#d32f2f", alpha.f = 0.12)

  plot(ret_dates, res_core6$Tu, type = "n",
       xlim = range(ret_dates), ylim = c(-ymax_t, ymax_t),
       xlab = "", ylab = expression(T[n](u)),
       main = "", # Title removed as requested
       las = 1, xaxt = "n")

  grid(col = "gray90", lty = 2, lwd = 0.6)
  axis.Date(1, at = year_dates, labels = year_labels, cex.axis = 0.85)

  # Rejection shading
  polygon(c(dates_v, rev(dates_v)), c(pmax(Tu_v, q95), rep(q95, length(dates_v))), col = fill_col, border = NA)
  polygon(c(dates_v, rev(dates_v)), c(pmin(Tu_v, -q95), rep(-q95, length(dates_v))), col = fill_col, border = NA)

  # Reference & Critical Value lines
  abline(h = 0, col = "gray45", lty = 2, lwd = 0.7)
  abline(h = c(q95, -q95), col = "#b71c1c", lty = 3, lwd = 0.9)
  abline(h = c(q90, -q90), col = "#e65100", lty = 2, lwd = 0.8)

  # Trajectory line
  lines(dates_v, Tu_v, col = "#1565c0", lwd = 1.1)

  # Break marker & text label
  abline(v = cp_date, col = "#b71c1c", lty = 2, lwd = 1.0)
  text(cp_date, ymax_t * 0.85,
       labels = sprintf(" Break: %s", format(cp_date, "%Y-%m-%d")),
       col = "#b71c1c", font = 2, adj = 0, cex = 0.78)

  # Only critical values in legend as requested
  legend("bottomleft", inset = c(0.02, 0.02),
         legend = c("95% Critical Value", "90% Critical Value"),
         col = c("#b71c1c", "#e65100"),
         lty = c(3, 2), lwd = c(0.9, 0.8),
         bty = "o", box.col = "gray85", bg = "white", cex = 0.75)
}

# Save Figures
pdf_core6 <- "../paper/img/sedcd_multivariate_core6.pdf"
png_core6 <- "../paper/img/sedcd_multivariate_core6.png"
png_loc   <- "sedcd_multivariate_core6.png"
png_art   <- "C:/Users/Antonio/.gemini/antigravity/brain/57f79b8c-2c63-47f6-9959-52420cf4e82f/sedcd_multivariate_core6.png"

plot_core6_sedcd(pdf_core6, is_pdf = TRUE)
plot_core6_sedcd(png_core6, is_pdf = FALSE)
plot_core6_sedcd(png_loc, is_pdf = FALSE)
if (dir.exists("C:/Users/Antonio/.gemini/antigravity/brain/57f79b8c-2c63-47f6-9959-52420cf4e82f")) {
  plot_core6_sedcd(png_art, is_pdf = FALSE)
}

cat("   [Saved] Publication PDF: ", pdf_core6, "\n")
cat("   [Saved] Publication PNG: ", png_core6, "\n")
cat("   [Saved] Local PNG:       ", png_loc, "\n")
cat("   [Saved] Brain Artifact:  ", png_art, "\n\n")
# ------------------------------------------------------------------------------
# 5. Additional Empirical Figures: Time Series & Returns + Bivariate Residuals
# ------------------------------------------------------------------------------
cat("5. Generating Additional Empirical Figures (Prices/Returns & Bivariate Residuals)...\n")

# Figure 1: 6x2 Prices & Daily Returns
plot_electricity_ts_returns(df, core_countries = core_countries,
                            filename = "../paper/img/electricity_ts_returns.pdf", is_pdf = TRUE)
plot_electricity_ts_returns(df, core_countries = core_countries,
                            filename = "../paper/img/electricity_ts_returns.png", is_pdf = FALSE)
plot_electricity_ts_returns(df, core_countries = core_countries,
                            filename = "electricity_ts_returns.png", is_pdf = FALSE)

# Figure 2: 3x5 Bivariate Residual Scatter
plot_residual_bivariate_scatter(df, core_countries = core_countries,
                                residuals_type = "mscale", point_alpha = 0.45,
                                filename = "../paper/img/electricity_residuals_bivariate.pdf", is_pdf = TRUE)
plot_residual_bivariate_scatter(df, core_countries = core_countries,
                                residuals_type = "mscale", point_alpha = 0.45,
                                filename = "../paper/img/electricity_residuals_bivariate.png", is_pdf = FALSE)
plot_residual_bivariate_scatter(df, core_countries = core_countries,
                                residuals_type = "mscale", point_alpha = 0.45,
                                filename = "electricity_residuals_bivariate.png", is_pdf = FALSE)

cat("   [Saved] electricity_ts_returns.pdf / .png\n")
cat("   [Saved] electricity_residuals_bivariate.pdf / .png\n\n")

cat("================================================================================\n")
cat("CORE 6 COUNTRIES SEDCD TEST & ALL EMPIRICAL PLOTS COMPLETE!\n")
cat("================================================================================\n")
