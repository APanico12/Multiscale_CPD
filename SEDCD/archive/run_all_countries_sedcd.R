# ==============================================================================
# run_all_countries_sedcd.R
# Purpose: Multivariate Robust M-Scale SEDCD CUSUM Change-Point Detection Test
#          for All European Wholesale Electricity Markets Together (2020 - 2026)
# ==============================================================================

cat("\n================================================================================\n")
cat("   MULTIVARIATE ROBUST M-SCALE SEDCD CUSUM TEST: ALL EUROPEAN MARKETS\n")
cat("================================================================================\n\n")

suppressPackageStartupMessages({
  library(zoo)
  library(MASS)
  library(numDeriv)
})

source("sedcd_mscale.R")

# 1. Load Price Data
file_path <- "C:/Users/Antonio/Desktop/prices.csv"
if (!file.exists(file_path)) stop("Data file not found at: ", file_path)

df <- read.csv(file_path)
calendar_dates <- as.Date(substr(df$Date, 1, 10))
all_countries <- c("DE_LU", "IT_NORD", "FR", "AT", "PL", "ES", "BE", "NL")
colnames(df) <- c("Date", all_countries)

country_names <- c(
  "DE_LU"   = "Germany / Luxembourg",
  "IT_NORD" = "Italy (North)",
  "FR"      = "France",
  "AT"      = "Austria",
  "PL"      = "Poland",
  "ES"      = "Spain",
  "BE"      = "Belgium",
  "NL"      = "Netherlands"
)

# Compute price returns
price_mat <- as.matrix(df[, all_countries])
ret_mat <- na.omit(diff(price_mat))
ret_dates <- calendar_dates[-1]
N <- nrow(ret_mat)
d_all <- ncol(ret_mat)

cat(sprintf("1. Data Loaded: %d observations (%s to %s)\n",
            N, format(ret_dates[1], "%Y-%m-%d"), format(tail(ret_dates, 1), "%Y-%m-%d")))
cat(sprintf("   All European Markets (d = %d): %s\n\n",
            d_all, paste(all_countries, collapse = ", ")))

# 2. Parameters
alpha_tail  <- 0.10                       # 10% lower tail quantile
tau_val     <- qnorm(alpha_tail)          # tau = -1.28155
gamma_val   <- 0.05                       # Logistic smoother bandwidth
c_scale     <- 2.985                      # Dennis-Welsch tuning parameter (95% efficiency)
k_val       <- 140                        # Window bandwidth M_n
L_n         <- max(1, floor(0.1 * (log(N))^2)) # Separation lag (6)
block_size  <- L_n                        # Block size (6)
cutoff      <- k_val + L_n + block_size   # Boundary cutoff (152)
MC_reps     <- 500                        # Multiplier bootstrap replications
dh_step_val <- 10                         # Jacobian step caching

cat("2. Theoretical Calibration:\n")
cat(sprintf("   - Downside Tail:        alpha_tail = %.2f -> tau = Phi^-1(alpha) = %.4f\n", alpha_tail, tau_val))
cat(sprintf("   - Smoothing Bandwidth:  gamma = %.2f\n", gamma_val))
cat(sprintf("   - Dennis-Welsch M-Scale: c = %.3f (delta = %.5f)\n", c_scale, welsch_delta_constant(c_scale)))
cat(sprintf("   - Window & Cutoff:      k = %d, L_n = %d, b = %d, cutoff = %d\n", k_val, L_n, block_size, cutoff))
cat(sprintf("   - Bootstrap:            MC = %d replications\n\n", MC_reps))

# 3. Execute Joint Multivariate Test (All 8 Countries)
cat("3. Executing Joint Multivariate SEDCD Test on All 8 Countries...\n")
t0 <- Sys.time()
res_all8 <- cusum_sedcd_multivariate(
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
t_elapsed_all8 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cp_idx <- res_all8$tau_idx
cp_date <- ret_dates[cp_idx]
cp_u <- res_all8$tau_u

sedcd_s <- res_all8$sedcd
pre_mean_all8  <- mean(sedcd_s[(cutoff + 1):cp_idx], na.rm = TRUE)
post_mean_all8 <- mean(sedcd_s[(cp_idx + 1):N], na.rm = TRUE)
overall_mean   <- mean(sedcd_s[(cutoff + 1):N], na.rm = TRUE)

# 4. Also Execute Core 6-Market System for direct comparison
core_countries <- c("DE_LU", "IT_NORD", "FR", "BE", "NL", "PL")
cat("\n4. Executing Core 6-Market System (CWE + Italy + Poland) for comparison...\n")
t0_core <- Sys.time()
res_core6 <- cusum_sedcd_multivariate(
  X        = ret_mat[, core_countries],
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
t_elapsed_core6 <- as.numeric(difftime(Sys.time(), t0_core, units = "secs"))

cp_idx_core <- res_core6$tau_idx
cp_date_core <- ret_dates[cp_idx_core]

# Print Comparative Summary
cat("\n================================================================================\n")
cat("      MULTIVARIATE ROBUST M-SCALE SEDCD TEST RESULTS COMPARISON                 \n")
cat("================================================================================\n")
cat(sprintf("%-28s %-25s %-25s\n", "Metric", "All 8 European Markets", "Core 6 European Markets"))
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("%-28s %-25s %-25s\n", "Dimensions (d)", "d = 8", "d = 6"))
cat(sprintf("%-28s %-25s %-25s\n", "Parameters (p)", "p = 89", "p = 55"))
cat(sprintf("%-28s %-25s %-25s\n", "Markets Included", "DE,IT,FR,AT,PL,ES,BE,NL", "DE,IT,FR,BE,NL,PL"))
cat(sprintf("%-28s %-25.4f %-25.4f\n", "Test Statistic (Z)", res_all8$test_stat, res_core6$test_stat))
cat(sprintf("%-28s %-25.4f %-25.4f\n", "90% Critical Value (q90)", res_all8$q90, res_core6$q90))
cat(sprintf("%-28s %-25.4f %-25.4f\n", "95% Critical Value (q95)", res_all8$q95, res_core6$q95))
cat(sprintf("%-28s %-25.4f %-25.4f\n", "Bootstrap P-Value", res_all8$p_value, res_core6$p_value))

dec8 <- if (res_all8$p_value < 0.05) "REJECT (5%)" else if (res_all8$p_value < 0.10) "REJECT (10%)" else "FAIL REJ"
dec6 <- if (res_core6$p_value < 0.05) "REJECT (5%)" else if (res_core6$p_value < 0.10) "REJECT (10%)" else "FAIL REJ"
cat(sprintf("%-28s %-25s %-25s\n", "Decision", dec8, dec6))
cat(sprintf("%-28s %-25s %-25s\n", "Detected Break Date", format(cp_date, "%Y-%m-%d"), format(cp_date_core, "%Y-%m-%d")))
cat(sprintf("%-28s %-25.4f %-25.4f\n", "Fractional Time (u)", cp_u, res_core6$tau_u))
cat(sprintf("%-28s %-25.4f %-25.4f\n", "Overall Mean SEDCD", overall_mean, mean(res_core6$sedcd[(cutoff+1):N], na.rm=TRUE)))
cat(sprintf("%-28s %-25.4f %-25.4f\n", "Pre-Break Mean SEDCD", pre_mean_all8, mean(res_core6$sedcd[(cutoff+1):cp_idx_core], na.rm=TRUE)))
cat(sprintf("%-28s %-25.4f %-25.4f\n", "Post-Break Mean SEDCD", post_mean_all8, mean(res_core6$sedcd[(cp_idx_core+1):N], na.rm=TRUE)))
cat(sprintf("%-28s %-25.2f s %-25.2f s\n", "Execution Time", t_elapsed_all8, t_elapsed_core6))
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# 5. Publication Plots
cat("5. Generating Publication Figures...\n")

plot_multivariate_comparison <- function(filename = NULL, is_pdf = FALSE) {
  if (!is.null(filename)) {
    if (is_pdf) {
      pdf(filename, width = 10, height = 8.5)
      on.exit(dev.off(), add = TRUE)
    } else {
      png(filename, width = 3000, height = 2550, res = 300)
      on.exit(dev.off(), add = TRUE)
    }
  }

  par(mfrow = c(2, 2), oma = c(1.5, 1.5, 2.8, 1), mar = c(3.8, 4.5, 2.2, 1.2))
  fill_col <- adjustcolor("#d32f2f", alpha.f = 0.16)

  dates_v <- ret_dates[(cutoff + 1):N]

  # --- Panel A: All 8 Markets CUSUM Process ---
  Tu8 <- res_all8$Tu[(cutoff + 1):N]
  q95_8 <- res_all8$q95
  q90_8 <- res_all8$q90
  ymax8 <- max(abs(c(Tu8, q95_8))) * 1.25

  plot(ret_dates, res_all8$Tu, type = "n", ylim = c(-ymax8, ymax8),
       xlab = "", ylab = expression(T[n](u)),
       main = sprintf("Panel A: All 8 Markets CUSUM Process (Z = %.2f, p = %.3f)", res_all8$test_stat, res_all8$p_value),
       cex.main = 0.95, font.main = 2, las = 1)
  grid(col = "gray88", lty = 2)

  polygon(c(dates_v, rev(dates_v)), c(pmax(Tu8, q95_8), rep(q95_8, length(dates_v))), col = fill_col, border = NA)
  polygon(c(dates_v, rev(dates_v)), c(pmin(Tu8, -q95_8), rep(-q95_8, length(dates_v))), col = fill_col, border = NA)

  abline(h = 0, col = "gray40", lty = 2, lwd = 0.9)
  abline(h = c(q95_8, -q95_8), col = "#b71c1c", lty = 3, lwd = 1.3)
  abline(h = c(q90_8, -q90_8), col = "#e65100", lty = 2, lwd = 1.0)
  lines(dates_v, Tu8, col = "#1565c0", lwd = 1.9)
  abline(v = cp_date, col = "#b71c1c", lty = 2, lwd = 1.6)
  text(cp_date, ymax8 * 0.85, labels = sprintf(" Peak: %s", format(cp_date, "%b %Y")), col = "#b71c1c", font = 2, cex = 0.82)

  legend("bottomleft", legend = c(expression(T[n](u)), "95% Crit", "90% Crit"),
         col = c("#1565c0", "#b71c1c", "#e65100"), lty = c(1, 3, 2), lwd = c(1.9, 1.3, 1.0),
         bty = "o", box.col = "gray85", bg = "white", cex = 0.72)

  # --- Panel B: Core 6 Markets CUSUM Process ---
  Tu6 <- res_core6$Tu[(cutoff + 1):N]
  q95_6 <- res_core6$q95
  q90_6 <- res_core6$q90
  ymax6 <- max(abs(c(Tu6, q95_6))) * 1.25

  plot(ret_dates, res_core6$Tu, type = "n", ylim = c(-ymax6, ymax6),
       xlab = "", ylab = expression(T[n](u)),
       main = sprintf("Panel B: Core 6 Markets CUSUM Process (Z = %.2f, p = %.3f)", res_core6$test_stat, res_core6$p_value),
       cex.main = 0.95, font.main = 2, las = 1)
  grid(col = "gray88", lty = 2)

  polygon(c(dates_v, rev(dates_v)), c(pmax(Tu6, q95_6), rep(q95_6, length(dates_v))), col = fill_col, border = NA)
  polygon(c(dates_v, rev(dates_v)), c(pmin(Tu6, -q95_6), rep(-q95_6, length(dates_v))), col = fill_col, border = NA)

  abline(h = 0, col = "gray40", lty = 2, lwd = 0.9)
  abline(h = c(q95_6, -q95_6), col = "#b71c1c", lty = 3, lwd = 1.3)
  abline(h = c(q90_6, -q90_6), col = "#e65100", lty = 2, lwd = 1.0)
  lines(dates_v, Tu6, col = "#1565c0", lwd = 1.9)
  abline(v = cp_date_core, col = "#b71c1c", lty = 2, lwd = 1.6)
  text(cp_date_core, ymax6 * 0.85, labels = sprintf(" Peak: %s", format(cp_date_core, "%b %Y")), col = "#b71c1c", font = 2, cex = 0.82)

  legend("bottomleft", legend = c(expression(T[n](u)), "95% Crit", "90% Crit"),
         col = c("#1565c0", "#b71c1c", "#e65100"), lty = c(1, 3, 2), lwd = c(1.9, 1.3, 1.0),
         bty = "o", box.col = "gray85", bg = "white", cex = 0.72)

  # --- Panel C: All 8 Markets Aggregate SEDCD(t) ---
  sedcd_8_s <- res_all8$sedcd
  sedcd_8_s[1:cutoff] <- NA
  ymax_s8 <- max(c(sedcd_8_s, 0.6), na.rm = TRUE) * 1.2

  plot(ret_dates, sedcd_8_s, type = "l", col = "#2e7d32", lwd = 1.8,
       ylim = c(0, ymax_s8), xlab = "Calendar Date", ylab = expression(widehat(SEDCD)[tau](t)),
       main = expression(bold("Panel C: All 8 Markets Downside Density ") * bolditalic(SEDCD[tau](t))),
       cex.main = 0.95, font.main = 2, las = 1)
  grid(col = "gray88", lty = 2)
  abline(v = cp_date, col = "#b71c1c", lty = 2, lwd = 1.6)
  segments(ret_dates[cutoff + 1], pre_mean_all8, cp_date, pre_mean_all8, col = "#1565c0", lwd = 2.2)
  segments(cp_date, post_mean_all8, tail(ret_dates, 1), post_mean_all8, col = "#c62828", lwd = 2.2)

  legend("topleft",
         legend = c(expression(widehat(SEDCD)[tau](t)),
                    sprintf("Pre-Crisis:  %.3f", pre_mean_all8),
                    sprintf("Post-Crisis: %.3f", post_mean_all8)),
         col = c("#2e7d32", "#1565c0", "#c62828"), lty = 1, lwd = c(1.8, 2.2, 2.2),
         bty = "o", box.col = "gray85", bg = "white", cex = 0.72)

  # --- Panel D: Core 6 Markets Aggregate SEDCD(t) ---
  sedcd_6_s <- res_core6$sedcd
  sedcd_6_s[1:cutoff] <- NA
  ymax_s6 <- max(c(sedcd_6_s, 0.6), na.rm = TRUE) * 1.2
  pre_mean_6 <- mean(res_core6$sedcd[(cutoff+1):cp_idx_core], na.rm=TRUE)
  post_mean_6 <- mean(res_core6$sedcd[(cp_idx_core+1):N], na.rm=TRUE)

  plot(ret_dates, sedcd_6_s, type = "l", col = "#2e7d32", lwd = 1.8,
       ylim = c(0, ymax_s6), xlab = "Calendar Date", ylab = expression(widehat(SEDCD)[tau](t)),
       main = expression(bold("Panel D: Core 6 Markets Downside Density ") * bolditalic(SEDCD[tau](t))),
       cex.main = 0.95, font.main = 2, las = 1)
  grid(col = "gray88", lty = 2)
  abline(v = cp_date_core, col = "#b71c1c", lty = 2, lwd = 1.6)
  segments(ret_dates[cutoff + 1], pre_mean_6, cp_date_core, pre_mean_6, col = "#1565c0", lwd = 2.2)
  segments(cp_date_core, post_mean_6, tail(ret_dates, 1), post_mean_6, col = "#c62828", lwd = 2.2)

  legend("topleft",
         legend = c(expression(widehat(SEDCD)[tau](t)),
                    sprintf("Pre-Crisis:  %.3f", pre_mean_6),
                    sprintf("Post-Crisis: %.3f", post_mean_6)),
         col = c("#2e7d32", "#1565c0", "#c62828"), lty = 1, lwd = c(1.8, 2.2, 2.2),
         bty = "o", box.col = "gray85", bg = "white", cex = 0.72)

  mtext(expression(bold("Multivariate Linearized SEDCD CUSUM Tests: European Wholesale Electricity Markets")),
        outer = TRUE, side = 3, line = 0.8, cex = 1.15)
}

# Save Figures
pdf_out <- "../paper/img/sedcd_multivariate_all_countries.pdf"
png_out <- "../paper/img/sedcd_multivariate_all_countries.png"
png_loc <- "sedcd_multivariate_all_countries.png"
png_art <- "C:/Users/Antonio/.gemini/antigravity/brain/57f79b8c-2c63-47f6-9959-52420cf4e82f/sedcd_multivariate_all_countries.png"

plot_multivariate_comparison(pdf_out, is_pdf = TRUE)
plot_multivariate_comparison(png_out, is_pdf = FALSE)
plot_multivariate_comparison(png_loc, is_pdf = FALSE)
if (dir.exists("C:/Users/Antonio/.gemini/antigravity/brain/57f79b8c-2c63-47f6-9959-52420cf4e82f")) {
  plot_multivariate_comparison(png_art, is_pdf = FALSE)
}

cat("   [Saved] Publication PDF: ", pdf_out, "\n")
cat("   [Saved] Publication PNG: ", png_out, "\n")
cat("   [Saved] Local PNG:       ", png_loc, "\n")
cat("   [Saved] Brain Artifact:  ", png_art, "\n\n")

cat("================================================================================\n")
cat("ALL MULTIVARIATE ANALYSES COMPLETE!\n")
cat("================================================================================\n")

