# ─────────────────────────────────────────────────────────────────────────────
# Plot_functions.R
# Consistent plotting style for CUSUM / local-density article figures.
#
# Design
#   * Full rectangle box around plot area  (bty = "o")
#   * Compact font sizes  (cex ~ 0.75–0.88)
#   * Horizontal legend placed BELOW the panel
#   * Font: Inter via showtext; falls back to "sans"
# ─────────────────────────────────────────────────────────────────────────────

# ── 0. Font ───────────────────────────────────────────────────────────────────
.FONT_FAMILY <- "sans"

if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts",  quietly = TRUE)) {
  tryCatch({
    sysfonts::font_add_google("Inter", "inter")
    showtext::showtext_auto()
    showtext::showtext_opts(dpi = 300)
    .FONT_FAMILY <- "inter"
    message("[Plot_functions] Font: Inter via showtext.")
  }, error = function(e) {
    message("[Plot_functions] Inter unavailable; using 'sans'.")
  })
} else {
  message("[Plot_functions] showtext not installed; using 'sans'.")
}

# ── 1. Palette ────────────────────────────────────────────────────────────────
.COL_PROCESS  <- "#2166ac"   # steel blue
.COL_MEAN     <- "#d73027"   # red
.COL_QUANTILE <- "#a50026"   # dark red
.COL_ALT      <- "#4dac26"   # green
.COL_AXIS     <- "grey10"

.LWD_MAIN  <- 1.5
.LWD_REF   <- 1.2
.LTY_MEAN  <- 2
.LTY_QUANT <- 3

# ── 2. Theme ──────────────────────────────────────────────────────────────────
.set_theme <- function() {
  par(
    family    = .FONT_FAMILY,
    cex.main  = 0.55,        # title
    cex.lab   = 0.45,        # axis labels
    cex.axis  = 0.45,        # axis tick labels
    font.main = 1,           # plain (not bold)
    col.axis  = .COL_AXIS,
    col.lab   = .COL_AXIS,
    col.main  = "grey10",
    tcl       = -0.28,
    mgp       = c(1.8, 0.40, 0),
    bty       = "o",         # full rectangle around plot area
    las       = 1
  )
}
.set_theme()

# ── 3. Legend helper ──────────────────────────────────────────────────────────
# Horizontal legend placed below the panel.
# Requires the calling function to have set mar[1] >= 6.5 first.
.below_legend <- function(labels, cols, ltys, lwds) {
  legend(
    x       = "bottom",
    inset   = c(0, -0.3),
    legend  = labels,
    col     = cols,
    lty     = ltys,
    lwd     = lwds,
    bty     = "n",
    horiz   = TRUE,
    xpd     = NA,
    cex     = 0.45,
    seg.len = 5
  )
}

# Widen bottom margin for the legend row and restore on exit.
.with_legend_margin <- function(expr) {
  old <- par(mar = c(7.5, 4.0, 3.0, 2.0))
  on.exit(par(old), add = TRUE)
  force(expr)
}



# ── 6. plot_cusum_test ────────────────────────────────────────────────────────
# Full CUSUM plot with bootstrap critical bands.
plot_cusum_test <- function(u_seq,
                            T_n_u,
                            q95,
                            q90        = NULL,
                            p_value    = NULL,
                            main_title = "CUSUM test for density stability") {
  
  ymax       <- max(abs(c(T_n_u, q95))) * 1.2
  leg_labels <- c("CUSUM process", "H\u2080 reference", "95% critical value")
  leg_cols   <- c(.COL_PROCESS, .COL_MEAN, .COL_QUANTILE)
  leg_ltys   <- c(1, .LTY_MEAN, .LTY_QUANT)
  leg_lwds   <- c(.LWD_MAIN, .LWD_REF, .LWD_REF)
  
  if (!is.null(q90)) {
    leg_labels <- c(leg_labels, "90% critical value")
    leg_cols   <- c(leg_cols, adjustcolor(.COL_QUANTILE, 0.55))
    leg_ltys   <- c(leg_ltys, .LTY_QUANT)
    leg_lwds   <- c(leg_lwds, .LWD_REF)
  }
  
  .with_legend_margin({
    # 1. Initialize Plot
    plot(u_seq, T_n_u,
         type = "n", # Create empty plot first to layer shading underneath
         main = main_title,
         xlab = "Time proportion (u)",
         ylab = expression(T[n](u)),
         ylim = c(-ymax, ymax))
    
    # 2. Add THE SHADING (FILL GAPS)
    # Define a light, transparent red
    fill_col <- adjustcolor("red", alpha.f = 0.2)
    
    # Upper Gap: Shade between q95 and T_n_u (only where T_n_u > q95)
    polygon(c(u_seq, rev(u_seq)), 
            c(pmax(T_n_u, q95), rep(q95, length(u_seq))), 
            col = fill_col, border = NA)
    
    # Lower Gap: Shade between -q95 and T_n_u (only where T_n_u < -q95)
    polygon(c(u_seq, rev(u_seq)), 
            c(pmin(T_n_u, -q95), rep(-q95, length(u_seq))), 
            col = fill_col, border = NA)
    
    # 3. Add Reference Lines
    abline(h = 0, col = .COL_MEAN, lty = .LTY_MEAN, lwd = .LWD_REF)
    abline(h = c(q95, -q95), col = .COL_QUANTILE, lty = .LTY_QUANT, lwd = .LWD_REF)
    
    if (!is.null(q90))
      abline(h = c(q90, -q90),
             col = adjustcolor(.COL_QUANTILE, 0.55),
             lty = .LTY_QUANT, lwd = .LWD_REF)
    
    # 4. Draw the actual process line on top
    lines(u_seq, T_n_u, lwd = .LWD_MAIN, col = .COL_PROCESS)
    
    if (!is.null(p_value))
      mtext(sprintf("p = %.3f", p_value),
            side = 3, adj = 1, cex = 0.72, col = "grey40")
    
    .below_legend(leg_labels, leg_cols, leg_ltys, leg_lwds)
  })
}


# ── 7. plot_stat_histograms_by_level ──────────────────────────────────────────
plot_pvalue_histograms_fabian <- function(sim_results,
                                          main_title_prefix = "P-value Distribution",
                                          n_breaks = 10,
                                          n_to_plot = NULL,
                                          dgp_to_plot = NULL,
                                          tau_to_plot = NULL,
                                          Ln_to_plot = NULL) {
  # Save original graphical parameters to restore them on exit
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  # If all parameters for a single plot are specified, plot only that one.
  is_single_plot <- !is.null(n_to_plot) && !is.null(dgp_to_plot) &&
                    !is.null(tau_to_plot) && !is.null(Ln_to_plot)

  if (is_single_plot) {
    sub_df <- sim_results %>%
      filter(n_values == n_to_plot, dgp_name == dgp_to_plot,
             tau == tau_to_plot, L_n_type == Ln_to_plot)
    
    if (nrow(sub_df) == 0) {
      warning("No data for the specified single combination.")
      return(invisible(NULL))
    }
    
    # Setup for a single plot
    par(pty = "s", mar = c(4.5, 4.5, 3, 1))
    breaks_fixed <- seq(0, 1, length.out = n_breaks + 1)
    pvals <- sub_df$stat[!is.na(sub_df$stat)]
    M_obs <- length(pvals)
    h <- hist(pvals, breaks = breaks_fixed, plot = FALSE)
    h$counts <- h$counts / sum(h$counts)
    rej_rate <- mean(pvals <= 0.05)
    
    plot(h, freq = TRUE, col = adjustcolor(.COL_PROCESS, 0.70), border = "white",
         xlim = c(0, 1), ylim = c(0, max(max(h$counts), 1/n_breaks) * 1.4),
         main = sprintf("%s: n=%d, tau=%s, Ln=%s", dgp_to_plot, n_to_plot, tau_to_plot, Ln_to_plot),
         xlab = "Bootstrap p-value", ylab = "Relative frequency", axes = FALSE)
    axis(1, at = seq(0, 1, by = 0.2)); axis(2, las = 1); box()
    abline(h = 1 / n_breaks, col = .COL_MEAN, lty = .LTY_MEAN, lwd = .LWD_REF)
    legend("topleft", legend = c(sprintf("Rate (5%%) = %.3f", rej_rate), sprintf("M = %d", M_obs)),
           col = c(.COL_MEAN, NA), lty = c(.LTY_MEAN, NA), lwd = c(.LWD_REF, NA),
           bty = "n", cex = 0.8, inset = c(-0.02, 0))
    return(invisible(NULL))
  }

  # --- Multi-panel plot logic ---
  n_vals <- sort(unique(sim_results$n_values))

  # Loop over each sample size to create a separate plot for each
  for (n_val in n_vals) {
    
    sub_df_n <- sim_results %>% filter(n_values == n_val)
    
    dgp_names <- sort(unique(sub_df_n$dgp_name))
    tau_vals  <- sort(unique(sub_df_n$tau))
    ln_types  <- sort(unique(sub_df_n$L_n_type))
    
    # Rows: DGP x tau, Columns: L_n_type
    rows <- length(dgp_names) * length(tau_vals)
    cols <- length(ln_types)

    par(mfrow = c(rows, cols), pty = "s", oma = c(4, 6, 5, 1), mar = c(2, 2, 2, 1))

    breaks_fixed <- seq(0, 1, length.out = n_breaks + 1)

    # Loop to create each panel
    for (dgp in dgp_names) {
      for (tau_val in tau_vals) {
        for (ln_type in ln_types) {
          
          pvals <- sub_df_n$stat[sub_df_n$dgp_name == dgp & sub_df_n$tau == tau_val & sub_df_n$L_n_type == ln_type]
          pvals <- pvals[!is.na(pvals)]
          
          M_obs <- length(pvals)
          
          if (M_obs == 0) {
            plot.new(); box(); text(0.5, 0.5, "No data", cex=0.8); next
          }
          
          h <- hist(pvals, breaks = breaks_fixed, plot = FALSE)
          h$counts <- h$counts / sum(h$counts)
          
          plot(h, freq = TRUE, col = adjustcolor(.COL_PROCESS, 0.70), border = "white",
               xlim = c(0, 1), ylim = c(0, max(max(h$counts), 1/n_breaks) * 1.4),
               main = "", xlab = "", ylab = "", axes = FALSE)
          
          axis(1, at = seq(0, 1, by = 0.5), cex.axis = 0.7, mgp = c(3, 0.1, 0))
          axis(2, las = 1, cex.axis = 0.7, mgp = c(3, 0.4, 0))
          box()
          abline(h = 1 / n_breaks, col = .COL_MEAN, lty = .LTY_MEAN, lwd = .LWD_REF)
          
          # Add column titles (L_n_type) to the top row of plots
          if (dgp == dgp_names[1] && tau_val == tau_vals[1]) {
            mtext(ln_type, side = 3, line = 0.5, cex = 0.8)
          }
          # Add row titles (DGP, tau) to the first column of plots
          if (ln_type == ln_types[1]) {
            mtext(paste(dgp, ", tau=", tau_val, sep=""), side = 2, line = 3, cex = 0.8, las = 1)
          }
        }
      }
    }

    mtext(paste(main_title_prefix, "| n =", n_val), outer = TRUE, font = 2, cex = 1.2, col = "grey10")
    mtext("Bootstrap p-value", side = 1, outer = TRUE, line = 2.5, cex = 0.9)
    mtext("Relative Frequency", side = 2, outer = TRUE, line = 3.5, cex = 0.9, las = 0)
  }
}


# ── 8. plot_prices_and_returns ────────────────────────────────────────────────
#' Plot Time Series and Returns
#'
#' Creates a multi-panel plot showing the time series of prices and their
#' corresponding returns for multiple assets. The layout will be N x 2, where N
#' is the number of assets, with prices on the left and returns on the right.
#'
#' @param df A data frame containing the time series data. Must have a date
#'   column and one or more numeric price columns.
#' @param date_col_name The name of the column containing dates or date-time
#'   objects. Defaults to "Date".
#' @param filename The name of the output PDF file.
#' @param pdf_width The width of the PDF in inches.
#' @param pdf_height The height of the PDF in inches. If `NULL` (the default),
#'   the height is calculated based on the number of assets.
#' @param units The unit for `pdf_width` and `pdf_height`. Can be "in" (inches, default) or "pt" (points).
#'
plot_prices_and_returns <- function(df, date_col_name = "Date", filename = "prices_and_returns.pdf", pdf_width = 7, pdf_height = NULL, units = "in") {
  
  # Ensure the date ticks are in English
  original_locale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "English")
  on.exit(Sys.setlocale("LC_TIME", original_locale), add = TRUE)

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  .set_theme() # Apply consistent theme

  price_cols <- setdiff(names(df), date_col_name)
  n_assets <- length(price_cols)

  if (n_assets == 0) {
    warning("No price columns found to plot.")
    return(invisible(NULL))
  }

  # Determine the base height in the specified units
  base_height_in_units <- if (is.null(pdf_height)) n_assets * 3.5 else pdf_height

  # Convert dimensions to inches for the pdf() function
  if (units == "pt") {
    conversion_factor <- 72 # 1 inch = 72 points
    actual_pdf_width_in_inches <- pdf_width / conversion_factor
    actual_pdf_height_in_inches <- base_height_in_units / conversion_factor
  } else if (units == "in") {
    actual_pdf_width_in_inches <- pdf_width
    actual_pdf_height_in_inches <- base_height_in_units
  } else {
    warning("Unsupported unit. Using inches. Please use 'in' or 'pt'.")
    actual_pdf_width_in_inches <- pdf_width
    actual_pdf_height_in_inches <- base_height_in_units
  }

  # --- Save to PDF ---
  # Determine height if not provided, maintaining the scaling logic
  pdf(filename, width = actual_pdf_width_in_inches, height = actual_pdf_height_in_inches)
  on.exit(dev.off(), add = TRUE)
  
  # Calculate simple returns
  returns_df <- as.data.frame(lapply(df[price_cols], function(x) c(NA, diff(x))))
  names(returns_df) <- paste0("ret_", price_cols)

  dates <- df[[date_col_name]]
  if (!inherits(dates, c("Date", "POSIXt"))) {
    dates <- as.Date(dates)
  }

  tick_locations <- seq(from = min(dates, na.rm = TRUE), to = max(dates, na.rm = TRUE), by = "3 months")

  par(
    mfrow = c(n_assets, 2),
    oma   = c(2, 2, 2, 1),
    mar   = c(3.0, 2, 2.5, 0),
    cex.main  = 0.75,
    cex.axis  = 0.75,
    mgp       = c(2.2, 0.7, 0)
  )

  for (asset in price_cols) {
    ret_col <- paste0("ret_", asset)

    # Plot Price Time Series
    plot(dates, df[[asset]], type = 'l', col = .COL_PROCESS, lwd = .LWD_MAIN,
         xlab = "", ylab = "", main = paste(asset, "Prices"), xaxt = "n")
    axis.Date(1, at = tick_locations, format = "%b-%y")

    # Plot Returns Time Series
    plot(dates, returns_df[[ret_col]], type = 'l', col = .COL_PROCESS, lwd = .LWD_MAIN,
         xlab = "", ylab = "", main = paste(asset, "Returns"), xaxt = "n")
    axis.Date(1, at = tick_locations, format = "%b-%y")
    abline(h = 0, col = .COL_MEAN, lty = .LTY_MEAN, lwd = .LWD_REF)
  }

}

# ── 9. plot_electricity_ts_returns ───────────────────────────────────────────
#' Plot Historical Electricity Prices and Daily Returns for Core Countries
#'
#' Creates a 6x2 multi-panel publication figure displaying historical prices
#' (left column) and daily returns / price differences (right column) for the
#' core 6 European electricity markets.
#' No titles or y-axis labels. X-axis labels only at the bottom row.
#'
#' @param df Data frame containing Date and country price series
#' @param core_countries Vector of country codes (default: c("DE_LU", "IT_NORD", "FR", "BE", "NL", "PL"))
#' @param filename Optional output filename (.pdf or .png)
#' @param is_pdf Logical, whether output is PDF (TRUE) or PNG (FALSE)
plot_electricity_ts_returns <- function(df,
                                        core_countries = c("DE_LU", "IT_NORD", "FR", "BE", "NL", "PL"),
                                        filename = NULL,
                                        is_pdf = FALSE) {
  
  if (!is.null(filename)) {
    if (is_pdf) {
      pdf(filename, width = 11.5, height = 11.0)
      on.exit(dev.off(), add = TRUE)
    } else {
      png(filename, width = 3450, height = 3300, res = 300)
      on.exit(dev.off(), add = TRUE)
    }
  }

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  dates <- as.Date(substr(df$Date, 1, 10))
  price_mat <- as.matrix(df[, core_countries])
  ret_mat <- diff(price_mat)
  ret_dates <- dates[-1]

  year_dates <- as.Date(c("2020-01-01", "2021-01-01", "2022-01-01", "2023-01-01", 
                          "2024-01-01", "2025-01-01", "2026-01-01"))
  year_labels <- c("2020", "2021", "2022", "2023", "2024", "2025", "2026")

  # Layout: 6 rows x 2 columns
  # No outer titles. Y-labels removed. Dates removed. X-labels only at the bottom 2 plots.
  par(mfrow = c(6, 2),
      oma   = c(1.5, 0.8, 1.0, 0.8),
      mar   = c(2.2, 3.6, 1.3, 1.0),
      mgp   = c(2.0, 0.6, 0),
      tcl   = -0.35)

  col_line <- "#1565c0"
  col_grid <- "gray90"
  n_rows <- length(core_countries)

  for (i in seq_along(core_countries)) {
    c_code <- core_countries[i]
    is_last_row <- (i == n_rows)
    
    # --- Column 1: Historical Prices ---
    p_series <- price_mat[, i]
    y_lim_p <- extendrange(p_series, f = 0.12)
    
    plot(dates, p_series, type = "n",
         xlim = range(dates), ylim = y_lim_p,
         xlab = if (is_last_row) "Date" else "",
         ylab = "",
         font.lab = 2, cex.lab = 1.0, las = 1, xaxt = "n")
    grid(col = col_grid, lty = 2, lwd = 0.8)
    
    if (is_last_row) {
      axis.Date(1, at = year_dates, labels = year_labels, cex.axis = 0.95)
    } else {
      axis.Date(1, at = year_dates, labels = FALSE, cex.axis = 0.95)
    }
    
    # Trajectory line
    lines(dates, p_series, col = col_line, lwd = 1.1)
    
    # Subpanel header (only panel identifier, no dates, no break marker)
    text(dates[1] + 15, y_lim_p[2] * 0.98, labels = paste0("Historical prices (", c_code, ")"),
         font = 2, cex = 1.0, adj = c(0, 1), col = "gray10")
    
    box(which = "plot", lty = "solid", lwd = 1.3, col = "black")
    
    # --- Column 2: Daily Returns ---
    r_series <- ret_mat[, i]
    y_lim_r <- extendrange(r_series, f = 0.12)
    
    plot(ret_dates, r_series, type = "n",
         xlim = range(dates), ylim = y_lim_r,
         xlab = if (is_last_row) "Date" else "",
         ylab = "",
         font.lab = 2, cex.lab = 1.0, las = 1, xaxt = "n")
    grid(col = col_grid, lty = 2, lwd = 0.8)
    
    if (is_last_row) {
      axis.Date(1, at = year_dates, labels = year_labels, cex.axis = 0.95)
    } else {
      axis.Date(1, at = year_dates, labels = FALSE, cex.axis = 0.95)
    }
    
    # Baseline & trajectory line
    abline(h = 0, col = "gray75", lty = 1, lwd = 0.8)
    lines(ret_dates, r_series, col = col_line, lwd = 0.85)
    
    # Subpanel header (only panel identifier, no dates, no break marker)
    text(dates[1] + 15, y_lim_r[2] * 0.98, labels = paste0("Daily returns (", c_code, ")"),
         font = 2, cex = 1.0, adj = c(0, 1), col = "gray10")
    
    box(which = "plot", lty = "solid", lwd = 1.3, col = "black")
  }

  invisible(NULL)
}

# ── 10. plot_residual_bivariate_scatter ───────────────────────────────────────
#' Plot Pairwise Bivariate Residual Scatter Plots with Time Progression Hue
#'
#' Generates a 3x5 multi-panel grid of bivariate scatter plots for all 15 pairs
#' of the 6 core electricity markets. Points are colored with a continuous time
#' hue (light sky blue in 2020 to deep navy in 2026) to illustrate dynamic co-dependence.
#'
#' @param df Data frame containing Date and price/return series
#' @param core_countries Vector of country codes (default: c("DE_LU", "IT_NORD", "FR", "BE", "NL", "PL"))
#' @param residuals_type "mscale" (robust Dennis-Welsch M-scale, default), "garch", or "returns"
#' @param point_alpha Point opacity (default 0.45)
#' @param filename Optional output filename (.pdf or .png)
#' @param is_pdf Logical, whether output is PDF or PNG
plot_residual_bivariate_scatter <- function(df,
                                            core_countries = c("DE_LU", "IT_NORD", "FR", "BE", "NL", "PL"),
                                            residuals_type = c("mscale", "garch", "returns"),
                                            point_alpha = 0.45,
                                            filename = NULL,
                                            is_pdf = FALSE) {
  
  residuals_type <- match.arg(residuals_type)
  
  if (!is.null(filename)) {
    if (is_pdf) {
      pdf(filename, width = 12.0, height = 7.8)
      on.exit(dev.off(), add = TRUE)
    } else {
      png(filename, width = 3600, height = 2340, res = 300)
      on.exit(dev.off(), add = TRUE)
    }
  }

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  dates <- as.Date(substr(df$Date, 1, 10))
  price_mat <- as.matrix(df[, core_countries])
  ret_mat <- diff(price_mat)
  ret_dates <- dates[-1]
  N <- nrow(ret_mat)
  d <- ncol(ret_mat)

  # Compute standardized residuals
  if (residuals_type == "mscale") {
    Z_mat <- matrix(0, nrow = N, ncol = d)
    colnames(Z_mat) <- core_countries
    for (j in 1:d) {
      diff_j <- ret_mat[, j] - median(ret_mat[, j])
      mad_j  <- median(abs(diff_j)) / 0.6745
      if (is.na(mad_j) || mad_j < 1e-4) mad_j <- sd(diff_j)
      sig_j  <- if (exists("solve_mscale_fast", mode = "function")) {
        solve_mscale_fast(diff_j, mad_j, c = 2.985)
      } else {
        mad_j
      }
      Z_mat[, j] <- diff_j / sig_j
    }
  } else if (residuals_type == "garch" && requireNamespace("rugarch", quietly = TRUE)) {
    Z_mat <- matrix(0, nrow = N, ncol = d)
    colnames(Z_mat) <- core_countries
    spec <- rugarch::ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
                                distribution.model = "std")
    for (j in 1:d) {
      fit <- rugarch::ugarchfit(spec = spec, data = ret_mat[, j], solver = "hybrid")
      Z_mat[, j] <- as.numeric(rugarch::residuals(fit, standardize = TRUE))
    }
  } else {
    Z_mat <- scale(ret_mat)
    colnames(Z_mat) <- core_countries
  }

  # Time hue gradient
  date_num <- as.numeric(ret_dates)
  date_norm <- (date_num - min(date_num)) / (max(date_num) - min(date_num))
  col_ramp <- colorRampPalette(c("#90caf9", "#42a5f5", "#1976d2", "#0d47a1", "#051d40"))
  palette_100 <- col_ramp(100)
  point_colors <- palette_100[pmin(100, pmax(1, floor(date_norm * 99) + 1))]
  point_colors_alpha <- adjustcolor(point_colors, alpha.f = point_alpha)

  # Layout: 15 pairwise plots (3x5) + colorbar panel
  pairs_list <- combn(core_countries, 2, simplify = FALSE)
  n_pairs <- length(pairs_list)

  layout_mat <- matrix(c(1:15, rep(16, 5)), nrow = 4, ncol = 5, byrow = TRUE)
  layout(layout_mat, heights = c(1, 1, 1, 0.28))

  z_range <- c(-10, 14)
  ticks_at <- c(-10, -5, 0, 5, 10)

  par(mar = c(3.2, 3.4, 1.2, 1.0),
      mgp = c(2.0, 0.6, 0),
      tcl = -0.35,
      pty = "s")

  for (k in 1:n_pairs) {
    c1 <- pairs_list[[k]][1]
    c2 <- pairs_list[[k]][2]
    
    plot(Z_mat[, c1], Z_mat[, c2],
         xlim = z_range, ylim = z_range,
         xlab = c1, ylab = c2,
         font.lab = 2, cex.lab = 1.05,
         las = 1, bty = "n",
         xaxt = "n", yaxt = "n",
         pch = 20, cex = 0.55,
         col = point_colors_alpha)
    
    grid(col = "gray92", lty = 2, lwd = 0.6)
    abline(h = 0, col = "gray55", lty = 2, lwd = 0.9)
    abline(v = 0, col = "gray55", lty = 2, lwd = 0.9)
    points(Z_mat[, c1], Z_mat[, c2], pch = 20, cex = 0.55, col = point_colors_alpha)
    
    axis(1, at = ticks_at, labels = ticks_at, cex.axis = 0.85)
    axis(2, at = ticks_at, labels = ticks_at, cex.axis = 0.85, las = 1)
    box(which = "plot", lty = "solid", lwd = 1.3, col = "black")
  }

  # --- Colorbar panel at bottom ---
  # Sleek horizontal color bar without title or break marker
  par(mar = c(1.4, 6.0, 0.6, 6.0), pty = "m")
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")

  x_left  <- 0.08
  x_right <- 0.92
  n_bars  <- 150
  x_steps <- seq(x_left, x_right, length.out = n_bars + 1)
  bar_cols <- col_ramp(n_bars)

  for (b in 1:n_bars) {
    rect(x_steps[b], 0.38, x_steps[b + 1], 0.78, col = bar_cols[b], border = NA)
  }
  rect(x_left, 0.38, x_right, 0.78, border = "gray30", lwd = 1.2)

  # Explicit English date labels
  month_abbrs <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  fmt_en <- function(d) paste(month_abbrs[as.integer(format(d, "%m"))], format(d, "%Y"))

  text(x_left, 0.16, labels = fmt_en(min(ret_dates)), font = 2, cex = 0.95, adj = c(0, 1), col = "gray20")
  text(x_right, 0.16, labels = fmt_en(tail(ret_dates, 1)), font = 2, cex = 0.95, adj = c(1, 1), col = "gray20")

  invisible(NULL)
}

# Alias for backwards compatibility
plot_return_correlations <- function(df, ...) {
  plot_residual_bivariate_scatter(df, ...)
}
