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
#'   objects.
#' @param main_title The main title for the entire plot grid.
plot_prices_and_returns <- function(df, date_col_name = "Date") {

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  .set_theme() # Apply consistent theme

  price_cols <- setdiff(names(df), date_col_name)
  n_assets <- length(price_cols)

  if (n_assets == 0) {
    warning("No price columns found to plot.")
    return(invisible(NULL))
  }

  # Calculate simple returns (price changes): P_t - P_{t-1}
  # This is used because electricity prices can be negative, making log returns invalid.
  # This also matches the analysis in Application_new.R which uses diff(price_matrix).
  returns_df <- as.data.frame(lapply(df[price_cols], function(x) c(NA, diff(x))))
  names(returns_df) <- paste0("ret_", price_cols)

  dates <- df[[date_col_name]]
  # Ensure the date column is a Date object for correct plotting on the x-axis
  if (!inherits(dates, c("Date", "POSIXt"))) {
    dates <- as.Date(dates)
  }

  # Generate tick locations at 3-month intervals
  tick_locations <- seq(from = min(dates, na.rm = TRUE), to = max(dates, na.rm = TRUE), by = "3 months")

  # Adjust margins for a compact layout with no outer titles or y-axis labels.
  par(
    mfrow = c(n_assets, 2),
    oma   = c(2, 2, 2, 1),        # Outer margins for spacing
    mar   = c(3.0, 3.5, 2.5, 1.0), # c(bottom, left, top, right)
    # Match font sizes from plot_return_correlations for consistency
    cex.main  = 0.75,
    cex.axis  = 0.75,
    # Pull axis labels closer to the plot, matching other functions
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



# ── 9. plot_return_correlations ───────────────────────────────────────────────
#' Create a Pairs Plot of Asset Returns Colored by Date
#'
#' This function generates a scatterplot matrix (pairs plot) for multivariate time
#' series data (e.g., standardized residuals) using base R graphics. It is
#' designed to visualize how correlations between pairs of assets change over
#' time.
#'
#' @param df A dataframe where one column is the Date and subsequent columns are returns.
#' @param date_col_name A string specifying the name of the date column.
#' @param point_alpha Numeric; controls the transparency of the scatter points (0 to 1).
plot_return_correlations <- function(df, date_col_name = "Date", point_alpha = 0.4) {

  # --- 1. Setup plotting environment ---
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  # --- 2. Prepare data and colors ---
  if (!inherits(df[[date_col_name]], c("Date", "POSIXt"))) {
    df[[date_col_name]] <- as.Date(df[[date_col_name]])
  }
  dates <- df[[date_col_name]]
  numeric_cols <- setdiff(names(df), date_col_name)
  n_assets <- length(numeric_cols)

  if (n_assets < 2) {
    warning("Need at least two numeric columns to create a pairs plot.")
    return(invisible(NULL))
  }
  
  # Get all unique pairs of countries
  country_pairs <- combn(numeric_cols, 2, simplify = FALSE)
  n_plots <- length(country_pairs)
  n_cols <- 5
  n_rows <- 3 # We have 15 pairs from 6 countries

  # Create a color gradient based on date
  date_numeric <- as.numeric(dates)
  date_norm <- (date_numeric - min(date_numeric, na.rm = TRUE)) /
               (max(date_numeric, na.rm = TRUE) - min(date_numeric, na.rm = TRUE))

  col_ramp <- colorRampPalette(c("white", .COL_PROCESS)) # White to Blue
  gradient_colors <- col_ramp(100)
  point_colors <- gradient_colors[floor(date_norm * 99) + 1]
  point_colors_alpha <- adjustcolor(point_colors, alpha.f = point_alpha)

  # --- 3. Setup multi-panel layout with space for a legend ---
  # Create a layout with 3x5 for plots and a row at the bottom for the legend
  layout_matrix <- matrix(c(1:n_plots, rep(n_plots + 1, n_cols)), nrow = n_rows + 1, ncol = n_cols, byrow = TRUE)
  # Increase relative height of the legend panel for better readability
  layout(layout_matrix, heights = c(rep(1, n_rows), 0.4))

  # Set graphical parameters for the scatterplots
  par(
    # c(bottom, left, top, right). Increased left margin for y-label visibility
    mar       = c(3.0, 3, 1, 1.0),
    family    = .FONT_FAMILY,
    # Use consistent and readable font sizes for axes and labels
    cex.lab   = 0.75,
    cex.axis  = 0.75,
    col.axis  = .COL_AXIS,
    col.lab   = .COL_AXIS,
    tcl       = -0.28,
    mgp       = c(2.2, 0.7, 0), # Adjust label positions for readability
    bty       = "o",
    las       = 1,
    pty       = "s" # Make plots square
  )


  # --- 4. Create the grid of scatterplots ---
  for (pair in country_pairs) {
    country1 <- pair[1]
    country2 <- pair[2]

    # Create the plot with country names as axis labels
    plot(df[[country1]], df[[country2]],
         xlab = country1,
         ylab = country2,
         main = "", # No individual plot titles
         pch = 20, cex = 0.5, col = point_colors_alpha)

    # Add more visible reference lines at zero to highlight the tails
    abline(h = 0, v = 0, col = "grey50", lty = "dashed", lwd = 0.8)
  }
  
  # --- 5. Draw the color scale legend ---
  # Reduce horizontal margins to make the legend wider.
  par(mar = c(0.5, 1, 0.5, 1))
  plot(c(0, 1), c(0, 1), type = 'n', axes = FALSE, xlab = '', ylab = '')
  
  # Draw a very wide color bar
  legend_colors <- col_ramp(100)
  # Define coordinates for the bar, spanning most of the horizontal space (e.g., 5% to 95%)
  x_coords <- seq(0.05, 0.95, length.out = 101)
  
  # Draw the gradient as a series of thin rectangles
  rect(x_coords[-101], 0.4, x_coords[-1], 0.7, col = legend_colors, border = NA)
  # Draw a border around the entire color bar
  rect(min(x_coords), 0.4, max(x_coords), 0.7, border = "grey30")
  
  # Add larger Date Labels below the bar
  min_date_label <- format(min(dates, na.rm = TRUE), "%b %Y")
  max_date_label <- format(max(dates, na.rm = TRUE), "%b %Y")
  
  # Align the date labels with the edges of the color bar
  text(x = min(x_coords), y = 0.15, labels = min_date_label, cex = 1.0, col = .COL_AXIS, adj = 0)
  text(x = max(x_coords), y = 0.15, labels = max_date_label, cex = 1.0, col = .COL_AXIS, adj = 1)
  
  # Add a larger, non-bold Legend Title above the bar
  text(x = 0.5, y = 0.85, "Date", cex = 1.2, font = 1, col = .COL_AXIS)
}