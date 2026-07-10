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
                                          n_breaks = 10) {
  # Save original graphical parameters to restore them on exit
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  # Identify all unique combinations of non-grid parameters (tau, L_n_type)
  param_combos <- sim_results %>%
    select(tau, L_n_type) %>%
    distinct()

  # Loop over each parameter combination to create a separate plot for each
  for (i in 1:nrow(param_combos)) {
    current_tau <- param_combos$tau[i]
    current_Ln <- param_combos$L_n_type[i]

    # Filter the main dataframe for the current combination
    sub_df <- sim_results %>%
      filter(tau == current_tau, L_n_type == current_Ln)

    # Get the DGPs and Sample sizes in a fixed order for this subset
    dgp_names <- sort(unique(sub_df$dgp_name))
    n_vals    <- sort(unique(sub_df$n_values))

    # Vertical DGPs (rows), Horizontal n values (columns)
    rows <- length(dgp_names)
    cols <- length(n_vals)

    # Setup multi-panel layout for the current plot
    par(
      mfrow = c(rows, cols),
      pty   = "s",               # Keep panels square
      oma   = c(2, 0, 4.5, 0),   # Outer margin for main title
      mar   = c(3.5, 3.5, 2, 1)  # Tighter margins
    )

    breaks_fixed <- seq(0, 1, length.out = n_breaks + 1)

    # Loop over DGPs and sample sizes to create each panel
    for (dgp in dgp_names) {
      for (n_val in n_vals) {
        
        pvals <- sub_df$stat[sub_df$dgp_name == dgp & sub_df$n_values == n_val]
        pvals <- pvals[!is.na(pvals)]
        
        M_obs <- length(pvals)
        
        if (M_obs == 0) {
           plot.new()
           title(main = sprintf("%s | n = %d", dgp, n_val))
           text(0.5, 0.5, "No data")
           next
        }
        
        h <- hist(pvals, breaks = breaks_fixed, plot = FALSE)
        h$counts <- h$counts / sum(h$counts) # Convert to relative frequency
        
        rej_rate <- mean(pvals <= 0.05)
        
        plot(h,
           freq    = TRUE,
           col     = adjustcolor(.COL_PROCESS, 0.70),
           border  = "white",
           xlim    = c(0, 1),
           ylim    = c(0, max(max(h$counts), 1/n_breaks) * 1.4),
           main    = sprintf("%s | n = %d", dgp, n_val),
           xlab    = "Bootstrap p-value",
           ylab    = "Relative frequency",
           axes    = FALSE)
      
      axis(1, at = seq(0, 1, by = 0.2))
      axis(2, las = 1)
      box()
      
      abline(h   = 1 / n_breaks,
             col = .COL_MEAN,
             lty = .LTY_MEAN,
             lwd = .LWD_REF)
      
      legend("topleft", 
             legend = c(sprintf("Rate (5%%) = %.3f", rej_rate),
                        sprintf("M = %d", M_obs)),
             col    = c(.COL_MEAN, NA),
             lty    = c(.LTY_MEAN, NA),
             lwd    = c(.LWD_REF, NA),
             bty    = "n",
             cex    = 0.45,            
             inset  = c(-0.02, 0))
      }
    }

    # Construct and add a dynamic main title for the entire plot
    main_title <- sprintf(
      "%s (Tau = %s, Ln = %s)",
      main_title_prefix,
      current_tau,
      current_Ln
    )
    mtext(main_title, outer = TRUE, font = 2, cex = 0.65, col = "grey10")
  }
}