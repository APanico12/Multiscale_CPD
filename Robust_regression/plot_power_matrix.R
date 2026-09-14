#!/usr/bin/env Rscript
# ==============================================================================
# File: plot_power_matrix.R
# Description: Generates publication-ready 2x4 minimal square power curve matrices:
#              - Twice bigger, prominent axis ticks (tcl = -0.6, lwd.ticks = 2.0)
#              - Rows: Abrupt Break (H1) & Gradual Break (H2)
#              - Columns: Clean, AO, IO, RO
#              - Separate figures for Gaussian and Student-t3 noise
#              - Panel [1, 1] displays axis labels: X label = delta, Y label = Power
#              - Panel [1, 1] displays sample size n legend
#              - Panel [2, 1] displays test procedure (Robust M-Test vs Classical L2) legend
#              - Alpha = 0.05 nominal size reference line
#              - Square panels (pty = "s")
# ==============================================================================

plot_power_matrix <- function(csv_file = "sim_summary_power_n01.csv",
                              output_prefix = "power_curves_matrix") {
  
  # Smart file resolution for both project root and Robust_regression directory
  find_file <- function(fname) {
    cands <- c(
      fname,
      file.path("Robust_regression", fname),
      file.path("Multiscale_CPD", "Robust_regression", fname)
    )
    for (f in cands) {
      if (file.exists(f)) return(f)
    }
    return(NULL)
  }
  
  resolved_csv <- find_file(csv_file)
  
  if (is.null(resolved_csv)) {
    raw_csv <- find_file("sim_results_power.csv")
    if (!is.null(raw_csv)) {
      cat(sprintf("Reading raw '%s' and aggregating...\n", raw_csv))
      raw_df <- read.csv(raw_csv, stringsAsFactors = FALSE)
      df <- aggregate(cbind(rej_lin, rej_l2) ~ model_type + n + innov_dist + contamination + epsilon + hp_scenario + delta,
                      data = raw_df, FUN = function(x) round(mean(x, na.rm = TRUE), 4))
      names(df)[names(df) == "rej_lin"] <- "rate_linearized"
      names(df)[names(df) == "rej_l2"]  <- "rate_l2"
    } else {
      stop(sprintf("Data file '%s' not found in current directory or Robust_regression!", csv_file))
    }
  } else {
    df <- read.csv(resolved_csv, stringsAsFactors = FALSE)
  }
  
  # Ensure standard column names
  col_names <- names(df)
  lin_col <- if ("rate_linearized" %in% col_names) "rate_linearized" else if ("rej_lin" %in% col_names) "rej_lin" else NULL
  l2_col  <- if ("rate_l2" %in% col_names) "rate_l2" else if ("rej_l2" %in% col_names) "rej_l2" else NULL
  
  if (is.null(lin_col) || is.null(l2_col)) {
    stop("Could not find rejection rate columns in summary data!")
  }
  
  # Design elements
  col_scenarios <- c("Clean", "AO", "IO", "RO")
  col_titles    <- c("Clean", "AO", "IO", "RO")
  row_scenarios <- c("H1", "H2")
  row_titles    <- c("Abrupt Break (H1)", "Gradual Break (H2)")
  
  # Colors for sample sizes matching reference image: red, green, blue
  n_vals <- sort(unique(df$n))
  n_colors <- c(
    "500"  = "#e53935", # Crimson Red
    "1000" = "#2e7d32", # Forest Green
    "5000" = "#1565c0"  # Deep Blue
  )
  palette_base <- c("#1565c0", "#2e7d32", "#e53935", "#8e24aa")
  for (i in seq_along(n_vals)) {
    n_str <- as.character(n_vals[i])
    if (!n_str %in% names(n_colors)) {
      n_colors[n_str] <- palette_base[(i - 1) %% length(palette_base) + 1]
    }
  }
  
  # Innovation distributions to generate
  target_innovs <- list(
    list(id = "normal", tag = "gaussian", name = "Gaussian Innovations", file_tag = "gaussian"),
    list(id = "t3",     tag = "t3",       name = "Student-t(3) Innovations", file_tag = "t3")
  )
  
  # Internal plotting function for one 2x4 matrix
  draw_matrix <- function(innov_id, innov_name) {
    inn_filter <- if (tolower(innov_id) %in% c("t3", "t_3")) c("t3", "t_3") else c("normal", "gaussian")
    inn_df <- df[tolower(df$innov_dist) %in% inn_filter, ]
    
    if (nrow(inn_df) == 0) {
      cat(sprintf("Warning: No data for innovation distribution '%s'\n", innov_id))
      return(NULL)
    }
    
    # Global graphical setup: square panels (pty = "s")
    # Twice bigger tick marks: tcl = -0.6 (twice length of -0.3)
    par(mfrow = c(2, 4), 
        pty   = "s",
        mar   = c(4.8, 4.8, 3.2, 0.8),
        oma   = c(1.2, 5.2, 1.2, 0.8),
        mgp   = c(3.6, 1.6, 0),
        tcl   = -0.6)
    
    y_ticks <- seq(0, 1, by = 0.25)
    
    for (row_idx in 1:2) {
      cur_hp <- row_scenarios[row_idx]
      
      for (col_idx in 1:4) {
        cur_contam <- col_scenarios[col_idx]
        
        # Filter sub-data for panel
        sub_p <- inn_df[toupper(inn_df$hp_scenario) == cur_hp & 
                        toupper(inn_df$contamination) == toupper(cur_contam), ]
        
        # Ticks and limits
        delta_vals <- sort(unique(sub_p$delta))
        if (length(delta_vals) == 0) delta_vals <- c(0, 0.25, 0.5, 0.75, 1, 1.25)
        x_labels   <- as.character(delta_vals)
        y_labels   <- as.character(y_ticks)
        
        # Axis labels: only delta on first plot (Row 1, Col 1); "Power" removed as requested
        is_first_panel <- (row_idx == 1 && col_idx == 1)
        cur_xlab <- if (is_first_panel) expression(delta) else ""
        cur_ylab <- ""
        
        # Base empty plot with square box (bty = "n" so final box() controls uniform frame)
        plot(NA, NA, 
             xlim = c(min(delta_vals), max(delta_vals)), 
             ylim = c(0, 1.0),
             xlab = cur_xlab, ylab = cur_ylab,
             cex.lab = 2.2, font.lab = 2, col.lab = "black",
             xaxt = "n", yaxt = "n",
             bty = "n",
             xaxs = "r", yaxs = "i")
        
        # Interior light dashed gridlines matching tick levels (strictly inside boundaries)
        x_grid <- delta_vals[delta_vals > min(delta_vals) & delta_vals < max(delta_vals)]
        y_grid <- y_ticks[y_ticks > 0 & y_ticks < 1.0]
        abline(v = x_grid, col = "gray88", lty = 2, lwd = 0.9)
        abline(h = y_grid, col = "gray88", lty = 2, lwd = 0.9)
        
        # Nominal alpha = 0.05 line
        abline(h = 0.05, lty = 2, col = "gray20", lwd = 2.0)
        
        # Axis ticks: formatted as (0, 0.25, 0.5, 0.75, 1, 1.25) without float on 0, 1 and 0.5
        axis(1, at = delta_vals, labels = x_labels,
             cex.axis = 2.2, gap.axis = 0, lwd = 0, lwd.ticks = 2.0)
        axis(2, at = y_ticks, labels = y_labels,
             las = 1, cex.axis = 2.2, gap.axis = 0, lwd = 0, lwd.ticks = 2.0)
        
        # Scenario titles on top row (Row 1) in bold black, scaled proportionally
        if (row_idx == 1) {
          title(main = col_titles[col_idx], col.main = "black", font.main = 2, cex.main = 3, line = 1.1)
        }
        
        # Plot curves for each sample size n
        for (ni in seq_along(n_vals)) {
          cur_n   <- n_vals[ni]
          row_sub <- sub_p[sub_p$n == cur_n, ]
          if (nrow(row_sub) == 0) next
          row_sub <- row_sub[order(row_sub$delta), ]
          
          c_col <- n_colors[as.character(cur_n)]
          
          # 1. Proposed Linearized Robust M-Test: Solid line, filled circle
          lines(row_sub$delta, row_sub[[lin_col]], 
                col = c_col, lty = 1, lwd = 2.2)
          points(row_sub$delta, row_sub[[lin_col]], 
                 col = c_col, pch = 16, cex = 1.8)
          
          # 2. Classical L2 Test: Dashed line, triangle
          lines(row_sub$delta, row_sub[[l2_col]], 
                col = c_col, lty = 2, lwd = 2.2)
          points(row_sub$delta, row_sub[[l2_col]], 
                 col = c_col, pch = 17, cex = 1.8)
        }
        
        # Legend 1: In Panel [1, 1] (Row 1, Col 1 -> H1 / Clean) display sample sizes n
        if (row_idx == 1 && col_idx == 1) {
          leg_labels <- sapply(n_vals, function(nv) as.expression(bquote(n == .(nv))))
          leg_cols   <- sapply(as.character(n_vals), function(nv) n_colors[nv])
          legend("bottomright", 
                 legend   = leg_labels,
                 col      = leg_cols,
                 lty      = 1,
                 lwd      = 2.2,
                 pch      = 16,
                 pt.cex   = 1.25,
                 bty      = "o",
                 box.col  = "gray35",
                 box.lwd  = 1.0,
                 bg       = "white",
                 cex      = 2.0,
                 inset    = c(0.03, 0.03))
        }
        
        # Legend 2: In Panel [2, 1] (Row 2, Col 1 -> H2 / Clean) display test procedures
        if (row_idx == 2 && col_idx == 1) {
          legend("topleft",
                 legend   = c(expression(psi), expression(L[2])),
                 col      = "gray20",
                 lty      = c(1, 2),
                 lwd      = c(2.2, 1.8),
                 pch      = c(16, 17),
                 pt.cex   = c(1.25, 1.15),
                 bty      = "o",
                 box.col  = "gray35",
                 box.lwd  = 1.0,
                 bg       = "white",
                 cex      = 2.0,
                 inset    = c(0.03, 0.03))
        }
        
        # Final bounding frame on top of all lines and curves (solid, uniform lwd = 2.0 on all 4 sides)
        box(which = "plot", lty = "solid", lwd = 2.0, col = "black")
      }
    }
    
    # Outer labels:
    # Left vertical labels for rows in bold black, scaled proportionally
    mtext("Abrupt Break (H1)", side = 2, line = 3.0, at = 0.75, outer = TRUE, 
          col = "black", font = 2, cex = 2.0)
    mtext("Gradual Break (H2)", side = 2, line = 3.0, at = 0.25, outer = TRUE, 
          col = "black", font = 2, cex = 2.0)
  }
  
  # Determine output directory (where csv is or current directory)
  out_dir <- if (!is.null(resolved_csv)) dirname(resolved_csv) else "."
  if (out_dir == "") out_dir <- "."
  
  # Generate both Gaussian and t3 figures
  for (t_inn in target_innovs) {
    pdf_path <- file.path(out_dir, sprintf("%s_%s.pdf", output_prefix, t_inn$file_tag))
    png_path <- file.path(out_dir, sprintf("%s_%s.png", output_prefix, t_inn$file_tag))
    
    # 1. Vector PDF (16x8.8 inches for square panels and spacious layout)
    pdf(pdf_path, width = 16.0, height = 8.8)
    draw_matrix(t_inn$id, t_inn$name)
    dev.off()
    cat(sprintf("Successfully saved vector PDF: %s\n", pdf_path))
    
    # 2. High-resolution PNG (300 DPI)
    png(png_path, width = 4800, height = 2640, res = 300)
    draw_matrix(t_inn$id, t_inn$name)
    dev.off()
    cat(sprintf("Successfully saved high-res PNG: %s\n", png_path))
  }
}

# Automatically execute when run interactively or via Rscript if data exists
has_summary <- file.exists("sim_summary_power_n01.csv") || file.exists(file.path("Robust_regression", "sim_summary_power_n01.csv"))
has_raw     <- file.exists("sim_results_power_n01.csv") || file.exists(file.path("Robust_regression", "sim_results_powern01.csv"))

if (has_summary || has_raw) {
  cat("Auto-running plot_power_matrix()...\n")
  plot_power_matrix()
}
