# ==============================================================================
# File: plot_cpd_comparison_by_sigma.R
# Description: Generates publication-ready matrix plots for change-point tests
#              separately for each variance scenario:
#                - Scenario (i): Smooth Trending Variance (sigma(u) = exp(u))
#                - Scenario (ii): Time-Varying Piecewise Variance (sigma(u) = 1 + |u - 0.5|)
#
# Layout for each scenario:
#   - 4 Rows: Contamination Scenarios (Clean, AO, IO, RO)
#   - 3 Columns: Hypotheses (H0 [Size], H1 [Abrupt Power], H2 [Gradual Power])
#   - Inside each plot: Divided by vertical dotted line into:
#       * Left half: Gaussian Innovations (X-axis: L, H, W, S)
#       * Right half: Student-t3 Innovations (X-axis: L, H, W, S)
#   - Single-letter test labels:
#       * L = Proposed Linearized Welsh CUSUM
#       * H = Huberized CUSUM
#       * W = Wilcoxon-Mann-Whitney rank test
#       * S = Schmidt (2021) Gini Heteroscedastic Test
#   - Sample sizes n in {500, 1000, 5000}:
#       * Distinct high-contrast colors (Crimson Red for 500, Forest Green for 1000, Deep Blue for 5000)
#       * Open circle markers scaling with n (cex = 1.7, 2.4, 3.2)
#       * Horizontal dodge offsets so sample sizes are clearly side-by-side without occlusion
#   - Dashed line at nominal alpha = 0.05
#   - High-resolution PNG (300 DPI) and vector PDF
# ==============================================================================

plot_one_sigma_scenario <- function(csv_summary_file,
                                    csv_raw_file = NULL,
                                    output_prefix,
                                    main_title_expr,
                                    dodge = 0.20,
                                    copy_to_paper = TRUE,
                                    copy_to_artifacts = TRUE,
                                    artifacts_dir = "C:/Users/Antonio/.gemini/antigravity/brain/57f79b8c-2c63-47f6-9959-52420cf4e82f") {

  # 1. Load summary data (or aggregate raw if summary not found)
  if (file.exists(csv_summary_file)) {
    cat(sprintf("Loading summary data from '%s'...\n", csv_summary_file))
    df <- read.csv(csv_summary_file, stringsAsFactors = FALSE)
  } else if (!is.null(csv_raw_file) && file.exists(csv_raw_file)) {
    cat(sprintf("Summary '%s' not found. Aggregating raw '%s'...\n", csv_summary_file, csv_raw_file))
    df_raw <- read.csv(csv_raw_file, stringsAsFactors = FALSE)
    rate_cols <- intersect(c("rej_our", "rej_hl", "rej_huber", "rej_wmw", "rej_schmidt"), names(df_raw))
    group_cols <- c("hp_scenario", "contamination", "innov_dist", "var_scenario", "n")
    df <- aggregate(
      df_raw[rate_cols],
      by = df_raw[group_cols],
      FUN = function(v) if (all(is.na(v))) NA_real_ else round(mean(v, na.rm = TRUE), 4)
    )
    names(df)[names(df) %in% rate_cols] <- paste0("rate_", sub("^rej_", "", rate_cols))
    write.csv(df, csv_summary_file, row.names = FALSE)
  } else {
    stop(sprintf("Neither summary file '%s' nor raw file '%s' exists!", csv_summary_file, csv_raw_file))
  }

  cat(sprintf("Loaded dataset with %d rows.\n", nrow(df)))

  # Standardize casing
  df$hp_scenario   <- toupper(df$hp_scenario)
  df$contamination <- tolower(df$contamination)
  df$innov_dist    <- tolower(df$innov_dist)

  # Row layout: Clean, AO, IO, RO
  all_row_scens <- c("clean", "ao", "io", "ro")
  all_row_labs  <- c("Clean", "AO", "IO", "RO")
  present_scens <- unique(df$contamination)
  row_scenarios <- all_row_scens[all_row_scens %in% present_scens]
  row_labels    <- all_row_labs[all_row_scens %in% present_scens]
  n_rows        <- length(row_scenarios)

  # Column layout: Hypotheses
  col_hypotheses <- c("H0", "H1", "H2")
  col_labels     <- c(expression(bold(paste(H[0], ": No Shift (Size)"))),
                      expression(bold(paste(H[1], ": Abrupt Shift (Power)"))),
                      expression(bold(paste(H[2], ": Gradual Shift (Power)"))))

  # Test columns to display:
  has_hl      <- ("rate_hl" %in% names(df)) && !all(is.na(df$rate_hl))
  has_schmidt <- ("rate_schmidt" %in% names(df)) && !all(is.na(df$rate_schmidt))

  if (has_hl && has_schmidt) {
    test_cols   <- c("rate_our", "rate_huber", "rate_wmw", "rate_hl", "rate_schmidt")
    test_labels <- c("L", "H", "W", "T", "S")
    x_gauss     <- 1:5
    mid_sep     <- 6.0
    x_t3        <- 7:11
    x_lim       <- c(0.4, 11.6)
    gauss_at    <- 3.0
    t3_at       <- 9.0
  } else if (!has_hl && has_schmidt) {
    test_cols   <- c("rate_our", "rate_huber", "rate_wmw", "rate_schmidt")
    test_labels <- c("L", "H", "W", "S")
    x_gauss     <- 1:4
    mid_sep     <- 5.0
    x_t3        <- 6:9
    x_lim       <- c(0.4, 9.6)
    gauss_at    <- 2.5
    t3_at       <- 7.5
  } else {
    test_cols   <- c("rate_our", "rate_huber", "rate_wmw")
    test_labels <- c("L", "H", "W")
    x_gauss     <- 1:3
    mid_sep     <- 4.0
    x_t3        <- 5:7
    x_lim       <- c(0.4, 7.6)
    gauss_at    <- 2.0
    t3_at       <- 6.0
  }

  x_ticks       <- c(x_gauss, x_t3)
  x_tick_labels <- c(test_labels, test_labels)

  # Sample sizes setup: colors, scaling, and dodge offsets
  n_vals <- sort(unique(df$n))
  num_n  <- length(n_vals)

  # Color palette matching regression figures: Red (500), Green (1000), Blue (5000)
  n_color_map <- c(
    "200"   = "#fb8c00", # Amber
    "500"   = "#e53935", # Crimson Red
    "700"   = "#8e24aa", # Purple
    "1000"  = "#2e7d32", # Forest Green
    "5000"  = "#1565c0", # Deep Blue
    "10000" = "#4a148c"  # Deep Indigo
  )
  fallback_palette <- c("#e53935", "#2e7d32", "#1565c0", "#8e24aa", "#fb8c00")

  n_colors <- character(num_n)
  names(n_colors) <- as.character(n_vals)
  for (idx in seq_along(n_vals)) {
    nv_str <- as.character(n_vals[idx])
    if (nv_str %in% names(n_color_map)) {
      n_colors[nv_str] <- n_color_map[nv_str]
    } else {
      n_colors[nv_str] <- fallback_palette[((idx - 1) %% length(fallback_palette)) + 1]
    }
  }

  # Point size scaling: larger n has slightly larger point
  base_cex_map <- c("200" = 1.3, "500" = 1.7, "700" = 2.1, "1000" = 2.4, "5000" = 3.2)
  n_cexs <- numeric(num_n)
  names(n_cexs) <- as.character(n_vals)
  for (idx in seq_along(n_vals)) {
    nv_str <- as.character(n_vals[idx])
    if (nv_str %in% names(base_cex_map)) {
      n_cexs[nv_str] <- base_cex_map[nv_str]
    } else {
      n_cexs[nv_str] <- 1.4 + (idx - 1) * 0.7
    }
  }

  # Compute horizontal dodge offsets for sample sizes
  if (num_n == 1 || dodge <= 0) {
    n_offsets <- setNames(0, as.character(n_vals))
  } else if (num_n == 2) {
    n_offsets <- setNames(c(-dodge * 0.6, dodge * 0.6), as.character(n_vals))
  } else if (num_n == 3) {
    n_offsets <- setNames(c(-dodge, 0, dodge), as.character(n_vals))
  } else {
    n_offsets <- setNames(seq(-dodge, dodge, length.out = num_n), as.character(n_vals))
  }

  # Dimensions: adjust canvas height based on presence of main title
  has_title <- !is.null(main_title_expr)
  canvas_height_png <- if (n_rows == 4) {
    if (has_title) 4400 else 4200
  } else {
    if (has_title) 3400 else 3200
  }
  canvas_height_pdf <- if (n_rows == 4) {
    if (has_title) 14.6 else 14.0
  } else {
    if (has_title) 11.2 else 10.6
  }
  top_oma <- if (has_title) 4.2 else 1.2

  # Generate both PNG and PDF
  for (dev_type in c("png", "pdf")) {
    out_file <- paste0(output_prefix, ".", dev_type)
    if (dev_type == "png") {
      png(out_file, width = 3400, height = canvas_height_png, res = 300)
    } else {
      pdf(out_file, width = 11.5, height = canvas_height_pdf)
    }

    # Plot parameters: square subplots, clean margins with distinct header bands
    par(mfrow = c(n_rows, 3),
        pty   = "s",
        mar   = c(2.6, 3.4, 3.6, 0.5),
        oma   = c(0.8, 3.6, top_oma, 0.6),
        mgp   = c(2.5, 0.8, 0),
        tcl   = -0.45)

    for (row_idx in 1:n_rows) {
      cur_contam <- row_scenarios[row_idx]

      for (col_idx in 1:3) {
        cur_hp <- col_hypotheses[col_idx]

        # Y-axis scaling
        is_h0    <- (cur_hp == "H0")
        y_lim    <- if (is_h0) c(0, 0.25) else c(0, 1.05)
        y_ticks  <- if (is_h0) seq(0, 0.25, by = 0.05) else seq(0, 1.0, by = 0.25)
        y_labels <- if (is_h0) sprintf("%.2f", y_ticks) else as.character(y_ticks)

        # Initialize blank square canvas
        plot(NA, NA,
             xlim = x_lim,
             ylim = y_lim,
             xlab = "", ylab = "",
             xaxt = "n", yaxt = "n",
             bty  = "n",
             xaxs = "i", yaxs = "i")

        # Interior light dashed horizontal gridlines
        y_grid <- y_ticks[y_ticks > 0 & y_ticks < max(y_ticks)]
        if (length(y_grid) > 0) {
          abline(h = y_grid, col = "gray88", lty = 2, lwd = 1.0)
        }

        # Vertical dividing dotted line between Gaussian and Student-t3
        abline(v = mid_sep, col = "gray40", lty = 3, lwd = 1.8)

        # Horizontal nominal alpha = 0.05 significance line
        abline(h = 0.05, col = "gray20", lty = 2, lwd = 2.0)

        # Top sub-headers for Gaussian and t3
        mtext("Gaussian", side = 3, line = 0.35, at = gauss_at, cex = 1.15, font = 2, col = "black")
        mtext(expression(bold(t[3])), side = 3, line = 0.35, at = t3_at, cex = 1.35, font = 2, col = "black")

        # Plot points for Gaussian (left) and t3 (right)
        for (inn_idx in 1:2) {
          inn_name <- if (inn_idx == 1) c("gaussian", "normal") else c("t3", "student-t")
          x_base   <- if (inn_idx == 1) x_gauss else x_t3

          for (i_nv in seq_along(n_vals)) {
            nv      <- n_vals[i_nv]
            nv_str  <- as.character(nv)
            nv_off  <- n_offsets[[nv_str]]
            nv_col  <- n_colors[[nv_str]]
            nv_cex  <- n_cexs[[nv_str]]

            x_pts_n <- x_base + nv_off

            sub_df <- df[df$hp_scenario == cur_hp &
                         df$contamination == cur_contam &
                         df$innov_dist %in% inn_name &
                         df$n == nv, ]

            if (nrow(sub_df) > 0) {
              y_pts <- as.numeric(sub_df[1, test_cols])

              points(x   = x_pts_n, y = y_pts,
                     pch = 1,          # Open circle
                     col = nv_col,     # Distinct color per n
                     cex = nv_cex,     # Scaling size per n
                     lwd = 2.0)
            }
          }
        }

        # Custom informative legend in Panel (Row 1, Col 1) [Clean, H0]
        if (row_idx == 1 && col_idx == 1) {
          legend("topright",
                 legend  = paste("n =", n_vals),
                 col     = n_colors[as.character(n_vals)],
                 pch     = 1,
                 pt.lwd  = 2.2,
                 pt.cex  = n_cexs[as.character(n_vals)],
                 bty     = "o",
                 box.col = "gray85",
                 bg      = "white",
                 cex     = 1.15,
                 inset   = c(0.02, 0.04))
        }

        # X-axis with single-letter test labels
        axis(1, at = x_ticks, labels = x_tick_labels,
             cex.axis = 1.45, font.axis = 2, lwd = 0, lwd.ticks = 1.8, gap.axis = -1)

        # Y-axis ticks and labels
        axis(2, at = y_ticks, labels = y_labels,
             cex.axis = 1.5, font.axis = 1, lwd = 0, lwd.ticks = 1.8, las = 1, gap.axis = -1)

        # Square border box
        box(which = "plot", lwd = 1.8)

        # Main Hypothesis column titles on top row
        if (row_idx == 1) {
          mtext(col_labels[[col_idx]], side = 3, line = 2.0, cex = 1.40, font = 2, xpd = NA)
        }

        # Contamination scenario labels on left margin of column 1
        if (col_idx == 1) {
          mtext(row_labels[row_idx], side = 2, line = 4.0, cex = 2.2, font = 2, las = 0, xpd = NA)
        }
      }
    }

    # Outer overarching Scenario title if requested
    if (has_title) {
      mtext(main_title_expr,
            outer = TRUE, side = 3, line = 1.6, cex = 1.75, font = 2)
    }

    dev.off()
    cat(sprintf("Generated: %s\n", out_file))
  }

  # Copy generated PDF/PNG to paper/img if directory exists
  if (copy_to_paper) {
    paper_img_dir <- file.path("..", "paper", "img")
    if (dir.exists(paper_img_dir)) {
      for (ext in c("pdf", "png")) {
        src <- paste0(output_prefix, ".", ext)
        if (file.exists(src)) {
          dst <- file.path(paper_img_dir, paste0(output_prefix, ".", ext))
          file.copy(src, dst, overwrite = TRUE)
          cat(sprintf("Copied to paper/img: %s\n", dst))
        }
      }
    }
  }

  # Copy generated PNG to artifacts directory for direct user viewing
  if (copy_to_artifacts && dir.exists(artifacts_dir)) {
    src_png <- paste0(output_prefix, ".png")
    if (file.exists(src_png)) {
      dst_png <- file.path(artifacts_dir, paste0(output_prefix, ".png"))
      file.copy(src_png, dst_png, overwrite = TRUE)
      cat(sprintf("Copied to artifacts: %s\n", dst_png))
    }
  }
}

# Run generation for both Scenario i and Scenario ii (without main title)
cat("\n=== GENERATING PLOT FOR SCENARIO (i) ===\n")
plot_one_sigma_scenario(
  csv_summary_file = "sim_summary_cpd_comparison_sigma_i.csv",
  csv_raw_file     = "sim_results_cpd_comparison_sigma_i.csv",
  output_prefix    = "cpd_comparison_matrix_sigma_i",
  main_title_expr  = NULL
)

cat("\n=== GENERATING PLOT FOR SCENARIO (ii) ===\n")
plot_one_sigma_scenario(
  csv_summary_file = "sim_summary_cpd_comparison_sigma_ii.csv",
  csv_raw_file     = "sim_results_cpd_comparison_sigma_ii.csv",
  output_prefix    = "cpd_comparison_matrix_sigma_ii",
  main_title_expr  = NULL
)

cat("\nAll plots generated successfully!\n")
