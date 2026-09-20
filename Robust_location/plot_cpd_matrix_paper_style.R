# ==============================================================================
# File: plot_cpd_matrix_paper_style.R
# Description: Generates publication-ready matrix plot for change-point tests:
#              - Rows: Contamination Scenarios (Clean, AO, IO, RO)
#              - 3 Columns: Hypotheses (H0 [Size], H1 [Abrupt Power], H2 [Gradual Power])
#              - Inside each plot: Split into 2 halves by a vertical dotted line:
#                  * Left half: Gaussian Innovations (X-axis: L, H, W, T, S)
#                  * Right half: Student-t3 Innovations (X-axis: L, H, W, T, S)
#              - X-axis single-letter test labels:
#                  * L = Proposed Linearized Welsh CUSUM
#                  * H = Huberized CUSUM
#                  * W = Wilcoxon-Mann-Whitney
#                  * T = Two-sample Hodges-Lehmann
#                  * S = Schmidt (2021) Gini Heteroscedastic Test
#              - Color represents Heteroskedasticity Scenarios:
#                  * Scenario (i) Smooth trending variance (Deep Blue)
#                  * Scenario (ii) Cyclic variance (Forest Green)
#                  * Scenario (iii) Abrupt break variance (Crimson Red)
#              - Points for different sample sizes n (open circles scaling with n)
#              - Horizontal dodge so variance scenario markers do not occlude each other
#              - Dashed line at nominal alpha = 0.05
#              - Style matching plot_power_matrix.R for journal publication
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Summarize Empirical Rejection Rates (Size and Power)
# ------------------------------------------------------------------------------
if (file.exists("sim_results_cpd_comparison.csv")) {
  df_raw <- read.csv("sim_results_cpd_comparison.csv", stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
  })

  if (!"var_scenario" %in% names(df_raw)) {
    df_raw$var_scenario <- "i"
  }

  has_schmidt_raw <- "rej_schmidt" %in% names(df_raw)

  summary_df <- df_raw %>%
    group_by(hp_scenario, contamination, innov_dist, var_scenario, n) %>%
    summarise(
      reps         = n(),
      rate_our     = mean(rej_our, na.rm = TRUE),
      rate_hl      = mean(rej_hl, na.rm = TRUE),
      rate_huber   = mean(rej_huber, na.rm = TRUE),
      rate_wmw     = mean(rej_wmw, na.rm = TRUE),
      rate_schmidt = if (has_schmidt_raw) mean(rej_schmidt, na.rm = TRUE) else NA_real_,
      .groups      = "drop"
    )

  if (!has_schmidt_raw) summary_df$rate_schmidt <- NULL

  write.csv(summary_df, "sim_summary_cpd_comparison.csv", row.names = FALSE)
  cat("Saved summary rates to sim_summary_cpd_comparison.csv\n")
}

# ------------------------------------------------------------------------------
# 2. Main Publication Plotting Function
# ------------------------------------------------------------------------------

plot_cpd_matrix_paper_style <- function(csv_file = "sim_summary_cpd_comparison.csv",
                                        output_prefix = "cpd_comparison_matrix",
                                        dodge = 0.16) {

  if (!file.exists(csv_file)) {
    raw_file <- "sim_results_cpd_comparison.csv"
    if (!file.exists(raw_file)) {
      stop(sprintf("Data file '%s' not found!", csv_file))
    }
    cat(sprintf("Aggregating raw '%s'...\n", raw_file))
    raw_df <- read.csv(raw_file, stringsAsFactors = FALSE)
    if (!"var_scenario" %in% names(raw_df)) raw_df$var_scenario <- "i"

    agg_cols <- c("rej_our", "rej_hl", "rej_huber", "rej_wmw")
    if ("rej_schmidt" %in% names(raw_df)) agg_cols <- c(agg_cols, "rej_schmidt")

    formula_str <- paste("cbind(", paste(agg_cols, collapse = ", "), ") ~ hp_scenario + contamination + innov_dist + var_scenario + n")
    df <- aggregate(as.formula(formula_str), data = raw_df, FUN = function(x) round(mean(x, na.rm = TRUE), 4))

    names(df)[names(df) == "rej_our"]     <- "rate_our"
    names(df)[names(df) == "rej_hl"]      <- "rate_hl"
    names(df)[names(df) == "rej_huber"]   <- "rate_huber"
    names(df)[names(df) == "rej_wmw"]     <- "rate_wmw"
    if ("rej_schmidt" %in% names(df)) names(df)[names(df) == "rej_schmidt"] <- "rate_schmidt"
  } else {
    df <- read.csv(csv_file, stringsAsFactors = FALSE)
    if (!"var_scenario" %in% names(df)) df$var_scenario <- "i"
  }

  cat(sprintf("Loaded data with %d rows.\n", nrow(df)))

  # Standardize casing & alias mapping
  df$hp_scenario   <- toupper(df$hp_scenario)
  df$contamination <- tolower(df$contamination)
  df$innov_dist    <- tolower(df$innov_dist)
  df$var_scenario  <- tolower(as.character(df$var_scenario))
  df$var_scenario  <- ifelse(df$var_scenario %in% c("i", "trending", "smooth", "exp"), "i",
                      ifelse(df$var_scenario %in% c("ii", "cyclic", "sin"), "ii",
                      ifelse(df$var_scenario %in% c("iii", "break", "abrupt"), "iii", df$var_scenario)))

  # Layout configurations: support Clean, AO, IO, and RO
  all_row_scens <- c("clean", "ao", "io", "ro")
  all_row_labs  <- c("Clean", "AO", "IO", "RO")
  present_scens <- unique(df$contamination)

  if ("ro" %in% present_scens) {
    row_scenarios <- all_row_scens
    row_labels    <- all_row_labs
  } else {
    row_scenarios <- c("clean", "ao", "io")
    row_labels    <- c("Clean", "AO", "IO")
  }
  n_rows <- length(row_scenarios)

  col_hypotheses <- c("H0", "H1", "H2")
  col_labels     <- c(expression(bold(paste(H[0], ": No Shift (Size)"))),
                      expression(bold(paste(H[1], ": Abrupt Shift (Power)"))),
                      expression(bold(paste(H[2], ": Gradual Shift (Power)"))))

  # Tests and single-letter codes: include Schmidt test S if present
  has_schmidt <- "rate_schmidt" %in% names(df) && !all(is.na(df$rate_schmidt))
  if (has_schmidt) {
    test_cols   <- c("rate_our", "rate_huber", "rate_wmw", "rate_hl", "rate_schmidt")
    test_labels <- c("L", "H", "W", "T", "S")
    x_gauss     <- 1:5
    mid_sep     <- 6.0
    x_t3        <- 7:11
    x_lim       <- c(0.4, 11.6)
    gauss_at    <- 3.0
    t3_at       <- 9.0
  } else {
    test_cols   <- c("rate_our", "rate_huber", "rate_wmw", "rate_hl")
    test_labels <- c("L", "H", "W", "T")
    x_gauss     <- 1:4
    mid_sep     <- 5.0
    x_t3        <- 6:9
    x_lim       <- c(0.4, 9.6)
    gauss_at    <- 2.5
    t3_at       <- 7.5
  }
  x_ticks       <- c(x_gauss, x_t3)
  x_tick_labels <- c(test_labels, test_labels)

  # Heteroskedasticity Scenarios setup: colors and dodging
  var_scenarios_order <- c("i", "ii", "iii")
  var_colors <- c(
    "i"   = "#1565c0", # Scenario (i) Smooth trending (Deep Blue)
    "ii"  = "#2e7d32", # Scenario (ii) Cyclic (Forest Green)
    "iii" = "#c62828"  # Scenario (iii) Abrupt break (Crimson Red)
  )

  present_vars <- var_scenarios_order[var_scenarios_order %in% unique(df$var_scenario)]
  if (length(present_vars) == 0) present_vars <- unique(df$var_scenario)

  # Compute horizontal offsets for dodge
  if (length(present_vars) <= 1 || dodge <= 0) {
    var_offsets <- setNames(rep(0, length(present_vars)), present_vars)
  } else if (length(present_vars) == 2) {
    var_offsets <- setNames(c(-dodge * 0.7, dodge * 0.7), present_vars)
  } else {
    var_offsets <- c("i" = -dodge, "ii" = 0, "iii" = dodge)
  }

  # Sample sizes and scaling
  n_vals <- sort(unique(df$n))
  base_cexs   <- c("200" = 1.3, "500" = 2.1, "700" = 2.6, "1000" = 3.0, "5000" = 4.0)
  cex_mapping <- setNames(seq(1.3, by = 0.8, length.out = length(n_vals)), as.character(n_vals))
  for (nv_str in names(base_cexs)) {
    if (nv_str %in% names(cex_mapping)) cex_mapping[nv_str] <- base_cexs[nv_str]
  }

  # Graphics driver setup
  canvas_height_png <- if (n_rows == 4) 4200 else 3200
  canvas_height_pdf <- if (n_rows == 4) 14.0 else 10.8

  for (dev_type in c("png", "pdf")) {
    out_name <- paste0(output_prefix, ".", dev_type)
    if (dev_type == "png") {
      png(out_name, width = 3400, height = canvas_height_png, res = 300)
    } else {
      pdf(out_name, width = 11.5, height = canvas_height_pdf)
    }

    par(mfrow = c(n_rows, 3),
        pty   = "s",
        mar   = c(2.8, 3.4, 1.8, 0.5),
        oma   = c(0.6, 3.4, 3.6, 0.6),
        mgp   = c(2.5, 0.8, 0),
        tcl   = -0.45)

    for (row_idx in 1:n_rows) {
      cur_contam <- row_scenarios[row_idx]

      for (col_idx in 1:3) {
        cur_hp <- col_hypotheses[col_idx]

        # Y-axis setup
        is_h0    <- (cur_hp == "H0")
        y_lim    <- if (is_h0) c(0, 0.25) else c(0, 1.05)
        y_ticks  <- if (is_h0) seq(0, 0.25, by = 0.05) else seq(0, 1.0, by = 0.25)
        y_labels <- if (is_h0) sprintf("%.2f", y_ticks) else as.character(y_ticks)

        # Base empty plot
        plot(NA, NA,
             xlim = x_lim,
             ylim = y_lim,
             xlab = "", ylab = "",
             xaxt = "n", yaxt = "n",
             bty  = "n",
             xaxs = "i", yaxs = "i")

        # Interior light dashed gridlines
        y_grid <- y_ticks[y_ticks > 0 & y_ticks < max(y_ticks)]
        if (length(y_grid) > 0) {
          abline(h = y_grid, col = "gray88", lty = 2, lwd = 1.0)
        }

        # Vertical dividing dotted line in the middle
        abline(v = mid_sep, col = "gray40", lty = 3, lwd = 1.8)

        # Horizontal nominal alpha = 0.05 line
        abline(h = 0.05, col = "gray20", lty = 2, lwd = 2.0)

        # Top labels for Gaussian and t3
        mtext("Gaussian", side = 3, line = 0.35, at = gauss_at, cex = 1.25, font = 2, col = "black")
        mtext(expression(bold(t[3])), side = 3, line = 0.35, at = t3_at, cex = 1.45, font = 2, col = "black")

        # Plot data points
        for (inn_idx in 1:2) {
          inn_name <- if (inn_idx == 1) c("gaussian", "normal") else c("t3", "student-t")
          x_base   <- if (inn_idx == 1) x_gauss else x_t3

          for (v_scen in present_vars) {
            v_offset   <- if (v_scen %in% names(var_offsets)) var_offsets[[v_scen]] else 0
            v_col      <- if (v_scen %in% names(var_colors)) var_colors[[v_scen]] else "blue4"
            x_pts_scen <- x_base + v_offset

            sub_df <- df[df$hp_scenario == cur_hp &
                         df$contamination == cur_contam &
                         df$innov_dist %in% inn_name &
                         df$var_scenario == v_scen, ]

            if (nrow(sub_df) > 0) {
              for (i in rev(seq_along(n_vals))) {
                nv    <- n_vals[i]
                row_n <- sub_df[sub_df$n == nv, ]

                if (nrow(row_n) > 0) {
                  y_pts <- as.numeric(row_n[1, test_cols])

                  points(x   = x_pts_scen, y = y_pts,
                         pch = 1,              # Open circle without fill
                         col = v_col,          # Color according to variance scenario
                         cex = cex_mapping[as.character(nv)],
                         lwd = 1.8)
                }
              }
            }
          }
        }

        # In Panel (Row 1, Col 1) [H0: Clean], add informative legends
        if (row_idx == 1 && col_idx == 1) {
          legend("topleft",
                 legend  = c(expression(sigma[(i)] ~ "(Trending)"),
                             expression(sigma[(ii)] ~ "(Cyclic)"),
                             expression(sigma[(iii)] ~ "(Break)")),
                 col     = c("#1565c0", "#2e7d32", "#c62828"),
                 pch     = 1,
                 pt.lwd  = 2.0,
                 pt.cex  = 1.6,
                 bty     = "o",
                 box.col = "gray85",
                 bg      = "white",
                 cex     = 1.05,
                 inset   = c(0.02, 0.04))

          legend("topright",
                 legend  = paste("n =", n_vals),
                 col     = "gray30",
                 pch     = 1,
                 pt.lwd  = 2.0,
                 pt.cex  = cex_mapping[as.character(n_vals)],
                 bty     = "o",
                 box.col = "gray85",
                 bg      = "white",
                 cex     = 1.05,
                 inset   = c(0.02, 0.04))
        }

        # X-axis with tick labels
        axis(1, at = x_ticks, labels = x_tick_labels,
             cex.axis = 1.4, font.axis = 2, lwd = 0, lwd.ticks = 1.8, gap.axis = -1)

        # Y-axis ticks
        axis(2, at = y_ticks, labels = y_labels,
             cex.axis = 1.5, font.axis = 1, lwd = 0, lwd.ticks = 1.8, las = 1, gap.axis = -1)

        # Uniform outer square box border
        box(which = "plot", lwd = 1.8)

        # Main Hypothesis titles on top row
        if (row_idx == 1) {
          title(main = col_labels[[col_idx]], line = 2.8, cex.main = 2.1, font.main = 2, xpd = NA)
        }

        # Scenario titles on left margin of column 1
        if (col_idx == 1) {
          mtext(row_labels[row_idx], side = 2, line = 3.8, cex = 2.2, font = 2, las = 0, xpd = NA)
        }
      }
    }

    dev.off()
    cat(sprintf("Saved: %s\n", out_name))
  }

  # Copy PDF to paper img directory if directory exists
  paper_img_dir <- file.path("..", "paper", "img")
  if (dir.exists(paper_img_dir) && file.exists(paste0(output_prefix, ".pdf"))) {
    dest <- file.path(paper_img_dir, paste0(output_prefix, ".pdf"))
    file.copy(paste0(output_prefix, ".pdf"), dest, overwrite = TRUE)
    cat(sprintf("Copied %s.pdf to %s\n", output_prefix, dest))
  }
}

if (!interactive()) {
  plot_cpd_matrix_paper_style()
}
