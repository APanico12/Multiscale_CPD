# =============================================================================
# WBS.R: Wild Binary Segmentation for Robust Linearized M-Estimators
#
# Implements:
#   - Linearized pseudo-observation generator \hat{Z}_t
#   - Vectorized CUSUM contrast evaluation over random intervals
#   - Full Wild Binary Segmentation (WBS) solution path
#   - Model Selection via Robust Strengthened Schwarz Information Criterion (sSIC)
#   - Segment level estimation: Robust M-estimator (default), Median, or Mean
#   - Diagnostic changepoint visualization with segment levels
# =============================================================================

if (file.exists("CUSUM.R")) {
  source("CUSUM.R")
}

# -----------------------------------------------------------------------------
# 1. Linearized Pseudo-Observation Generator
# -----------------------------------------------------------------------------

#' Compute linearized pseudo-observations Z_hat[t] across the full time series
#'
#' Evaluates: \hat{Z}_t = pilot + (H(X_t, pilot) / DH_pilot)
#'
#' @param x     numeric vector of observations
#' @param k     rolling window size for pilot estimator (integer or exponent < 1)
#' @param lag   decoupling lag L_n
#' @param c     Welsch/Tukey tuning constant (default: 2.985)
#' @param loss  "Welsh", "Tukey", or "L2"
#' @return numeric vector of linearized observations of length(x)
get_linearized_series <- function(x, k = 0.45, lag = NULL, c = 2.985,
                                  loss = c("Welsh", "Tukey", "L2"),
                                  teta_pilot = NULL) {
  loss <- match.arg(loss)
  N <- length(x)
  if (k < 1) k_win <- floor(N^k) else k_win <- floor(k)
  k_win <- max(3, k_win)

  if (is.null(lag)) lag <- max(1, floor(0.1 * (log(N))^2))
  if (lag < 1) lag <- max(1, ceiling(N^lag))

  is_l2 <- (loss == "L2")
  psi_fn   <- if (loss == "Welsh") Welsh.psi else if (loss == "Tukey") tukey_loss_derivative else function(r, c) r
  psi_p_fn <- if (loss == "Welsh") Welsh.psi.prime else if (loss == "Tukey") tukey_loss_2nd_derivative else function(r, c) rep(1, length(r))

  # Pilot estimates from rolling window
  if (is.null(teta_pilot)) {
    teta_pilot <- get_teta(x, k = k_win, c = c, loss = loss)
  }

  Z_hat <- numeric(N)
  ws <- k_win + lag

  # Initial burn-in initialization using first stable pilot
  init_pilot <- teta_pilot[k_win]
  r_init <- x[1:k_win] - init_pilot
  sig_init <- median(abs(r_init)) / 0.6745
  if (is.na(sig_init) || sig_init < 1e-5) sig_init <- 1.0

  DH_init <- if (is_l2) 1.0 else mean(psi_p_fn(r_init / sig_init, c)) / sig_init
  if (is.na(DH_init) || DH_init < 1e-4) DH_init <- 1e-4

  for (t in 1:(ws - 1)) {
    if (is_l2) {
      Z_hat[t] <- x[t]
    } else {
      H_val <- psi_fn((x[t] - init_pilot) / sig_init, c)
      Z_hat[t] <- init_pilot + (H_val / DH_init)
    }
  }

  # Linearized updates for t >= ws
  for (t in ws:N) {
    pilot <- teta_pilot[t - lag]
    start_idx <- max(1, t - k_win - lag + 1)
    end_idx   <- t - lag
    r_hist    <- x[start_idx:end_idx] - pilot

    if (is_l2) {
      Z_hat[t] <- x[t]
    } else {
      sigma_hat <- median(abs(r_hist)) / 0.6745
      if (is.na(sigma_hat) || sigma_hat < 1e-5) sigma_hat <- 1.0

      DH <- mean(psi_p_fn(r_hist / sigma_hat, c)) / sigma_hat
      if (is.na(DH) || DH < 1e-4) DH <- 1e-4

      H_future <- psi_fn((x[t] - pilot) / sigma_hat, c)
      Z_hat[t] <- pilot + (H_future / DH)
    }
  }

  return(Z_hat)
}

# -----------------------------------------------------------------------------
# 2. Vectorized CUSUM Contrast Calculation over an Interval
# -----------------------------------------------------------------------------

#' Fast CUSUM contrast calculation using precomputed cumulative sums
#'
#' @param cum_y     numeric vector c(0, cumsum(y))
#' @param s         start index
#' @param e         end index
#' @param min_dist  minimum distance from boundaries
#' @return list(max_contrast, best_b)
fast_cusum_contrast <- function(cum_y, s, e, min_dist = 5) {
  n <- e - s + 1
  if (n < 2 * min_dist) {
    return(list(max_contrast = 0, best_b = NA))
  }

  b_candidates <- (s + min_dist - 1):(e - min_dist)
  n1 <- b_candidates - s + 1
  n2 <- e - b_candidates

  # Local partial sums
  S1 <- cum_y[b_candidates + 1] - cum_y[s]
  S2 <- cum_y[e + 1] - cum_y[b_candidates + 1]

  contrasts <- abs(sqrt(n2 / (n * n1)) * S1 - sqrt(n1 / (n * n2)) * S2)

  max_idx <- which.max(contrasts)
  return(list(
    max_contrast = contrasts[max_idx],
    best_b       = b_candidates[max_idx]
  ))
}

# -----------------------------------------------------------------------------
# 3. Wild Binary Segmentation Core Recursion
# -----------------------------------------------------------------------------

#' Recursive WBS path explorer
wbs_recursive <- function(s, e, cum_y, intervals, min_dist = 5) {
  # Find drawn intervals entirely contained in [s, e]
  sub_idx <- which(intervals[, 1] >= s & intervals[, 2] <= e)

  best_contrast <- -1
  best_b <- NA

  if (length(sub_idx) > 0) {
    for (m in sub_idx) {
      sm <- intervals[m, 1]
      em <- intervals[m, 2]
      res_m <- fast_cusum_contrast(cum_y, sm, em, min_dist = min_dist)
      if (!is.na(res_m$best_b) && res_m$max_contrast > best_contrast) {
        best_contrast <- res_m$max_contrast
        best_b <- res_m$best_b
      }
    }
  }

  # Augment with evaluation on [s, e] directly
  res_full <- fast_cusum_contrast(cum_y, s, e, min_dist = min_dist)
  if (!is.na(res_full$best_b) && res_full$max_contrast > best_contrast) {
    best_contrast <- res_full$max_contrast
    best_b <- res_full$best_b
  }

  if (is.na(best_b) || best_contrast <= 0) {
    return(NULL)
  }

  node <- data.frame(
    b        = best_b,
    contrast = best_contrast,
    s        = s,
    e        = e
  )

  left_tree  <- if (best_b - s >= 2 * min_dist) wbs_recursive(s, best_b, cum_y, intervals, min_dist) else NULL
  right_tree <- if (e - best_b >= 2 * min_dist) wbs_recursive(best_b + 1, e, cum_y, intervals, min_dist) else NULL

  return(rbind(node, left_tree, right_tree))
}

# -----------------------------------------------------------------------------
# 4. Robust sSIC Model Selection
# -----------------------------------------------------------------------------

#' Select optimal number of changepoints using Robust sSIC
#'
#' @param x            raw time series vector
#' @param candidates   ordered candidate change-points
#' @param Kmax         maximum change-points to consider
#' @param alpha        sSIC penalty exponent (default: 1.01)
#' @param c            tuning constant
#' @param loss         loss function ("Welsh", "Tukey", or "L2")
#' @param segment_level estimator for segment level ("robust", "median", "mean")
#' @return list(k_opt, cpts, ssic_table, segment_levels)
model_selection_ssic <- function(x, candidates, Kmax = 20, alpha = 1.01,
                                 c = 2.985, loss = "Welsh",
                                 segment_level = c("robust", "median", "mean")) {
  segment_level <- match.arg(segment_level)
  N <- length(x)
  K_eval <- min(Kmax, length(candidates))

  # Global robust scale estimate
  sigma_global <- median(abs(diff(x))) / (sqrt(2) * 0.6745)
  if (is.na(sigma_global) || sigma_global < 1e-5) sigma_global <- 1.0

  rho_fn <- if (loss == "Welsh") Welsh.rho else function(u, c) pmin(u^2, c^2)

  ssic_values <- numeric(K_eval + 1)
  loss_values <- numeric(K_eval + 1)

  # Function to estimate level on segment
  est_level <- function(seg) {
    if (segment_level == "robust") {
      solve_teta(seg, c = c, loss = loss, method = "IRLS")
    } else if (segment_level == "median") {
      median(seg)
    } else {
      mean(seg)
    }
  }

  for (k in 0:K_eval) {
    cpts_k <- if (k == 0) integer(0) else sort(candidates[1:k])
    bounds <- c(0, cpts_k, N)

    total_loss <- 0
    for (j in 1:(length(bounds) - 1)) {
      seg_idx <- (bounds[j] + 1):bounds[j + 1]
      seg_data <- x[seg_idx]
      level_hat <- est_level(seg_data)

      r_norm <- (seg_data - level_hat) / sigma_global
      if (loss == "L2") {
        total_loss <- total_loss + sum((seg_data - level_hat)^2)
      } else {
        total_loss <- total_loss + sum(rho_fn(r_norm, c = c))
      }
    }

    if (loss == "L2") {
      sig_k_sq <- total_loss / N
    } else {
      # Standardized empirical robust variance: 2 * sigma^2 * mean(rho)
      sig_k_sq <- (2 * (sigma_global^2) * total_loss) / N
    }

    sig_k_sq <- max(sig_k_sq, 1e-8)
    loss_values[k + 1] <- sig_k_sq

    # sSIC(k) = (N / 2) * log(sigma_k^2) + k * (log(N))^alpha
    penalty <- k * (log(N)^alpha)
    ssic_values[k + 1] <- (N / 2) * log(sig_k_sq) + penalty
  }

  best_k <- which.min(ssic_values) - 1
  best_cpts <- if (best_k == 0) integer(0) else sort(candidates[1:best_k])

  # Compute final segment levels
  final_bounds <- c(0, best_cpts, N)
  levels <- numeric(length(final_bounds) - 1)
  for (j in 1:(length(final_bounds) - 1)) {
    seg_data <- x[(final_bounds[j] + 1):final_bounds[j + 1]]
    levels[j] <- est_level(seg_data)
  }

  ssic_table <- data.frame(
    k           = 0:K_eval,
    rob_var     = round(loss_values, 5),
    sSIC        = round(ssic_values, 3)
  )

  return(list(
    k_opt          = best_k,
    cpts           = best_cpts,
    ssic_table     = ssic_table,
    segment_levels = levels
  ))
}

# -----------------------------------------------------------------------------
# 5. Main User-Facing Function: WBS.mean
# -----------------------------------------------------------------------------

#' Wild Binary Segmentation for Robust Location Changepoints
#'
#' @param x             numeric data vector
#' @param M             number of random intervals (default: 5000)
#' @param Kmax          maximum candidate change-points for sSIC path (default: 20)
#' @param selection     "sSIC" (default) or "threshold"
#' @param threshold_C   threshold multiplier if selection = "threshold" (default: 1.0)
#' @param alpha         sSIC penalty exponent (default: 1.01)
#' @param c             loss tuning parameter (default: 2.985)
#' @param loss          "Welsh", "Tukey", or "L2"
#' @param k             rolling window size for pilot linearization (default: 0.45)
#' @param lag           decoupling lag L_n
#' @param linearized    use linearized pseudo-observations if TRUE (default: TRUE)
#' @param min_dist      minimum segment spacing
#' @param segment_level method for segment levels: "robust" (M-estimator), "median", or "mean"
#' @param plot          plot the segmentation results if TRUE (default: TRUE)
#' @param main          plot title
#' @param dates         optional vector of dates (Date, character, or POSIXct)
#' @param plot          plot the segmentation results if TRUE (default: TRUE)
#' @param main          plot title
#' @param ylab          y-axis label
#' @param xlab          x-axis label (defaults to "Date" if dates provided, else "Time (t)")
#' @param flag_outliers whether to detect and highlight severe isolated outliers (default: TRUE)
#' @param outlier_thresh threshold on Hampel filter score for outliers (default: 4.0)
#' @return an object of class 'wbs' containing changepoints, levels, fitted values, dates, and outliers
WBS.mean <- function(x, M = 5000, Kmax = 20,
                     selection = c("sSIC", "threshold"),
                     threshold_C = 1.0, alpha = 1.01,
                     c = 2.985, loss = c("Welsh", "Tukey", "L2"),
                     k = 0.45, lag = NULL, linearized = TRUE,
                     use_cv = FALSE, cv_grid = c(0.45, 0.65),
                     min_dist = NULL,
                     segment_level = c("robust", "median", "mean"),
                     dates = NULL,
                     plot = TRUE, main = NULL,
                     ylab = expression(bold(paste("Observed ", italic(X[t])))),
                     xlab = NULL,
                     flag_outliers = TRUE, outlier_thresh = 4.0,
                     cex_break = 1.15, date_format = "%Y-%m-%d", min_label_gap = 0.15) {

  loss <- match.arg(loss)
  selection <- match.arg(selection)
  segment_level <- match.arg(segment_level)
  x <- as.numeric(x)
  N <- length(x)

  # Parse dates if provided
  dates_vec <- if (!is.null(dates)) as.Date(dates) else NULL
  if (!is.null(dates_vec) && length(dates_vec) != N) {
    warning("Length of 'dates' does not match length of 'x'. Ignoring dates.")
    dates_vec <- NULL
  }

  # Predictive cross-validation for bandwidth selection if requested
  cv_info  <- NULL
  teta_opt <- NULL
  if (linearized && (isTRUE(use_cv) || (is.character(k) && tolower(k) == "cv"))) {
    cv_info  <- cv_optimal_bandwidth_location(x = x, k_grid = cv_grid, lag = lag, loss = loss, c = c)
    k        <- cv_info$k_opt
    teta_opt <- cv_info$teta_opt
  }

  if (is.null(min_dist)) min_dist <- max(5, floor(N^0.15))

  # 1. Transform via linearization if selected
  if (linearized) {
    y_series <- get_linearized_series(x, k = k, lag = lag, c = c, loss = loss, teta_pilot = teta_opt)
  } else {
    y_series <- x
  }

  cum_y <- c(0, cumsum(y_series))

  # 2. Draw M random intervals [s_m, e_m]
  s_draws <- sample.int(N - 2 * min_dist, size = M, replace = TRUE)
  len_max <- N - s_draws
  len_draws <- mapply(function(lm) sample(2 * min_dist:lm, 1), len_max)
  e_draws <- s_draws + len_draws - 1
  intervals <- cbind(s_draws, e_draws)

  # 3. Explore tree of candidate change-points
  tree_res <- wbs_recursive(1, N, cum_y = cum_y, intervals = intervals, min_dist = min_dist)

  if (is.null(tree_res) || nrow(tree_res) == 0) {
    message("No candidate changepoints found satisfying minimum spacing.")
    cpts_final <- integer(0)
    levels_final <- if (segment_level == "robust") solve_teta(x, c = c, loss = loss) else median(x)
    res_obj <- list(
      x = x, dates = dates_vec, cpts = cpts_final, cpts_dates = NULL,
      k_opt = 0, segment_levels = levels_final, ssic_table = NULL,
      linearized_series = y_series, fitted_values = rep(levels_final, N),
      segment_level = segment_level, loss = loss
    )
    class(res_obj) <- c("wbs", "list")
    return(res_obj)
  }

  # Sort candidates by descending contrast magnitude
  tree_res <- tree_res[order(tree_res$contrast, decreasing = TRUE), ]
  candidates <- unique(tree_res$b)

  # 4. Model selection
  if (selection == "sSIC") {
    sel_res <- model_selection_ssic(
      x             = x,
      candidates    = candidates,
      Kmax          = Kmax,
      alpha         = alpha,
      c             = c,
      loss          = loss,
      segment_level = segment_level
    )
    cpts_final   <- sel_res$cpts
    k_opt        <- sel_res$k_opt
    levels_final <- sel_res$segment_levels
    ssic_table   <- sel_res$ssic_table

  } else {
    # Thresholding selection
    sig_mad <- median(abs(diff(y_series))) / (sqrt(2) * 0.6745)
    zeta <- threshold_C * sig_mad * sqrt(2 * log(N))

    sig_nodes <- tree_res[tree_res$contrast > zeta, ]
    cpts_final <- if (nrow(sig_nodes) > 0) sort(unique(sig_nodes$b)) else integer(0)
    k_opt <- length(cpts_final)

    bounds <- c(0, cpts_final, N)
    levels_final <- numeric(length(bounds) - 1)
    for (j in 1:(length(bounds) - 1)) {
      seg_data <- x[(bounds[j] + 1):bounds[j + 1]]
      levels_final[j] <- if (segment_level == "robust") {
        solve_teta(seg_data, c = c, loss = loss)
      } else if (segment_level == "median") {
        median(seg_data)
      } else {
        mean(seg_data)
      }
    }
    ssic_table <- NULL
  }

  # 5. Build fitted curve
  bounds <- c(0, cpts_final, N)
  fit_curve <- numeric(N)
  for (j in 1:(length(bounds) - 1)) {
    fit_curve[(bounds[j] + 1):bounds[j + 1]] <- levels_final[j]
  }

  cpts_dates <- if (!is.null(dates_vec) && length(cpts_final) > 0) dates_vec[cpts_final] else NULL

  res_obj <- list(
    x                 = x,
    dates             = dates_vec,
    cpts              = cpts_final,
    cpts_dates        = cpts_dates,
    k_opt             = k_opt,
    segment_levels    = levels_final,
    ssic_table        = ssic_table,
    linearized_series = y_series,
    fitted_values     = fit_curve,
    segment_level     = segment_level,
    loss              = loss,
    bandwidth_k       = k,
    cv_info           = cv_info
  )
  class(res_obj) <- c("wbs", "list")

  # 6. Plotting
  if (plot) {
    plot_wbs(res_obj, dates = dates_vec, main = main, ylab = ylab, xlab = xlab,
             flag_outliers = flag_outliers, outlier_thresh = outlier_thresh,
             cex_break = cex_break, date_format = date_format, min_label_gap = min_label_gap)
  }

  return(res_obj)
}

#' Plot WBS Segmentation Results (Publication Quality matching DGP.R Style)
#'
#' @param wbs_res        result object from WBS.mean()
#' @param dates          vector of dates (optional, defaults to wbs_res$dates if available)
#' @param main           plot title
#' @param ylab           y-axis label
#' @param xlab           x-axis label (defaults to "Date" if dates present, else "Time (t)")
#' @param flag_outliers  whether to detect and highlight severe outliers with red open circles
#' @param outlier_thresh Hampel filter score threshold for severe outliers (default: 4.0)
#' @param level_col      color for estimated segment levels (default: "#1565c0", deep blue)
#' @param level_lwd      line width for segment levels (default: 2.8)
#' @param break_col      color for changepoint vertical lines (default: "#c62828", crimson red)
#' @param break_lwd      line width for changepoint lines (default: 2.0)
#' @param break_lty      line type for changepoints (default: 2, dashed)
#' @param outlier_col    color for severe outliers (default: "#c62828")
#' @param series_col     color for observed series (default: "gray25")
#' @param series_lwd     line width for observed series (default: 1.3)
#' @param cex_axis       axis font size (default: 1.6)
#' @param cex_lab        axis labels font size (default: 1.7)
#' @param cex_main       title font size (default: 1.8)
#' @param cex_legend     legend font size (default: 1.3)
#' @param show_legend    logical, display boxed legend (default: TRUE)
#' @param save_files     logical, whether to save PDF and PNG (default: FALSE)
#' @param output_prefix  filename prefix for saved files
plot_wbs <- function(wbs_res,
                     dates          = NULL,
                     main           = NULL,
                     ylab           = expression(bold(paste("Observed ", italic(X[t])))),
                     xlab           = NULL,
                     flag_outliers  = TRUE,
                     outlier_thresh = 4,
                     level_col      = "#1565c0",
                     level_lwd      = 2.8,
                     break_col      = "#c62828",
                     break_lwd      = 2.0,
                     break_lty      = 2,
                     outlier_col    = "#c62828",
                     series_col     = "gray25",
                     series_lwd     = 1.3,
                     cex_axis       = 1.6,
                     cex_lab        = 1.7,
                     cex_main       = 1.8,
                     cex_legend     = 1.3,
                     cex_break      = 1.15,
                     date_format    = "%Y-%m-%d",
                     min_label_gap  = 0.15,
                     show_legend    = TRUE,
                     save_files     = FALSE,
                     output_prefix  = "wbs_segmentation") {

  x <- wbs_res$x
  N <- length(x)

  # Resolve dates
  dates_vec <- if (!is.null(dates)) {
    as.Date(dates)
  } else if (!is.null(wbs_res$dates)) {
    as.Date(wbs_res$dates)
  } else {
    NULL
  }

  has_dates <- !is.null(dates_vec) && length(dates_vec) == N
  x_axis <- if (has_dates) dates_vec else 1:N
  x_lims <- range(x_axis)
  y_lims <- extendrange(x, f = 0.08)

  if (is.null(xlab)) {
    xlab <- if (has_dates) expression(bold("Date")) else expression(bold(paste("Time ", italic(t))))
  }
  if (is.null(main)) {
    main <- "Robust Linearized WBS Changepoint Detection"
  }

  # Identify severe isolated outliers via two-sided Hampel filter
  out_idx <- integer(0)
  if (isTRUE(flag_outliers)) {
    k_hamp <- 5
    local_med <- sapply(1:N, function(i) median(x[max(1, i - k_hamp):min(N, i + k_hamp)]))
    local_mad <- sapply(1:N, function(i) median(abs(x[max(1, i - k_hamp):min(N, i + k_hamp)] - local_med[i])) / 0.6745)
    sig_glob  <- median(abs(diff(x))) / (sqrt(2) * 0.6745)
    local_scale <- pmax(local_mad, 0.5 * sig_glob)
    r_score   <- abs(x - local_med) / local_scale
    out_idx   <- which(r_score > outlier_thresh)
  }
  draw_canvas <- function() {
    top_mar <- if (length(wbs_res$cpts) > 0) 7.0 else 3.2
    par(mar = c(4.8, 5.4, top_mar, 1.2),
        mgp = c(3.3, 1.1, 0),
        tcl = -0.5)

    plot(NA, NA,
         xlim = x_lims,
         ylim = y_lims,
         xlab = xlab,
         ylab = ylab,
         cex.lab = cex_lab,
         font.lab = 2,
         col.lab = "black",
         xaxt = "n", yaxt = "n",
         bty = "n")

    # Interior light dashed gridlines matching tick levels
    grid(col = "gray88", lty = 2, lwd = 0.9)

    # Observed series path (dark charcoal / gray25, lwd = 1.3)
    lines(x_axis, x, col = series_col, lwd = series_lwd)

    # Piecewise segment levels
    bounds <- c(0, wbs_res$cpts, N)
    for (j in 1:(length(bounds) - 1)) {
      x_start <- bounds[j] + 1
      x_end   <- bounds[j + 1]
      lines(c(x_axis[x_start], x_axis[x_end]), rep(wbs_res$segment_levels[j], 2),
            col = level_col, lwd = level_lwd)
      if (j < (length(bounds) - 1)) {
        lines(rep(x_axis[x_end], 2), c(wbs_res$segment_levels[j], wbs_res$segment_levels[j + 1]),
              col = level_col, lty = 3, lwd = 1.5)
      }
    }

    # Detected changepoints: crimson dashed vertical lines & tiered horizontal date labels
    max_tier <- 1
    if (length(wbs_res$cpts) > 0) {
      break_x <- x_axis[wbs_res$cpts]
      abline(v = break_x, col = break_col, lty = break_lty, lwd = break_lwd)

      # Format horizontal break labels
      top_labels <- if (has_dates) format(break_x, date_format) else paste0("t=", wbs_res$cpts)

      # Measure string widths in user coordinates
      str_w <- sapply(top_labels, function(s) strwidth(s, units = "user", cex = cex_break))
      w_safe <- str_w * (1 + min_label_gap)

      # Assign tiers greedily to avoid horizontal overlap
      tiers <- integer(length(break_x))
      tier_rights <- numeric(0)
      for (i in seq_along(break_x)) {
        xi <- as.numeric(break_x[i])
        wi <- w_safe[i]
        left_i  <- xi - 0.5 * wi
        right_i <- xi + 0.5 * wi

        assigned <- FALSE
        if (length(tier_rights) > 0) {
          for (t in seq_along(tier_rights)) {
            if (left_i >= tier_rights[t]) {
              tiers[i] <- t
              tier_rights[t] <- right_i
              assigned <- TRUE
              break
            }
          }
        }
        if (!assigned) {
          new_t <- length(tier_rights) + 1
          tiers[i] <- new_t
          tier_rights <- c(tier_rights, right_i)
        }
      }
      max_tier <- max(tiers, 1)

      # Margin drawing coordinates
      usr <- par("usr")
      pin <- par("pin")
      csi <- par("csi")
      dy_line <- csi * (usr[4] - usr[3]) / pin[2]

      # Draw guide stems and bold horizontal labels
      for (i in seq_along(break_x)) {
        t_i <- tiers[i]
        stem_end_line <- 0.35 + (t_i - 1) * 1.80
        label_line    <- 0.55 + (t_i - 1) * 1.80

        y_stem_end <- usr[4] + stem_end_line * dy_line
        y_label    <- usr[4] + label_line * dy_line

        # Guide stem connecting plot border to the label
        segments(x0 = break_x[i], y0 = usr[4], x1 = break_x[i], y1 = y_stem_end,
                 col = break_col, lwd = 1.4, xpd = NA)

        # Bold horizontal label (larger cex)
        text(x = break_x[i], y = y_label, labels = top_labels[i],
             col = break_col, font = 2, cex = cex_break, adj = c(0.5, 0), xpd = NA)
      }
    }

    # Severe outliers: prominent red open circles matching DGP.R style
    if (length(out_idx) > 0) {
      points(x_axis[out_idx], x[out_idx],
             col = outlier_col, pch = 1, lwd = 2.0, cex = 1.3)
    }

    # Prominent axis ticks
    if (has_dates) {
      d_ticks <- pretty(dates_vec, n = 8)
      axis(1, at = d_ticks, labels = format(d_ticks, "%Y"),
           cex.axis = cex_axis, font.axis = 1, lwd = 0, lwd.ticks = 2.0)
    } else {
      axis(1, cex.axis = cex_axis, font.axis = 1, lwd = 0, lwd.ticks = 2.0)
    }
    axis(2, cex.axis = cex_axis, font.axis = 1, lwd = 0, lwd.ticks = 2.0, las = 1)

    # Sharp outer bounding frame
    box(which = "plot", lty = "solid", lwd = 2.0, col = "black")

    # Prominent bold title positioned above the highest label tier
    title_line <- if (length(wbs_res$cpts) > 0) 2.9 + (max_tier - 1) * 1.80 else 1.8
    title(main = main, col.main = "black", font.main = 2, cex.main = cex_main, line = title_line)

    # Boxed legend matching DGP.R
    if (show_legend) {
      seg_lbl <- if (!is.null(wbs_res$segment_level)) sprintf("Segment Level (%s)", wbs_res$segment_level) else "Estimated Level"
      leg_text <- c("Observed Series",
                    seg_lbl,
                    paste0("Detected Breaks (K = ", length(wbs_res$cpts), ")"))
      leg_cols <- c(series_col, level_col, break_col)
      leg_lty  <- c(1, 1, break_lty)
      leg_lwd  <- c(series_lwd, level_lwd, break_lwd)
      leg_pch  <- c(NA, NA, NA)
      leg_pt_lwd <- c(NA, NA, NA)
      leg_pt_cex <- c(NA, NA, NA)

      if (length(out_idx) > 0) {
        leg_text   <- c(leg_text, paste0("Severe Outliers (N = ", length(out_idx), ")"))
        leg_cols   <- c(leg_cols, outlier_col)
        leg_lty    <- c(leg_lty, NA)
        leg_lwd    <- c(leg_lwd, NA)
        leg_pch    <- c(leg_pch, 1)
        leg_pt_lwd <- c(leg_pt_lwd, 2.0)
        leg_pt_cex <- c(leg_pt_cex, 1.3)
      }

      legend("topleft",
             legend  = leg_text,
             col     = leg_cols,
             lty     = leg_lty,
             lwd     = leg_lwd,
             pch     = leg_pch,
             pt.lwd  = leg_pt_lwd,
             pt.cex  = leg_pt_cex,
             bty     = "o",
             box.col = "gray80",
             box.lwd = 1.0,
             bg      = "white",
             cex     = cex_legend,
             inset   = c(0.02, 0.03))
    }
  }

  draw_canvas()

  if (save_files) {
    png_file <- paste0(output_prefix, ".png")
    png(png_file, width = 16, height = 7.5, units = "in", res = 300)
    draw_canvas()
    dev.off()
    cat(sprintf("[SUCCESS] Saved PNG: %s\n", png_file))

    pdf_file <- paste0(output_prefix, ".pdf")
    pdf(pdf_file, width = 16, height = 7.5)
    draw_canvas()
    dev.off()
    cat(sprintf("[SUCCESS] Saved PDF: %s\n", pdf_file))
  }

  invisible(list(outliers = out_idx, break_locs = if (has_dates) dates_vec[wbs_res$cpts] else wbs_res$cpts))
}

#' S3 plot method for class 'wbs'
plot.wbs <- function(x, ...) {
  plot_wbs(x, ...)
}

############################
#Example
############################

# # Load implementations
# source("CUSUM.R")
# source("RobWBS.R")

# set.seed(42)
# n <- 1200

# # 1. Construct a piecewise-constant signal with 3 structural breaks
# # Changes occur at t = 300, 700, and 950
# signal <- rep(c(0, 3.0, -1.5, 1.5), c(300, 400, 250, 250))

# # 2. Add ARMA(2,1) stationary dependent noise
# noise <- arima.sim(model = list(ar = c(0.2, -0.1), ma = 0.2), n = n, sd = 0.8)

# # 3. Additive Outliers (AO): 5% contamination rate of magnitude K = 10
# contam_idx <- sample(1:n, size = floor(0.05 * n))
# outlier_shocks <- sample(c(-10, 10), size = length(contam_idx), replace = TRUE)

# x <- signal + noise
# x[contam_idx] <- x[contam_idx] + outlier_shocks

# # 4. Run Linearized Robust WBS with default sSIC selection and Robust M-estimator levels
# res <- WBS.mean(
#   x             = x,
#   M             = 5000,
#   Kmax          = 10,
#   selection     = "sSIC",
#   loss          = "Welsh",
#   c             = 2.985,
#   linearized    = TRUE,
#   segment_level = "robust",
#   plot          = TRUE,
#   main          = "Robust Linearized WBS (Energy Outlier Contamination)"
# )

# # Inspect detected breaks and sSIC progression
# print(res$ssic_table)
# cat("True breaks at: 300, 700, 950\n")
# cat("Detected breaks at:", res$cpts, "\n")
# cat("Robust Segment Levels:", round(res$segment_levels, 3), "\n")