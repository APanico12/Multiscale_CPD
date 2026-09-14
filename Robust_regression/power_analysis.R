# ==============================================================================
# File: power_analysis.R
# Description: Monte Carlo Empirical Power Analysis for Robust CUSUM Regression.
#              Evaluates rejection rates under H1 across varying break magnitudes (delta),
#              sample sizes (N), contrast matrices (C_mat), and contamination levels.
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(foreach)
  library(ggplot2)
})
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
}

# Source dependencies
if (file.exists("DGP.R")) source("DGP.R")
if (file.exists("CUSUM.R")) source("CUSUM.R")

#' Run Monte Carlo Power Analysis across a grid of jump sizes (delta)
#'
#' @param delta_grid     numeric vector of break magnitudes (e.g. seq(0.0, 1.5, by = 0.25))
#' @param n_values       numeric vector of sample sizes (e.g. c(250, 500, 1000))
#' @param MC_reps        number of Monte Carlo replications per grid point (e.g. 100 - 500)
#' @param B_boot         number of bootstrap iterations for critical values (default 100)
#' @param C_mat          contrast matrix (default: diag(d) for joint, or c(1, 0) for beta_1)
#' @param model_type     "LMHC" or "LMAT"
#' @param error_dist     "normal" or "t3"
#' @param contamination  "Clean", "AO", "IO", or "RO"
#' @param epsilon        outlier contamination proportion (default 0.05)
#' @param loss           "Welsh" or "Tukey"
#' @param n_cores        number of CPU cores for parallel execution
#' @param seed           master random seed for reproducibility
#' @return data.frame with columns: delta, n, power, se_power, ci_lower, ci_upper, mean_break_u, rmse_break_u
run_power_analysis <- function(delta_grid    = seq(0.0, 1.5, by = 0.25),
                               n_values      = c(300, 600),
                               MC_reps       = 1000,
                               B_boot        = 100,
                               C_mat         = NULL,
                               model_type    = "LMHC",
                               error_dist    = "normal",
                               contamination = "Clean",
                               epsilon       = 0.05,
                               loss          = "Welsh",
                               n_cores       = max(1, parallel::detectCores() - 2),
                               seed          = 123) {
  
  set.seed(seed)
  cat("======================================================================\n")
  cat(" Starting Monte Carlo Power Analysis for Robust CUSUM Regression\n")
  cat("======================================================================\n")
  cat(sprintf("Delta Grid    : %s\n", paste(delta_grid, collapse = ", ")))
  cat(sprintf("Sample Sizes  : %s\n", paste(n_values, collapse = ", ")))
  cat(sprintf("MC Reps       : %d per design point\n", MC_reps))
  cat(sprintf("Bootstrap B   : %d\n", B_boot))
  cat(sprintf("Model Type    : %s | Error: %s | Contamination: %s (eps = %.2f)\n", 
              model_type, error_dist, contamination, epsilon))
  cat(sprintf("Parallel Cores: %d\n", n_cores))
  cat("----------------------------------------------------------------------\n")
  
  # Setup parallel cluster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit({
    stopCluster(cl)
  })
  
  # Export required environment and functions to workers
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(zoo)
      library(robustbase)
    })
    source("DGP.R")
    source("CUSUM.R")
  })
  
  results_list <- list()
  
  for (n_cur in n_values) {
    cat(sprintf("\n--> Running for Sample Size N = %d ...\n", n_cur))
    
    for (delta_cur in delta_grid) {
      t_start <- Sys.time()
      
      # Determine whether this is H0 (delta = 0) or H1 (delta > 0)
      hp_scen <- if (abs(delta_cur) < 1e-6) "H0" else "H1"
      
      # Run MC_reps in parallel
      mc_res <- foreach(m = 1:MC_reps, .combine = rbind) %dopar% {
        # Unique seed per worker & iteration
        rep_seed <- seed + m * 1000 + round(delta_cur * 100) + n_cur
        
        # Generate DGP
        dgp <- generate_tv_regression_dgp(
          n             = n_cur,
          model_type    = model_type,
          error_dist    = error_dist,
          contamination = contamination,
          epsilon       = epsilon,
          hp_scenario   = hp_scen,
          delta         = delta_cur,
          seed          = rep_seed
        )
        
        # Run CUSUM test
        test_out <- tryCatch({
          CUSUM.regression(
            Y          = dgp$Y,
            X          = dgp$X,
            C_mat      = C_mat,
            B          = B_boot,
            loss       = loss,
            linearized = TRUE,
            plotting   = FALSE
          )
        }, error = function(e) NULL)
        
        if (is.null(test_out)) {
          return(c(reject = NA, break_u = NA))
        }
        
        c(reject = as.numeric(test_out$reject_95),
          break_u = as.numeric(test_out$break_u))
      }
      
      # Aggregate MC statistics
      valid_reps <- na.omit(mc_res)
      n_valid    <- nrow(valid_reps)
      
      if (n_valid > 0) {
        rej_vec   <- valid_reps[, 1]
        loc_vec   <- valid_reps[, 2]
        
        pow_hat   <- mean(rej_vec)
        se_pow    <- sqrt(pow_hat * (1 - pow_hat) / n_valid)
        ci_low    <- max(0, pow_hat - 1.96 * se_pow)
        ci_high   <- min(1, pow_hat + 1.96 * se_pow)
        
        # Localization accuracy (true change point is tau* = 0.5)
        mean_loc  <- mean(loc_vec)
        rmse_loc  <- sqrt(mean((loc_vec - 0.5)^2))
      } else {
        pow_hat <- se_pow <- ci_low <- ci_high <- mean_loc <- rmse_loc <- NA
      }
      
      elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
      cat(sprintf("  delta = %4.2f | Power = %5.1f%% (SE = %4.1f%%) | Loc u* = %.3f | Elapsed: %4.1fs\n",
                  delta_cur, pow_hat * 100, se_pow * 100, mean_loc, elapsed))
      
      results_list[[length(results_list) + 1]] <- data.frame(
        N             = n_cur,
        delta         = delta_cur,
        power         = pow_hat,
        se_power      = se_pow,
        ci_lower      = ci_low,
        ci_upper      = ci_high,
        mean_break_u  = mean_loc,
        rmse_break_u  = rmse_loc,
        contamination = contamination,
        model_type    = model_type,
        stringsAsFactors = FALSE
      )
    }
  }
  
  df_res <- do.call(rbind, results_list)
  df_res$N <- factor(paste0("N = ", df_res$N), levels = paste0("N = ", n_values))
  
  return(df_res)
}

#' Plot Power Curves with Confidence Ribbons
#'
#' @param df_res data.frame returned by run_power_analysis
#' @param title  plot title
#' @return ggplot object

plot_power_curves <- function(df_res, title = "Empirical Power Curves under H1 (Abrupt Break in beta_1)") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("Notice: 'ggplot2' is not installed. Plotting using base R graphics...\n")
    deltas <- sort(unique(df_res$delta))
    n_vals <- unique(df_res$N)
    cols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
    
    plot(deltas, seq(0, 1, length.out = length(deltas)), type = "n", ylim = c(0, 1.05),
         xlab = expression(paste("Break Magnitude (", delta, ")")),
         ylab = "Empirical Power P(Reject H0)",
         main = title, font.main = 2)
    abline(h = 0.05, lty = 2, col = "gray40")
    grid(col = "gray85")
    
    for (i in seq_along(n_vals)) {
      sub_n <- df_res[df_res$N == n_vals[i], ]
      sub_n <- sub_n[order(sub_n$delta), ]
      col_cur <- cols[(i - 1) %% length(cols) + 1]
      lines(sub_n$delta, sub_n$power, col = col_cur, lwd = 2)
      points(sub_n$delta, sub_n$power, col = col_cur, pch = 16, cex = 1.1)
    }
    legend("bottomright", legend = as.character(n_vals), col = cols[seq_along(n_vals)],
           lwd = 2, pch = 16, bty = "n")
    return(invisible(NULL))
  }

  p <- ggplot(df_res, aes(x = delta, y = power, color = N, group = N, fill = N)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkgray", linewidth = 0.8) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.5) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = unique(df_res$delta)) +
    labs(
      title    = title,
      subtitle = expression(paste("True break at ", tau^"*", " = 0.50 | Nominal Size ", alpha, " = 0.05 (dashed line)")),
      x        = expression(paste("Break Magnitude (", delta, ")")),
      y        = "Empirical Power P(Reject H0)",
      color    = "Sample Size",
      fill     = "Sample Size"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle   = element_text(hjust = 0.5, color = "dimgray", size = 11),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  return(p)
}

