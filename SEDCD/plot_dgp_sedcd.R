# ==============================================================================
# File: plot_dgp_sedcd.R
# Description: Generates publication-ready diagnostic figures for the
#              Smoothed Extremal Downside Correlation Density (SEDCD) DGP
#              matching exact aesthetic styling (bold annotations,
#              prominent tick marks, interior dashed grids, sharp borders).
# ==============================================================================

suppressPackageStartupMessages({
  library(grDevices)
  library(graphics)
})

# Source DGP if not already loaded
if (!exists("simulate_sedcd_dgp")) {
  cands <- c("dgp_sedcd.R", file.path("SEDCD", "dgp_sedcd.R"), file.path("..", "SEDCD", "dgp_sedcd.R"))
  found <- FALSE
  for (cand in cands) {
    if (file.exists(cand)) {
      source(cand)
      found <- TRUE
      break
    }
  }
  if (!found) stop("Could not locate dgp_sedcd.R!")
}

#' Plot Publication Figure for SEDCD DGP
#'
#' @param n Sample size (default 1000)
#' @param seed Random seed for reproducibility (default 42)
#' @param output_prefix Output filename prefix (default "sedcd_dgp")
#' @param save_files Logical, whether to save PDF and PNG (default TRUE)
#' @param cex_axis Axis tick font size (default 1.8)
#' @param cex_lab Axis label font size (default 1.9)
#' @param cex_main Main title font size (default 2.1)
plot_sedcd_dgp <- function(n = 1000,
                           seed = 42,
                           output_prefix = "sedcd_dgp",
                           save_files = TRUE,
                           cex_axis = 1.6,
                           cex_lab = 1.75,
                           cex_main = 1.9) {
  
  cat("Generating SEDCD DGP simulation data for visualization...\n")
  set.seed(seed)
  
  # 1. Generate Clayton scenarios: H0, H1, H2
  dgp_c_h0 <- simulate_sedcd_dgp(n = n, d = 2, scenario = "H0", copula_type = "clayton", seed = seed)
  dgp_c_h1 <- simulate_sedcd_dgp(n = n, d = 2, scenario = "H1", copula_type = "clayton", seed = seed + 1)
  dgp_c_h2 <- simulate_sedcd_dgp(n = n, d = 2, scenario = "H2", copula_type = "clayton", seed = seed + 2)
  
  # 2. Generate Student-t scenario H1
  dgp_t_h1 <- simulate_sedcd_dgp(n = n, d = 2, scenario = "H1", copula_type = "t", seed = seed + 3)
  
  u_grid <- (1:n) / n
  
  # Colors
  col_blue   <- "#1565c0"
  col_green  <- "#2e7d32"
  col_red    <- "#c62828"
  col_purple <- "#6a1b9a"
  col_gray   <- "gray30"
  col_grid   <- "gray88"
  
  draw_panels <- function() {
    # 2x2 layout with wide panels
    par(mfrow = c(2, 2),
        mar   = c(4.8, 5.8, 3.2, 1.4),
        mgp   = c(3.6, 1.1, 0),
        tcl   = -0.5)
    
    # --------------------------------------------------------------------------
    # Panel (a): Bivariate Returns & GJR-GARCH Volatility Bands
    # --------------------------------------------------------------------------
    y_lim_a <- extendrange(dgp_c_h1$X[, 1:2], f = 0.10)
    
    plot(NA, NA,
         xlim = c(0, 1.0),
         ylim = y_lim_a,
         xlab = expression(bold(paste("Rescaled Time ", italic(u == t/n)))),
         ylab = expression(bold(paste("Observed Returns ", italic(X[t])))),
         cex.lab = cex_lab, font.lab = 2, col.lab = "black",
         xaxt = "n", yaxt = "n", bty = "n")
    
    grid(col = col_grid, lty = 2, lwd = 0.9)
    
    # Volatility envelope for series 1 (+- 2 sigma)
    sig1 <- dgp_c_h1$sigma_matrix[, 1]
    mu1  <- dgp_c_h1$mu_matrix[, 1]
    polygon(c(u_grid, rev(u_grid)),
            c(mu1 + 2.0 * sig1, rev(mu1 - 2.0 * sig1)),
            col = rgb(0.13, 0.40, 0.75, 0.15), border = NA)
    
    # Lines for X1 and X2
    lines(u_grid, dgp_c_h1$X[, 1], col = col_blue, lwd = 1.2)
    lines(u_grid, dgp_c_h1$X[, 2], col = col_green, lwd = 1.0)
    
    # Boundary and axis ticks
    axis(1, at = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"),
         cex.axis = cex_axis, lwd = 0, lwd.ticks = 2.0)
    axis(2, cex.axis = cex_axis, lwd = 0, lwd.ticks = 2.0, las = 1)
    box(which = "plot", lty = "solid", lwd = 2.0, col = "black")
    
    title(main = "(a) Bivariate Returns & GJR-GARCH Volatility",
          font.main = 2, cex.main = cex_main, line = 1.1)
    
    legend("topright",
           legend = c(expression(italic(X[1]) ~ "(Series 1)"),
                      expression(italic(X[2]) ~ "(Series 2)"),
                      expression(paste("\u00B1", 2, italic(sigma[1])) ~ "Band")),
           col = c(col_blue, col_green, rgb(0.13, 0.40, 0.75, 0.4)),
           lty = c(1, 1, 1),
           lwd = c(1.8, 1.8, 6.0),
           bty = "o", box.col = "gray80", bg = "white",
           cex = 1.25, inset = c(0.02, 0.03))
    
    # --------------------------------------------------------------------------
    # Panel (b): Joint Copula Realizations & Downside Tail Clustering
    # --------------------------------------------------------------------------
    u1 <- dgp_c_h1$U[, 1]
    u2 <- dgp_c_h1$U[, 2]
    alpha_tail <- 0.10
    
    plot(NA, NA,
         xlim = c(0, 1.0),
         ylim = c(0, 1.0),
         xlab = expression(bold(paste("Copula Margin ", italic(U[1])))),
         ylab = expression(bold(paste("Copula Margin ", italic(U[2])))),
         cex.lab = cex_lab, font.lab = 2, col.lab = "black",
         xaxt = "n", yaxt = "n", bty = "n")
    
    grid(col = col_grid, lty = 2, lwd = 0.9)
    
    # Shaded lower tail downside quadrant: [0, alpha_tail] x [0, alpha_tail]
    rect(0, 0, alpha_tail, alpha_tail, col = rgb(0.85, 0.20, 0.20, 0.20), border = col_red, lty = 2, lwd = 1.5)
    
    # Scatter points: regular vs joint downside
    is_downside <- (u1 <= alpha_tail) & (u2 <= alpha_tail)
    points(u1[!is_downside], u2[!is_downside], col = "gray45", pch = 16, cex = 0.75)
    points(u1[is_downside], u2[is_downside], col = col_red, pch = 16, cex = 1.4)
    
    # Threshold dashed lines
    abline(v = alpha_tail, col = col_red, lty = 3, lwd = 1.5)
    abline(h = alpha_tail, col = col_red, lty = 3, lwd = 1.5)
    
    axis(1, at = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"),
         cex.axis = cex_axis, lwd = 0, lwd.ticks = 2.0)
    axis(2, at = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"),
         cex.axis = cex_axis, lwd = 0, lwd.ticks = 2.0, las = 1)
    box(which = "plot", lty = "solid", lwd = 2.0, col = "black")
    
    title(main = "(b) Joint Clayton Copula Realizations",
          font.main = 2, cex.main = cex_main, line = 1.1)
    
    legend("topright",
           legend = c("Copula Draws", expression(paste("Joint Downside (", italic(U[j] <= 0.10), ")"))),
           col = c("gray45", col_red),
           pch = c(16, 16),
           pt.cex = c(1.0, 1.5),
           bty = "o", box.col = "gray80", bg = "white",
           cex = 1.25, inset = c(0.02, 0.03))
    
    # --------------------------------------------------------------------------
    # Panel (c): Tail Dependence Trajectories across Hypotheses (lambda_L(t))
    # --------------------------------------------------------------------------
    y_lim_c <- c(0.40, 0.95)
    
    plot(NA, NA,
         xlim = c(0, 1.0),
         ylim = y_lim_c,
         xlab = expression(bold(paste("Rescaled Time ", italic(u == t/n)))),
         ylab = expression(bold(paste("Lower Tail Dependence ", italic(lambda[L](t))))),
         cex.lab = cex_lab, font.lab = 2, col.lab = "black",
         xaxt = "n", yaxt = "n", bty = "n")
    
    grid(col = col_grid, lty = 2, lwd = 0.9)
    
    # H0 path (stable)
    lines(u_grid, dgp_c_h0$lambda_L_path, col = col_gray, lwd = 1.6, lty = 1)
    
    # H2 path (gradual drift)
    lines(u_grid, dgp_c_h2$lambda_L_path, col = col_green, lwd = 2.0, lty = 1)
    
    # H1 path (abrupt break at u = 0.5)
    lines(u_grid, dgp_c_h1$lambda_L_path, col = col_red, lwd = 2.2, lty = 1)
    
    # Vertical line indicating break point u* = 0.5
    abline(v = 0.50, col = col_red, lty = 2, lwd = 1.8)
    
    axis(1, at = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"),
         cex.axis = cex_axis, lwd = 0, lwd.ticks = 2.0)
    axis(2, at = seq(0.4, 0.9, by = 0.1), cex.axis = cex_axis, lwd = 0, lwd.ticks = 2.0, las = 1)
    box(which = "plot", lty = "solid", lwd = 2.0, col = "black")
    
    title(main = expression(bold(paste("(c) Tail Dependence Trajectories ", italic(lambda[L](t))))),
          font.main = 2, cex.main = cex_main, line = 1.1)
    
    legend("topleft",
           legend = c(expression(italic(H[0]) ~ "(Stable Baseline)"),
                      expression(italic(H[1]) ~ "(Abrupt Shift at " * italic(u) * " = 0.5)"),
                      expression(italic(H[2]) ~ "(Gradual Drift)")),
           col = c(col_gray, col_red, col_green),
           lty = c(1, 1, 1),
           lwd = c(2.0, 2.4, 2.2),
           bty = "o", box.col = "gray80", bg = "white",
           cex = 1.25, inset = c(0.02, 0.03))
    
    # --------------------------------------------------------------------------
    # Panel (d): Latent GAS Score Dynamics psi_t & Baseline Intercept omega_t
    # --------------------------------------------------------------------------
    y_lim_d <- extendrange(c(dgp_c_h1$psi_path, dgp_c_h1$omega_path * 10), f = 0.15)
    
    plot(NA, NA,
         xlim = c(0, 1.0),
         ylim = y_lim_d,
         xlab = expression(bold(paste("Rescaled Time ", italic(u == t/n)))),
         ylab = expression(bold(paste("Latent Parameter ", italic(psi[t])))),
         cex.lab = cex_lab, font.lab = 2, col.lab = "black",
         xaxt = "n", yaxt = "n", bty = "n")
    
    grid(col = col_grid, lty = 2, lwd = 0.9)
    
    # Score innovation spikes in background
    # Scale score to fit pleasantly
    score_scaled <- dgp_c_h1$score_path * 0.1 + mean(dgp_c_h1$psi_path)
    segments(u_grid, mean(dgp_c_h1$psi_path), u_grid, score_scaled,
             col = rgb(0.6, 0.6, 0.6, 0.35), lwd = 0.8)
    
    # Latent state psi_t
    lines(u_grid, dgp_c_h1$psi_path, col = col_blue, lwd = 1.8)
    
    # Intercept omega_t scaled to unconstrained stationary level: omega_t / (1 - beta)
    omega_stationary <- dgp_c_h1$omega_path / (1.0 - 0.90)
    lines(u_grid, omega_stationary, col = col_red, lwd = 2.2, lty = 2)
    
    abline(v = 0.50, col = col_red, lty = 2, lwd = 1.5)
    
    axis(1, at = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"),
         cex.axis = cex_axis, lwd = 0, lwd.ticks = 2.0)
    axis(2, cex.axis = cex_axis, lwd = 0, lwd.ticks = 2.0, las = 1)
    box(which = "plot", lty = "solid", lwd = 2.0, col = "black")
    
    title(main = expression(bold(paste("(d) Latent GAS Parameter Dynamics ", italic(psi[t]), " & Intercept ", italic(omega[t])))),
          font.main = 2, cex.main = cex_main, line = 1.1)
    
    legend("topleft",
           legend = c(expression(italic(psi[t]) ~ "(Latent Path)"),
                      expression(italic(bar(psi)[t]) == italic(omega[t]) / (1 - beta) ~ "(Mean Level)"),
                      expression(italic(s[t]) ~ "(Scaled Score)")),
           col = c(col_blue, col_red, "gray60"),
           lty = c(1, 2, 1),
           lwd = c(2.0, 2.2, 1.2),
           bty = "o", box.col = "gray80", bg = "white",
           cex = 1.25, inset = c(0.02, 0.03))
  }
  
  # Save vector PDF & high-res PNG
  if (save_files) {
    pdf_file <- paste0(output_prefix, ".pdf")
    png_file <- paste0(output_prefix, ".png")
    
    # Sized for 1/3-1/2 of LaTeX page across full linewidth
    pdf(pdf_file, width = 16.0, height = 10.5)
    draw_panels()
    dev.off()
    cat(sprintf("Saved vector PDF: %s\n", pdf_file))
    
    png(png_file, width = 4800, height = 3150, res = 300)
    draw_panels()
    dev.off()
    cat(sprintf("Saved high-res PNG (300 DPI): %s\n", png_file))
    
    # Sync to paper/img/ directory
    img_dirs <- c(
      file.path("..", "paper", "img"),
      file.path("Multiscale_CPD", "paper", "img"),
      file.path("paper", "img")
    )
    for (idir in img_dirs) {
      if (dir.exists(idir)) {
        target_png <- file.path(idir, paste0(output_prefix, ".png"))
        target_pdf <- file.path(idir, paste0(output_prefix, ".pdf"))
        file.copy(png_file, target_png, overwrite = TRUE)
        file.copy(pdf_file, target_pdf, overwrite = TRUE)
        cat(sprintf("Synced figure to paper directory: %s\n", target_png))
        break
      }
    }
  }
  
  return(invisible(list(
    dgp_c_h0 = dgp_c_h0,
    dgp_c_h1 = dgp_c_h1,
    dgp_c_h2 = dgp_c_h2,
    dgp_t_h1 = dgp_t_h1
  )))
}

# Execute script directly if run via Rscript
if (!interactive()) {
  plot_sedcd_dgp(n = 1000, seed = 42, output_prefix = "sedcd_dgp", save_files = TRUE)
}
