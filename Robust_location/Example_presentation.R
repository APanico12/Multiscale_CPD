# ==============================================================================
# File: Example_presentation.R
# Description: Generates and saves the time-varying location series (X_t)
#              overlaid with the true mean trajectory (mu(u)) under H1.
#              Saves publication figures in both Robust_location and paper/img.
# ==============================================================================

# Dynamic path resolution to source DGP.R
dgp_cands <- c("DGP.R", file.path("Robust_location", "DGP.R"), file.path("..", "Robust_location", "DGP.R"))
found_dgp <- FALSE
for (cand in dgp_cands) {
  if (file.exists(cand)) {
    source(cand)
    found_dgp <- TRUE
    break
  }
}
if (!found_dgp) stop("Could not locate DGP.R")

# 1. Simulate data: clean series with abrupt location break (H1)
set.seed(42)
X <- ARMA_mu(n = 500, var_scenario = "none", mu_scenario = "H1", delta = 0.5, seed = 42)

# 2. Plotting function
draw_location_plot <- function() {
  y_lims <- extendrange(c(X$Xt, X$mu_u), f = 0.10)
  
  par(mar = c(4.8, 5.2, 3.2, 1.2), mgp = c(3.2, 1.1, 0), tcl = -0.4)
  plot(NA, NA,
       xlim = c(0, 1),
       ylim = y_lims,
       xlab = expression(bold(paste("Rescaled Time ", italic(u == t/n)))),
       ylab = expression(bold(paste("Observed ", italic(X[t])))),
       cex.lab = 1.35,
       font.lab = 2,
       las = 1,
       bty = "n")
  
  # Light interior grid
  grid(col = "gray88", lty = 2, lwd = 0.9)
  
  # Draw observed series path
  lines(X$u, X$Xt, col = "gray30", lwd = 1.3)
  
  # Overlay true mean trajectory mu(u)
  lines(X$u, X$mu_u, col = "#c62828", lwd = 2.6)
  
  # Break location dashed line at u* = 0.5
  abline(v = 0.50, col = "#c62828", lty = 2, lwd = 1.8)
  
  # Outer boundary box
  box(which = "plot", lty = "solid", lwd = 2.0, col = "black")
  
  # Title
  title(main = expression(bold(paste("Observed Series ", italic(X[t]), " & True Mean ", italic(mu(u))))),
        cex.main = 1.45, font.main = 2, line = 1.1)
  
  # Legend
  legend("topleft",
         legend = c(expression(italic(X[t]) ~ "(Observed Series)"),
                    expression(italic(mu(u)) ~ "(True Mean Shift)"),
                    expression("Break Point (" * italic(u) * " = 0.5)")),
         col = c("gray30", "#c62828", "#c62828"),
         lty = c(1, 1, 2),
         lwd = c(1.3, 2.6, 1.8),
         bg = "white",
         box.col = "gray80",
         box.lwd = 1.0,
         cex = 1.15,
         inset = c(0.02, 0.03))
}

# 3. Output filenames
base_out <- "example_presentation"
pdf_file <- paste0(base_out, ".pdf")
png_file <- paste0(base_out, ".png")

# Generate vector PDF
pdf(pdf_file, width = 9.0, height = 5.5)
draw_location_plot()
dev.off()

# Generate high-resolution PNG (300 DPI)
png(png_file, width = 2700, height = 1650, res = 300)
draw_location_plot()
dev.off()

cat(sprintf("Saved %s and %s in current directory.\n", pdf_file, png_file))

# 4. Sync to paper/img/ directory
target_img_dirs <- c(
  file.path("..", "paper", "img"),
  file.path("paper", "img"),
  file.path("Multiscale_CPD", "paper", "img"),
  file.path("..", "..", "paper", "img")
)

synced <- FALSE
for (pdir in target_img_dirs) {
  if (dir.exists(pdir)) {
    file.copy(pdf_file, file.path(pdir, pdf_file), overwrite = TRUE)
    file.copy(png_file, file.path(pdir, png_file), overwrite = TRUE)
    cat(sprintf("Successfully synced to paper directory: %s/%s\n", pdir, pdf_file))
    cat(sprintf("Successfully synced to paper directory: %s/%s\n", pdir, png_file))
    synced <- TRUE
    break
  }
}

if (!synced) {
  warning("Could not find paper/img directory to sync files.")
}