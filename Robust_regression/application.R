source("CUSUM.R")

suppressPackageStartupMessages({
  library(ggplot2)
})
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
if (has_ggplot) {
  suppressPackageStartupMessages(library(ggplot2))
}

# Load prepared dataset
data <- read.csv("wind_scada_10min.csv")
data$datetime <- as.POSIXct(data$datetime)

# Extract arrays
Y <- data$Y_MW
X <- as.matrix(data[, c("X_intercept", "X_v1", "X_v2", "X_v3")])

# ------------------------------------------------------------------------------
# Plot Response Y and Regressors (No Intercept)
# ------------------------------------------------------------------------------

# 1. Faceted Time Series Plot (Response + Variables over time)
plot_df <- rbind(
  data.frame(datetime = data$datetime, Value = data$Y_MW, Variable = "Response: Y (Active Power, MW)"),
  data.frame(datetime = data$datetime, Value = data$X_v1, Variable = "Variable: X_v1 (Wind Speed v / 10)"),
  data.frame(datetime = data$datetime, Value = data$X_v2, Variable = "Variable: X_v2 (v / 10)^2"),
  data.frame(datetime = data$datetime, Value = data$X_v3, Variable = "Variable: X_v3 (v / 10)^3")
)
if (has_ggplot) {
  # 1. Faceted Time Series Plot (Response + Variables over time)
  plot_df <- rbind(
    data.frame(datetime = data$datetime, Value = data$Y_MW, Variable = "Response: Y (Active Power, MW)"),
    data.frame(datetime = data$datetime, Value = data$X_v1, Variable = "Variable: X_v1 (Wind Speed v / 10)"),
    data.frame(datetime = data$datetime, Value = data$X_v2, Variable = "Variable: X_v2 (v / 10)^2"),
    data.frame(datetime = data$datetime, Value = data$X_v3, Variable = "Variable: X_v3 (v / 10)^3")
  )

plot_df$Variable <- factor(plot_df$Variable, levels = c(
  "Response: Y (Active Power, MW)",
  "Variable: X_v1 (Wind Speed v / 10)",
  "Variable: X_v2 (v / 10)^2",
  "Variable: X_v3 (v / 10)^3"
))
  plot_df$Variable <- factor(plot_df$Variable, levels = c(
    "Response: Y (Active Power, MW)",
    "Variable: X_v1 (Wind Speed v / 10)",
    "Variable: X_v2 (v / 10)^2",
    "Variable: X_v3 (v / 10)^3"
  ))

p_timeseries <- ggplot(plot_df, aes(x = datetime, y = Value, color = Variable)) +
  geom_line(linewidth = 0.35, alpha = 0.85) +
  facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728")) +
  labs(
    title = "Wind Turbine SCADA: Response Y and Regressors (No Intercept)",
    subtitle = "10-Minute Resolution Time Series across Active Operating Range (3 - 18 m/s)",
    x = "Date / Time",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "dimgray", size = 10),
    panel.grid.minor = element_blank()
  )
  p_timeseries <- ggplot(plot_df, aes(x = datetime, y = Value, color = Variable)) +
    geom_line(linewidth = 0.35, alpha = 0.85) +
    facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728")) +
    labs(
      title = "Wind Turbine SCADA: Response Y and Regressors (No Intercept)",
      subtitle = "10-Minute Resolution Time Series across Active Operating Range (3 - 18 m/s)",
      x = "Date / Time",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray95", color = "gray80"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "dimgray", size = 10),
      panel.grid.minor = element_blank()
    )

print(p_timeseries)
ggsave("wind_scada_variables_time_series.png", plot = p_timeseries, width = 11, height = 8, dpi = 300)
ggsave("wind_scada_variables_time_series.pdf", plot = p_timeseries, width = 11, height = 8)
cat("Saved: wind_scada_variables_time_series.png and .pdf\n")
  print(p_timeseries)
  ggsave("wind_scada_variables_time_series.png", plot = p_timeseries, width = 11, height = 8, dpi = 300)
  ggsave("wind_scada_variables_time_series.pdf", plot = p_timeseries, width = 11, height = 8)
  cat("Saved: wind_scada_variables_time_series.png and .pdf\n")

# 2. Empirical Power Curve: Response Y vs. Wind Speed Regressor (with curtailment outliers)
p_power_curve <- ggplot(data, aes(x = wind_speed, y = Y_MW)) +
  geom_point(aes(color = as.logical(is_curtailment_outlier)), alpha = 0.35, size = 0.8) +
  scale_color_manual(
    name = "Observation",
    values = c("FALSE" = "#1f77b4", "TRUE" = "#d62728"),
    labels = c("FALSE" = "Normal Operation", "TRUE" = "Curtailment Outlier")
  ) +
  labs(
    title = "Empirical Wind Turbine Power Curve: Y (MW) vs Wind Speed (m/s)",
    subtitle = "Nonlinear relationship modeled via cubic regressors [X_v1, X_v2, X_v3]",
    x = "Wind Speed (m/s)",
    y = "Active Power Output Y (MW)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "dimgray", size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
  # 2. Empirical Power Curve: Response Y vs. Wind Speed Regressor (with curtailment outliers)
  p_power_curve <- ggplot(data, aes(x = wind_speed, y = Y_MW)) +
    geom_point(aes(color = as.logical(is_curtailment_outlier)), alpha = 0.35, size = 0.8) +
    scale_color_manual(
      name = "Observation",
      values = c("FALSE" = "#1f77b4", "TRUE" = "#d62728"),
      labels = c("FALSE" = "Normal Operation", "TRUE" = "Curtailment Outlier")
    ) +
    labs(
      title = "Empirical Wind Turbine Power Curve: Y (MW) vs Wind Speed (m/s)",
      subtitle = "Nonlinear relationship modeled via cubic regressors [X_v1, X_v2, X_v3]",
      x = "Wind Speed (m/s)",
      y = "Active Power Output Y (MW)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "dimgray", size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

print(p_power_curve)
ggsave("wind_scada_power_curve.png", plot = p_power_curve, width = 8.5, height = 6, dpi = 300)
ggsave("wind_scada_power_curve.pdf", plot = p_power_curve, width = 8.5, height = 6)
cat("Saved: wind_scada_power_curve.png and .pdf\n")
  print(p_power_curve)
  ggsave("wind_scada_power_curve.png", plot = p_power_curve, width = 8.5, height = 6, dpi = 300)
  ggsave("wind_scada_power_curve.pdf", plot = p_power_curve, width = 8.5, height = 6)
  cat("Saved: wind_scada_power_curve.png and .pdf\n")

} else {
  cat("Notice: 'ggplot2' is not installed. Generating figures using base R graphics...\n")
  
  # Function to render base R time series
  draw_timeseries_base <- function() {
    par(mfrow = c(4, 1), mar = c(2.5, 4.5, 2, 1), oma = c(2, 0, 2.5, 0))
    plot(data$datetime, data$Y_MW, type = "l", col = "#1f77b4", lwd = 0.8,
         xlab = "", ylab = "Power (MW)", main = "Response: Y (Active Power, MW)", font.main = 2, cex.main = 1.0)
    grid(col = "gray85")
    plot(data$datetime, data$X_v1, type = "l", col = "#2ca02c", lwd = 0.8,
         xlab = "", ylab = "v / 10", main = "Variable: X_v1 (Wind Speed v / 10)", font.main = 2, cex.main = 1.0)
    grid(col = "gray85")
    plot(data$datetime, data$X_v2, type = "l", col = "#ff7f0e", lwd = 0.8,
         xlab = "", ylab = "(v / 10)^2", main = "Variable: X_v2 (v / 10)^2", font.main = 2, cex.main = 1.0)
    grid(col = "gray85")
    plot(data$datetime, data$X_v3, type = "l", col = "#d62728", lwd = 0.8,
         xlab = "Date / Time", ylab = "(v / 10)^3", main = "Variable: X_v3 (v / 10)^3", font.main = 2, cex.main = 1.0)
    grid(col = "gray85")
    mtext("Wind Turbine SCADA: Response Y and Regressors (No Intercept)", outer = TRUE, cex = 1.2, font = 2)
  }

  pdf("wind_scada_variables_time_series.pdf", width = 11, height = 8)
  draw_timeseries_base()
  dev.off()

  if (capabilities("png")) {
    tryCatch({
      png("wind_scada_variables_time_series.png", width = 11, height = 8, units = "in", res = 300)
      draw_timeseries_base()
      dev.off()
    }, error = function(e) invisible(NULL))
  }
  cat("Saved: wind_scada_variables_time_series.pdf (and .png)\n")

  # Function to render base R power curve
  draw_power_curve_base <- function() {
    par(mar = c(4.5, 4.5, 3.5, 1.5))
    is_outlier <- as.logical(data$is_curtailment_outlier)
    plot(data$wind_speed[!is_outlier], data$Y_MW[!is_outlier], pch = 16, cex = 0.6,
         col = rgb(31/255, 119/255, 180/255, 0.4),
         xlab = "Wind Speed (m/s)", ylab = "Active Power Output Y (MW)",
         main = "Empirical Wind Turbine Power Curve: Y (MW) vs Wind Speed (m/s)")
    points(data$wind_speed[is_outlier], data$Y_MW[is_outlier], pch = 16, cex = 0.7,
           col = rgb(214/255, 39/255, 40/255, 0.6))
    grid(col = "gray85")
    legend("topleft", legend = c("Normal Operation", "Curtailment Outlier"),
           col = c("#1f77b4", "#d62728"), pch = 16, bty = "n", cex = 0.95)
  }

  pdf("wind_scada_power_curve.pdf", width = 8.5, height = 6)
  draw_power_curve_base()
  dev.off()

  if (capabilities("png")) {
    tryCatch({
      png("wind_scada_power_curve.png", width = 8.5, height = 6, units = "in", res = 300)
      draw_power_curve_base()
      dev.off()
    }, error = function(e) invisible(NULL))
  }
  cat("Saved: wind_scada_power_curve.pdf (and .png)\n")
}




# Test for change in any of the power curve parameters:
# H0: beta(t) is constant vs H1: beta(t) has structural change
cusum_res <- CUSUM.regression(
  Y = Y, 
  X = X, 
  k = 0.65, 
  loss = "Welsh", 
  c = 2.985, 
  MC = 500, 
  linearized = TRUE, 
  plotting = TRUE
)

cat("Omnibus Test Statistic Z:", cusum_res$stat, "\n")
cat("Critical Value (95%):", cusum_res$crit_value, "\n")
cat("P-value:", cusum_res$p_value, "\n")
cat("Estimated Break Location u:", cusum_res$break_u, "\n")
