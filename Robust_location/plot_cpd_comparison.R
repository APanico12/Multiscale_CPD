# ==============================================================================
# Visualization: Robust Change-Point Detection Comparison in Location
# Generates a 3x3 multi-panel plot:
#   - 3 Rows: Contamination Scenarios (Clean, AO, IO)
#   - 3 Columns: Hypotheses (H0 [Size], H1 [Abrupt Power], H2 [Gradual Power])
#   - Each cell divided into 2 parts: Left = Gaussian Noise, Right = Student-t3 Noise
#   - X-axis: Test Methods (Our Welsh, Hodges-Lehmann, Huberized CUSUM, Wilcoxon)
#   - Y-axis: Empirical Size / Power (Rejection Rate)
#   - Point size and text labels grow larger as sample size n increases
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

csv_file <- "sim_summary_cpd_comparison.csv"

if (!file.exists(csv_file)) {
  # Fallback to raw results if summary not found
  raw_file <- "sim_results_cpd_comparison.csv"
  if (!file.exists(raw_file)) {
    stop("Neither sim_summary_cpd_comparison.csv nor sim_results_cpd_comparison.csv found.")
  }
  raw_df <- read.csv(raw_file)
  df <- raw_df %>%
    group_by(hp_scenario, contamination, innov_dist, n) %>%
    summarise(
      reps       = n(),
      rate_our   = mean(rej_our, na.rm = TRUE),
      rate_hl    = mean(rej_hl, na.rm = TRUE),
      rate_huber = mean(rej_huber, na.rm = TRUE),
      rate_wmw   = mean(rej_wmw, na.rm = TRUE),
      .groups    = "drop"
    )
} else {
  df <- read.csv(csv_file)
}

cat(sprintf("Loaded summary data with %d rows.\n", nrow(df)))

# ------------------------------------------------------------------------------
# 1. Reshape Data to Long Format for Plotting
# ------------------------------------------------------------------------------

df_long <- df %>%
  pivot_longer(
    cols = c(rate_our, rate_hl, rate_huber, rate_wmw),
    names_to = "test_raw",
    values_to = "rejection_rate"
  ) %>%
  mutate(
    # Clean method names
    method = case_when(
      test_raw == "rate_our"   ~ "Our (Welsh)",
      test_raw == "rate_hl"    ~ "HL",
      test_raw == "rate_huber" ~ "Huber",
      test_raw == "rate_wmw"   ~ "Wilcoxon",
      TRUE                     ~ test_raw
    ),
    method = factor(method, levels = c("Our (Welsh)", "HL", "Huber", "Wilcoxon")),
    
    # Clean contamination scenario names (Rows)
    contamination_clean = case_when(
      tolower(contamination) == "clean" ~ "Clean",
      toupper(contamination) == "AO"    ~ "Additive Outliers (AO)",
      toupper(contamination) == "IO"    ~ "Innovation Outliers (IO)",
      TRUE                              ~ contamination
    ),
    contamination_clean = factor(contamination_clean, levels = c("Clean", "Additive Outliers (AO)", "Innovation Outliers (IO)")),
    
    # Clean hypothesis labels (Columns)
    hp_clean = case_when(
      hp_scenario == "H0" ~ "H[0]:~No~Shift~(Size)",
      hp_scenario == "H1" ~ "H[1]:~Abrupt~Shift~(Power)",
      hp_scenario == "H2" ~ "H[2]:~Gradual~Shift~(Power)",
      TRUE                ~ hp_scenario
    ),
    hp_clean = factor(hp_clean, levels = c("H[0]:~No~Shift~(Size)", "H[1]:~Abrupt~Shift~(Power)", "H[2]:~Gradual~Shift~(Power)")),
    
    # Clean innovation distribution labels (Left vs Right sub-columns)
    innov_clean = case_when(
      tolower(innov_dist) %in% c("gaussian", "normal") ~ "Gaussian~Noise",
      tolower(innov_dist) %in% c("t3", "student-t")    ~ "t[3]~Noise",
      TRUE                                             ~ innov_dist
    ),
    innov_clean = factor(innov_clean, levels = c("Gaussian~Noise", "t[3]~Noise")),
    
    # Factor for sample size n
    n_fac = factor(paste0("n = ", n), levels = paste0("n = ", sort(unique(n))))
  )

# Compute dynamic point size and text size scaling based on n
unique_n <- sort(unique(df$n))
size_mapping <- setNames(seq(2.5, by = 2.0, length.out = length(unique_n)), paste0("n = ", unique_n))
text_size_mapping <- setNames(seq(2.2, by = 0.9, length.out = length(unique_n)), paste0("n = ", unique_n))

# ------------------------------------------------------------------------------
# 2. Build Multi-Panel Figure
# ------------------------------------------------------------------------------

p <- ggplot(df_long, aes(x = method, y = rejection_rate, color = n_fac, shape = n_fac)) +
  # Reference line at nominal alpha = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  # Points dodged by sample size n
  geom_point(
    position = position_dodge(width = 0.65),
    aes(size = n_fac),
    alpha = 0.9
  ) +
  # Value labels formatted with label size growing with n
  geom_text(
    aes(
      label = sprintf("%.2f", rejection_rate),
      size = n_fac
    ),
    position = position_dodge(width = 0.65),
    vjust = -0.85,
    show.legend = FALSE,
    fontface = "bold"
  ) +
  # Scale sizes dynamically with n
  scale_size_manual(
    name = "Sample Size",
    values = size_mapping
  ) +
  scale_shape_manual(
    name = "Sample Size",
    values = c(16, 17, 15, 18)[1:length(unique_n)]
  ) +
  scale_color_manual(
    name = "Sample Size",
    values = c("#1f77b4", "#2ca02c", "#d62728", "#9467bd")[1:length(unique_n)]
  ) +
  # 3 Rows (Contamination) x 3 Columns (Hypothesis) split by Left (Gaussian) / Right (t3)
  facet_grid(
    contamination_clean ~ hp_clean + innov_clean,
    labeller = labeller(
      hp_clean    = label_parsed,
      innov_clean = label_parsed,
      contamination_clean = label_value
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.12),
    breaks = seq(0, 1, by = 0.25)
  ) +
  labs(
    title = "Empirical Size and Power Comparison of Robust Location Change-Point Tests",
    subtitle = "Comparing Proposed Welsh Linearized CUSUM against Hodges-Lehmann, Huberized CUSUM, and Wilcoxon tests",
    x = "Change-Point Test Method",
    y = "Rejection Rate (Empirical Size / Power)",
    caption = "Dashed horizontal line denotes nominal significance level alpha = 0.05. Larger markers correspond to larger sample sizes n."
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f0f2f5", color = "gray60", linewidth = 0.6),
    strip.text       = element_text(face = "bold", size = 9),
    axis.text.x      = element_text(angle = 35, hjust = 1, vjust = 1, face = "bold", size = 8.5),
    axis.title       = element_text(face = "bold", size = 10.5),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle    = element_text(size = 10, hjust = 0.5, color = "gray30", margin = margin(b = 10)),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linetype = "dotted"),
    panel.spacing    = unit(0.5, "lines")
  )

# ------------------------------------------------------------------------------
# 3. Save Outputs
# ------------------------------------------------------------------------------

output_png <- "plot_cpd_comparison.png"
output_pdf <- "plot_cpd_comparison.pdf"

ggsave(output_png, plot = p, width = 14, height = 9, dpi = 300)
ggsave(output_pdf, plot = p, width = 14, height = 9)

cat(sprintf("Successfully generated plots:\n  - %s\n  - %s\n", output_png, output_pdf))
