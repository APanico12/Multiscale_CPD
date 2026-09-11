# ==============================================================================
# File: read_sim_res_power.R
# Description: Generates publication-ready LaTeX Tables and Power Curve Plots
#              for Empirical Power across all contamination regimes, models,
#              transition types (Abrupt H1 vs Gradual H2), and test procedures.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

#' Format numeric value to 3 decimal places
fmt3 <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("-")
  sprintf("%.3f", as.numeric(x))
}

#' Generate LaTeX Table for Empirical Power for a Specific Contamination Regime
generate_latex_power_table <- function(summary_df, 
                                       output_file   = "table_power_output.tex",
                                       contam_filter = "Clean",
                                       eps_filter    = 0.00) {
  col_names <- names(summary_df)
  lin_col   <- if ("rate_linearized" %in% col_names) "rate_linearized" else if ("rej_lin" %in% col_names) "rej_lin" else "rate_proposed_lin"
  l2_col    <- if ("rate_l2" %in% col_names) "rate_l2" else if ("rej_l2" %in% col_names) "rej_l2" else if ("rate_plugin" %in% col_names) "rate_plugin" else if ("rej_plug" %in% col_names) "rej_plug" else NULL
  has_l2    <- !is.null(l2_col)
  
  # Filter to specified contamination regime if available
  sub_df <- summary_df
  if (!is.null(contam_filter) && ("contamination" %in% names(sub_df))) {
    match_rows <- which(toupper(sub_df$contamination) == toupper(contam_filter) &
                        abs(sub_df$epsilon - eps_filter) < 1e-4)
    if (length(match_rows) > 0) {
      sub_df <- sub_df[match_rows, ]
    } else {
      first_c <- sub_df$contamination[1]
      first_e <- sub_df$epsilon[1]
      sub_df  <- sub_df[sub_df$contamination == first_c & sub_df$epsilon == first_e, ]
      contam_filter <- first_c
      eps_filter    <- first_e
    }
  }
  
  n_vals     <- sort(unique(sub_df$n))
  delta_vals <- sort(unique(sub_df$delta))
  models     <- c("LMHC" = "Model I (LMHC)", "LMAT" = "Model II (LMAT)")
  hyp_types  <- c("H1" = "Abrupt ($\\mathcal{H}_1$)", "H2" = "Gradual ($\\mathcal{H}_2$)")
  
  lines <- c()
  lines <- c(lines, "% ==============================================================================")
  lines <- c(lines, sprintf("%% Empirical Power: %s (epsilon = %.2f)", contam_filter, eps_filter))
  lines <- c(lines, "% ==============================================================================")
  lines <- c(lines, "\\begin{table}[!htbp]")
  lines <- c(lines, "\\centering")
  lines <- c(lines, sprintf("\\caption{Empirical power rejection frequencies ($\\alpha = 0.05$) under %s contamination ($\\varepsilon = %.2f$) across shift magnitudes $\\delta \\in \\{0.00, 0.25, 0.50, 0.75, 1.00, 1.25\\}$ for sample sizes $n \\in \\{500, 1{,}000, 5{,}000\\}$. Compares the Proposed Linearized Robust M-Test with the Classical $L_2$/OLS Test under Abrupt ($\\mathcal{H}_1$) and Gradual ($\\mathcal{H}_2$) structural transitions.}",
                            contam_filter, eps_filter))
  lines <- c(lines, sprintf("\\label{tab:power_%s_eps%s}", tolower(contam_filter), sub("\\.", "", sprintf("%.2f", eps_filter))))
  lines <- c(lines, "\\vspace{0.2cm}")
  lines <- c(lines, "\\resizebox{\\textwidth}{!}{%")
  
  num_deltas <- length(delta_vals)
  
  if (has_l2) {
    col_spec <- paste0("lllc", paste(rep("c", num_deltas), collapse = ""))
    lines <- c(lines, sprintf("\\begin{tabular}{%s}", col_spec))
    lines <- c(lines, "\\toprule")
    delta_headers <- paste(sprintf("$\\delta = %.2f$", delta_vals), collapse = " & ")
    lines <- c(lines, sprintf("\\textbf{Model} & \\textbf{Change Type} & \\textbf{Sample Size} & \\textbf{Procedure} & %s \\\\", delta_headers))
    lines <- c(lines, "\\midrule")
    
    for (m_idx in seq_along(models)) {
      m_key  <- names(models)[m_idx]
      m_name <- models[[m_key]]
      
      lines <- c(lines, sprintf("\\multirow{%d}{*}{\\textbf{%s}}", length(hyp_types) * length(n_vals) * 2, m_name))
      
      for (h_idx in seq_along(hyp_types)) {
        h_key  <- names(hyp_types)[h_idx]
        h_name <- hyp_types[[h_key]]
        
        for (nv_idx in seq_along(n_vals)) {
          nv <- n_vals[nv_idx]
          
          row_match <- sub_df[toupper(sub_df$model_type) == toupper(m_key) &
                              toupper(sub_df$hp_scenario) == toupper(h_key) &
                              sub_df$n == nv, ]
          
          pow_vals_lin <- sapply(delta_vals, function(d_val) {
            d_match <- row_match[abs(row_match$delta - d_val) < 1e-4, ]
            if (nrow(d_match) > 0) fmt3(d_match[[lin_col]][1]) else "-"
          })
          pow_vals_l2 <- sapply(delta_vals, function(d_val) {
            d_match <- row_match[abs(row_match$delta - d_val) < 1e-4, ]
            if (nrow(d_match) > 0) fmt3(d_match[[l2_col]][1]) else "-"
          })
          
          type_prefix <- if (nv_idx == 1) sprintf("\\multirow{%d}{*}{%s}", length(n_vals) * 2, h_name) else ""
          n_prefix    <- sprintf("\\multirow{2}{*}{$n = %-4d$}", nv)
          
          lines <- c(lines, sprintf("& %-25s & %-20s & Proposed M     & %s \\\\",
                                   type_prefix, n_prefix, paste(pow_vals_lin, collapse = " & ")))
          lines <- c(lines, sprintf("& %-25s & %-20s & Classical $L_2$ & %s \\\\",
                                   "", "", paste(pow_vals_l2, collapse = " & ")))
          
          if (nv_idx < length(n_vals)) {
            lines <- c(lines, sprintf("\\cmidrule{3-%d}", 4 + num_deltas))
          }
        }
        
        if (h_idx < length(hyp_types)) {
          lines <- c(lines, sprintf("\\cmidrule{2-%d}", 4 + num_deltas))
        }
      }
      
      if (m_idx < length(models)) {
        lines <- c(lines, "\\midrule")
      }
    }
  } else {
    col_spec <- paste0("llc", paste(rep("c", num_deltas), collapse = ""))
    lines <- c(lines, sprintf("\\begin{tabular}{%s}", col_spec))
    lines <- c(lines, "\\toprule")
    delta_headers <- paste(sprintf("$\\delta = %.2f$", delta_vals), collapse = " & ")
    lines <- c(lines, sprintf("\\textbf{Model} & \\textbf{Change Type} & \\textbf{Sample Size} & %s \\\\", delta_headers))
    lines <- c(lines, "\\midrule")
    
    for (m_idx in seq_along(models)) {
      m_key  <- names(models)[m_idx]
      m_name <- models[[m_key]]
      
      lines <- c(lines, sprintf("\\multirow{%d}{*}{\\textbf{%s}}", length(hyp_types) * length(n_vals), m_name))
      
      for (h_idx in seq_along(hyp_types)) {
        h_key  <- names(hyp_types)[h_idx]
        h_name <- hyp_types[[h_key]]
        
        for (nv_idx in seq_along(n_vals)) {
          nv <- n_vals[nv_idx]
          
          pow_vals <- sapply(delta_vals, function(d_val) {
            row_match <- sub_df[toupper(sub_df$model_type) == toupper(m_key) &
                                toupper(sub_df$hp_scenario) == toupper(h_key) &
                                sub_df$n == nv &
                                abs(sub_df$delta - d_val) < 1e-4, ]
            if (nrow(row_match) > 0) fmt3(row_match[[lin_col]][1]) else "-"
          })
          
          type_prefix <- if (nv_idx == 1) sprintf("\\multirow{%d}{*}{%s}", length(n_vals), h_name) else ""
          row_str <- sprintf("& %-25s & $n = %-4d$ & %s \\\\",
                             type_prefix, nv, paste(pow_vals, collapse = " & "))
          lines <- c(lines, row_str)
        }
        
        if (h_idx < length(hyp_types)) {
          lines <- c(lines, sprintf("\\cmidrule{2-%d}", 3 + num_deltas))
        }
      }
      
      if (m_idx < length(models)) {
        lines <- c(lines, "\\midrule")
      }
    }
  }
  
  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular}%")
  lines <- c(lines, "}")
  lines <- c(lines, "\\vspace{0.15cm}")
  lines <- c(lines, sprintf("\\parbox{\\textwidth}{\\footnotesize \\textit{Note:} Regime: %s ($\\varepsilon = %.2f$). Nominal size $\\alpha = 0.05$ evaluated at $\\delta = 0.00$. Under $\\mathcal{H}_1$, break occurs at $\\tau^* = 0.50$ via $\\beta_1(u) = 1.0 + \\delta \\cdot \\mathbf{1}_{\\{u > 0.5\\}}$ (LMHC) and $\\beta_1(u) = 1.0 + \\delta \\cdot \\mathbf{1}_{\\{u > 0.5\\}}$ (LMAT). Under $\\mathcal{H}_2$, transition is gradual via $\\beta_1(u) = 1.0 + \\delta u$. Comparison between Proposed Linearized Robust M-Test and Classical $L_2$/OLS Test based on $N_{\\text{sim}} = 1{,}000$ Monte Carlo replications with $B = 100$ multiplier bootstrap draws.}",
                            contam_filter, eps_filter))
  lines <- c(lines, "\\end{table}")
  
  writeLines(lines, output_file)
  cat(sprintf("LaTeX Power Table successfully written to: %s\n", output_file))
}

#' Generate Comprehensive Appendix LaTeX Tables for ALL Contamination Regimes
generate_latex_power_appendix_table <- function(summary_df, output_file = "table_power_appendix_all.tex") {
  if (!("contamination" %in% names(summary_df))) {
    generate_latex_power_table(summary_df, output_file)
    return()
  }
  
  regimes <- unique(summary_df[, c("contamination", "epsilon")])
  regimes <- regimes[order(regimes$contamination, regimes$epsilon), ]
  
  all_lines <- c()
  all_lines <- c(all_lines, "% ==============================================================================")
  all_lines <- c(all_lines, "% APPENDIX: Empirical Power Across All Contamination Scenarios")
  all_lines <- c(all_lines, "% ==============================================================================\n")
  
  for (r_idx in seq_len(nrow(regimes))) {
    c_name <- regimes$contamination[r_idx]
    c_eps  <- regimes$epsilon[r_idx]
    
    temp_file <- tempfile(fileext = ".tex")
    generate_latex_power_table(summary_df, output_file = temp_file, contam_filter = c_name, eps_filter = c_eps)
    table_content <- readLines(temp_file)
    unlink(temp_file)
    
    single_tex <- sprintf("table_power_%s_eps%.2f.tex", c_name, c_eps)
    writeLines(table_content, single_tex)
    
    all_lines <- c(all_lines, table_content, "\n\\clearpage\n")
  }
  
  writeLines(all_lines, output_file)
  cat(sprintf("Comprehensive Appendix LaTeX Table saved to: %s (contains %d regimes)\n",
              output_file, nrow(regimes)))
}

#' Plot Faceted Empirical Power Curves Separate by Contamination Regime
generate_power_plots <- function(summary_df, 
                                 output_png   = "power_curves_comparison.png", 
                                 output_pdf   = "power_curves_comparison.pdf",
                                 appendix_pdf = "power_curves_appendix_all.pdf") {
  col_names <- names(summary_df)
  has_lin   <- "rate_linearized" %in% col_names || "rej_lin" %in% col_names
  has_l2    <- "rate_l2" %in% col_names || "rej_l2" %in% col_names || "rate_plugin" %in% col_names || "rej_plug" %in% col_names
  
  lin_col <- if ("rate_linearized" %in% col_names) "rate_linearized" else "rej_lin"
  l2_col  <- if ("rate_l2" %in% col_names) "rate_l2" else if ("rej_l2" %in% col_names) "rej_l2" else if ("rate_plugin" %in% col_names) "rate_plugin" else "rej_plug"
  
  if (has_lin && has_l2) {
    df_lin <- summary_df
    df_lin$power  <- as.numeric(df_lin[[lin_col]])
    df_lin$Method <- "Linearized Robust M-Test (Proposed)"
    
    df_l2 <- summary_df
    df_l2$power  <- as.numeric(df_l2[[l2_col]])
    df_l2$Method <- "Classical L2 Test"
    
    df_plot <- rbind(df_lin, df_l2)
    df_plot$Method <- factor(df_plot$Method, levels = c("Linearized Robust M-Test (Proposed)", "Classical L2 Test"))
  } else {
    df_plot <- summary_df
    df_plot$power  <- as.numeric(df_plot[[lin_col]])
    df_plot$Method <- "Linearized Robust M-Test (Proposed)"
  }
  
  df_plot$N_label <- factor(paste0("n = ", df_plot$n), levels = paste0("n = ", sort(unique(df_plot$n))))
  
  df_plot$Hypothesis <- ifelse(toupper(df_plot$hp_scenario) == "H1",
                               "Abrupt Shift (H1, tau* = 0.50)",
                               "Gradual Shift (H2, beta0 + delta*u)")
  
  df_plot$Model <- ifelse(toupper(df_plot$model_type) == "LMHC",
                          "Model I: LMHC (Heteroscedastic)",
                          "Model II: LMAT (Autoregressive)")
  
  if ("contamination" %in% names(df_plot) && "epsilon" %in% names(df_plot)) {
    regimes <- unique(df_plot[, c("contamination", "epsilon")])
    regimes <- regimes[order(regimes$contamination, regimes$epsilon), ]
  } else {
    regimes <- data.frame(contamination = "Clean", epsilon = 0.00, stringsAsFactors = FALSE)
  }
  
  plot_list <- list()
  
  for (r_idx in seq_len(nrow(regimes))) {
    cur_contam <- regimes$contamination[r_idx]
    cur_eps    <- regimes$epsilon[r_idx]
    
    if ("contamination" %in% names(df_plot)) {
      sub_data <- df_plot[df_plot$contamination == cur_contam &
                          abs(df_plot$epsilon - cur_eps) < 1e-4, ]
    } else {
      sub_data <- df_plot
    }
    
    regime_title <- if (toupper(cur_contam) == "CLEAN") {
      "Clean Innovations (No Outliers)"
    } else if (toupper(cur_contam) == "AO") {
      sprintf("Additive Outliers (AO, epsilon = %.2f)", cur_eps)
    } else if (toupper(cur_contam) == "IO") {
      sprintf("Innovation Outliers (IO, epsilon = %.2f)", cur_eps)
    } else if (toupper(cur_contam) == "RO") {
      sprintf("Replacement Outliers (RO, epsilon = %.2f)", cur_eps)
    } else {
      sprintf("%s (epsilon = %.2f)", cur_contam, cur_eps)
    }
    
    p <- ggplot(sub_data, aes(x = delta, y = power, 
                              color = N_label, 
                              linetype = Method, 
                              shape = Method, 
                              group = interaction(N_label, Method))) +
      geom_hline(yintercept = 0.05, linetype = "dotted", color = "gray30", linewidth = 0.7) +
      geom_line(linewidth = 1.05) +
      geom_point(size = 2.4) +
      facet_grid(Model ~ Hypothesis) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                         limits = c(0, 1.02), 
                         breaks = seq(0, 1, 0.2)) +
      scale_x_continuous(breaks = sort(unique(sub_data$delta))) +
      scale_color_brewer(palette = "Set1") +
      scale_linetype_manual(values = c("Linearized Robust M-Test (Proposed)" = "solid", 
                                       "Classical L2 Test" = "dashed")) +
      scale_shape_manual(values = c("Linearized Robust M-Test (Proposed)" = 16, 
                                    "Classical L2 Test" = 17)) +
      labs(
        title    = sprintf("Empirical Power Curves: %s", regime_title),
        subtitle = expression(paste("Evaluation across shift magnitudes ", delta, " and sample sizes ", n, " | Nominal size ", alpha, " = 0.05")),
        x        = expression(paste("Shift Magnitude (", delta, ")")),
        y        = "Empirical Power P(Reject H0)",
        color    = "Sample Size",
        linetype = "Test Procedure",
        shape    = "Test Procedure"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle    = element_text(hjust = 0.5, color = "dimgray", size = 10),
        strip.background = element_rect(fill = "gray94", color = "gray80"),
        strip.text       = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.box       = "horizontal"
      )
    
    clean_tag <- sprintf("power_curves_%s_eps%.2f", cur_contam, cur_eps)
    reg_png   <- paste0(clean_tag, ".png")
    reg_pdf   <- paste0(clean_tag, ".pdf")
    
    ggsave(reg_png, plot = p, width = 10.5, height = 6.8, dpi = 300)
    tryCatch(ggsave(reg_pdf, plot = p, width = 10.5, height = 6.8), error = function(e) NULL)
    cat(sprintf("  --> Saved regime plot: %s and %s\n", reg_png, reg_pdf))
    
    if (toupper(cur_contam) == "CLEAN" || r_idx == 1) {
      ggsave(output_png, plot = p, width = 10.5, height = 6.8, dpi = 300)
      tryCatch(ggsave(output_pdf, plot = p, width = 10.5, height = 6.8), error = function(e) NULL)
    }
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  if (length(plot_list) > 0 && !is.null(appendix_pdf)) {
    tryCatch({
      pdf(appendix_pdf, width = 11, height = 7.5)
      for (plt in plot_list) {
        print(plt)
      }
      dev.off()
      cat(sprintf("Multi-Page Appendix PDF successfully saved: %s (%d contamination regimes)\n",
                  appendix_pdf, length(plot_list)))
    }, error = function(e) {
      cat(sprintf("Note: Could not compile multi-page appendix PDF (%s)\n", e$message))
    })
  }
}

# Standalone execution
if (sys.nframe() == 0) {
  input_csv <- "sim_summary_power.csv"
  if (!file.exists(input_csv)) input_csv <- "sim_results_power.csv"
  if (file.exists(input_csv)) {
    df <- read.csv(input_csv)
    if ("rej_lin" %in% names(df) && !("rate_linearized" %in% names(df))) {
      df <- aggregate(cbind(rej_lin, rej_plug) ~ model_type + n + innov_dist + contamination + epsilon + hp_scenario + delta,
                      data = df, FUN = mean)
      names(df)[names(df) == "rej_lin"]  <- "rate_linearized"
      names(df)[names(df) == "rej_plug"] <- "rate_plugin"
    }
    generate_latex_power_table(df, "table_power_output.tex", contam_filter = "Clean", eps_filter = 0.00)
    generate_latex_power_appendix_table(df, "table_power_appendix_all.tex")
    generate_power_plots(df, "power_curves_comparison.png", "power_curves_comparison.pdf", "power_curves_appendix_all.pdf")
  } else {
    cat(sprintf("Input CSV not found: %s\n", input_csv))
  }
}


