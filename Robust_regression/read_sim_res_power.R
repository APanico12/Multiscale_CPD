# ==============================================================================
# File: read_sim_res_power.R
# Description: Generates publication-ready LaTeX Tables and Power Curve Plots
#              for Empirical Power across all contamination regimes, models,
#              transition types (Abrupt H1 vs Gradual H2), and test procedures.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(ggplot2)
  })
}

#' Format numeric value to 3 decimal places
fmt3 <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("-")
  sprintf("%.3f", as.numeric(x))
}

#' Generate LaTeX Table for Empirical Power for a Specific Contamination Regime and Noise
generate_latex_power_table <- function(summary_df, 
                                       output_file   = "table_power_output.tex",
                                       contam_filter = "Clean",
                                       eps_filter    = 0.00,
                                       innov_filter  = "normal") {
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

  # Filter to specified innovation distribution if available
  if (!is.null(innov_filter) && ("innov_dist" %in% names(sub_df))) {
    t_targets <- if (tolower(innov_filter) %in% c("t3", "t_3")) c("t3", "t_3") else c("normal", "gaussian")
    inn_match <- which(tolower(sub_df$innov_dist) %in% t_targets)
    if (length(inn_match) > 0) {
      sub_df <- sub_df[inn_match, ]
    }
  }
  
  n_vals     <- sort(unique(sub_df$n))
  delta_vals <- sort(unique(sub_df$delta))
  all_models <- c("LMHC" = "Model I (LMHC)", "LMAT" = "Model II (LMAT)")
  present_models <- unique(toupper(sub_df$model_type))
  models     <- all_models[names(all_models) %in% present_models]
  if (length(models) == 0) models <- all_models["LMHC"]
  hyp_types  <- c("H1" = "Abrupt ($\\mathcal{H}_1$)", "H2" = "Gradual ($\\mathcal{H}_2$)")
  
  noise_tex <- if (is.null(innov_filter) || tolower(innov_filter) %in% c("normal", "gaussian")) {
    "Gaussian $\\mathcal{N}(0, 1)$"
  } else if (tolower(innov_filter) %in% c("t3", "t_3")) {
    "$t_3$"
  } else {
    innov_filter
  }
  
  lines <- c()
  lines <- c(lines, "% ==============================================================================")
  lines <- c(lines, sprintf("%% Empirical Power: %s (epsilon = %.2f, noise = %s)", contam_filter, eps_filter, noise_tex))
  lines <- c(lines, "% ==============================================================================")
  lines <- c(lines, "\\begin{table}[!htbp]")
  lines <- c(lines, "\\centering")
  lines <- c(lines, sprintf("\\caption{Empirical power rejection frequencies ($\\alpha = 0.05$) under %s contamination ($\\varepsilon = %.2f$) with %s innovations across shift magnitudes $\\delta \\in \\{0.00, 0.25, 0.50, 0.75, 1.00, 1.25\\}$ for sample sizes $n \\in \\{500, 1{,}000, 5{,}000\\}$. Compares the Proposed Linearized Robust M-Test with the Classical $L_2$/OLS Test under Abrupt ($\\mathcal{H}_1$) and Gradual ($\\mathcal{H}_2$) structural transitions.}",
                            contam_filter, eps_filter, noise_tex))
  inn_label_tag <- if (tolower(innov_filter) %in% c("t3", "t_3")) "t3" else "normal"
  lines <- c(lines, sprintf("\\label{tab:power_%s_%s_eps%s}", inn_label_tag, tolower(contam_filter), sub("\\.", "", sprintf("%.2f", eps_filter))))
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

#' Generate Comprehensive Appendix LaTeX Tables for ALL Contamination Regimes and Noises
generate_latex_power_appendix_table <- function(summary_df, output_file = "table_power_appendix_all.tex") {
  if (!("contamination" %in% names(summary_df))) {
    generate_latex_power_table(summary_df, output_file)
    return()
  }
  
  if ("innov_dist" %in% names(summary_df)) {
    regimes <- unique(summary_df[, c("innov_dist", "contamination", "epsilon")])
    regimes <- regimes[order(regimes$innov_dist, regimes$contamination, regimes$epsilon), ]
  } else {
    regimes <- unique(summary_df[, c("contamination", "epsilon")])
    regimes$innov_dist <- "normal"
    regimes <- regimes[order(regimes$contamination, regimes$epsilon), ]
  }
  
  all_lines <- c()
  all_lines <- c(all_lines, "% ==============================================================================")
  all_lines <- c(all_lines, "% APPENDIX: Empirical Power Across All Contamination Scenarios and Noises")
  all_lines <- c(all_lines, "% ==============================================================================\n")
  
  for (r_idx in seq_len(nrow(regimes))) {
    c_inn  <- regimes$innov_dist[r_idx]
    c_name <- regimes$contamination[r_idx]
    c_eps  <- regimes$epsilon[r_idx]
    
    temp_file <- tempfile(fileext = ".tex")
    generate_latex_power_table(summary_df, output_file = temp_file, contam_filter = c_name, eps_filter = c_eps, innov_filter = c_inn)
    table_content <- readLines(temp_file)
    unlink(temp_file)
    
    inn_tag <- if (tolower(c_inn) %in% c("t3", "t_3")) "t3" else "normal"
    single_tex <- sprintf("table_power_%s_%s_eps%.2f.tex", inn_tag, tolower(c_name), c_eps)
    writeLines(table_content, single_tex)
    
    all_lines <- c(all_lines, table_content, "\n\\clearpage\n")
  }
  
  writeLines(all_lines, output_file)
  cat(sprintf("Comprehensive Appendix LaTeX Table saved to: %s (contains %d regimes)\n",
              output_file, nrow(regimes)))
}

#' Plot Faceted Empirical Power Curves Separate by Contamination Regime and Noise
generate_power_plots <- function(summary_df, 
                                 output_png   = "power_curves_comparison.png", 
                                 output_pdf   = "power_curves_comparison.pdf",
                                 appendix_pdf = "power_curves_appendix_all.pdf") {
  col_names <- names(summary_df)
  has_lin   <- "rate_linearized" %in% col_names || "rej_lin" %in% col_names
  has_l2    <- "rate_l2" %in% col_names || "rej_l2" %in% col_names || "rate_plugin" %in% col_names || "rej_plug" %in% col_names
  
  lin_col <- if ("rate_linearized" %in% col_names) "rate_linearized" else "rej_lin"
  l2_col  <- if ("rate_l2" %in% col_names) "rate_l2" else if ("rej_l2" %in% col_names) "rej_l2" else if ("rate_plugin" %in% col_names) "rate_plugin" else "rej_plug"

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("Notice: 'ggplot2' is not installed in this R environment. Generating base R PDF power curves...\n")
    tryCatch({
      pdf_file <- if (!is.null(appendix_pdf)) appendix_pdf else output_pdf
      pdf(pdf_file, width = 11, height = 7)
      
      if ("innov_dist" %in% names(summary_df) && "contamination" %in% names(summary_df)) {
        groups <- unique(summary_df[, c("innov_dist", "contamination", "epsilon")])
        groups <- groups[order(groups$innov_dist, groups$contamination, groups$epsilon), ]
      } else if ("contamination" %in% names(summary_df)) {
        groups <- unique(summary_df[, c("contamination", "epsilon")])
        groups$innov_dist <- "normal"
      } else {
        groups <- data.frame(innov_dist = "normal", contamination = "Clean", epsilon = 0, stringsAsFactors = FALSE)
      }
      
      par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))
      for (g_idx in seq_len(nrow(groups))) {
        c_inn <- groups$innov_dist[g_idx]
        c_con <- groups$contamination[g_idx]
        c_eps <- groups$epsilon[g_idx]
        
        for (hp in c("H1", "H2")) {
          sub_d <- summary_df
          if ("contamination" %in% names(sub_d)) {
            sub_d <- sub_d[sub_d$contamination == c_con & abs(sub_d$epsilon - c_eps) < 1e-4, ]
          }
          if ("innov_dist" %in% names(sub_d)) {
            sub_d <- sub_d[tolower(sub_d$innov_dist) == tolower(c_inn), ]
          }
          sub_d <- sub_d[toupper(sub_d$hp_scenario) == hp, ]
          if (nrow(sub_d) == 0) next
          
          hp_title <- if (hp == "H1") "Abrupt Shift (H1, tau*=0.5)" else "Gradual Shift (H2)"
          main_title <- sprintf("%s | %s (eps=%.2f, %s)", hp_title, c_con, c_eps, c_inn)
          
          deltas <- sort(unique(sub_d$delta))
          plot(deltas, seq(0, 1, length.out = length(deltas)), type = "n", ylim = c(0, 1.05),
               xlab = expression(paste("Shift Magnitude ", delta)),
               ylab = expression(paste("Empirical Power ", P(Reject~H[0]))),
               main = main_title)
          abline(h = 0.05, lty = 3, col = "gray50")
          
          n_vals <- sort(unique(sub_d$n))
          cols <- c("blue", "forestgreen", "firebrick", "purple")
          
          for (ni in seq_along(n_vals)) {
            nv <- n_vals[ni]
            row_n <- sub_d[sub_d$n == nv, ]
            row_n <- row_n[order(row_n$delta), ]
            cur_col <- cols[(ni - 1) %% length(cols) + 1]
            
            lines(row_n$delta, row_n[[lin_col]], col = cur_col, lwd = 2, lty = 1)
            points(row_n$delta, row_n[[lin_col]], col = cur_col, pch = 16, cex = 1.1)
            
            if (has_l2) {
              lines(row_n$delta, row_n[[l2_col]], col = cur_col, lwd = 1.5, lty = 2)
              points(row_n$delta, row_n[[l2_col]], col = cur_col, pch = 17, cex = 0.9)
            }
          }
          leg_labels <- c(paste0("n = ", n_vals, " (Proposed)"), if (has_l2) paste0("n = ", n_vals, " (L2)") else NULL)
          leg_cols   <- c(cols[seq_along(n_vals)], if (has_l2) cols[seq_along(n_vals)] else NULL)
          leg_ltys   <- c(rep(1, length(n_vals)), if (has_l2) rep(2, length(n_vals)) else NULL)
          leg_pchs   <- c(rep(16, length(n_vals)), if (has_l2) rep(17, length(n_vals)) else NULL)
          legend("bottomright", legend = leg_labels, col = leg_cols, lty = leg_ltys, pch = leg_pchs, bty = "n", cex = 0.75)
        }
      }
      dev.off()
      cat(sprintf("  --> Base R PDF power curves successfully saved to: %s\n", pdf_file))
    }, error = function(e) {
      cat(sprintf("Note: Could not generate base R plot (%s)\n", e$message))
    })
    return(invisible(NULL))
  }
  
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
  
  if ("innov_dist" %in% names(df_plot) && "contamination" %in% names(df_plot)) {
    regimes <- unique(df_plot[, c("innov_dist", "contamination", "epsilon")])
    regimes <- regimes[order(regimes$innov_dist, regimes$contamination, regimes$epsilon), ]
  } else if ("contamination" %in% names(df_plot) && "epsilon" %in% names(df_plot)) {
    regimes <- unique(df_plot[, c("contamination", "epsilon")])
    regimes$innov_dist <- "normal"
    regimes <- regimes[order(regimes$contamination, regimes$epsilon), ]
  } else {
    regimes <- data.frame(innov_dist = "normal", contamination = "Clean", epsilon = 0.00, stringsAsFactors = FALSE)
  }
  
  plot_list <- list()
  
  for (r_idx in seq_len(nrow(regimes))) {
    cur_inn    <- regimes$innov_dist[r_idx]
    cur_contam <- regimes$contamination[r_idx]
    cur_eps    <- regimes$epsilon[r_idx]
    
    if ("contamination" %in% names(df_plot)) {
      sub_data <- df_plot[df_plot$contamination == cur_contam &
                          abs(df_plot$epsilon - cur_eps) < 1e-4, ]
    } else {
      sub_data <- df_plot
    }
    
    if ("innov_dist" %in% names(df_plot)) {
      inn_targets <- if (tolower(cur_inn) %in% c("t3", "t_3")) c("t3", "t_3") else c("normal", "gaussian")
      sub_data <- sub_data[tolower(sub_data$innov_dist) %in% inn_targets, ]
    }
    
    noise_label <- if (tolower(cur_inn) %in% c("t3", "t_3")) "t3" else "Gaussian"
    
    regime_title <- if (toupper(cur_contam) == "CLEAN") {
      sprintf("Clean Innovations (%s, No Outliers)", noise_label)
    } else if (toupper(cur_contam) == "AO") {
      sprintf("Additive Outliers (%s | AO, epsilon = %.2f)", noise_label, cur_eps)
    } else if (toupper(cur_contam) == "IO") {
      sprintf("Innovation Outliers (%s | IO, epsilon = %.2f)", noise_label, cur_eps)
    } else if (toupper(cur_contam) == "RO") {
      sprintf("Replacement Outliers (%s | RO, epsilon = %.2f)", noise_label, cur_eps)
    } else {
      sprintf("%s (%s | epsilon = %.2f)", cur_contam, noise_label, cur_eps)
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
      scale_y_continuous(limits = c(0, 1.0), 
                         breaks = seq(0, 1, 0.2), 
                         labels = function(y) sprintf("%.1f", y)) +
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
    
    inn_tag   <- if (tolower(cur_inn) %in% c("t3", "t_3")) "t3" else "normal"
    clean_tag <- sprintf("power_curves_%s_%s_eps%.2f", inn_tag, tolower(cur_contam), cur_eps)
    reg_png   <- paste0(clean_tag, ".png")
    reg_pdf   <- paste0(clean_tag, ".pdf")
    
    ggsave(reg_png, plot = p, width = 10.5, height = 6.8, dpi = 300)
    tryCatch(ggsave(reg_pdf, plot = p, width = 10.5, height = 6.8), error = function(e) NULL)
    cat(sprintf("  --> Saved regime plot: %s and %s\n", reg_png, reg_pdf))
    
    if ((toupper(cur_contam) == "CLEAN" && inn_tag == "normal") || r_idx == 1) {
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

# Standalone execution when executed directly via Rscript (not when sourced)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
is_main_script <- length(file_arg) > 0 && grepl("read_sim_res_power\\.R$", file_arg)

if (is_main_script) {
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


