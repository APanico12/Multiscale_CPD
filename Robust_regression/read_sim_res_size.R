# ==============================================================================
# File: read_sim_res_size.R
# Description: Generates publication-ready LaTeX Table 2 for Empirical Size
# Description: Generates publication-ready LaTeX Table 2 for Empirical Size under H0
# ==============================================================================

generate_latex_size_table <- function(summary_df, output_file = "table_size_output.tex") {
generate_latex_size_table <- function(summary_df, 
                                      output_file    = "table_size_output.tex",
                                      models         = c("LMHC" = "Model I (LMHC)", "LMAT" = "Model II (LMAT)"),
                                      epsilon_levels = c(0.00, 0.05, 0.10),
                                      models         = c("LMHC" = "Model I (LMHC)"),
                                      epsilon_levels = c(0.05),
                                      include_clean  = TRUE) {
  # Normalize column names
  col_names <- names(summary_df)
  lin_col   <- if ("rate_linearized" %in% col_names) "rate_linearized" else if ("rej_lin" %in% col_names) "rej_lin" else "rate_proposed_lin"
  plug_col  <- if ("rate_l2" %in% col_names) "rate_l2" else if ("rej_l2" %in% col_names) "rej_l2" else if ("rate_plugin" %in% col_names) "rate_plugin" else if ("rej_plug" %in% col_names) "rej_plug" else NULL
  
  has_plug <- !is.null(plug_col) && plug_col %in% col_names
  has_plug <- !is.null(plug_col) && (plug_col %in% col_names)
  
  n_vals <- sort(unique(summary_df$n))
  fmt <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return("-")
    sprintf("%.3f", as.numeric(x))
  }
  
  models <- c("LMHC" = "Model I (LMHC)", "LMAT" = "Model II (LMAT)")
  noises <- list(
    list(key = c("normal", "gaussian"), label = "$\\mathcal{N}(0, 1)$"),
    list(key = c("t3"),                 label = "$t_3$")
  )
  regimes <- list(
  
  all_regimes <- list(
    list(regime = "Clean", eps = 0.00),
    list(regime = "AO",    eps = 0.05),
    list(regime = "AO",    eps = 0.10),
    list(regime = "IO",    eps = 0.05),
    list(regime = "IO",    eps = 0.10),
    list(regime = "RO",    eps = 0.05),
    list(regime = "RO",    eps = 0.10)
  )
  
  regimes <- Filter(function(r) {
    if (r$regime == "Clean") return(isTRUE(include_clean))
    return(any(abs(r$eps - epsilon_levels) < 1e-4))
  }, all_regimes)
  
  model_caption_desc <- if (length(models) == 1) sprintf("under %s", models[1]) else "across Model~I (LMHC) and Model~II (LMAT)"
  eps_caption_desc   <- if (length(epsilon_levels) == 1 && abs(epsilon_levels[1] - 0.05) < 1e-4) {
    "across Clean ($\\varepsilon = 0.00$) and contaminated ($\\varepsilon = 0.05$) regimes"
  } else {
    "across various contamination regimes"
  }
  
  lines <- c()
  lines <- c(lines, "% ==============================================================================")
  lines <- c(lines, "% Table 2: Empirical Size Rejection Frequencies under H0 (alpha = 0.05)")
  lines <- c(lines, "% Table: Empirical Size Rejection Frequencies under H0 (alpha = 0.05)")
  lines <- c(lines, "% ==============================================================================")
  lines <- c(lines, "\\begin{table}[!htbp]")
  lines <- c(lines, "\\centering")
  lines <- c(lines, "\\caption{Empirical size rejection frequencies under the null hypothesis $\\mathcal{H}_0$ ($\\alpha = 0.05$) across $N_{\\text{sim}} = 1{,}000$ Monte Carlo replications for sample sizes $n \\in \\{500, 1{,}000, 5{,}000\\}$. Compares the Proposed Linearized Robust M-Test with the Classical $L_2$/OLS Test.}")
  lines <- c(lines, "\\label{tab:empirical_size}")
  lines <- c(lines, sprintf("\\caption{Empirical size rejection frequencies under the null hypothesis $\\mathcal{H}_0$ (nominal level $\\alpha = 0.05$) across $N_{\\text{sim}} = 1{,}000$ Monte Carlo replications for sample sizes $n \\in \\{500, 1{,}000, 5{,}000\\}$ %s. Compares the Proposed Linearized Robust M-Test with the Classical $L_2$/OLS Test %s.}",
  lines <- c(lines, sprintf("\\caption{Empirical size rejection frequencies under the null hypothesis $\\mathcal{H}_0$ (nominal level $\\alpha = 0.05$) across $N_{\\text{sim}} = 1{,}000$ Monte Carlo replications for sample sizes $n \\in \\{%s\\}$ %s. Compares the Proposed Linearized Robust M-Test with the Classical $L_2$/OLS Test %s.}",
                            paste(sprintf("%d", n_vals), collapse = ", "),
                            model_caption_desc, eps_caption_desc))
  lines <- c(lines, sprintf("\\label{tab:empirical_size%s}", if (length(models) == 1) paste0("_", tolower(names(models)[1])) else ""))
  lines <- c(lines, "\\vspace{0.2cm}")
  lines <- c(lines, "\\resizebox{\\textwidth}{!}{%")
  
  tot_model_cols <- length(n_vals)
  if (has_plug) {
    lines <- c(lines, "\\begin{tabular}{lllcccccccc}")
    lines <- c(lines, "\\toprule")
    lines <- c(lines, sprintf("& & & & \\multicolumn{%d}{c}{\\textbf{Proposed Linearized M-Test}} & & \\multicolumn{%d}{c}{\\textbf{Classical $L_2$/OLS Test}} \\\\", length(n_vals), length(n_vals)))
    lines <- c(lines, sprintf("\\cmidrule{5-%d} \\cmidrule{%d-%d}", 4 + length(n_vals), 6 + length(n_vals), 5 + 2 * length(n_vals)))
    lines <- c(lines, sprintf("& & & & \\multicolumn{%d}{c}{\\textbf{Proposed Linearized M-Test}} & & \\multicolumn{%d}{c}{\\textbf{Classical $L_2$/OLS Test}} \\\\", tot_model_cols, tot_model_cols))
    lines <- c(lines, sprintf("\\cmidrule{5-%d} \\cmidrule{%d-%d}", 4 + tot_model_cols, 6 + tot_model_cols, 5 + 2 * tot_model_cols))
    header_n <- paste(sprintf("$n = %d$", n_vals), collapse = " & ")
    lines <- c(lines, sprintf("\\textbf{Model} & \\textbf{Noise} & \\textbf{Regime} & $\\boldsymbol{\\varepsilon}$ & %s & & %s \\\\", header_n, header_n))
  } else {
    lines <- c(lines, sprintf("\\begin{tabular}{lllcccc}"))
    lines <- c(lines, "\\begin{tabular}{lllcccc}")
    lines <- c(lines, "\\toprule")
    lines <- c(lines, sprintf("& & & & \\multicolumn{%d}{c}{\\textbf{Proposed Linearized M-Test}} \\\\", length(n_vals)))
    lines <- c(lines, sprintf("\\cmidrule{5-%d}", 4 + length(n_vals)))
    lines <- c(lines, sprintf("& & & & \\multicolumn{%d}{c}{\\textbf{Proposed Linearized M-Test}} \\\\", tot_model_cols))
    lines <- c(lines, sprintf("\\cmidrule{5-%d}", 4 + tot_model_cols))
    header_n <- paste(sprintf("$n = %d$", n_vals), collapse = " & ")
    lines <- c(lines, sprintf("\\textbf{Model} & \\textbf{Noise} & \\textbf{Regime} & $\\boldsymbol{\\varepsilon}$ & %s \\\\", header_n))
  }
  lines <- c(lines, "\\midrule")
  
  n_rows_per_model <- length(noises) * length(regimes)
  
  for (m_idx in seq_along(models)) {
    m_key  <- names(models)[m_idx]
    m_name <- models[[m_key]]
    
    lines <- c(lines, sprintf("\\multirow{14}{*}{\\textbf{%s}}", m_name))
    lines <- c(lines, sprintf("\\multirow{%d}{*}{\\textbf{%s}}", n_rows_per_model, m_name))
    
    for (inn_idx in seq_along(noises)) {
      inn_keys <- noises[[inn_idx]]$key
      inn_name <- noises[[inn_idx]]$label
      
      for (r_idx in seq_along(regimes)) {
        reg  <- regimes[[r_idx]]$regime
        eps  <- regimes[[r_idx]]$eps
        
        # Filter matching rows
        sub_df <- summary_df[toupper(summary_df$model_type) == toupper(m_key) &
                             tolower(summary_df$innov_dist) %in% inn_keys &
                             toupper(summary_df$contamination) == toupper(reg) &
                             abs(summary_df$epsilon - eps) < 1e-4, ]
        
        lin_vals <- sapply(n_vals, function(nv) {
          row_n <- sub_df[sub_df$n == nv, ]
          if (nrow(row_n) > 0) fmt(row_n[[lin_col]][1]) else "-"
        })
        
        plug_vals <- if (has_plug) {
          sapply(n_vals, function(nv) {
            row_n <- sub_df[sub_df$n == nv, ]
            if (nrow(row_n) > 0) fmt(row_n[[plug_col]][1]) else "-"
          })
        } else NULL
        
        # Format line prefix
        noise_str <- if (r_idx == 1) inn_name else ""
        
        row_str <- if (has_plug) {
          sprintf("& %-20s & %-5s & %-4.2f & %s & & %s \\\\",
                  noise_str, reg, eps,
                  paste(lin_vals, collapse = " & "),
                  paste(plug_vals, collapse = " & "))
        } else {
          sprintf("& %-20s & %-5s & %-4.2f & %s \\\\",
                  noise_str, reg, eps,
                  paste(lin_vals, collapse = " & "))
        }
        lines <- c(lines, row_str)
      }
      
      if (inn_idx < length(noises)) {
        tot_cols <- if (has_plug) 5 + 2 * length(n_vals) else 4 + length(n_vals)
        lines <- c(lines, sprintf("\\cmidrule{2-%d}", tot_cols))
      }
    }
    
    if (m_idx < length(models)) {
      lines <- c(lines, "\\midrule")
    }
  }
  
  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular}%")
  lines <- c(lines, "}")
  lines <- c(lines, "\\vspace{0.15cm}")
  lines <- c(lines, "\\parbox{\\textwidth}{\\footnotesize \\textit{Note:} Nominal significance level is set to $\\alpha = 0.05$ across $N_{\\text{sim}} = 1{,}000$ Monte Carlo replications. Evaluated using multiplier bootstrap with $B = 100$ iterations.}")
  lines <- c(lines, sprintf("\\parbox{\\textwidth}{\\footnotesize \\textit{Note:} Nominal significance level is set to $\\alpha = 0.05$ across $N_{\\text{sim}} = 1{,}000$ Monte Carlo replications with $B = 100$ multiplier bootstrap draws. Evaluated under %s for Gaussian $\\mathcal{N}(0, 1)$ and Student-$t_3$ innovations across Clean ($\\varepsilon = 0.00$) and contamination regimes (AO, IO, RO) at $\\varepsilon = %s$.}",
  lines <- c(lines, sprintf("\\parbox{\\textwidth}{\\footnotesize \\textit{Note:} Nominal significance level is set to $\\alpha = 0.05$ across $N_{\\text{sim}} = 1{,}000$ Monte Carlo replications with $B = 100$ multiplier bootstrap draws. Evaluated under %s for Gaussian $\\mathcal{N}(0, 1)$ and Student-$t_3$ innovations across Clean ($\\varepsilon = 0.00$) and contamination regimes (AO, IO, RO) at $\\varepsilon \\in \\{%s\\}$.}",
                            paste(models, collapse = " and "),
                            paste(sprintf("%.2f", epsilon_levels[epsilon_levels > 0]), collapse = ", ")))
                            paste(sprintf("%.2f", sort(unique(c(if (include_clean) 0.00 else NULL, epsilon_levels)))), collapse = ", ")))
  lines <- c(lines, "\\end{table}")
  
  writeLines(lines, output_file)
  cat(sprintf("LaTeX Size Table successfully generated: %s\n", output_file))
}

# Standalone execution
if (sys.nframe() == 0) {
if (!interactive()) {
  input_csv <- "sim_summary_size.csv"
  if (!file.exists(input_csv)) input_csv <- "sim_results_size.csv"
  if (file.exists(input_csv)) {
  
  # Fallback to power results under H0 (delta == 0)
  if (!file.exists(input_csv) && file.exists("sim_summary_power.csv")) {
    df_power <- read.csv("sim_summary_power.csv")
    df <- subset(df_power, delta == 0 & hp_scenario == "H1")
    df_power <- read.csv("sim_summary_power.csv", stringsAsFactors = FALSE)
    df <- subset(df_power, delta == 0 & toupper(hp_scenario) == "H1")
    input_csv <- "sim_summary_power.csv (delta = 0, H1)"
  } else if (!file.exists(input_csv) && file.exists("sim_summary_power_n01.csv")) {
    df_power <- read.csv("sim_summary_power_n01.csv")
    df <- subset(df_power, delta == 0 & hp_scenario == "H1")
    df_power <- read.csv("sim_summary_power_n01.csv", stringsAsFactors = FALSE)
    df <- subset(df_power, delta == 0 & toupper(hp_scenario) == "H1")
    input_csv <- "sim_summary_power_n01.csv (delta = 0, H1)"
  } else if (file.exists(input_csv)) {
    df <- read.csv(input_csv)
    df <- read.csv(input_csv, stringsAsFactors = FALSE)
    if ("rej_lin" %in% names(df) && !("rate_linearized" %in% names(df))) {
      df <- aggregate(cbind(rej_lin, rej_plug) ~ model_type + n + innov_dist + contamination + epsilon,
      df <- aggregate(cbind(rej_lin, rej_l2) ~ model_type + n + innov_dist + contamination + epsilon,
                      data = df, FUN = mean)
      names(df)[names(df) == "rej_lin"]  <- "rate_linearized"
      names(df)[names(df) == "rej_plug"] <- "rate_plugin"
      names(df)[names(df) == "rej_l2"]   <- "rate_l2"
    }
  } else {
    df <- NULL
  }
  
  if (!is.null(df) && nrow(df) > 0) {
    cat(sprintf("Loaded data from: %s (%d rows)\n", input_csv, nrow(df)))
    cat(sprintf("Loaded size data from: %s (%d rows)\n", input_csv, nrow(df)))
    
    # 1. Generate full Table 2 (if multiple models available)
    generate_latex_size_table(df, "table_size_output.tex")
    # 1. Generate full Table 2
    avail_models <- unique(df$model_type)
    model_map    <- c("LMHC" = "Model I (LMHC)", "LMAT" = "Model II (LMAT)")
    use_models   <- model_map[avail_models[avail_models %in% names(model_map)]]
    if (length(use_models) == 0) use_models <- setNames(avail_models, avail_models)
    
    # 2. Generate specialized Table for Model I (LMHC) and eps = 0.05
    generate_latex_size_table(df, "table_size_lmhc_eps0.05.tex",
                              models = c("LMHC" = "Model I (LMHC)"),
                              epsilon_levels = c(0.05),
                              include_clean = TRUE)
    avail_eps <- setdiff(unique(df$epsilon), 0.00)
    if (length(avail_eps) == 0) avail_eps <- c(0.05)
    
    generate_latex_size_table(df, "table_size_output.tex",
                              models         = use_models,
                              epsilon_levels = avail_eps,
                              include_clean  = TRUE)
  } else {
    cat(sprintf("Input CSV not found: %s\n", input_csv))
    cat("Input CSV not found for size analysis.\n")
  }
}
