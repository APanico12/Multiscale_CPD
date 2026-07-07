#This file contains the functions for the DGP used into New_Simulation

# ----- Package dependencies --------------------------------------------------
if (!requireNamespace("copula", quietly = TRUE)) {
  install.packages("copula")
}
if (!requireNamespace("evd", quietly = TRUE)) {
  install.packages("evd")
}
library(copula)
library(evd)# for frechet marginals


# ==========================================
# t-Copula
# ==========================================

Gen.from.t <- function(n = 2000, d = 5, rho = 0.5, nu = 5, margins_df = 5, time_varying = FALSE) {
  
  # 1. Define the t-copula 
  t_cop <- tCopula(param = rho, dim = d, df = nu, dispstr = "ex")
  
  # 2. Apply Student's t margins
  margins_dist <- rep("t", d)
  # The paramMargins argument requires a list of lists.
  margins_params <- rep(list(list(df = margins_df)), d)
  my_model <- mvdc(copula = t_cop, margins = margins_dist, paramMargins = margins_params)
  sim_data <- rMvdc(n, my_model)
  
  if(time_varying){
    sim_data = Apply.sine(sim_data)}
  
  return(sim_data)
}
# ==========================================
# Clayton Copula
# ==========================================

Gen.from.clayton <- function(n = 2000, d = 5, theta = 2, margins_df = 5, time_varying = FALSE) {
  
  # Define the Clayton blueprint
  clayton_cop <- claytonCopula(param = theta, dim = d)
  
  # Apply Student's t margins
  margins_dist <- rep("t", d)
  margins_params <- rep(list(list(df = margins_df)), d)
  
  my_model <- mvdc(copula = clayton_cop, 
                   margins = margins_dist, 
                   paramMargins = margins_params)
  sim_data <- rMvdc(n, my_model)

  if(time_varying){sim_data = Apply.sine(sim_data)}
  return(sim_data)
}

# ==========================================
# Frank Copula 
# ==========================================
Gen.from.frank <- function(n = 2000, d = 5, theta = 5, margins_df = 5, time_varying = FALSE) {
  
  # Define the Frank blueprint
  frank_cop <- frankCopula(param = theta, dim = d)
  
  # Apply Student's t margins
  margins_dist <- rep("t", d)
  margins_params <- rep(list(list(df = margins_df)), d)
  
  my_model <- mvdc(copula = frank_cop, 
                   margins = margins_dist, 
                   paramMargins = margins_params)
  sim_data <- rMvdc(n, my_model)
  if(time_varying){sim_data = Apply.sine(sim_data)}
  return(sim_data)
}

# ==========================================
# Time varying factor
# ==========================================

Apply.sine <- function(X) {
  n <- nrow(X)
  # Generate the time-varying factor
  t_grid <- seq_len(n) / n
  c_t <- 1 + sin(2 * pi * t_grid) / 2
  #apply to each column
  X_new <- apply(X, MARGIN = 2, function(col) {
    col * c_t
  })
  return(X_new)
}
