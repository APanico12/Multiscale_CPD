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

Gen.from.t <- function(n = 2000, d = 5, rho = 0.5, nu = 5, alpha_param = NULL, time_varying = FALSE) {
  
  # 1. Define the t-copula 
  t_cop <- tCopula(param = rho, dim = d, df = nu, dispstr = "ex")
  
  # 2. Check if there is an alpha_param
  if(!is.null(alpha_param)){
    
    # Define margins and bind them to the copula
    margins_dist <- rep("frechet", d)
    margins_params <- rep(list(list(shape = alpha_param)), d)
    
    my_model <- mvdc(copula = t_cop, 
                     margins = margins_dist, 
                     paramMargins = margins_params)
    
    # Generate the simulated data with Frechet margins
    sim_data <- rMvdc(n, my_model)
    
  } else {
    # If no alpha_param is provided, just generate raw uniform copula data
    sim_data <- rCopula(n, t_cop)
  }
  
  if(time_varying){sim_data = Apply.sine(sim_data)}
  
  # Return the actual simulated matrix, not just the blueprint
  return(sim_data)
}
# ==========================================
# Clayton Copula
# ==========================================

Gen.from.clayton <- function(n = 2000, d = 5, theta = 2, alpha_param = NULL, time_varying = FALSE) {
  
  # Define the Clayton blueprint
  clayton_cop <- claytonCopula(param = theta, dim = d)
  
  if(!is.null(alpha_param)){
    margins_dist <- rep("frechet", d)
    margins_params <- rep(list(list(shape = alpha_param)), d)
    
    my_model <- mvdc(copula = clayton_cop, 
                     margins = margins_dist, 
                     paramMargins = margins_params)
    sim_data <- rMvdc(n, my_model)
  } else {
    sim_data <- rCopula(n, clayton_cop)
  }

  if(time_varying){sim_data = Apply.sine(sim_data)}
  return(sim_data)
}

# ==========================================
# Frank Copula 
# ==========================================
Gen.from.frank <- function(n = 2000, d = 5, theta = 5, alpha_param = NULL, time_varying = FALSE) {
  
  # Define the Frank blueprint
  frank_cop <- frankCopula(param = theta, dim = d)
  
  if(!is.null(alpha_param)){
    margins_dist <- rep("frechet", d)
    margins_params <- rep(list(list(shape = alpha_param)), d)
    
    my_model <- mvdc(copula = frank_cop, 
                     margins = margins_dist, 
                     paramMargins = margins_params)
    sim_data <- rMvdc(n, my_model)
  } else {
    sim_data <- rCopula(n, frank_cop)
  }
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
