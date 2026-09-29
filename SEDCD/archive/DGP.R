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

Gen.from.t <- function(n = 2000, d = 5, rho = 0, nu = 5, margins_df = 5, time_varying = FALSE) {
  
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

# ==========================================
# Alternative DGP: Clayton Copula with linear trend in theta
# ==========================================
Gen.from.clayton.alt <- function(n, d, teta = 10 , marginal_df = 5) {
  U <- matrix(0, nrow = n, ncol = d)
  for (i in 1:n) {
    u <- i / n
    # Linear trend for theta: theta(u) = 2 + 4u
    theta <- teta * u
    cop <- claytonCopula(theta, dim = d)
    U[i, ] <- rCopula(1, cop)
  }
  X <- qt(U, df = marginal_df)
  return(X)
}

# ==========================================
# Alternative DGP: t-Copula with linear trend or jump in rho
# ==========================================
Gen.from.t.alt <- function(n, d, scenario = "linear", final_rho = 0.95, marginal_df = 5, copula_df = 5) {
  U <- matrix(0, nrow = n, ncol = d)
  rho_vec <- numeric(n)
   
  if (scenario == "linear") {
    # use linspace to create a linear sequence of rho values from 0 to final_rho
    rho_vec <- seq(0, final_rho, length.out = n)
  } else if (scenario == "jump") {
    break_point <- floor(n * 0.5)
    rho_vec[1:break_point] <- 0
    rho_vec[(break_point + 1):n] <- final_rho
  } else {
    stop("Scenario must be 'linear' or 'jump'")
  }

  for (i in 1:n) {
    # Create a vector of the single correlation parameter, repeated as many times as required
    # by the 'unstructured' dispersion structure for a given dimension 'd'.
    num_params <- d * (d - 1) / 2
    if (num_params > 0) {
      param_vector <- rep(rho_vec[i], num_params)
    } else {
      param_vector <- numeric(0) # For d=1, no correlation parameters
    }
    # Note: tCopula requires df > 2 for the correlation to be defined.
    # The default copula_df = 4 is safe.
    cop <- tCopula(param = param_vector, df = copula_df, dim = d, dispstr = "un")
    U[i, ] <- rCopula(1, cop)
  }

  X <- qt(U, df = marginal_df)
  return(X)
}
