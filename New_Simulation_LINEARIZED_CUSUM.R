# Code for: 
# Author: Antonio Panico

library(MASS)
library(zoo)
library(cumulcalib)
require(doParallel)
require(foreach)
require(SymTS)

# ─────────────────────────────────────────────
# Source Dependences
# ─────────────────────────────────────────────
source("LinearLocalFunctions.R")
source("DGP.R")
source("Plot_functions.R")
source("Linearized_CUSUM.R")

# ─────────────────────────────────────────────
# Simulation settings
# ─────────────────────────────────────────────

MC.Bootstrap <- 100
MC.levels <- c(1000) 

#Set up DGP
n_values = c(100, 500, 2000, 5000)
d = 5
rho = 0.5
nu = 5
theta.Clayton = 2
alpha_param = 4 #frechet 

#tau and L_n types
tau_values = c(-2, -2.5, -3)
L_n_types = c("log_n", "n_0.15", "n_0.20")

#Implement the DGP
dgp_functions <- list(
  
  # Base Copulas (Constant Marginals)
  "T_Copula" = function(n_val) {
    Gen.from.t(n = n_val, d = d, rho = rho, nu = nu, 
               alpha_param = alpha_param)
  },
  "Clayton" = function(n_val) {
    Gen.from.clayton(n = n_val, d = d, theta = theta.Clayton,
                     alpha_param = alpha_param)
  },
  
  # Time-Varying Copulas (Sine Wave Marginals)
  "T_Copula_Sine" = function(n_val) {
    Gen.from.t(n = n_val, d = d, rho = rho, nu = nu,
               alpha_param = alpha_param, time_varying = TRUE)
  },
  "Clayton_Sine" = function(n_val) {
    Gen.from.clayton(n = n_val, d = d, theta = theta.Clayton,
                     alpha_param = alpha_param, time_varying = TRUE)
   }
)

# ─────────────────────────────────────────────
# Build the Grid
# ─────────────────────────────────────────────

grid_list <- lapply(MC.levels, function(m) {
  expand.grid(
    iteration = 1:m,         
    MC_level  = m,           
    n         = n_values,
    dgp_name  = names(dgp_functions),
    tau       = tau_values,
    L_n_type  = L_n_types,
    stringsAsFactors = FALSE
  )
})

# Combine the list of grids into one 
sim_grid <- do.call(rbind, grid_list)

# ─────────────────────────────────────────────
# Optimal Bandwidth Selection function
# ─────────────────────────────────────────────
find_k_opt <- function(X, L_n, tau, gamma) {
  N <- nrow(X)
  k_min <- max(1, floor(N^0.35))
  k_max <- floor(N^0.55)
  # evaluate on a small grid of 5 points to save computation time
  k_grid <- unique(floor(seq(k_min, k_max, length.out = 5)))
  
  best_k <- k_max
  best_val <- Inf
  
  for (k in k_grid) {
    theta_hat <- get_theta(X, k, tau, gamma)
    
    val <- 0
    for (t in k:(N - L_n)) {
      h_val <- H_eqn(X[t + L_n, , drop = FALSE], theta_hat[t, ], tau, gamma)
      val <- val + sum(h_val^2)
    }
    
    if (val < best_val) {
      best_val <- val
      best_k <- k
    }
  }
  return(best_k)
}

# ─────────────────────────────────────────────
# Single replication function
# ─────────────────────────────────────────────
run_one <- function(n, dgp_name, dgp_functions, tau, L_n_type) {
  
  X_sim <- dgp_functions[[dgp_name]](n)
  
  if (L_n_type == "log_n") {
    lag <- max(1, floor(log(n)))
  } else if (L_n_type == "n_0.15") {
    lag <- max(1, floor(n^0.15))
  } else if (L_n_type == "n_0.20") {
    lag <- max(1, floor(n^0.20))
  }
  
  gamma <- 0.1
  
  k_opt <- find_k_opt(X_sim, lag, tau, gamma)
  
  b      <- max(1, floor(n^0.15))
  cutoff <- max(1, floor(n^0.7))
  
  res <- CUSUM_SEDCD(X_sim, k_opt, cutoff, lag, b, tau, gamma, MC=MC.Bootstrap)
  return(res$p_value)
  
}

# ─────────────────────────────────────────────
# Parallel backend
# ─────────────────────────────────────────────

set.seed(123, kind = "L'Ecuyer-CMRG")

cores <- detectCores() 
cl    <- makeCluster(cores)

registerDoParallel(cl)
clusterSetRNGStream(cl, 123)
clusterExport(cl, c("dgp_functions", "run_one", "find_k_opt", "d", "rho", "nu",
                    "theta.Clayton", "alpha_param", "MC.Bootstrap"))


clusterEvalQ(cl, {
  source("DGP.R")
  source("LinearLocalFunctions.R")
  source("Linearized_CUSUM.R")
  
  library(zoo)
  library(MASS)
  library(numDeriv)
  library(copula)
  library(evd)
})

# ─────────────────────────────────────────────
# Parallel simulation
# ─────────────────────────────────────────────

cat("Running", nrow(sim_grid), "simulations across", cores, "cores...\n")

sim_results <- foreach(
  row      = iter(sim_grid, by = "row"),
  .combine = rbind,
  .packages = c("zoo", "MASS", "numDeriv", "copula", "evd")
) %dopar% {
  stat <- run_one(
    n         = row$n,
    dgp_name  = row$dgp_name,
    dgp_functions = dgp_functions,
    tau       = row$tau,
    L_n_type  = row$L_n_type
  )
  data.frame(
    iteration = row$iteration,
    level     = row$MC_level,
    dgp_name  = row$dgp_name,
    n_values  = row$n,
    tau       = row$tau,
    L_n_type  = row$L_n_type,
    stat      = stat
  )
}

stopCluster(cl)


# ─────────────────────────────────────────────
# Generate PDF Plot 
# ─────────────────────────────────────────────
cat("Generating plot report...\n")

#pdf("Simulation_Results.pdf", width = 10, height = 8)
plot_pvalue_histograms_fabian(sim_results)

cat("Plot saved successfully to Simulation_Results.pdf\n")
