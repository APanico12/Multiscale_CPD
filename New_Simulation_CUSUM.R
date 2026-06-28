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
source("LocalFunctions.R")
source("DGP.R")
source("Plot_functions.R")
source("CUMSUM.R")

# ─────────────────────────────────────────────
# Simulation settings
# ─────────────────────────────────────────────

MC.Bootstrap <-1e2
MC.levels <-c(1e3) #<- c(100,500)#,1000, 2000, 5000)

#Set up DGP
n_values = 2000; d = 2; rho = 0.5; nu = 5; 
theta.Clayton = 2 ; theta.Frank = 5
alpha_param = 4 #frechet 

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
  "Frank" = function(n_val) {
    Gen.from.frank(n = n_val, d = d, theta = theta.Frank,
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
  },
  "Frank_Sine" = function(n_val) {
    Gen.from.frank(n = n_val, d = d, theta = theta.Frank,
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
    stringsAsFactors = FALSE
  )
})

# Combine the list of grids into one 
sim_grid <- do.call(rbind, grid_list)

# ─────────────────────────────────────────────
# Single replication function
# ─────────────────────────────────────────────
run_one <- function(n, dgp_name, dgp_functions) {
  
  X_sim <- dgp_functions[[dgp_name]](n)
  p_value <- CUSUM_TEST(X_sim)
  return(p_value)
  
}

# ─────────────────────────────────────────────
# Parallel backend
# ─────────────────────────────────────────────

set.seed(123, kind = "L'Ecuyer-CMRG")

cores <- detectCores() 
cl    <- makeCluster(cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 123)
clusterExport(cl, c("dgp_functions", "run_one", "d", "rho", "nu", 
                    "theta.Clayton", "theta.Frank", "alpha_param"))


clusterEvalQ(cl, {
  source("DGP.R")
  source("CUMSUM.R")
  source("Plot_functions.R")
  source("LocalFunctions.R")
  library(zoo)
  library(copula)
  library(evd)
})

# ─────────────────────────────────────────────
# Parallel simulation
# ─────────────────────────────────────────────

cat("Running", nrow(sim_grid), "simulations across", cores, "cores...\n")

sim_results <- foreach(
  row      = iter(sim_grid, by = "row"),
  .combine = rbind
) %dopar% {
  stat <- run_one(
    n         = row$n,
    dgp_name  = row$dgp_name,
    dgp_functions = dgp_functions
  )
  data.frame(
    iteration = row$iteration,
    level     = row$MC_level,
    dgp_name  = row$dgp_name,
    n_values  = row$n,
    stat      = stat
  )
}

stopCluster(cl)


# ─────────────────────────────────────────────
# Generate PDF Plot 
# ─────────────────────────────────────────────
cat("Generating plot report...\n")

  
plot_pvalue_histograms_fabian(sim_results)
  


