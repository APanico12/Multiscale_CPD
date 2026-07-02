
# ==============================================================================
# Run_Analysis_SEDCD.R
# Purpose: Execution script for the Linearized SEDCD CUSUM Test
# ==============================================================================

cat("Loading Linearized SEDCD CUSUM framework...\n")

############################################################
# 1. LOAD LIBRARIES
############################################################

library(readxl)
library(copula)
library(zoo)

# Source the file containing the SEDCD functions
# The old source was "Linearized_CUSUM.R"
source("sedcd_try.R")


############################################################
# 2. IMPORT DATA
############################################################

file_path <- "C:/Users/Antonio/Desktop/euro_prices.xlsx"

if(!file.exists(file_path)){
    stop("Data file not found: ", file_path)
}

df <- read_excel(file_path)

str(df)


############################################################
# 3. BUILD MATRIX
############################################################

X_matrix <- as.matrix(df[, -1])

# Calculate arithmetic returns (price differences: p_t - p_{t-1}).
# This is the most suitable type of return for data with negative values.
# It also makes the time series stationary, which is a prerequisite
# for the dependence analysis that follows.
X_matrix <- diff(X_matrix)

# Remove any rows with missing values (NAs) that might result from the diff or were in the original data.
X_matrix <- na.omit(X_matrix)




# Select the first 5 countries (columns) for the analysis
X_matrix <- X_matrix[, 1:5]

N <- nrow(X_matrix)

cat(
    "Data loaded successfully:",
    N,
    "rows x",
    ncol(X_matrix),
    "columns\n"
)

############################################################
# 3.5 PREPROCESS DATA FOR EXTREMAL ANALYSIS
############################################################

# The SEDCD test is designed for data with heavy tails.
# To make the return data suitable for the test, we transform its marginal distributions
# to be heavy-tailed (Fréchet) while preserving the dependence structure (copula).

cat("Preprocessing data: Transforming marginals to be heavy-tailed (Fréchet)...\n")

# 1. Transform each column to pseudo-observations (Uniform on [0, 1])
# The pobs() function converts each observation to its empirical cumulative probability.
uniform_data <- pobs(X_matrix)
plot.zoo(uniform_data)
# 2. Transform from Uniform to Fréchet marginals
# The alpha parameter controls the heaviness of the tail. alpha=3 is used in the paper's simulations.
alpha <- 3
frechet_quantile_function <- function(u, a) {
  u_safe <- pmin(u, 1 - 1e-9) # Avoid log(0) which can happen with pobs
  return((-log(u_safe))^(-1/a))
}
X_matrix_transformed <- apply(uniform_data, 2, function(col) frechet_quantile_function(col, alpha))
plot.zoo(X_matrix_transformed)
############################################################
# 4. PARAMETER SETUP
############################################################

# Cross-validation for optimal bandwidth k
cat("Starting cross-validation for optimal bandwidth k...\n")

#' Cross-validation loss function to find the optimal bandwidth k.
#' This function implements the formula from the paper:
#' Lambda(k) = sum || H(X_{t+L}, \hat{theta}_t(k)) ||^2
cv_loss <- function(x, k, lag, tau, gamma) {
    N <- nrow(x)
    # Get all rolling parameter estimates for the given bandwidth k
    teta_k <- get.teta.sedcd(x, k, tau = tau, gamma = gamma)
    
    loss <- 0
    # Sum the squared norm of H over the validation range
    # The loop starts at k and ends at N-lag to ensure all indices are valid.
    for (t in k:(N - lag)) {
        # Parameter estimated at time t using window (t-k+1):t
        theta_t_k <- teta_k[t, ]
        # Data point from the future at t+lag
        X_future <- x[t + lag, ]
        # Calculate H for the future data point using parameters from the past
        loss <- loss + sum(Heq.sedcd(X_future, theta_t_k, tau = tau, gamma = gamma)^2)
    }
    return(loss)
}
# Define search space for k
k_grid <- floor(seq(N^0.35, N^0.65, length.out = 5))

# Define lag for cross-validation
L_n <- floor(0.1*log(N)^2)

cv_results <- sapply(k_grid, function(k_val) {
  cat("Testing k =", k_val, "\n")
  cv_loss(x = X_matrix_transformed, k = k_val, lag = L_n, tau = -3, gamma = 0.1)
})

# Find optimal k
optimal_k_index <- which.min(cv_results)
k <- k_grid[optimal_k_index]

cat("Cross-validation complete. Optimal k selected:", k, "\n")


# Initial observations discarded
cutoff <- k

# Separation between estimation and evaluation
lag <- L_n

# SEDCD threshold
# Changed from 2 to -3 to match the paper's focus on "downside" correlation.
tau <- -3

#block size
b <- L_n 
# Logistic smoothing parameter
gamma <- 1

# Bootstrap replications
MC <- 1000


############################################################
# 5. RUN LINEARIZED CUSUM
############################################################

cat("Executing Linearized SEDCD CUSUM...\n")

# Set a random seed for reproducibility of the bootstrap p-value
# Using the same seed will guarantee the exact same results on every run.
set.seed(123)

# The function name is CUSUM.sedcd (with a period) in sedcd_try.R
# The argument for block size is `block`, not `b`.
out <- CUSUM.sedcd(
  x = X_matrix_transformed, 
  k = k , 
  cutoff = cutoff, 
  lag = lag,
  block = b,
  tau = tau, 
  gamma = gamma,
  MC = 1e3, 
  plotting = TRUE) 
