
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
library(rugarch)

# Source the file containing the SEDCD functions
# The old source was "Linearized_CUSUM.R"
source("sedcd_try.R")
source("Plot_functions.R")


############################################################
# 2. IMPORT DATA
############################################################

file_path <- "C:/Users/Antonio/Desktop/prices.csv"

if(!file.exists(file_path)){
    stop("Data file not found: ", file_path)
}

df <- read.csv(file_path)

str(df)
 
 

############################################################
# 3. BUILD MATRIX
############################################################
countries <- c("DE_LU", "IT_NORD", "FR", "AT", "PL", "ES", "BE", "NL")
colnames(df)<-c("Date",countries)
sel_countries <- c("DE_LU", "IT_NORD", "BE", "PL", "NL","FR")

# --- Plot prices and returns before analysis ---
cat("Generating and saving prices and returns plot...\n")

df_selected <- df[, c("Date", sel_countries)]

plot_prices_and_returns(df_selected, date_col_name = "Date",
  filename = "../paper/img/prices_returns.pdf",
  pdf_width = 455,  # Now interpreted as 455 points
  pdf_height = 711, # Now interpreted as 711 points
  units = "pt"      # Specify that the dimensions are in points
  )

cat("Plot saved to paper/img/prices_returns.pdf\n")

price_matrix <- as.matrix(df[, -1])

X_matrix <- diff(price_matrix) 

# Remove any rows with missing values (NAs) that might result from the diff or were in the original data.
X_matrix <- na.omit(X_matrix)
X_matrix <- X_matrix[, sel_countries]

cat("Applying GARCH(1,1) filter to remove volatility clustering...\n")

# Define a standard GARCH(1,1) model specification.
# We assume the mean of the price differences is zero (or close to it).
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm" # Use "std" for Student's t-distribution if tails are very fat
)

# Create a matrix to store the standardized residuals
standardized_residuals <- matrix(NA, nrow = nrow(X_matrix), ncol = ncol(X_matrix))
colnames(standardized_residuals) <- colnames(X_matrix)

# Loop over each country's series, fit the GARCH model, and extract standardized residuals
for (i in 1:ncol(X_matrix)) {
  cat("  - Fitting GARCH for:", colnames(X_matrix)[i], "\n")
  fit <- ugarchfit(spec = spec, data = X_matrix[, i], solver = 'hybrid')
  standardized_residuals[, i] <- residuals(fit, standardize = TRUE)
}

# Replace the raw differences with the GARCH-filtered standardized residuals
X_matrix <- na.omit(standardized_residuals)

plot.zoo(X_matrix)
N <- nrow(X_matrix)

cat(
  "Data loaded successfully:",
  N,
  "rows x",
  ncol(X_matrix),
  "columns\n"
)

############################################################
# 4. PARAMETER SETUP
############################################################

#' Cross-validation loss function to find the optimal bandwidth k.
#' This function implements the formula from the paper:
#' Lambda(k) = sum || H(X_{t+L}, \hat{theta}_t(k)) ||^2

# Since the data is now standardized residuals, we can set tau based on
# standard normal quantiles to isolate the extreme tail.
# We use the 5% quantile, as suggested.
tau = -1.5 # qnorm(0.05)
gamma <- 0.1
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

# Cross-validation for optimal bandwidth k
cat("Starting cross-validation for optimal bandwidth k...\n")

# Define search space for k
k_grid <- floor(seq(N^0.35, N^0.65, length.out = 5))

# Define lag for cross-validation
L_n <- floor(0.1*log(N)^2)

cv_results <- sapply(k_grid, function(k_val) {
  cat("Testing k =", k_val, "\n")
  cv_loss(x = X_matrix, k = k_val, lag = L_n, tau = tau, gamma = gamma)
})

# Find optimal k
optimal_k_index <- which.min(cv_results)
k <- k_grid[optimal_k_index]
#k = 140 we also test it an it yields the same restuls we will procede with the optimal k 
cat("Cross-validation complete. Optimal k selected:", k, "\n")


# Initial observations discarded
cutoff <- k 

# Separation between estimation and evaluation
lag <- L_n

#block size
b <- L_n

# Bootstrap replications
MC <- 1000


############################################################
# 5. RUN LINEARIZED CUSUM
############################################################

cat("Executing Linearized SEDCD CUSUM...\n")

# Set a random seed for reproducibility of the bootstrap p-value
# Using the same seed will guarantee the exact same results on every run.
set.seed(123)


out <- CUSUM.sedcd(
  x = X_matrix, 
  k = k , 
  cutoff = cutoff, 
  lag = lag,
  block = b,
  tau = tau, 
  gamma = gamma,
  MC = 1e3, 
  plotting = TRUE) 
