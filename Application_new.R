
# ==============================================================================
# Run_Analysis_SEDCD.R
# Purpose: Execution script for the Linearized SEDCD CUSUM Test
# ==============================================================================

cat("Loading Linearized SEDCD CUSUM framework...\n")

############################################################
# 1. LOAD LIBRARIES
############################################################

library(readxl)
library(zoo)

source("Linearized_CUSUM.R")


############################################################
# 2. IMPORT DATA
############################################################

file_path <- "C:/Users/Antonio/Desktop/log_yields.csv"

if(!file.exists(file_path)){
    stop("Data file not found: ", file_path)
}

df <- read.csv(file_path)

str(df)


############################################################
# 3. BUILD MATRIX
############################################################

X_matrix <- as.matrix(df[, -1])

rownames(X_matrix) <- df[,1]

# Select first five markets
#X_matrix <- X_matrix[,1:5]

# Remove missing values
X_matrix <- na.omit(X_matrix)

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

# Rolling window size
k <- ceiling(N^0.65)

# Initial observations discarded
cutoff <- ceiling(N^0.65)

# Separation between estimation and evaluation
lag <- 1

# SEDCD threshold
tau <- 2

#block size
b <- 1 
# Logistic smoothing parameter
gamma <- 0.1

# Bootstrap replications
MC <- 1000


############################################################
# 5. RUN LINEARIZED CUSUM
############################################################

cat("Executing Linearized SEDCD CUSUM...\n")

out <- CUSUM_SEDCD(
  x = X_matrix, 
  k = k , 
  cutoff = cutoff, 
  lag = lag,
  b = b,
  tau = tau, 
  gamma = gamma,
  MC = 1e3, 
  plotting = TRUE) 
