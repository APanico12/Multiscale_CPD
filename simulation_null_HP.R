#Author Antotnio Panico

# Load the required library
library(MASS)
library(zoo)
source("network_functions.R")

# Simulation of the null distribution of the CUMSUM test for Network Density 
# Set initial variables 

M <- 1000 # number of simulations
n <- c(100,200,500,1000,2000,5000) # number of observations
d <- 5 # number of dimensions (node of the networks)
tau <- 2 # we need this in EDC 


# Data generating process
mu <- rep(0, d) # mean vector under the null hypothesis
sigma <- diag(5) # covariance matrix under the null hypothesis


#create data under risk scenario
p <- 0.8 # correlation coefficient
sigma_risk <- matrix(p, nrow = d, ncol = d)
diag(sigma_risk) <- 1 # set the diagonal elements to 1


# Create a matrix to store the results
results <- matrix(NA, nrow = M, ncol = length(n))

# Simulation loop
  # Generate M samples of multivariate normal data under the null hypothesis

compute_local_density <- function(X_slice) {
  local_EDC <- EDC(X_slice, tau = 0.10)
  return(D_n(local_EDC))
}


i <-1 

  W <- floor((n[i])^(2/3)) # set the window length for the local estimator
  
  for (j in 1:M) {
    
    X_sim <- mvrnorm(n = n[i], mu = mu, Sigma = sigma)
  
    #for w from 1 to n-1 compute the local estimator and the CUMSUM test statistic
    
    D_array <- rollapply(data = X_sim, 
                         width = W, 
                         FUN = compute_local_density, 
                         by.column = FALSE, 
                         align = "right", 
                         fill = NA)
    
    # compute t statitics  
    m_n_full <- cumsum(D_array) / n[i]
    total_m_n <- tail(m_n_full, 1)
    steps <- length(m_n_full)
    u_seq <- seq(1, steps) / steps
    
    t_n_u <- m_n_full - (u_seq * total_m_n)
    
    # Final Test Statistic 
    test_stat <- max(abs(t_n_u))
    results[j,i] <- test_stat
  
  }

  
  
  


