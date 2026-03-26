#Author Antotnio Panico

# Load the required library
library(MASS)
source("network_functions.R")

# Simulation of the null distribution of the CUMSUM test for Network Density 
# Set initial varaibele 

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
i <-1 
  for (j in 1:M) {
    
    #DGP = mu + L %*% rnorm(d) for the null hypothesis
    
    X <- mvrnorm(n = n[i], mu = mu, Sigma = sigma)
  
    W <- floor((n[i])^(2/3)) # set the window length for the local estimator
    
    #for w from 1 to n-1 compute the local estimator and the CUMSUM test statistic
  
    for (w in 1:(n[i] - W)) {
      
    #filter out 
    }
    
    
    
    
    
  }
  


