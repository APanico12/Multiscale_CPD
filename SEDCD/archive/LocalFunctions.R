# Here are the function to compute EDC and Dn and Variance 

sigmoid <- function(z, gamma =0.1){
  1/(1+exp(-z*gamma))
}


edc <- function(x, tau=2){
  # x should be a n by d matrix
  n = nrow(x)
  d = ncol(x)
  for(i in 1:d){
    x[,i] = (x[,i] - mean(x[,i]))/sd(x[,i]) #standardize columns
    x[,i] = x[,i] * sigmoid(-x[,i] - tau) #smooth threshold, only consider data below tau
  }
  
  res = list()
  res$matrix  = cor(x)
  res$density = mean(res$matrix[upper.tri(res$matrix)])
  return(res)
}

# Local density wrapper for rollapply
compute_local_density <- function(x_window, tau= 2) {
  edc_result <- edc(x_window,tau)
  return(edc_result$density)
}

# Local variance wrapper for rollapply
loc_var <- function(D){
  v = (tail(D,1)-D[1])^2/2*length(D)
  return(v)
}





