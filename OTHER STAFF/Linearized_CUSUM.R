# ─────────────────────────────────────────────
# 3. Final CUSUM SEDCD with Plotting
# ─────────────────────────────────────────────

#' @name cv_loss
#' @description Computes the cross-validation loss for a given bandwidth k.
#' This function implements the predictive loss function defined in the paper:
#' L(k) = sum_{t=k}^{N-L_n} ||H(X_{t+L_n}, \hat{\theta}_t(k))||^2
#' @param x The input data matrix (N x d).
#' @param k The rolling window size (bandwidth) to test.
#' @param lag The predictive lag L_n.
#' @param tau The threshold for SEDCD calculation.
#' @param gamma The smoothing parameter for the logistic function.
#' @return The scalar cross-validation loss for the given k.
cv_loss <- function(x, k, lag, tau, gamma) {
  N <- nrow(x)
  
  # 1. Get local parameter estimates \hat{\theta}_t(k) for all t
  # This returns a matrix where each row is \hat{\theta}_t
  theta_hat_t <- get_theta(x, k, tau, gamma)
  
  total_loss <- 0
  
  # 2. Sum the squared norm of the predictive error
  # Loop from t=k to N - lag
  for (t in k:(N - lag)) {
    # Get the parameter estimate at time t
    theta_t <- theta_hat_t[t, ]
    
    # Get the future data point at t + lag
    X_future <- x[t + lag, , drop = FALSE]
    
    # Calculate H(X_{t+L_n}, \hat{\theta}_t(k))
    H_val <- H_eqn(X_future, theta_t, tau, gamma)
    
    # Add the squared L2 norm to the total loss
    total_loss <- total_loss + sum(H_val^2)
  }
  
  return(total_loss)
}

source("LinearLocalFunctions.R")
source("Plot_functions.R")

# ─────────────────────────────────────────────
# 3. Final CUSUM SEDCD with Plotting
# ─────────────────────────────────────────────
CUSUM_SEDCD <- function(x, k=0.65, cutoff=1, lag=1, b=1, tau=-3, gamma=0.1, MC=1e3, plotting=FALSE) {
  N <- nrow(x)
  if(k < 1) k <- ceiling(k * N)
  lag <- ceiling(lag)
  cutoff <- ceiling(cutoff)
  if (cutoff < k+ lag + b) cutoff <-  k + lag + b + 1 
  b <- ceiling(b)
  
  # 1. Estimate time-varying parameters
  thetahat <- get_theta(x, k, tau, gamma)
  P <- ncol(thetahat) # The last column is SEDCD
  
  # 2. Get Linearized Process (Mn) and its Integrated Variance (Qn)
  Mn_all <- int.par.implicit.multi(x, thetahat, lag, cutoff, k, tau, gamma)
  Qn_all <- estimate_var(x, thetahat, b, lag, cutoff, k, tau, gamma)
  Mn <- Mn_all[, P]
  Qn <- Qn_all[, P]
  qn <- pmax(diff(c(0, Qn)), 0) # Non-negative increments for variance
  
  # 3. Calculate Tu Process
  time_weights <- (1:(N-cutoff)) / (N-cutoff)
  Tu <- sqrt(N) * (Mn[(cutoff+1):N] - time_weights * Mn[N])
  Tu <- c(rep(0, cutoff), Tu)
  Z <- max(abs(Tu))
  
  # 4. Bootstrap MC for Critical Values
  Z.mc <- rep(0, MC)
  for(i in 1:MC) {
    BM <- cumsum(rnorm(N) * sqrt(qn))
    Z.mc[i] <- max(abs(BM[(cutoff+1):N] - BM[cutoff] - time_weights * BM[N]))
  }
  
  q10 <- quantile(Z.mc, 0.90)
  q5  <- quantile(Z.mc, 0.95)
  p_val <- mean(Z.mc > Z)
  if(plotting) {
    u_axis <- (1:N) / N
    # Call the function from Plot_functions
    plot_cusum_test(u_axis, Tu, q5 )
  }
 
  
  # Return the results
  return(list(
    test_stat = Z,
    p_value = p_val,
    SEDCD = thetahat[, P],
    Mn = Mn,
    Tu_process = Tu
  ))
}
