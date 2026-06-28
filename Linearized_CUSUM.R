# ─────────────────────────────────────────────
# 3. Final CUSUM SEDCD with Plotting
# ─────────────────────────────────────────────
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
