library(numDeriv)
library(zoo)

sigmoid <- function(z, gamma = 1.0){
  1/(1+exp(-z/gamma))
}

solve.teta.sedcd <- function(X_window, tau = -3, gamma = 1.0){
  # This function implements the parameter estimation for the SEDCD example.
  d <- ncol(X_window)
  # 1. Estimate marginal means and variances (mu_i, sigma_i^2)
  mu_vec <- colMeans(X_window)
  var_vec <- apply(X_window, 2, var)
  
  # 2. Standardize the data to get Z_i = (X_i - mu_i) / sigma_i
  Z <- scale(X_window, center = mu_vec, scale = sqrt(var_vec))
  
  # 3. Calculate the smoothed Y matrix using the formula from the paper:
  # Y_i = Z_i * S_gamma(tau - Z_i)
  Y_smoothed <- Z * sigmoid(tau - Z, gamma = gamma)
  
  # 4. Estimate moments of the smoothed variables
  nu_vec <- colMeans(Y_smoothed)      # E(Y_i)
  gamma_mat <- cov(Y_smoothed)       # Cov(Y_i, Y_j)
  
  # 5. Estimate the SEDCD value
  cor_mat <- cor(Y_smoothed)
  sedcd_val <- mean(abs(cor_mat[upper.tri(cor_mat)]))
  # 6. Combine all parameters into a single vector 'teta' matching the paper's order:
  # (SEDCD, means, variances, nu's, flattened covariance matrix)
  teta <- c(sedcd_val, mu_vec, var_vec, nu_vec, as.vector(gamma_mat))
  
  return(teta)
}


#sedcd uses MTS
get.teta.sedcd <-function(X, k, tau = -3, gamma = 1.0){ # X MTS
  N = nrow(X)
  d = ncol(X)
  # The number of columns must match the length of the 'teta' vector.
  # Length = 1 (SEDCD) + d (means) + d (variances) + d (nus) + d*d (cov matrix)
  p = 1 + 3*d + d^2
  res <- matrix(0, nrow=N, ncol=p)
  for(t in k:N) {
    res[t, ] <- solve.teta.sedcd(X[(t-k+1):t,], tau = tau, gamma = gamma)
  }
  return(res)

}


Heq.sedcd <- function(X_window, teta, tau = -3, gamma = 1.0){
  # Ensure X_window is a matrix, even for a single row
  if (is.null(dim(X_window))) {
    X_window <- matrix(X_window, nrow = 1)
  }
  d = ncol(X_window)
  n_obs <- nrow(X_window)
  # 1. Extract parameters from teta vector
  sedcd_param <- teta[1]
  mu_vec      <- teta[2:(1+d)]
  var_vec     <- teta[(2+d):(1+2*d)]
  nu_vec      <- teta[(2+2*d):(1+3*d)]
  gamma_mat   <- matrix(teta[(2+3*d):length(teta)], nrow=d, ncol=d)

  # 2. Compute Y_smoothed using local parameters
  Z <- scale(X_window, center = mu_vec, scale = sqrt(var_vec)) #? 
  Y_smoothed <- Z * sigmoid(tau - Z, gamma = gamma)

  # h_sedcd: 1 equation (a constraint on parameters, independent of X_window)
  pairs <- combn(d, 2)
  gamma_diag <- diag(gamma_mat)
  denom <- sqrt(gamma_diag[pairs[1,]] * gamma_diag[pairs[2,]])
  safe_denom <- ifelse(denom > 1e-9, denom, 1) # Avoid division by zero
  correlations <- gamma_mat[lower.tri(gamma_mat)] / safe_denom
  sedcd_from_gamma <- mean(abs(correlations))
  h_sedcd <- matrix(sedcd_param - sedcd_from_gamma, nrow = n_obs, ncol = 1) #? 
  # h_means: d equations
  h_means <- sweep(X_window, 2, mu_vec, FUN = "-")
  # h_vars: d equations
  h_vars <- sweep(X_window^2, 2, mu_vec^2, FUN = "-")
  h_vars <- sweep(h_vars, 2, var_vec, FUN = "-")
  # h_nus: d equations
  h_nus <- sweep(Y_smoothed, 2, nu_vec, FUN = "-")
  # h_covs: d*d equations
  h_covs <- matrix(0, nrow = n_obs, ncol = d*d)
  col_idx <- 1
  for (j in 1:d) {
    for (i in 1:d) {
      # H_ij = Y_i*Y_j - E[Y_i*Y_j] = Y_i*Y_j - (nu_i*nu_j + gamma_ij)
      h_covs[, col_idx] <- Y_smoothed[, i] * Y_smoothed[, j] - (nu_vec[i] * nu_vec[j] + gamma_mat[i, j])
      col_idx <- col_idx + 1
    }
  }
  # 4. Combine all h vectors in the correct order
  return(cbind(h_sedcd, h_means, h_vars, h_nus, h_covs))
}

#jacobian from numDeriv
DH.sedcd <- function(X_window, teta, tau, gamma){
  mean_H <- function(th) colMeans(Heq.sedcd(X_window, th, tau = tau, gamma = gamma))
  return(jacobian(mean_H, teta, method="simple"))
}

int.par.sedcd <- function(x, teta, lag=1, cutoff=1, k, tau = -3, gamma = 1.0){
#check cutoff > k+lag
if(cutoff< k + lag) cutoff = k + lag 
N = nrow(teta)
P= ncol(teta)
Lin.teta = matrix(0, nrow=N, ncol=P)
for(t in (cutoff:N)){
            teta_t_lag  = teta[t-lag,]
            ws = max(1, t - lag - k + 1)
            X_window = x[ws:(t - lag),]
            DH = DH.sedcd(X_window, teta_t_lag, tau = tau, gamma = gamma)
            DH_inv <- solve(DH)
            H_val <- Heq.sedcd(x[t,], teta_t_lag, tau = tau, gamma = gamma)
            Lin.teta[t, ] <- teta_t_lag - as.vector(DH_inv %*% t(H_val))
  }
  # The parameter of interest (SEDCD) is the first one
  teta.sedcd = Lin.teta[,1]
  return(cumsum(teta.sedcd)/N)
}

var.est.sedcd <- function(x, teta, block = 1, lag = 1, cutoff = 1, k, tau = -3, gamma = 1.0){
  # This function estimates the integrated variance Q_n(u) from Theorem 3.2.

  N <- nrow(teta)
  P <- ncol(teta)
  
  # Ensure cutoff is large enough to accommodate lags, blocks, and bandwidth
  min_cutoff <- k + lag + block
  if(cutoff < min_cutoff) {
    warning(paste("Cutoff is too small for the given k, lag, and block size. Setting cutoff to", min_cutoff))
    cutoff <- min_cutoff
  }

  q_sq <- matrix(0, nrow = N, ncol = P) # To store the squared values

  # --- Main Loop ---
  for(t in (cutoff:N)){
    # 1. Parameter estimate from the past to ensure (near) independence
    param_idx <- t - lag - block
    teta_t_lag  <- teta[param_idx, ]
    ws <- max(1, param_idx - k + 1)
    X_window <- x[ws:param_idx, ]
    DH <- DH.sedcd(X_window, teta_t_lag, tau = tau, gamma = gamma)
    DH_inv <- solve(DH)
    H_val_block <- Heq.sedcd(x[(t - block + 1):t, ], teta_t_lag, tau = tau, gamma = gamma)
    q_par <- -DH_inv %*% colSums(H_val_block) # Result is a Px1 vector
    q_sq[t, ] <- as.vector(q_par^2)
  }
  
  # 6. Calculate the integrated variance estimate Q_n(u)
  q_scaled <- q_sq / block
  Q.est_matrix <- apply(q_scaled, 2, cumsum) / N
  
  # Return the variance process for the parameter of interest (SEDCD, the first parameter)
  return(Q.est_matrix[, 1])

  }


CUSUM.sedcd <- function(x, teta = NULL, lag = 1, block = 1, cutoff = 1, k, tau = -3, gamma = 1.0, MC = 1000, plotting = FALSE){
check.cutoff <- k+lag+block
if(cutoff<check.cutoff) cutoff = check.cutoff
# If teta is not provided, calculate it. This is the main parameter estimation step.
if(is.null(teta)) teta = get.teta.sedcd(x, k, tau = tau, gamma = gamma)

# N and P must be calculated *after* teta is guaranteed to exist.
N = nrow(teta)
P = ncol(teta)
Mn = int.par.sedcd(x, teta, lag=lag, cutoff=cutoff, k=k, tau = tau, gamma = gamma)
Qn = var.est.sedcd(x, teta, block = block, lag = lag, cutoff = cutoff, k = k, tau = tau, gamma = gamma)

qn = diff(c(0,Qn))
Tu = sqrt(N) * (Mn[(cutoff+1):N] - (1:(N-cutoff))/(N-cutoff) * Mn[length(Mn)])
Tu = c(rep(0, cutoff), Tu)
Z = max(abs(Tu))

Z.mc = rep(0, MC)
for(i in 1:MC){
    BM = cumsum(rnorm(N)*sqrt(qn))
    Z.mc[i] = max(abs( BM[(cutoff+1):N]-BM[cutoff] - (1:(N-cutoff))/(N-cutoff)*BM[N]))
  }
  q10 = quantile(Z.mc, 0.9)
  q5  = quantile(Z.mc, 0.95)
  if(plotting){
    plot((1:N)/N, Tu, xlab="u", ylab="T(u)", type="l", ylim = c(-max(q5, Z)*1.1, max(q5,Z)*1.1))
    abline(h=q10,  lty=2)
    abline(h=-q10, lty=2)
    abline(h=q5,   lty=3)
    abline(h=-q5,  lty=3)
    legend("bottomleft", c("10% threshold", "5% threshold"), lty=c(2,3))
  }
  
p = mean(Z.mc > Z)
return(p)
}

# ###########################
# #Example of usage
# ###########################

source("DGP.R")
set.seed(123)
alpha_param = 3

X = Gen.from.clayton(n = 1000, d = 5, theta = 1, alpha_param = alpha_param, time_varying = FALSE)
X = rbind(X,Gen.from.clayton(n = 1000, d = 5, theta = 3, alpha_param = alpha_param, time_varying = FALSE))
#standardize the data using Variance stabilizing transformation

n = nrow(X)
d =ncol(X)
k = floor(n^0.65)
lag = floor(log(n))

tau_val <- -3
gamma_val <- 1.0

teta = get.teta.sedcd(X, k, tau = tau_val, gamma = gamma_val)
# int.sedcd.antonio = int.par.sedcd(X, teta, lag=lag, cutoff=k+lag, k, tau = tau_val, gamma = gamma_val)

# var.est.antonio = var.est.sedcd(X, teta, block = lag, lag = lag, cutoff = 1, k, tau = tau_val, gamma = gamma_val)

p_value_sedcd = CUSUM.sedcd(X, teta = teta, lag = lag, block = lag,
                            cutoff = k + lag + 10, k = k,
                            tau = tau_val, gamma = gamma_val,
                            MC = 1000, plotting = TRUE)
