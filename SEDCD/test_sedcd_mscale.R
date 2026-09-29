# ==============================================================================
# File: test_sedcd_mscale.R
# Description: Automated Unit Test & Validation Suite for SEDCD Refactoring
#              (Option B: Robust M-Scale & Section 5.3 GAS Copula DGP)
# ==============================================================================

cat("======================================================================\n")
cat("STARTING AUTOMATED UNIT TESTS: SEDCD OPTION B & GAS COPULA DGP\n")
cat("======================================================================\n\n")

# Robust path resolution for sourcing required scripts
find_and_source <- function(fname) {
  cands <- c(
    fname,
    file.path("SEDCD", fname),
    file.path("..", "SEDCD", fname),
    file.path("Multiscale_CPD", "SEDCD", fname)
  )
  for (cand in cands) {
    if (file.exists(cand)) {
      source(cand)
      return(invisible(TRUE))
    }
  }
  stop(sprintf("Could not locate required file: %s", fname))
}

find_and_source("sedcd_mscale.R")
find_and_source("dgp_sedcd.R")

test_failures <- 0

assert_test <- function(test_name, condition, details = "") {
  if (isTRUE(condition)) {
    cat(sprintf("[PASS] %s\n", test_name))
  } else {
    cat(sprintf("[FAIL] %s - %s\n", test_name, details))
    test_failures <<- test_failures + 1
  }
}

# ------------------------------------------------------------------------------
# Test 1: Fisher Consistency Constant delta & Score Boundedness
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 1: M-Scale Score & Consistency Constant delta ---\n")

c_val <- 2.985
delta_theoretical <- 1.0 - (1.0 + 2.0 / (c_val^2))^(-0.5)

# Numerical integration: E[chi_c(Z)] for Z ~ N(0, 1)
integrand <- function(z) {
  (1.0 - exp(-(z / c_val)^2) - delta_theoretical) * dnorm(z)
}
int_res <- integrate(integrand, lower = -Inf, upper = Inf)$value

assert_test(
  "Theoretical delta matches Fisher consistency integral",
  abs(int_res) < 1e-6,
  sprintf("Integral error: %e", abs(int_res))
)

# Boundedness: sup |chi_c(u)| <= 1
u_grid <- seq(-100, 100, length.out = 1000)
chi_vals <- score_welsch_chi(u_grid, c = c_val)
assert_test(
  "Score chi_c(u) is uniformly bounded in [-1, 1]",
  all(abs(chi_vals) <= 1.0),
  sprintf("Max absolute value: %f", max(abs(chi_vals)))
)

# ------------------------------------------------------------------------------
# Test 2: 1D Robust M-Scale Root Finder
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 2: Pilot 1D Robust M-Scale Solver ---\n")

set.seed(42)
n_test <- 5000
true_sigma <- 2.5
x_gaussian <- rnorm(n_test, mean = 1.5, sd = true_sigma)

sigma_hat <- solve_pilot_mscale(x_gaussian, mu = 1.5, c = c_val)
assert_test(
  "M-scale estimates true Gaussian standard deviation accurately",
  abs(sigma_hat - true_sigma) / true_sigma < 0.05,
  sprintf("Estimated: %.4f, True: %.4f", sigma_hat, true_sigma)
)

# Check zero-root condition: mean(chi_c((x - mu)/sigma_hat)) ~ 0
root_val <- mean(score_welsch_chi((x_gaussian - 1.5) / sigma_hat, c = c_val))
assert_test(
  "M-scale satisfies root condition |mean(chi_c)| < 1e-5",
  abs(root_val) < 1e-5,
  sprintf("Residual score: %e", abs(root_val))
)

# Robustness check under heavy-tailed Student-t5
x_t5 <- rt(n_test, df = 5) * true_sigma
sigma_hat_t5 <- solve_pilot_mscale(x_t5, mu = 0, c = c_val)
assert_test(
  "M-scale converges to a stable positive root for Student-t5",
  is.finite(sigma_hat_t5) && sigma_hat_t5 > 0,
  sprintf("Estimated t5 scale: %.4f", sigma_hat_t5)
)

# ------------------------------------------------------------------------------
# Test 3: System Solvability & Estimating Function Zero-Mean
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 3: Local Estimator & 11D Estimating Function H(X, theta) ---\n")

set.seed(123)
M_n <- 400
X_win <- cbind(
  X1 = rnorm(M_n, mean = 0.5, sd = 1.2),
  X2 = rnorm(M_n, mean = -0.3, sd = 1.8)
)

tau_val <- -1.28
gamma_val <- 0.05

theta_hat <- estimate_local_theta(X_win, tau = tau_val, gamma = gamma_val, c = c_val)

assert_test(
  "estimate_local_theta produces parameter vector of length 11",
  length(theta_hat) == 11,
  sprintf("Length: %d", length(theta_hat))
)

assert_test(
  "Estimated SEDCD is in [0, 1]",
  theta_hat[11] >= 0 && theta_hat[11] <= 1,
  sprintf("SEDCD value: %.4f", theta_hat[11])
)

# Check that mean(H(X_window, theta_hat)) is identically 0 across all 11 equations
H_matrix <- compute_H_vector(X_win, theta_hat, tau = tau_val, gamma = gamma_val, c = c_val)
mean_H <- colMeans(H_matrix)
max_err_H <- max(abs(mean_H))

assert_test(
  "Sample average of H(X, theta_hat) is zero: ||mean(H)||_inf < 1e-5",
  max_err_H < 1e-5,
  sprintf("Max absolute error across 11 equations: %e", max_err_H)
)

# ------------------------------------------------------------------------------
# Test 4: Analytical vs Numerical Jacobian DH in R^(11 x 11)
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 4: Empirical Jacobian DH (Analytic vs numDeriv) ---\n")

DH_analytic <- compute_jacobian_DH(X_win, theta_hat, tau = tau_val, gamma = gamma_val,
                                  c = c_val, method = "analytic")
DH_numeric  <- compute_jacobian_DH(X_win, theta_hat, tau = tau_val, gamma = gamma_val,
                                  c = c_val, method = "numDeriv")

diff_DH <- abs(DH_analytic - DH_numeric)
max_diff_DH <- max(diff_DH)

assert_test(
  "Analytical Jacobian matches numDeriv within numerical tolerance",
  max_diff_DH < 1e-4,
  sprintf("Max absolute difference: %e", max_diff_DH)
)

# Check non-singularity and condition number
cond_DH <- kappa(DH_analytic)
assert_test(
  "Jacobian DH is non-singular and well-conditioned",
  is.finite(cond_DH) && cond_DH < 1e6,
  sprintf("Condition number kappa(DH): %.2f", cond_DH)
)

# ------------------------------------------------------------------------------
# Test 5: Dynamic GAS Copula Simulation (Clayton Copula, d = 2)
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 5: Dynamic Clayton GAS Copula Simulation (d = 2) ---\n")

dgp_h0 <- simulate_sedcd_dgp(n = 2000, d = 2, scenario = "H0", copula_type = "clayton", seed = 101)

assert_test(
  "DGP generates expected matrix dimensions (2000 x 2)",
  nrow(dgp_h0$X) == 2000 && ncol(dgp_h0$X) == 2,
  sprintf("Dimensions: %d x %d", nrow(dgp_h0$X), ncol(dgp_h0$X))
)

# Uniformity of copula margins: U1 and U2 should be Uniform(0, 1)
ks_u1 <- ks.test(dgp_h0$U[, 1], "punif")$p.value
ks_u2 <- ks.test(dgp_h0$U[, 2], "punif")$p.value
assert_test(
  "Copula marginals U_1 and U_2 are uniformly distributed (KS p > 0.01)",
  ks_u1 > 0.01 && ks_u2 > 0.01,
  sprintf("KS p-values: U1 = %.4f, U2 = %.4f", ks_u1, ks_u2)
)

# Adaptive threshold check: 10th percentile of X
emp_p_tau <- mean(dgp_h0$X <= dgp_h0$tau)
assert_test(
  "Adaptive threshold tau targets empirical lower 10th percentile",
  abs(emp_p_tau - 0.10) < 0.02,
  sprintf("Proportion below tau: %.4f (target 0.10), tau = %.3f", emp_p_tau, dgp_h0$tau)
)

# Alternatives generation check: H1 (abrupt break) & H2 (gradual trend)
dgp_h1 <- simulate_sedcd_dgp(n = 1000, d = 2, scenario = "H1", copula_type = "clayton", seed = 102)
dgp_h2 <- simulate_sedcd_dgp(n = 1000, d = 2, scenario = "H2", copula_type = "clayton", seed = 103)

assert_test(
  "DGP generates H1 with break in tail dependence at u = 0.5",
  mean(dgp_h1$lambda_L_path[1:500]) < mean(dgp_h1$lambda_L_path[501:1000]),
  sprintf("Pre-break lambda_L = %.3f, Post-break lambda_L = %.3f",
          mean(dgp_h1$lambda_L_path[1:500]), mean(dgp_h1$lambda_L_path[501:1000]))
)

assert_test(
  "DGP generates H2 with gradual drift in Clayton theta",
  mean(dgp_h2$theta_path[1:200]) < mean(dgp_h2$theta_path[801:1000]),
  sprintf("Theta trajectory: %.2f -> %.2f",
          mean(dgp_h2$theta_path[1:200]), mean(dgp_h2$theta_path[801:1000]))
)

# ------------------------------------------------------------------------------
# Test 6: Dynamic Student-t GAS Copula Simulation (d = 5)
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 6: Dynamic Student-t GAS Copula Simulation (d = 5) ---\n")

dgp_t5_h0 <- simulate_sedcd_dgp(n = 1000, d = 5, scenario = "H0", copula_type = "t", seed = 104)
dgp_t5_h1 <- simulate_sedcd_dgp(n = 1000, d = 5, scenario = "H1", copula_type = "t", seed = 105)

assert_test(
  "t-copula generates expected matrix dimensions (1000 x 5)",
  nrow(dgp_t5_h0$X) == 1000 && ncol(dgp_t5_h0$X) == 5,
  sprintf("Dimensions: %d x %d", nrow(dgp_t5_h0$X), ncol(dgp_t5_h0$X))
)

assert_test(
  "t-copula generates H1 with shift in dynamic correlation rho and tail dependence",
  mean(dgp_t5_h1$theta_path[1:500]) < mean(dgp_t5_h1$theta_path[501:1000]),
  sprintf("Pre rho = %.3f, Post rho = %.3f | Pre lambda_L = %.3f, Post lambda_L = %.3f",
          mean(dgp_t5_h1$theta_path[1:500]), mean(dgp_t5_h1$theta_path[501:1000]),
          mean(dgp_t5_h1$lambda_L_path[1:500]), mean(dgp_t5_h1$lambda_L_path[501:1000]))
)

# ------------------------------------------------------------------------------
# Test 7: End-to-End CUSUM Test Execution with New DGP
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 7: Full CUSUM Change-Point Test Pipeline ---\n")

t0 <- proc.time()
cusum_res <- cusum_sedcd_test(
  X = dgp_h0$X[1:600, ],
  tau = dgp_h0$tau,
  k = 25,
  lag = 4,
  block = 4,
  gamma = 0.05,
  c = c_val,
  MC = 50,
  plotting = FALSE
)
t_elapsed <- (proc.time() - t0)[3]

assert_test(
  "CUSUM test executes smoothly and returns valid test statistic and p-value",
  is.finite(cusum_res$test_stat) && is.finite(cusum_res$p_value) &&
    cusum_res$p_value >= 0 && cusum_res$p_value <= 1,
  sprintf("Test Stat = %.4f, p-value = %.4f, Elapsed = %.2fs",
          cusum_res$test_stat, cusum_res$p_value, t_elapsed)
)

cat("\n======================================================================\n")
if (test_failures == 0) {
  cat("ALL TESTS PASSED SUCCESSFULLY! (0 Failures)\n")
} else {
  cat(sprintf("TESTS FINISHED WITH %d FAILURE(S).\n", test_failures))
}
cat("======================================================================\n")
