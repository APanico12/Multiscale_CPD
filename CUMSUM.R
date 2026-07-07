#This file contains the CUMSUM functions
library(zoo)
source("Plot_functions.R")
source("LocalFunctions.R")

# ─────────────────────────────────────────────
# Single replication function
# ─────────────────────────────────────────────

CUSUM_TEST <- function(x,MC=1e3,ratio= 1/3,plotting=FALSE) {
  #x is a MTS of dimension n x d
  n = nrow(x)
  d= ncol(x)
  
  kn = floor(n^(ratio))
  D_array_raw <- rollapply(
    data      = x,
    width     = kn,
    FUN       = compute_local_density,
    by.column = FALSE,
    align     = "right",
    fill      = NA
  )
  
  D_array <- D_array_raw[!is.na(D_array_raw)]
  
  # Replace the Qn rollapply block entirely:
  Qn <- rollapply(
    data      = D_array,
    width     = kn,
    FUN       = loc_var,
    by.column = FALSE,
    align     = "right",
    fill      = NA
  )
  Cutoff <- sum(is.na(Qn))
  Qn     <- Qn[!is.na(Qn)]
  D_array <- D_array[-(1:Cutoff)]
  
  steps     <- length(D_array)
  m_n_full  <- cumsum(D_array) / steps
  total_m_n <- tail(m_n_full, 1)
  u_seq     <- seq_len(steps) / steps
  T_n_u     <- sqrt(steps) * (m_n_full - u_seq * total_m_n)
  Z         <- max(abs(T_n_u))
  Qn = cumsum(Qn)/steps

  qn= diff(c(0,Qn))
  # Bootstrap — 
  Z.mc <- replicate(MC, {
    BM <- cumsum(rnorm(steps) * sqrt(qn))
    max(abs(BM - u_seq * tail(BM, 1)))
  })
  
  
  q10 = quantile(Z.mc, 0.9)
  q5  = quantile(Z.mc, 0.95)
  
  if(plotting) {
    # Call the function from Plot_functions
    plot_cusum_test(u_seq, T_n_u, q5 )
  }
  
  p = mean(Z.mc > Z)
  return(p)
}
