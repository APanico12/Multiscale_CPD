library(numDeriv)
library(zoo)

solve.teta.autocorr <-function(X_window){
Zt =  X_window[,1]
Zth =  X_window[,2]

m1 = mean(Zt)

m2 = mean(Zth)

mu12 = mean(Zt*Zth)

var1 = var(Zt)

var2 = var(Zth)

rho = cov(Zt,Zth)/sqrt(var1*var2)

teta = c(m1,m2,mu12,var1,var2,rho)

return(teta)

}


get.teta.autocorr <-function(X,k){ # X = Zt,Zth
#check dim X = n*2

if(ncol(X)!=2) warning("The dimension of X shold be n*2")
  N = nrow(X)
  res <- matrix(0, nrow=N, ncol=6)
  for(t in k:N) {
    res[t, ] <- solve.teta.autocorr(X[(t-k+1):t,])
  }
  return(res)
  
}


Heq.autocorr <- function(X_window,teta){

#Xwindow = Zt,Zth
# if Xt is a matrix or vector with 2 columns, extract the columns, otherwise treat it as a vector of length 2
  vars <- tryCatch(
    expr = { list(Zt = X_window[, 1], Zth = X_window[, 2]) },
    error = function(e) { list(Zt = X_window[1], Zth = X_window[2]) }
  )

  Zt <- vars$Zt
  Zth <- vars$Zth
#extract the moments
m1 = teta[1]
m2 = teta[2]
mu12 = teta[3]
var1 = teta[4]
var2 = teta[5]
rho = teta[6]
#compute h_eq as in the paper
h1 = Zt-m1
h2 = Zth-m2
h3 = Zt*Zth-mu12
h4 = Zt^2-m1^2-var1
h5 = Zth^2-m2^2-var2
h6 = mu12 - m1*m2 - rho*sqrt(var1*var2)
res = cbind(h1,h2,h3,h4,h5,h6)
return(res)

}

#jacobian from numDeriv
DH.autocorr <- function(X_window,teta){
mean_H <- function(th) colMeans(Heq.autocorr(X_window, th))
  return(jacobian(mean_H, teta, method="simple"))
}

int.par.autocorr <- function(x, teta, lag=1, cutoff=1, k){
#check cutoff > k+lag
if(cutoff<k+lag) cutoff = k + lag 
#x = Zt,Zth
N = nrow(teta)
P= ncol(teta)
Lin.teta = matrix(0, nrow=N, ncol=P)

for(t in (cutoff:N)){
            teta_t_lag  = teta[t-lag,]
            ws = max(1, t - lag - k + 1)
            X_window = x[ws:(t - lag),]
            DH = DH.autocorr(X_window, teta_t_lag)
            DH_inv <- solve(DH)
            H_val <- Heq.autocorr(x[t,], teta_t_lag)
            Lin.teta[t, ] <- teta_t_lag - as.vector(DH_inv %*% t(H_val))
  }
  teta.autocorr = Lin.teta[,P]
return(cumsum(teta.autocorr)/length(teta.autocorr))

}

var.est.autocorr <- function(x, teta, block = 1, lag = 1, cutoff = 1, k){

  # This function estimates the integrated variance Q_n(u) from Theorem 3.2.
  # It has been corrected to follow the mathematical formula and the logic
  # of the reference `var.est` function.

  
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
    DH <- DH.autocorr(X_window, teta_t_lag)
    DH_inv <- solve(DH)
    H_val_block <- Heq.autocorr(x[(t - block + 1):t, ], teta_t_lag)
    q_par <- -DH_inv %*% colSums(H_val_block) # Result is a 6x1 vector
    q_sq[t, ] <- as.vector(q_par^2)
  }
  
  # 6. Calculate the integrated variance estimate Q_n(u)
  q_scaled <- q_sq / block
  Q.est_matrix <- apply(q_scaled, 2, cumsum) / N
  
  # Return the variance process for the parameter of interest (autocorrelation)
  return(Q.est_matrix[, P])

  }


CUSUM.autocorr <- function(x, teta = NULL, lag = 1, block = 1, cutoff = 1, k, MC = 1000,plotting = FALSE){

check.cutoff <- k+lag+block
if(cutoff<check.cutoff) cutoff = check.cutoff

# If teta is not provided, calculate it.
if(is.null(teta)) teta = get.teta.autocorr(x,k)

# N and P must be calculated *after* teta is guaranteed to exist.
N = nrow(teta)
P = ncol(teta)
Mn = int.par.autocorr(x, teta, lag=lag, cutoff=cutoff, k=k)
Qn = var.est.autocorr(x, teta, block = block, lag = lag, cutoff = cutoff, k = k)

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



###############
#Example of usage
###############
set.seed(123)
rSymgamma = function(n, alpha){#sample from symmetrized gamma distribution
  r = rgamma(n, rate=1, shape=alpha)*(2*rbinom(n,1,0.5)-1)/sqrt(alpha+alpha^2) #standardized to unit variance
  r
}

rAR1 = function(n, a, sigma, alpha){#sample from tvAR(1)-process
  X = rep(0,n)
  X[1]=rSymgamma(1, alpha(0))/sqrt(1-a(0)^2)*sigma(0)
  for(t in 2:n){
    X[t] = a(t/n)*X[t-1] + sigma(t/n)*rSymgamma(1, alpha(t/n))
  }
  X
}

sigma = function(u) abs(sin(u*2*pi)^2)+0.5
alpha = Vectorize(function(u) 1+1*(u>0.7))

a = function(u) 0.2
n = 10000
k = floor(n^0.65)
X = rAR1(n, a, sigma, alpha)
h =1
lag = floor(log(n))



X.autocorr = cbind(  X[(h+1):length(X)], X[1:(length(X)-h)])

teta = get.teta.autocorr(X.autocorr,k)

#heq = Heq.autocorr(X.autocorr,teta[500,])
int.autocorr.antonio = int.par.autocorr(X.autocorr, teta, lag=lag, cutoff=k+lag, k)
# try estimating the variance process Q_n(u)
var.autocorr.antonio = var.est.autocorr(X.autocorr, teta, block = 10, lag = lag, cutoff = k + lag + 10, k = k)
#try the CUSUM test
p_value.antonio = CUSUM.autocorr(X.autocorr, teta = teta, lag = lag, block = 10, cutoff = k + lag + 10, k = k, MC = 1000, plotting = TRUE)