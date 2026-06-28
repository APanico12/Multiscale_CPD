library(zoo)

create_moment_matrix = function(X) {
  n = nrow(X)
  d = ncol(X)
  
  part1 = X
  part2 = X^2
  
  pairs = list()
  for(i in 1:(d-1)) {
    for(j in (i+1):d) {
      pairs[[length(pairs) + 1]] = X[, i] * X[, j]
    }
  }
  part3 = do.call(cbind, pairs)
  
  Y = cbind(part1, part2, part3)
  return(Y)
}

EDC = function(m,d =5) {
  # m is an p vector of mean values of (z_i, z_i^2, z_i*z_j) for i,j=1,...,d
  
  EDC_mat = matrix(1, nrow = d, ncol = d)
  
  # 2. Extract standard deviations first for efficiency
  # SD = sqrt( E[z^2] - (E[z])^2 )
  means = m[1:d]
  sq_moments = m[(d + 1):(2 * d)]
  sds = sqrt(sq_moments - means^2)
  
  # 3. Compute correlations for off-diagonals
  cross_idx = 2 * d + 1
  
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      # Covariance = E[z_i * z_j] - E[z_i] * E[z_j]
      cov_ij = m[cross_idx] - (means[i] * means[j])
      
      # Correlation = Covariance / 
      rho_ij = cov_ij / (sds[i] * sds[j])
      
      EDC_mat[i, j] = rho_ij
      EDC_mat[j, i] = rho_ij
      
      cross_idx = cross_idx + 1
    }
  }
  
  # 4. Return density and matrix
  density = mean(EDC_mat)
  
  return(density)
}


mu.NW <- function(x,k){#local smoothing moment estimator
  #x: n*d time series
  #k: window size for local averaging
  res = 0*x
  k=ceiling(k)
  res[k:nrow(res),] = as.numeric(rollmean(x, k)) #using package 'zoo'
  #rollmean does not allow for changing window size at the beginning...
  for(i in 1:ncol(res)){
    res[1:(k-1),i] = cumsum(x[1:(k-1),i])/(1:(k-1))
  }
  res
}

int.par <- function(x,muhat,f, lag=0, cutoff=1){
  #x: Nxp matrix of data
  #muhat: Nxp matrix of local estimators
  #f: R^d to R
  
  N = nrow(x)
  lag = ceiling(lag)
  cutoff = ceiling(cutoff)
  if(!all.equal(dim(x),dim(muhat))){
    warning('dimensions of data x and estimator muhat do not match')
  }
  if(cutoff < lag){
    warning('cutoff may not be smaller than the lag')
  }
  
  #evaluate f
  vals    = matrix(0, nrow=N, ncol=1)
  for(t in (cutoff-lag+1):N) {
    Fmom      = f(muhat[t,])
    vals[t]   = Fmom
  }
  
  
  #partial sum process
  Mn = 1/N*cumsum(vals)
  return(Mn)
}

var.est <- function(x,muhat,f,block=1, lag=0, cutoff=1){#estimate integrated variance
  #x: data
  #muhat: local moment estimator
  #f: nonlinear moment transformation
  #block: block size for bootstrap
  #cutoff: 
  
  cutoff = ceiling(cutoff)
  lag    = ceiling(lag)
  block  = ceiling(block)
  
  if(!all.equal(dim(x),dim(muhat))){
    warning('dimensions of data x and estimator muhat do not match')
  }
  if(cutoff <= block + lag){
    warning('cutoff may not be smaller than the block+lag')
    cutoff = block + lag +1
  }
  
  N = nrow(x)
  
  #evaluate f
  vals    = matrix(0, nrow=N, ncol=1)
  for(t in (cutoff-lag+1):N) {
    Fmom      = EDC(muhat[t,])
    vals[t]   = Fmom
  }
  q = rep(0, N)
  for(t in (cutoff+1):N){
    blocksum = vals[(t-block), ] - vals[t,]
    q[t] = blocksum^2
  }
  Q = cumsum(q)/N
  return(Q)
}

CUSUM <- function(x, muhat, f, cutoff=1, lag=1, block=1, MC=1e4, avar=NULL, mu=NULL, plotting=FALSE){#returns bootstrap-based p-value of CUSUM test
  #mu: true moment
  #avar: true asymptotic variance
  
  N = nrow(x)
  lag = ceiling(lag)
  cutoff = ceiling(cutoff)
  
  if(!is.null(mu)){
    muhat[1:(N-lag),]=mu[(lag+1):N,] #this sets muhat_{t-lag} = mu_t
  }
  
  Mn = int.par(x,muhat,f, lag, cutoff)
  Qn = var.est(x,muhat,f,block,lag, cutoff)
  if(!is.null(avar)){
    Qn = avar
  }
  qn = diff(c(0,Qn))
  Tu = sqrt(N) * (Mn[(cutoff+1):N] - (1:(N-cutoff))/(N-cutoff) * Mn[length(Mn)])
  Tu = c(rep(0,cutoff), Tu)
  Z = max(abs(Tu)) #Test statistic
  
  Z.mc = rep(0, MC)
  for(i in 1:MC){
    BM = cumsum(rnorm(N)*sqrt(qn)) #multiplier Gaussian bootsrap
    Z.mc[i] = max(abs(BM[(cutoff+1):N]-BM[cutoff] - (1:(N-cutoff))/(N-cutoff)*BM[N]))
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


CUSUM.EDC <- function(X, k, cutoff=1, lag=0, b=0, MC=1e4, avar=NULL, mu=NULL){
  #X: matrix(!) of length N 
  Y = create_moment_matrix(X)
  muhat = mu.NW(Y, k) #k bandwidth for local smoothing
  CUSUM(Y, muhat, EDC, cutoff, lag, b, MC, avar, mu)
}








