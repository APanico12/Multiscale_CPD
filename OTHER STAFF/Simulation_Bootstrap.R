#Code supplement for the manuscript "Functional estimation and change detection for non-stationary time series"
#2020-04-14
source("network_functions.R")
source("DGP.R")
source("Plot_functions.R")
source("NonLinearChange.R")
require(stabledist)
require(doParallel)
require(foreach)
require(SymTS)


#  ---- simulation setup ----

MC.bootstrap = 1e3
MC = 1e4
cores = 1 #total runtime approx. 5h, using 80 CPU cores
Ns = c(1e2, 5e2, 1e3, 5e3, 1e4) #sample sizes to be simulated

# EDC with known moments and asymptotic variance ----
d=5
VarCov = matrix(0.5, d, d); diag(VarCov) = 1
A = chol(VarCov)
n=1e6
X = gen_iid(A,n)
varX = var(X)
Y = create_moment_matrix(X)
means_vec = rep(0, d)                       # d means (all 0)
vars_vec  = diag(VarCov)                    # d variances
covs_vec  = VarCov[upper.tri(VarCov)]       # d(d-1)/2 covariances 
m_row = c(means_vec, vars_vec, covs_vec)
mu = matrix(m_row, nrow = n, ncol = length(m_row), byrow = TRUE)

# avar = var.est(Y, muhat=mu, f=autocor, block=n^0.15, lag=n^0.4, cutoff=n^0.6+1)
# plot(avar, type="l")

#perform MC simulation
pvals.EDC = c()

for(N in Ns){
  grid = (2:(N+1))/N
  Q = avar[ceiling((1:N)/N*n)]
  mu = 
  k = N^0.65
  b = N^0.15
  lag=N^0.45
  cutoff = N^0.7
  
  registerDoParallel(cores)
  p = foreach(i=1:MC, .combine=rbind, .packages=c("zoo")) %dopar% {
    #null hypothesis
    x = gen_iid(A,n)
    p = CUSUM.EDC(x, k, cutoff, lag, b, MC=MC.bootstrap, avar=NULL, mu=NULL)

  }
  stopImplicitCluster()
  pvals.autocor = cbind(pvals.EDC, p)
  par(mfrow=c(2,2))
  hist(p[,1], xlab="", main="Q known, mu known")
  hist(p[,2], xlab="", main="Q unknown, mu known")
  hist(p[,3], xlab="", main="Q known, mu unknown")
  hist(p[,4], xlab="", main="Q unknown, mu unknown")
  par(mfrow=c(1,1))
}

save("pvals.autocor", file="pvals-autocor.dat")

#evaluation for the paper

png("../img/hist-autocor.png", height=600, width=3600, res=72*4) #standard resolution is 72
par(mar=c(4,4,1,0))
par(mfrow=c(1,5))
hist(pvals.autocor[,4], breaks=19, freq=FALSE, ylab="relative frequency", xlab="bootstrap p-value", main="n = 100")
hist(pvals.autocor[,10], breaks=19, freq=FALSE, ylab="", xlab="", main="n = 500")
hist(pvals.autocor[,16], breaks=19, freq=FALSE, ylab="", xlab="", main="n = 1000")
hist(pvals.autocor[,22], breaks=19, freq=FALSE, ylab="", xlab="", main="n = 5000")
hist(pvals.autocor[,28], breaks=19, freq=FALSE, ylab="", xlab="", main="n = 10000")
par(mfrow=c(1,1))
dev.off()

sizpow = matrix(0, nrow=6, ncol=5)
for(i in 1:5){
  sizpow[1,i] = mean(pvals.autocor[,i*6-2]<0.1)  #size,  nominal 10%
  sizpow[2,i] = mean(pvals.autocor[,i*6-1]<0.1)  #power, nominal 10%
  sizpow[3,i] = mean(pvals.autocor[,i*6]<0.1)    #power, nominal 10%
  sizpow[4,i] = mean(pvals.autocor[,i*6-2]<0.05) #size,  nominal 5%
  sizpow[5,i] = mean(pvals.autocor[,i*6-1]<0.05) #power, nominal 5%
  sizpow[6,i] = mean(pvals.autocor[,i*6]<0.05)   #power, nominal 5%
}
require(xtable)
table = xtable(sizpow)


