#' Comparing the performances of this version of HMC with respect to HybridMC::hybridMC
#'
#' This script contains the example from the function hybridMC (in the package HybridMC), that calls C code.
#' At the end an overview of the performances of the two functions is presented.

library(HybridMC)
# ## log density function for double exponential
logDens <- function(x)   -sum(abs(x))
dlogDens <- function(x)   -sign(x)
L.func <- function(x){
  return(list( - sum(abs(x)),- sign(x)))
}

L <- function(x) {
  return(list(-sum(abs(x)), -sign(x)))
}

samples.hyb <- hybridMC(y.start=c(1,1), n.samp=2000, logDens=logDens, dLogDens=dlogDens, epsilon=.2, LFsteps=10)
samples.test <- hmc.nograd(theta.start = c(1,1), epsilon = 0.2, L = 10, logDensity = logDens, dlogDensity=dlogDens, M = 2000)
<<<<<<< HEAD
samples.nuts <- NutsDual(c(1,1),0.6, L, 2000, 500)
=======
samples.hmcdual <- HmcDual(theta0=c(1,1), delta = 0.65, lambda = 1.5, L = L.func, M = 2000, Madapt = 1000)
>>>>>>> 23a65b9a11cc6fe54823b618e66e83547d8a6902

# For NUTS::hmc.nograd plot the MCMC chains and densities
ts.plot(samples.test[,1])
ts.plot(samples.test[,2])

ts.plot(samples.nuts[,1])
ts.plot(samples.nuts[,2])

# Plot a histogram of the first variable, with true density
hist(samples.test[,2],freq=FALSE,breaks=50)
x = seq(-5,5,len=100)
lines(x,0.5*dexp(abs(x)),col="red") #true density in red
# Autocorrelation function
acf(samples.test)

hist(samples.nuts[,2],freq=FALSE,breaks=50)
x = seq(-5,5,len=100)
lines(x,0.5*dexp(abs(x)),col="red") #true density in red
# Autocorrelation function
acf(samples.nuts)


microbenchmark(hybridMC(y.start = c(1,2,3,4), n.samp=2000, logDens=logDens, dLogDens=dlogDens, epsilon=.2, LFsteps=10),
               #hmc.nograd(theta.start = c(1,2,3,4), M = 2000, logDensity = logDens, dlogDensity=dlogDens, epsilon = 0.2, L = 10),
               HmcDual(c(1,2,3,4), delta = 0.65, lambda = 1.5, L = L.func, M = 2000, Madapt = 1000))

# For HmcDual
ts.plot(samples.hmcdual[,1])
ts.plot(samples.hmcdual[,2])

# Plot a histogram of the first variable, with true density
hist(samples.hmcdual[,2],freq=FALSE,breaks=50)
x = seq(-5,5,len=100)
lines(x,0.5*dexp(abs(x)),col="red") #true density in red
# Autocorrelation function
acf(samples.hmcdual)
