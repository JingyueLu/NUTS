#' ExampleDoubleExponential

### log density function for double exponential
logDens <- function(x)   -sum(abs(x))
dlogDens <- function(x)   -sign(x)

L <- function(x) {
  return(list(- sum(abs(x)), - sign(x)))
}

### For Hmc plot the MCMC chains and densities
samples.hmc <- Hmc(theta0 = c(1,1), epsilon = 0.2, leap.nsteps = 10, L = L, M = 2000)
ts.plot(samples.hmc[,1])
ts.plot(samples.hmc[,2])
# Plot a histogram of the second variable, with true density
hist(samples.hmc[,2],freq=FALSE,breaks=50)
x = seq(-5,5,len=100)
lines(x,0.5*dexp(abs(x)),col="red") #true density in red
# Autocorrelation function
acf(samples.hmc)


### For HmcDual plot the MCMC chains and densities
samples.hmcdual <- HmcDual(theta0=c(1,1), delta = 0.65, lambda = 1.5, L = L.func, M = 2000, Madapt = 1000)
ts.plot(samples.hmcdual$samples[,1])
ts.plot(samples.hmcdual$samples[,2])
# Plot a histogram of the second variable, with true density
hist(samples.hmcdual$samples[,2],freq=FALSE,breaks=50)
x = seq(-5,5,len=100)
lines(x,0.5*dexp(abs(x)),col="red") #true density in red
# Autocorrelation function
acf(samples.hmcdual$samples)



### For NUTS plot the MCMC chains and densities
samples.nuts <- NutsDual(theta0 = c(1,1), delta = 0.6, L = L, M = 2000, Madapt = 500)
ts.plot(samples.nuts$samples[,1])
ts.plot(samples.nuts$samples[,2])
# Plot a histogram of the second variable, with true density
hist(samples.nuts$samples[,2],freq=FALSE,breaks=50)
x = seq(-5,5,len=100)
lines(x,0.5*dexp(abs(x)),col="red") #true density in red
# Autocorrelation function
acf(samples.nuts$samples)

################################################################
################################################################
################################################################

#' Comparing the performances of this version of HMC with respect to HybridMC::hybridMC
#'
#' This script contains the example from the function hybridMC (in the package HybridMC), that calls C code.
#' At the end an overview of the performances of the two functions is presented.


microbenchmark(hybridMC(y.start = c(1,1), n.samp=2000, logDens=logDens, dLogDens=dlogDens, epsilon=.2, LFsteps=10),
               Hmc(theta0 = c(1,1), epsilon = 0.2, leap.nsteps = 10, L = L, M = 2000))

