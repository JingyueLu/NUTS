### Replicating the example from the post

#' analysis
#'
#' @param samples the sample points computed using some sampling method
#' @param B the number of batches
#' @param N the number of iterations within each batch
#'
#' @return a list: the first two elements are the mean and standard deviaiton of the sample points the last element is a matrix contains the variance of the batch means and batch variances
analysis <- function(samples, B,N){
  mu <- apply(samples,2, mean)
  stddev <- apply(samples,2,sd)

  ## Initialise a matrix that stores variances of means of batches (first column)
  ## and variances of variances of batches(second column)
  SubMatrix <- matrix(0, nrow=30, ncol=2 )
  for (j in 1:length(theta0)){
    meanBatch <- c()
    varBatch <- c()
    for (i in 1:B){
      samples[,j]
      a <-samples[((i-1)*N+1):((i-1)*N+N), j]
      meanBatch <- c(meanBatch, mean(a))
      varBatch <- c(varBatch, var(a))
    }
    SubMatrix[j,] <- c(var(meanBatch), var(varBatch))
  }
  return(list(mu, stddev, SubMatrix))
}


### Constructing a 30-dimensional multivariate normal distribution, in which the 30 components
### are independent, with standard deviations of 110, 100, 1.1, 1.0, and 26 equally space between 8 and 16
### Since the operation of HMC is invariant to rotation, behaviour on the distribution is the same as on any
### 30-dimensional multivariate normal for which the square roots of the eigenvalues of the covariance matrix
### have these values

sdevs <- c (110, 100, seq(16,8,length=26), 1.1, 1.0)
V <- sdevs^2

# Number of batchs
B = 20
# Number of iterations in a batch
N = 625

M = B*N
Madapt = 2000
delta =0.6

####################################################################################

### Construct sample points by NUTS

# L returns log posterior density and its gradient
L <- function(theta){
  grad0 = - theta/V
  log0 = 0.5*t(theta) %*%grad0

  return(list(log0, grad0))
}

theta0 <- rnorm(30)

g(sampleNuts, acceprate, epsilonNor, height) %=% NutsDual(theta0, delta, L, M, Madapt)

### Construct independent sample points using build-in normal functions
sampleInd = matrix(0, nrow = M, ncol = length(theta0));
for (i in 1:M){
  sampleInd[i,] = rnorm(length(theta0)) * sdevs
}


### Construct sample point by HMC

## trajectory lengths is 100
HMD1st <- Sys.time()
g(sampleHMD1, alpha1) %=% HmcDual(theta0, delta, 100, L, M, Madapt)
HMD1ed <- Sys.time()
HMD1time<- difftime(HMD1ed,HMD1st)
## trajectory lengths is 170
HMD2st <- Sys.time()
g(sampleHMD2, alpha2) %=% HmcDual(theta0, delta, 170, L, M, Madapt)
HMD2ed <- Sys.time()
HMD2time<- difftime(HMD2ed,HMD2st)
## trajectory lengths is 200
HMD3st <- Sys.time()
g(sampleHMD3, alpha3) %=% HmcDual(theta0, delta, 200, L, M, Madapt)
HMD3ed <- Sys.time()
HMD3time<- difftime(HMD3ed,HMD3st)

### Construct sample point by HMC Modified

## trajectory lengths is 100
HMD1stM <- Sys.time()
g(sampleHMD1mod, alpha1mod) %=% HmcDualMod(theta0, delta, 100, L, M, Madapt)
HMD1edM <- Sys.time()
HMD1timeMod<- difftime(HMD1edM,HMD1stM)
## trajectory lengths is 170
HMD2stM <- Sys.time()
g(sampleHMD2mod, alpha2mod) %=% HmcDualMod(theta0, delta, 170, L, M, Madapt)
HMD2edM <- Sys.time()
HMD2timeMod<- difftime(HMD2edM,HMD2stM)
## trajectory lengths is 200
HMD3stM <- Sys.time()
g(sampleHMD3mod, alpha3mod) %=% HmcDualMod(theta0, delta, 200, L, M, Madapt)
HMD3edM <- Sys.time()
HMD3timeMod<- difftime(HMD3edM,HMD3stM)

######################################################################
###Analysing Results
######################################################################

############# Nuts
g(muNuts, sdNuts, SubMatrixNuts) %=% analysis(sampleNuts, B,N)

############# Independent
g(muInd, sdInd, SubMatrixInd) %=% analysis(sampleInd, B,N)

############# HMC with Dual Averaging
## HMD with trajectory 100
g(muHMD1, sdHMD1, SubMatrixHMD1) %=% analysis(sampleHMD1, B,N)
## HMD with trajectory 170
g(muHMD2, sdHMD2, SubMatrixHMD2) %=% analysis(sampleHMD2, B,N)
## ## HMD with trajectory 200
g(muHMD3, sdHMD3, SubMatrixHMD3) %=% analysis(sampleHMD3, B,N)

############## HMC with Dual Averaging and random trajectory
## HMD with trajectory 100
g(muHMD1mod, sdHMD1mod, SubMatrixHMD1mod) %=% analysis(sampleHMD1mod, B,N)
## HMD with trajectory 170
g(muHMD2mod, sdHMD2mod, SubMatrixHMD2mod) %=% analysis(sampleHMD2mod, B,N)
## ## HMD with trajectory 200
g(muHMD3mod, sdHMD3mod, SubMatrixHMD3mod) %=% analysis(sampleHMD3mod, B,N)


######################################################################
### Plot the results
######################################################################

x <- seq(1,30,by=1)


########### NUTS & INDEPENDENT

# variance of mean estimate (NUTS & Independent)
plot(x, SubMatrixNuts[,1],log="y", ylim = c(1e-03,1e+03), xlab ="NUTS(red) / Independent(black)", ylab="variance of mean estimate", col="red",pch=16,type="b")
points(x,SubMatrixInd[,1], log="y", col="black", pch=16,type="b")

# variance of variance estimate (NUTS & Independent)
plot(x, SubMatrixNuts[,2], ylim = c(1e-03,1e+07),log="y", xlab ="NUTS(red) / Independent(black)", ylab="variance of variance estimate", col="red",pch=16, type="b")
points(x,SubMatrixInd[,2], log="y", col="black", pch=16, type="b")

########### HMC Dual (trajectory = 100, 170, 200)

## variance of mean estimate ()
plot(x, SubMatrixHMD1[,1],log="y", type="b",xlab ="trajectory length 100(red); 170(green); 200(blue)", ylab="variance of mean estimate", col="red",pch=16)
points(x,SubMatrixHMD2[,1], log="y", type="b",col="green", pch=16)
points(x,SubMatrixHMD3[,1], log="y",type="b", col="blue", pch=16)

## variance of variance estimate (NUTS & Independent)
plot(x, SubMatrixHMD1[,2],log="y", type="b",ylim = c(1e-03,1e+07),xlab ="trajectory length 100(red); 170(green); 200(blue)", ylab="variance of mean estimate", col="red",pch=16)
points(x,SubMatrixHMD2[,2], log="y",type="b", col="green", pch=16)
points(x,SubMatrixHMD3[,2], log="y", type="b",col="blue", pch=16)

########### HMC Dual with random trajectory (trajectory = 100, 170, 200)
## variance of mean estimate ()
plot(x, SubMatrixHMD1mod[,1],log="y",ylim = c(1e-03,1e+03), type="b",xlab ="trajectory length 100(red); 170(green); 200(blue)", ylab="variance of mean estimate", col="red",pch=16)
points(x,SubMatrixHMD2mod[,1], log="y",type="b", col="green", pch=16)
points(x,SubMatrixHMD3mod[,1], log="y", type="b",col="blue", pch=16)

## variance of variance estimate (NUTS & Independent)
plot(x, SubMatrixHMD1mod[,2],log="y",ylim = c(1e-03,1e+07), type="b",xlab ="trajectory length 100(red); 170(green); 200(blue)", ylab="variance of mean estimate", col="red",pch=16)
points(x,SubMatrixHMD2mod[,2], log="y",type="b", col="green", pch=16)
points(x,SubMatrixHMD3mod[,2], log="y", type="b",col="blue", pch=16)
