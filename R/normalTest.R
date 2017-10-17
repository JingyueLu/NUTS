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
## trajectory lengths is 170
## trajectory lengths is 200

g(muNuts, sdNuts, SubMatrixNuts) %=% analysis(sampleNuts, B,N)
g(muInd, sdInd, SubMatrixInd) %=% analysis(sampleInd, B,N)

### Plot the results
x <- seq(1,30,by=1)

## variance of mean estimate (NUTS & Independent)
plot(x, SubMatrixNuts[,1],log="y", xlab ="NUTS(red) / Independent(black)", ylab="variance of mean estimate", col="red",pch=16)
points(x,SubMatrixInd[,1], log="y", col="black", pch=16)

## variance of variance estimate (NUTS & Independent)
plot(x, SubMatrixNuts[,2],log="y", xlab ="NUTS(red) / Independent(black)", ylab="variance of variance estimate", col="red",pch=16)
points(x,SubMatrixInd[,2], log="y", col="black", pch=16)



