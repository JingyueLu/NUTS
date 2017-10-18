#' Bayesian Logistic Regression on German credit data (Hoffman & Gelman, 2011)
#'
#' This script contains the logarithm of the posterior distribution of the Bayesian logistic regression model fit to
#' the German credit data set (1000 observations) available online in the UCI repository. The gradient of the logarithm
#' of the posterior distribution is also reported. Following Hoffman & Gelman (2011),  all predictors are normalized to
#' have zero mean #' and unit variance. Each element of the vector of parameters \theta=c(\alpha,\beta) is drawn
#' indipendently from a normal distribution with zero mean and variance \sigma^2=100.

german <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/statlog/german/german.data-numeric",sep = "")
y <- german[, 25]
X <- scale(german[,1:24])
sigma.sq <- 100   ### because sigma.sq <- 100
alpha.start <- rnorm(1, mean = 0, sd = 10)
beta.start <- rnorm(24, 0, sd = 10)

# Target <- function(theta){
#   X <- X
#   y <- y
#   sigma.sq <- 100
#   exp( - sum(log(1 + exp(- y * (theta[1] + X%*%theta[-1])))) - (theta[1]^2 - crossprod(theta[-1])) / (2 * 100))
# }

# LogTarget2 <- function(theta){
#   X <- X
#   y <- y
#   sigma.sq <- 100
#   - sum(log(1 + exp(- y * (theta[1] + X%*%theta[-1])))) - (theta[1]^2 - crossprod(theta[-1])) / (2 * 100)
# }

LogTarget <- function(theta){
  X.star <- cbind(1,X)
  y <- y
  sigma.sq <- 100
  temp <- X.star*y
  return(as.numeric(- sum(log(1 + exp(- temp %*% theta))) - crossprod(theta) / (2 * 100)))
}

GradLogTarget <- function(theta){
  X.star <- cbind(1,X)
  y <- y
  sigma.sq <- 100
  temp <- X.star*y
  return(as.vector(t(temp) %*% ((exp(- temp %*% theta)) / (1 + exp(- temp %*% theta))) - theta / sigma.sq))
}

L<-function(theta){
  X.star <- cbind(1,X)
  y <- y
  sigma.sq <- 100
  temp <- X.star*y
  dlog <-as.numeric(- sum(log(1 + exp(- temp %*% theta))) - crossprod(theta) / (2 * 100))
  gdlog <-as.vector(t(temp) %*% ((exp(- temp %*% theta)) / (1 + exp(- temp %*% theta))) - theta / sigma.sq)
  return(list(dlog, gdlog))
}
M = 1000
Madapt = 500
delta = 0.6
theta0 <- c(alpha.start, beta.start)
out <- NutsDual(theta0, 0.6, L, M,Madapt)

#all.equal(pracma::grad(LogTarget,aaa),GradLogTarget(aaa)) #Check that the function returns effectively the gradient

# theta.hmc<-hmc(theta.start = c(alpha.start, beta.start), epsilon = 0.06, L = 10, logDensity = LogTarget, M = 10000)
# ts.plot(theta.hmc[,1])
# ts.plot(theta.hmc[,2])
#
# # Plot a histogram of the first variable, with true density
# choose.var<-14
#
# hist(theta.hmc[, choose.var],freq=FALSE,breaks=50)
# x = seq(range(theta.hmc[,choose.var])[1],range(theta.hmc[,choose.var])[2],by=0.01)
# vec.theta<-rep(0,25)
#
#
# target.eval<-function(choose.var,x){
#   vec.theta<-rnorm(25,sd = 10)
#   dens<-numeric(length(x))
#   for (i in 1:length(x)) {
#     vec.theta[choose.var] <- x[i]
#     dens[i] <- target(vec.theta)
#   }
#   return(dens)
# }
