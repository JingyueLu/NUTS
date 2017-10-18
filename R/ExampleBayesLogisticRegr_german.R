#' Bayesian Logistic Regression on German credit data (Hoffman & Gelman, 2011)
#'
#' This script contains the logarithm of the posterior distribution of the Bayesian logistic regression model fit to
#' the German credit data set (1000 observations) available online in the UCI repository
#' (https://archive.ics.uci.edu/ml/machine-learning-databases/statlog/german/). The gradient of the logarithm
#' of the posterior distribution is also reported. Following Hoffman & Gelman (2011),  all predictors are normalized to
#' have zero mean and unit variance. Each element of the vector of parameters \theta=c(\alpha,\beta) is drawn
#' indipendently from a normal distribution with zero mean and variance \sigma^2=100.

german <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/statlog/german/german.data-numeric",sep = "")
y <- german[, 25]
X <- scale(german[,1:24])
sigma.sq <- 100   ### because sigma.sq <- 100
alpha.start <- rnorm(1, mean = 0, sd = 10)
beta.start <- rnorm(24, 0, sd = 10)

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

L <- function(theta){
  return(list(LogTarget(theta),GradLogTarget(theta)))
}
