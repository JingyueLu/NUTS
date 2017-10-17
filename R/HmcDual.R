#' %=%
#'
#' The function is copied from "https://strugglingthroughproblems.wordpress.com/2010/08/27/matlab-style-multiple-assignment-in%C2%A0r/" and is used for multiple assignment
#' @param l		the list on the left-hand side
#' @param r 	the list on the right-hand side
#'
#' @example
#' g(a,b) = list(c(3,5,6),4)
#'
'%=%' = function(l, r, ...) UseMethod('%=%')

#' %=%,lbunch
#'
#' The function is copied from "https://strugglingthroughproblems.wordpress.com/2010/08/27/matlab-style-multiple-assignment-in%C2%A0r/"
#' It is a specific implementation for the generic %=%
#' @param l		the list on the left-hand side
#' @param r 	the list on the right-hand side
#'
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)

  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")

  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }

  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}


#' g
#'
#' The function is copied from "https://strugglingthroughproblems.wordpress.com/2010/08/27/matlab-style-multiple-assignment-in%C2%A0r/".
#' This function grabs '...' part without evaluating it (so pass by name)
#'
#'
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}


#' Leapfrog
#'
#' This function perform a leapfrog step. This function is a modified version of Leapfrog in the paper.
#' It returns etra values: log posterior value and gradient of log posterior at the new position theta.tilde
#'
#' @param theta		starting position
#' @param r 	starting momentum
#' @param epsilon 	step size
#' @param L 	callable function: returns the value of log posterior and the gradient of log posterior probability at given input
#'
#' @return the list of updated theta, r and the log posterior value at the updated point


# Leapfrog222 <- function(theta, r, epsilon, logDens, dLogDensity) {
#   r.tilde <- r + dLogDensity(theta) * epsilon / 2
#   theta.tilde <- theta + epsilon * r.tilde
#   r.tilde <- r.tilde + dLogDensity(theta.tilde) * epsilon / 2
#   return(list(theta.tilde = theta.tilde, r.tilde = r.tilde, log.tilde = logDens(theta.tilde) ))
# }




Leapfrog <- function(theta, r, epsilon, L){
  g(trivial,grad.theta) %=% L(theta)
  r.tilde <- r+ 0.5*epsilon*grad.theta
  theta.tilde <- theta + epsilon*r.tilde
  g(log.tilde, grad.tilde) %=% L(theta.tilde)
  r.tilde <- r.tilde + 0.5*epsilon*grad.tilde
  return(list(theta.tilde, r.tilde, log.tilde))
}

####################################################################################
####################################################################################
####################################################################################

#' FindReasonableEpsilon
#'
#' Heuristic for choosing an initial value of epsilon
#'
#' @param theta
#' @param log.start the log posterior value at initial state
#' @param grad.start the gradient value at initial state
#' @param L callable function needed in Leapfrog
#'
#' @return initial epsilon
#'
#'

FindReasonableEpsilon <- function(theta,log.start, L){
  epsilon <- 1
  r = rnorm(length(theta))
  g(theta.prime, r.prime, log.prime) %=% Leapfrog(theta,r,epsilon, L)

  #Computing the ratio of p(theta.prime, r.prime) over p(theta, r)
  tempratio <- exp(log.prime-0.5*(crossprod(r.prime, r.prime)) - log.start + 0.5*(crossprod(r, r)))
  a <- 2* as.numeric(tempratio > 0.5)-1
  #print(c(tempratio, a))

  while (tempratio^a > 2^(-a) ){
    epsilon <-  epsilon * (2^(a))
    g(theta.prime, r.prime, log.prime) %=% Leapfrog(theta,r,epsilon, L)
    tempratio <- exp(log.prime - 0.5*(crossprod(r.prime, r.prime)) - log.start + 0.5*(crossprod(r, r)))
    #print(tempratio)
  }
  return(epsilon)

}


####################################################################################
####################################################################################
####################################################################################

#' Hamiltonian Monte Carlo with Dual Averaging
#'
#' @param theta0 a p-dimensional vector with the initial value for each parameter.
#' @param delta 	the desired average acceptance probability
#' @param lambda trajectory length for leapfrog
#' @param L a callable function that returns a 2-dimensional vector whose first element is the
#' logarithm of the density of the input and the second is the gradient of the logarithm
#' of the density evaluated at the input.
#' @param M an integer specifying the number of iterations.
#' @param Madapt an integer specifying the number of iterations of the warmup phase.
#'
#' @return This function returns a matrix which m-th row is a sample from the joint density.

HmcDual <- function(theta0, delta, lambda, L, M, Madapt) {
  len <- length(theta0)
  outcome <-  matrix(theta0 , nrow = Madapt + M, ncol = len , byrow = T)

  # Find a reasonable initial value of epsilon using heuristic
  g(log.temp, grad0) %=% L(theta0)
  epsilon =  FindReasonableEpsilon(theta0, log.temp, L)

  #L.func<-function(x) return(c(logDensity(x), dLogDensity(x)))
  #epsilon <-  FindReasonableEpsilon(theta0, logDensity(theta0), L=L.func)

  # Setting up other parameters for dual averaging
  mu <- log(10*epsilon)
  epsilon.bar <- 1
  H.bar <- 0
  gamma <- 0.05
  t0 <- 10
  kappa <- 0.75

  for (m in 2:(Madapt+M)) {
    r0 <- rnorm(len,  0, sd = 1)
    outcome[m, ] <- outcome[m-1, ]
    #theta.tilde <- outcome[m-1, ]
    proposal <- list(theta.tilde = outcome[m-1, ], r.tilde = r0)
    Lm<-max(1,round(lambda/epsilon))

    for (i in 1:Lm)  {
      #proposal <- leapfrog(proposal$theta.tilde, proposal$r.tilde, epsilon, dLogDensity)
      g(proposal$theta.tilde, proposal$r.tilde, trivial)%=%Leapfrog(outcome[m-1, ], r0, epsilon, L)
    }
    #alpha <- min(1, within(proposal,exp(L(theta.tilde)[[1]] - 1/2 * crossprod(r.tilde) - L(outcome[m-1, ])[[1]] + 1/2 * crossprod(r0))))
    alpha <- min(1,exp(L(proposal$theta.tilde)[[1]] - 1/2 * crossprod(proposal$r.tilde) - L(outcome[m-1, ])[[1]] + 1/2 * crossprod(r0)))
    if (runif(1) <= alpha){
      #outcome[m,] <- proposal$theta.tilde
      outcome[m,]<-proposal$theta.tilde
    }

    ###Dual Averaging
    if (m <= Madapt) {
      temp <- 1/(m+t0)
      H.bar <- (1 - temp)*H.bar + temp*(delta - alpha)
      epsilon <- exp(mu - sqrt(m)/gamma*H.bar)
      temp <- m^{-kappa}
      epsilon.bar <- exp(temp*log(epsilon)+(1-temp)*log(epsilon.bar))
    }
    else{
      epsilon <- epsilon.bar
    }

  }
  return(outcome[(Madapt+1):(M+Madapt),])
}
