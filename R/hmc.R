#' %=%
#'
#' The function is copied from "https://strugglingthroughproblems.wordpress.com/2010/08/27/matlab-style-multiple-assignment-in%C2%A0r/" and is used for multiple assignment
#' @param l		the list on the left-hand side
#' @param r 	the list on the right-hand side
#'
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
#' This function performs a leapfrog step. This function is a modified version of Leapfrog in the paper.
#' It returns extra values: log posterior value and gradient of log posterior at the new position theta.tilde
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

######################################################################
######################################################################
######################################################################

#' Hmc
#'
#' This function performs the Hamiltonian Monte Carlo.
#'
#' @param theta0 a p-dimensional vector with the initial value for each parameter.
#' @param epsilon a non-negative value specifying the leapfrog step size.
#' @param leap.nsteps an integer number specifying the number of leapfrog steps.
#' @param L callable function: returns the value of log posterior and the gradient of log posterior probability at given input.
#' @param M an integer specifying the number of iterations.
#'
#' @return This function returns a matrix whose m-th row is a sample from the joint density.

Hmc <- function(theta0, epsilon, leap.nsteps, L, M) {
  len <- length(theta0)
  outcome <- matrix(theta0, nrow = M+1, ncol = len, byrow = T)
  for (m in 2:(M+1)) {
    r0 <- rnorm(len,  0, sd = 1)
    outcome[m, ] <- outcome[m-1, ]
    proposal <- list(theta.tilde = outcome[m-1, ], r.tilde = r0)
    for (i in 1:leap.nsteps)  {
      g(proposal$theta.tilde, proposal$r.tilde, trivial) %=% Leapfrog(proposal$theta.tilde, proposal$r.tilde, epsilon, L)
    }
    alpha <- min(1, exp(L(proposal$theta.tilde)[[1]] - 1/2 * crossprod(proposal$r.tilde) - L(outcome[m-1, ])[[1]] + 1/2 * crossprod(r0)))
    if (runif(1) <= alpha){
      outcome[m,] <- proposal$theta.tilde
    }
  }
  HmcResult <- list(samples=outcome[-1,])
  class(HmcResult ) <- "Hmc"
  return(HmcResult)
}

print.Hmc <- function(obj){
  D <- length(obj[[1]][1,])
  m <- length(obj[[1]][,1])
  cat("Dimension of the parameter:", D, "\n")
  cat("Sample size generated: M = ", m, "\n")
  len <- min(m,6)
  dots <- ifelse(m>6, "...","\n")
  for (i in 1:D){
    cat("theta",i,": ", format(round(obj[[1]][(1:len),i], digits = 6), nsmall=6), dots, "\n")
  }
  invisible(obj)
}







