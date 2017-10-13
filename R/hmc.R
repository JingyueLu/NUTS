#' Hamiltonian Monte Carlo
#'
#' @param theta.start a p-dimensional vector with the initial value for each parameter.
#' @param epsilon a non-negative value specifying the leapfrog step size.
#' @param L an integer number specifying the length of the leapfrog step.
#' @param logDensity a function that computes the logarithm of the joint density of the variables.
#' @param M an integer specifying the number of iterations.
#'
#' @return This function returns a matrix which m-th row is a sample from the joint density.

hmc.nograd <- function(theta.start, epsilon, L, logDensity, dlogDensity, M) {

  leapfrog <- function(theta, r, epsilon, dlogDensity) {
    r.tilde <- r + dlogDensity(theta) * epsilon / 2
    theta.tilde <- theta + epsilon * r.tilde
    r.tilde <- r.tilde + dlogDensity(theta.tilde) * epsilon / 2
    return(list(theta.tilde = theta.tilde, r.tilde = r.tilde))
  }

  p <- length(theta.start)
  outcome <- list(theta = matrix(theta.start, nrow = M+1, ncol = p, byrow = T), r = matrix(NA, nrow = M+1, ncol = p))##, dimnames[[2]] = paste("Iteration ",1:M))
  for (m in 2:(M+1)) {
    r0 <- rnorm(p,  0, sd = 1)
    outcome$theta[m, ] <- outcome$theta[m-1, ]
    proposal <- list(theta.tilde = outcome$theta[m-1, ], r.tilde = r0)
    for (i in 1:L)  {
      proposal <- leapfrog(proposal$theta.tilde, proposal$r.tilde, epsilon, dlogDensity)
    }
    alpha <- min(1, with(proposal,exp(logDensity(theta.tilde) - 1/2 * crossprod(r.tilde) - logDensity(outcome$theta[m-1, ]) + 1/2 * crossprod(r0))))
    if (runif(1) <= alpha){
      outcome$theta[m,] <- proposal$theta.tilde
      outcome$r[m, ] <- - proposal$r.tilde
    }
  }
  return(outcome$theta[-1,])
}






hmc <- function(theta.start, epsilon, L=5, logDensity, M) {

  leapfrog <- function(theta, r, epsilon, logDensity) {
    require(pracma)
    r.tilde <- r + grad(logDensity,theta) * epsilon / 2
    theta.tilde <- theta + epsilon * r.tilde
    r.tilde <- r.tilde + grad(logDensity,theta.tilde) * epsilon / 2
    return(list(theta.tilde = theta.tilde, r.tilde = r.tilde, grad.logDens.tilde = grad(logDensity,theta.tilde)))
  }

  p <- length(theta.start)
  outcome <- list(theta = matrix(theta.start, nrow = M+1, ncol = p, byrow = T), r = matrix(NA, nrow = M+1, ncol = p))##, dimnames[[2]] = paste("Iteration ",1:M))
  for (m in 2:(M+1)) {
    r0 <- rnorm(p,  0, sd = 1)
    outcome$theta[m, ] <- outcome$theta[m-1, ]
    proposal <- list(theta.tilde = outcome$theta[m-1, ], r.tilde = r0)
    for (i in 1:L)  {
      proposal <- leapfrog(proposal$theta.tilde, proposal$r.tilde, epsilon, logDensity)[1:2]
    }
    alpha <- min(1, with(proposal,exp(logDensity(theta.tilde) - 1/2 * crossprod(r.tilde) - logDensity(outcome$theta[m-1, ]) + 1/2 * crossprod(r0))))
    if (runif(1) <= alpha){
      outcome$theta[m,] <- proposal$theta.tilde
      outcome$r[m, ] <- - proposal$r.tilde
    }
  }
  return(outcome$theta[-1,])
}

