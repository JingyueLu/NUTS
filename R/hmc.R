# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Hamiltonian Monte Carlo
#'
#' @param theta.start Initial value.
#' @param epsilon Leapfrog step size.
#' @param L Length of the leapfrog step.
#' @param logDensity Logarithm of the joint density of the variables theta.
#' @param M Number of iterations.
#'
#' @return This function returns an object of type ... .  Each row is a sample from the joint density.
#'
#'

hmc.nograd <- function(theta.start, epsilon, L, logDensity, dlogDensity, M) {

  leapfrog <- function(theta, r, epsilon, dlogDensity) {
    r.tilde <- r + dlogDensity(theta) * epsilon / 2
    theta.tilde <- theta + epsilon * r.tilde
    r.tilde <- r.tilde + dlogDensity(theta.tilde) * epsilon / 2
    return(list(theta.tilde = theta.tilde, r.tilde = r.tilde))
  }

  #require(mvtnorm)
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

  require(mvtnorm)
  p <- length(theta.start)
  outcome <- list(theta = matrix(theta.start, nrow = M, ncol = p, byrow = T), r = matrix(NA, nrow = M, ncol = p))##, dimnames[[2]] = paste("Iteration ",1:M))
  for (m in 2:M) {
    r0 <- rmvnorm(1, mean = rep(0,p), sigma = diag(p))
    outcome$theta[m, ] <- outcome$theta[m-1, ]
    proposal <- list(theta.tilde = outcome$theta[m-1, ], r.tilde = r0)
    for (i in 1:L)  proposal <- leapfrog(proposal$theta.tilde, proposal$r.tilde, epsilon, logDensity)[1:2]
    alpha <- min(1, with(proposal,exp(logDensity(theta.tilde) - 1/2 * crossprod(r.tilde) - logDensity(outcome$theta[m-1, ]) - 1/2 * crossprod(r0))))
    if (runif(1) <= alpha){
      outcome$theta[m,] <- proposal$theta.tilde
      outcome$r[m, ] <- - proposal$r.tilde
      }
    }
  return(outcome$theta)
}

