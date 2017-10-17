#' Hamiltonian Monte Carlo with Dual Averaging and random trajectory length
#'
#' The function is modified version of Hamiltonian Monte Carlo with Dual Averaging. Based on the post by Radford Neal,
#' Hamiltonian dynamics can be periodic for variables. It is possible for the period for a variabl to match its trajectory
#' length, resulting the ending point close to the starting point. To address this issue, Neal proposes to randomly vary the
#' trajectory length over some moderate range. Therefore, one modification is made:
#' the previous code: Lm<-max(1,round(lambda/epsilon)) is changed to
#' the new code: Lm<-max(1,round((0.9 + runif(1)/5)*lambda/epsilon))
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
#'
HmcDualMod <- function(theta0, delta, lambda, L, M, Madapt) {
  len <- length(theta0)
  outcome <-  matrix(theta0 , nrow = Madapt + M, ncol = len , byrow = T)
  alpha.vec <- rep(NA,Madapt + M)

  # Find a reasonable initial value of epsilon using heuristic
  g(log.temp, grad0) %=% L(theta0)
  epsilon =  FindReasonableEpsilon(theta0, log.temp, L)

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
    proposal <- list(theta.tilde = outcome[m-1, ], r.tilde = r0)
    Lm<-max(1,round((0.9 + runif(1)/5)*lambda/epsilon))
    for (i in 1:Lm)  {
      g(proposal$theta.tilde, proposal$r.tilde, trivial)%=%Leapfrog(proposal$theta.tilde, proposal$r.tilde, epsilon, L)
    }
    alpha.vec[m] <- min(1,exp(L(proposal$theta.tilde)[[1]] - 1/2 * crossprod(proposal$r.tilde) - L(outcome[m-1, ])[[1]] + 1/2 * crossprod(r0)))
    if (runif(1) <= alpha.vec[m]){
      outcome[m,]<-proposal$theta.tilde
    }

    ###Dual Averaging
    if (m <= Madapt) {
      temp <- 1/(m-1+t0)
      H.bar <- (1 - temp)*H.bar + temp*(delta - alpha.vec[m])
      epsilon <- exp(mu - sqrt(m-1)/gamma*H.bar)
      temp <- (m-1)^(-kappa)
      epsilon.bar <- exp(temp*log(epsilon)+(1-temp)*log(epsilon.bar))
    } else{
      epsilon <- epsilon.bar
    }
    #print(m)
  }
  return(list(samples = outcome[(Madapt + 1):(M + Madapt),],alpha = alpha.vec))
}
