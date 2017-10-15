# GIVE REFERENCES
# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
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


# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}


#' Leapfrog
#'
#' @param theta		starting position
#' @param r 	starting momentum
#' @param epsilon 	step size
#' @param L 	callable function: returns the value of log posterior and the gradient of log posterior prbability at given input
#'
#' @output theta.tilde
#'
#' This function perform a leapfrog step. This function is a modified version of Leapfrog in the paper. It returns etra values: log posterior value and gradient of log posterior at the new position theta.tilde
#'
Leapfrog <- function(theta, r, epsilon, L){
  g(trivial,grad.theta) %=% L(theta)
  #grad.theta is a vector
  r.tilde <- r+ 0.5*epsilon*grad.theta   #r.tilde is a vector
  theta.tilde <- theta + epsilon*r.tilde  #theta.tilde is a vector
  g(log.tilde, grad.tilde) %=% L(theta.tilde)
  r.tilde <- r.tilde + 0.5*epsilon*grad.tilde

  return(list(theta.tilde, r.tilde, log.tilde))

}


#' Build trees
#'
#' @param theta
#' @param r
#' @param u
#' @param v
#' @param j the height of a tree
#' @param epsilon step size
#' @param joint0 	to avoid computing the same value repeatedly, we compute joint0=log0- 0.5*(r0 %*% r0) beforehand and take this value as an input for the BUildTree function
#'
#' This function builds the tree for NUTS
#'
BuildTree <- function(theta, r,u,v, j, epsilon, L, joint0) {
  if (j==0){
    # Base case -- take one leapfrog step in the direction v
    g(theta.prime, r.prime, log.prime) %=% Leapfrog(theta, r, v*epsilon, L)
    temp <-log.prime - 0.5*(crossprod(r.prime, r.prime))
    # Checking whether the newly visted state is in the slice
    n.prime <- as.numeric(I(log(u)<temp))
    # Checking whether the simulation is moderately accurate
    s.prime <- as.numeric(I((log(u)-1000)< temp))

    alpha.prime <- min(1, exp(temp-joint0))
    n.alphaprime <- 1
    return(list(theta.prime,r.prime,theta.prime,r.prime,theta.prime, n.prime, s.prime, alpha.prime, n.alphaprime, log.prime ))
  }

  else {
    # Recursion -- implicitly build the left and right subtrees
    g(theta.minus, r.minus, theta.plus, r.plus, theta.prime, n.prime, s.prime, alpha.prime, n.alphaprime, log.prime) %=% BuildTree(theta, r,u,v,j-1, epsilon, L, joint0)
    # if s.prime =0, the stopping criteria were met
    if (s.prime==1){
      if(v == -1){
        g(theta.minus, r.minus, trivial, trivial, theta.prime2, n.prime2, s.prime2, alpha.prime2, n.alphaprime2, log.prime2) %=% BuildTree(theta.minus, r.minus, u,v,j-1, epsilon, L, joint0)
      }
      else{
        g(trivial, trivial, theta.plus, r.plus, theta.prime2, n.prime2, s.prime2, alpha.prime2, n.alphaprime2, log.prime2) %=% BuildTree(theta.plus, r.plus, u,v,j-1, epsilon)
      }
      # DO WE NEED TO CHECK WHETHER THE DENOMINATER IS GREATER THAN 0?
      if (runif(1)< n.prime2/(n.prime+n.prime2)){
        theta.prime = theta.prime2
        log.prime = log.prime2
      }
      n.prime <- n.prime + n.prime2
      s.prime <-  s.prime2*StopCon(theta.minus, theta.plus, r.minus,r.plus)

      alpha.prime <- alpha.prime + alpha.prime2
      n.alphaprime <- n.alphaprime + n.alphaprime2
    }
    return(list(theta.minus, r.minus, theta.plus, r.plus, theta.prime, n.prime, s.prime, alpha.prime, n.alphaprime, log.prime ))

  }
}

#' StopCon
#'
#' @param theta.minus
#' @param theta.plus
#' @param r.minus
#' @param r.plus
#'
#' This function computes the U-Turn stopping condition, which is used in NUTS and BuildTree functions
#'
StopCon <-function(theta.minus, theta.plus, r.minus,r.plus){
  theta.diff <- theta.plus - theta.minus
  temp1 <- as.numeric(crossprod(theta.diff,r.minus))
  temp2 <- as.numeric(crossprod(theta.diff, r.plus))
  return (temp1*temp2)
}


#' FindReasonableEpsilon
#'
#' @param theta
#' @param log.start the log posterior value at initial state
#' @param grad.start the gradient value at initial state
#' @param L callable function needed in Leapfrog
#'
#' @return epsilon
#'
#' Heuristic for choosing an initial value of epsilon
#'
FindReasonableEpsilon <- function(theta,log.start, grad.start, L){
  epsilon <- 0.5
  r = rnorm(length(theta))
  g(theta.prime, r.prime, log.prime) %=% Leapfrog(theta,r,epsilon, L)

  #Modifications may be needed. Initial epsilon is fairly large (Speed up possible?).

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


#' NutsDual
#'
#' @param theta0	inital state of theta
#' @param delta 	the desired average acceptance probability
#' @param L 	a callable function that returns log probability
#' @param M 	number of samples to generate
#' @param Madapt	number of iterations
#'
#' This function implements No-U-Turn Sampler with Dual Averaging.
#' The code is a modified version of algorithm6 in paper. Instead of including samples in burn-in period (m < Madapt), the function generates a sample where samples in the burn-in period are discarded.
#'
#' @return samples

NutsDual <- function(theta0, delta, L, M, Madapt){

  # Set up output structure: each theta is an array
  # samples is a matrix of size M x length(theta)
  out <- matrix(0 , nrow = Madapt + M, ncol = length(theta0) )
  out[1, ] <- theta0

  # Find a reasonable initial value of epsilon using heuristic
  g(log.temp, grad0) %=% L(theta0)
  epsilon =  FindReasonableEpsilon(theta0, log.temp, grad0, L)

  # Setting up other parameters for dual averaging
  mu <- log(10*epsilon)
  epsilon.bar <- 1
  H.bar <- 0
  gamma = 0.05
  t0 = 10
  kappa = 0.75

  for(m in 2:Madapt+M){
    r0 <- rnorm(length(theta0))
    #resample u (logtemp updated each iteration)
    ratio.temp <- log.temp - 0.5*crossprod(r0, r0)
    u = runif(1, min=0, max= exp(ratio.temp))
    joint0 <- ratio.temp

    #Setting up initial parameters for BuildTree
    theta.minus <- out[m-1, ]
    theta.plus <- out[m-1, ]
    r.minus <- r0
    r.plus <- r0
    out[m, ] <- out[m-1, ]

    j <- 0
    n <- 1
    s <- 1


    print(m)
    # Build Subtrees
    while(s==1){
      #Choose a direction : forward direction (+1), backward (-1)
      v  = 2*as.numeric(runif(1)>0.5)-1

      if (v==-1){
        g(theta.minus, r.minus, trivial, trivial, theta.prime, n.prime, s.prime, alpha, n.alpha, log.prime) %=% BuildTree(theta.minus,r.minus, u, v, j-1, epsilon, L, joint0  )
      }
      else{
        #### Correct this!!
        g(trivial,trivial, theta.plus, r.plus, theta.prime, n.prime, s.prime, alpha, n.alpha, log.prime) %=% BuildTree(theta.minus,r.minus, u, v, j-1, epsilon, L, joint0 )
      }

      if (s.prime==1 & runif(1)< min(1, n.prime/n) ) {
        out[m, ] = theta.prime
        log.temp = log.prime
      }
      n <- n+ n.prime
      s <- s.prime * StopCon(theta.minus, theta.plus, r.minus,r.plus )
      j <- j+1
    }

    # Using Dual Averaging to adapt epsilon (Burn-in period -- modify sample size?
    # discard the burn-in samples?)
    if (m < Madapt) {
      temp <- 1/(m+t0)
      Hbar <- (1 - temp)*Hbar + temp*(delta - alpha/n.alpha)
      epsilon <- exp(mu - sqrt(m)/gamma)*Hbar
      temp <- m^{-kappa}
      epsilonbar <- exp(temp*log(epsilon)+(1-temp)*log(epsilonbar))
    }
    else{
      epsilon <- epsilonbar
    }

  }
   return(out[Madapt+1:M+Madapt , ])

}


#testing two-dimensional multivariate normal

#Setting up callable function L








