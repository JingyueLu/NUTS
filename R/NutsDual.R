#' %=%
#'
#' The function is copied from "https://strugglingthroughproblems.wordpress.com/2010/08/27/matlab-style-multiple-assignment-in%C2%A0r/" and is used for multiple assignment
#' @param l		the list on the left-hand side
#' @param r 	the list on the right-hand side
#'
#' @example
#' g(a,b) %=% list(c(3,5,6),4)
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
#' This function perform a leapfrog step. This function is a modified version of Leapfrog in the paper. It returns etra values: log posterior value and gradient of log posterior at the new position theta.tilde
#'
#' @param theta		starting position
#' @param r 	starting momentum
#' @param epsilon 	step size
#' @param L 	callable function: returns the value of log posterior and the gradient of log posterior prbability at given input
#'
#' @return the list of updated theta, r and the log posterior value at the updated point
#'
#' @example

Leapfrog <- function(theta, r, epsilon, L){
  g(trivial,grad.theta) %=% L(theta)
  r.tilde <- r+ 0.5*epsilon*grad.theta
  theta.tilde <- theta + epsilon*r.tilde
  g(log.tilde, grad.tilde) %=% L(theta.tilde)
  r.tilde <- r.tilde + 0.5*epsilon*grad.tilde

  return(list(theta.tilde, r.tilde, log.tilde))

}


#' Build trees
#'
#' This function builds the tree for NUTS
#'
#' @param theta starting position
#' @param r starting momentum
#' @param u a slice variable, which is used to simplify the implementation
#' @param v a randomly generated value, which is used to determine the direction of movement
#' @param j the height of a tree
#' @param epsilon step size
#' @param joint0 	to avoid computing the same value repeatedly, we compute joint0=log0- 0.5*(r0 %*% r0) beforehand and take this value as an input for the BUildTree function
#'
#' @return a list of variables. It includes updated leftmost and rightmost states. It also returns updated values which are used to improve memory efficience.
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
        g(trivial, trivial, theta.plus, r.plus, theta.prime2, n.prime2, s.prime2, alpha.prime2, n.alphaprime2, log.prime2) %=% BuildTree(theta.plus, r.plus, u,v,j-1, epsilon, L, joint0)
      }
      # Make sure the DENOMINATER IS GREATER THAN 0
      if (runif(1)< n.prime2/max(n.prime+n.prime2, 1)){
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
#'#' This function computes the U-Turn stopping condition, which is used in NUTS and BuildTree functions
#'
#' @param theta.minus the leftmost position of a subtree
#' @param theta.plus the rightmost position of a subtree
#' @param r.minus the leftmost momentum of a subtree
#' @param r.plus the leftmost position of a subtree
#'
#' @return  1 if the stopping criterion is met by the subtree; 0 if the stopping criterion is not by the subtree
#'
StopCon <-function(theta.minus, theta.plus, r.minus,r.plus){
  theta.diff <- theta.plus - theta.minus
  temp1 <- as.numeric(I(crossprod(theta.diff,r.minus)>0))
  temp2 <- as.numeric(I(crossprod(theta.diff, r.plus)>0))
  return (temp1*temp2)
}


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
FindReasonableEpsilon <- function(theta0,log.start, L){
  epsilon <- 1
  r = rnorm(length(theta0))
  g(theta.prime, r.prime, log.prime) %=% Leapfrog(theta0,r,epsilon, L)

  #Check whether log density is infinite. If it is, half the epsilong and the process continues until log density is no longer infinite
  k = 0.5
  while (is.infinite(log.prime)){
    k <- k* 0.5
    g(trivial, r.prime, log.prime) %=% Leapfrog(theta0, r, epsilon * k, L)}

  epsilon = 0.5 * k * epsilon

  #Computing the ratio of p(theta.prime, r.prime) over p(theta, r)
  tempratio <- exp(log.prime-0.5*(crossprod(r.prime, r.prime)) - log.start + 0.5*(crossprod(r, r)))
  a <- 2* as.numeric(tempratio > 0.5)-1
  #print(c(tempratio, a))

  while (tempratio^a > 2^(-a) ){
    epsilon <-  epsilon * (2^(a))
    g(theta.prime, r.prime, log.prime) %=% Leapfrog(theta0,r,epsilon, L)
    tempratio <- exp(log.prime - 0.5*(crossprod(r.prime, r.prime)) - log.start + 0.5*(crossprod(r, r)))
    #print(tempratio)
  }
  return(epsilon)

}


#' NutsDual
#'
#' This function implements No-U-Turn Sampler with Dual Averaging.
#' The code is a modified version of algorithm6 in paper. Instead of including samples in burn-in period (m < Madapt), the function generates a sample where samples in the burn-in period are discarded.
#'
#' @param theta0	inital state of theta
#' @param delta 	the desired average acceptance probability
#' @param L 	a callable function that returns log probability
#' @param M 	number of samples to generate
#' @param Madapt	number of iterations
#'
#'
#' @return samples generated by NUTS

NutsDual <- function(theta0, delta, L, M, Madapt){

  # Set up output structure: each theta is an array
  # samples is a matrix of size M x length(theta)
  len <- length(theta0)
  out <- matrix(0 , nrow = Madapt + M, ncol = len )
  out[1, ] <- theta0

  # Find a reasonable initial value of epsilon using heuristic
  g(log.temp, grad0) %=% L(theta0)
  epsilon =  FindReasonableEpsilon(theta0, log.temp, L)

  # Setting up other parameters for dual averaging
  mu <- log(10*epsilon)
  epsilon.bar <- 1
  H.bar <- 0
  gamma = 0.05
  t0 = 10
  kappa = 0.75

  # Initialise a vector for storing average acceptance probability statistic
  accproV <- c()
  # Initialise a vector for storing epsilon bar
  epsibarV <-c()
  # Initialise a vector for storing the number of doublings j for each sample point (height of the tree)
  heightTV <-c()
  #stepsV <-c()

  for(m in 2:(Madapt+M)){
    r0 <- rnorm(len)

    #resample u (logtemp updated each iteration)
    joint <- log.temp - 0.5*crossprod(r0, r0)
    u = runif(1, min=0, max= exp(joint))


    #Setting up initial parameters for BuildTree
    theta.minus <- out[m-1, ]
    theta.plus <- out[m-1, ]
    r.minus <- r0
    r.plus <- r0

    out[m, ] <- out[m-1, ]

    j <- 0
    n <- 1
    s <- 1

    #print(m)
    # Build Subtrees
    while(s==1){
      #Choose a direction : forward direction (+1), backward (-1)

      v  = 2*as.numeric(runif(1)>0.5)-1

      if (v==-1){
        g(theta.minus, r.minus, trivial, trivial, theta.prime, n.prime, s.prime, alpha, n.alpha, log.prime) %=% BuildTree(theta.minus,r.minus, u, v, j, epsilon, L, joint)
      }
      else{
        g(trivial,trivial, theta.plus, r.plus, theta.prime, n.prime, s.prime, alpha, n.alpha, log.prime) %=% BuildTree(theta.plus,r.plus, u, v, j, epsilon, L, joint)
      }

      if (s.prime==1 & runif(1)< min(1, n.prime/n) ) {
        out[m, ] = theta.prime
        log.temp = log.prime
      }
      n <- n+ n.prime
      s <- s.prime * StopCon(theta.minus, theta.plus, r.minus,r.plus )
      j <- j+1
    }

    heightTV <-c(heightTV, j)
    #stepsV <-c(stepsV, steps)

    # Using Dual Averaging to adapt epsilon (Burn-in period -- modify sample size?
    # discard the burn-in samples?)

    acc<-alpha/n.alpha
    accproV <- c(accproV, acc)

    if (m < Madapt) {
      temp <- 1/(m-1+t0)
      H.bar <- (1 - temp)*H.bar + temp*(delta - acc)
      epsilon <- exp(mu - (sqrt(m-1)/gamma)*H.bar)
      temp <- (m-1)^{-kappa}
      epsilon.bar <- exp(temp*log(epsilon)+(1-temp)*log(epsilon.bar))
      epsibarV <- c(epsibarV, epsilon.bar)
    }
    else{
      epsilon <- epsilon.bar
    }

  }
  return(list(samples=out[(Madapt+1):(M+Madapt),], acceptrate= accproV[(Madapt+1): length(accproV)], epsilonNor = epsibarV/epsilon, Length=heightTV))

}












