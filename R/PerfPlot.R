#' PerfPlotHmc
#'
#' Create two matrices:
#' - one with the discrepancies between the desired average acceptance rate
#' \delta and the mean acceptance rate \alpha in Hamiltonian Monte Carlo with
#' Dual Averaging (Figure 3 in Hoffman & Gelman, 2014);
#' - the other with the effective sample size (ESS).
#'
#' @param theta0 inizialization
#' @param lambda the trajectory length
#' @param L callable function needed in Leapfrog
#' @param M the number of iterations of the algorithm
#' @param Madapt the number of iterations of the warmup phase
#'
#' @return A list whose first element is the matrix DualAvMatrix of 10 discrepancies for
#' each value of \delta and the second element is the matrix ESSMatrix
#' of 10 ESS for each value of \delta.



PerfPlotHmc<-function(theta0, lambda, L, M = 2000, Madapt = 1000){
  deltaSeq <-seq(0.25, 0.95, length = 8)
  seed <- 0:9
  DualAvMatrix <- expand.grid(Seed = seed,Delta = deltaSeq, Discrepancy = NA, Model = lambda)
  EssMatrix <- expand.grid(Seed = seed, Delta = deltaSeq, ESS = NA, Model = lambda)
  for(i in 1:dim(DualAvMatrix)[1]){
    set.seed(i %% 10)
    g(samples, acceptrate) %=% HmcDual(theta0, delta = DualAvMatrix$Delta[i], lambda = DualAvMatrix$Model[i], L, M, Madapt)
    DualAvMatrix$Discrepancy[i]<-mean(acceptrate,na.rm=T)-DualAvMatrix$Delta[i]
    EssMatrix$ESS[i] <- multiESS(samples)
    print(round(100*i/dim(DualAvMatrix)[1],2))
  }

  DualAvMatrix$Model<-paste("HMC",round(DualAvMatrix$Model,3))
  EssMatrix$Model<-paste("HMC",round(EssMatrix$Model,3))
  return(list(DualAvMatrix = DualAvMatrix, EssMatrix = EssMatrix))
}

########################################################################################
########################################################################################
########################################################################################

#' PerfPlotNuts
#'
#' Create two matrices:
#' - one with the discrepancies between the desired average acceptance rate
#' \delta and the mean acceptance rate \alpha in NUTS with
#' Dual Averaging (Figure 3 in Hoffman & Gelman, 2014);
#' - the other with the effective sample size (ESS).
#'
#' @param theta0 inizialization
#' @param L callable function needed in Leapfrog
#' @param M the number of iterations of the algorithm
#' @param Madapt the number of iterations of the warmup phase
#'
#' @return A list whose first element is the matrix DualAvMatrix of 10 discrepancies for
#' each value of \delta and the second element is the matrix ESSMatrix
#' of 10 ESS for each value of \delta.


PerfPlotNuts<-function(theta0, L, M = 2000, Madapt = 1000){
  require(mcmcse)
  deltaSeq <-seq(0.25, 0.95, length = 15)
  seed <- 0:9
  DualAvMatrix <- expand.grid(Seed = seed, Delta = deltaSeq, Discrepancy = NA, Model = "NUTS")
  EssMatrix <- expand.grid(Seed = seed, Delta = deltaSeq, ESS = NA, Model = "NUTS")
  for(i in 1:dim(DualAvMatrix)[1]){
    set.seed(i %% 10)
    g(samples, acceptrate) %=% NutsDual(theta0, delta = DualAvMatrix$Delta[i], L, M, Madapt)
    DualAvMatrix$Discrepancy[i]<-mean(acceptrate,na.rm=T)-DualAvMatrix$Delta[i]
    EssMatrix$ESS[i] <- multiESS(samples)
    print(round(100*i/dim(DualAvMatrix)[1],2))
  }
  return(list(DualAvMatrix = DualAvMatrix, EssMatrix = EssMatrix))
}



##################################################################
#Example 1 - Bivariate Normal Distribution #######################
##################################################################
#' Bivariate Normal - Gradient
#'
#' @param theta a 2-dimensional vector of parameters
#'
#' @return A list containing the evaluation of the log distribution
#' and its gradient evaluated at \theta


L <- function(theta){
  cov =c(1,4)
  grad0 = - theta/cov
  log0 = 0.5*t(theta) %*%grad0
  return(list(log0, grad0))
}

#' Sequence logarithmically scaled
#'
#' @param from the starting point of the desired sequence
#' @param scale the ratio between the first and last element of the sequence
#' @param len the length of the sequence
#'
#' @return A sequence of length len containing the sequence of logarithmically scaled values

lseq <- function(from=0.5, scale = 40, len = 9) {
  return(from*scale^(0:(len-1)/len))
}

### Run the functions PerfPlot

plotmatHMC<-PerfPlotHmc(theta0 = c(rnorm(1),rnorm(1,sd=2)), lambda=lseq(0.03318664)[c(1,3,5,7,9)], L = L, M = 2000, Madapt = 1000)
plotmatNUTS<-PerfPlotNuts(theta0 = c(rnorm(1),rnorm(1,sd=2)), L = L, M = 2000, Madapt = 1000)
finalplot.dualav <- as.data.frame(rbind(plotmatNUTS$DualAvMatrix,plotmatHMC$DualAvMatrix))
finalplot.ESS <- as.data.frame(rbind(plotmatNUTS$EssMatrix,plotmatHMC$EssMatrix))



### Plot the discrepancies as in Figure 3 in Hoffman & Gelman (2014)

library(ggplot2)
sp <- ggplot(finalplot.dualav, aes(x=Delta, y=Discrepancy)) + ylim(-0.1,0.2) + geom_point(shape=1)
sp <- sp + facet_grid( ~ Model)
sp <- sp + stat_summary(fun.y=mean, colour="red", geom="line")
sp <- sp + theme_bw()
sp <- sp + ggtitle("30D Normal") + theme(plot.title = element_text(hjust = 0.5))
sp <- sp + labs(x = expression("Target acceptance rate statistics "*delta), y = expression("h - "*delta))
sp

### Plot the discrepancies (not exactly) like in Figure 6 in Hoffman & Gelman (2014)

library(ggplot2)
spESS <- ggplot(finalplot.ESS, aes(x=Delta, y=ESS)) + geom_point(shape=1)
spESS <- spESS + facet_grid( ~ Model)
spESS <- spESS + stat_summary(fun.y=mean, colour="red", geom="line")
spESS <- spESS + theme_bw()
spESS <- spESS + ggtitle("30D Normal") + theme(plot.title = element_text(hjust = 0.5))
spESS <- spESS + labs(x = expression("Target acceptance rate statistics "*delta), y = "ESS")
spESS
