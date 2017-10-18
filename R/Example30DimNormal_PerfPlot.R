#' Simulate a 30-dimensional multivariate normal distribution
#' and reproduce Figure 3 and Figure 6 from Hoffman & Gelman (2014)

sdevs <- c (110, 100, seq(16,8,length=26), 1.1, 1.0)
V <- sdevs^2
L <- function(theta){
  grad0 = - theta/V
  log0 = 0.5*t(theta) %*%grad0
  return(list(log0, grad0))
}

theta0 <- rnorm(30)
time.hmc <- proc.time()
plotmatHMC.norm30 <- PerfPlotHmc(theta0 = theta0, lambda=lseq(0.5)[c(1,3,5,7,9)], L = L, M = 2000, Madapt = 1000)
plotmatHMC.smallnorm30 <- PerfPlotHmc(theta0 = theta0, lambda=sort(0.5*40^-(0:(9-1)/9)[c(3,5,7)]), L = L, M = 2000, Madapt = 1000)
time.nuts<-proc.time()
plotmatNUTS.norm30<-PerfPlotNuts(theta0 = theta0, L = L, M = 2000, Madapt = 1000)
finalplot.dualav <- as.data.frame(rbind(plotmatNUTS.norm30$DualAvMatrix,plotmatHMC.norm30$DualAvMatrix))
finalplot.ESS <- as.data.frame(rbind(plotmatNUTS.norm30$EssMatrix,plotmatHMC.norm30$EssMatrix))


###########################################################
### 250-dim Multivariate Normal ###########################
###########################################################
# A.mat <- rWishart(30, 30, diag(30))
# L <- function(theta){
#   grad0 = - theta %*% A.mat
#   log0 = 0.5*t(theta) %*%grad0
#   return(list(log0, grad0))
# }
