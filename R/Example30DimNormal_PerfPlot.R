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
samples.hmc <- Hmc(theta0 = theta0, epsilon = 0.2, leap.nsteps = 10, L = L, M = 2000)
samples.hmcdual <- HmcDual(theta0 = theta0, delta = 0.65, lambda = 1.5, L = L, M = 2000, Madapt = 1000)
samples.nuts <- NutsDual(theta0 = theta0, delta = 0.6, L = L, M = 2000, Madapt = 1000)
hist(samples.hmc[,2],freq=FALSE,breaks=50,col=rgb(1,0,0,0.5),ylim=c(0,0.25))
hist(samples.hmcdual$samples[,2],freq=FALSE,breaks=50,add=T,col=rgb(0,0,1,0.5))
hist(samples.nuts$samples[,2],freq=FALSE,breaks=50,add=T,col=rgb(0,1,0,0.5))
x = seq(-5,5,len=100)
lines(x,0.5*dexp(abs(x)),col="black") #true density in black
















theta0 <- rnorm(30)
time.hmc <- proc.time()
plotmatHMC.norm30 <- PerfPlotHmc(theta0 = theta0, lambda=lseq(0.5)[c(1,3,5,7,9)], L = L, M = 2000, Madapt = 1000)
plotmatHMC.smallnorm30 <- PerfPlotHmc(theta0 = theta0, lambda=sort(0.5*40^-(0:(9-1)/9)[c(3,5,7)]), L = L, M = 2000, Madapt = 1000)
time.nuts<-proc.time()
plotmatNUTS.norm30<-PerfPlotNuts(theta0 = theta0, L = L, M = 2000, Madapt = 1000)
finalplot.dualav <- as.data.frame(rbind(plotmatNUTS.norm30$DualAvMatrix,plotmatHMC.norm30$DualAvMatrix))
finalplot.ESS <- as.data.frame(rbind(plotmatNUTS.norm30$EssMatrix,plotmatHMC.norm30$EssMatrix))


plotmatHMCDualMod.norm30 <- PerfPlotHmcDualMod(theta0 = theta0, lambda=lseq(0.5)[c(1,3,5,7,9)], L = L, M = 2000, Madapt = 1000)



finalplot.dualav <- as.data.frame(rbind(plotmatNUTS.norm30$DualAvMatrix,plotmatHMCDualMod.norm30$DualAvMatrix))
finalplot.ESS <- as.data.frame(rbind(plotmatNUTS.norm30$EssMatrix,plotmatHMCDualMod.norm30$EssMatrix))

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
spESS <- ggplot(finalplot.ESS, aes(x=Delta, y=ESS)) + ylim(0,2000) + geom_point(shape=1)
spESS <- spESS + facet_grid( ~ Model)
spESS <- spESS + stat_summary(fun.y=mean, colour="red", geom="line")
spESS <- spESS + theme_bw()
spESS <- spESS + ggtitle("30D Normal") + theme(plot.title = element_text(hjust = 0.5))
spESS <- spESS + labs(x = expression("Target acceptance rate statistics "*delta), y = "ESS")
spESS


###########################################################
### 250-dim Multivariate Normal ###########################
###########################################################
# A.mat <- rWishart(30, 30, diag(30))
# L <- function(theta){
#   grad0 = - theta %*% A.mat
#   log0 = 0.5*t(theta) %*%grad0
#   return(list(log0, grad0))
# }


time.hmc <- proc.time()
plotmatHMC.norm30 <- PerfPlotHmc(theta0 = theta0, lambda=lseq(0.5)[c(1,3,5,7,9)], L = L, M = 2000, Madapt = 1000)
plotmatHMC.smallnorm30 <- PerfPlotHmc(theta0 = theta0, lambda=sort(0.5*40^-(0:(9-1)/9)[c(3,5,7)]), L = L, M = 2000, Madapt = 1000)
time.nuts<-proc.time()
plotmatNUTS.norm30<-PerfPlotNuts(theta0 = theta0, L = L, M = 2000, Madapt = 1000)
finalplot.dualav <- as.data.frame(rbind(plotmatNUTS.norm30$DualAvMatrix, plotmatHMC.smallnorm30$DualAvMatrix, plotmatHMC.norm30$DualAvMatrix))
finalplot.ESS <- as.data.frame(rbind(plotmatNUTS.norm30$EssMatrix, plotmatHMC.smallnorm30$EssMatrix, plotmatHMC.norm30$EssMatrix))

###########################################################
### 30-dim normal, NormalTestData.RData ###################
###########################################################

DualAvMatrix$Discrepancy[i]<-mean(acceptrate,na.rm=T)-DualAvMatrix$Delta[i]
ESS <- multiESS(sampleHMD1)   ### replaced multiESS with ess


