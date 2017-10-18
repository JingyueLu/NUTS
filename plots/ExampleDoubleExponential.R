#' #' ExampleDoubleExponential
#'
#' ### log density function for double exponential
#' logDens <- function(x)   -sum(abs(x))
#' dlogDens <- function(x)   -sign(x)
#'
#' L <- function(x) {
#'   return(list(- sum(abs(x)), - sign(x)))
#' }
#'
#'
#' ### For Hmc plot the MCMC chains and densities
#' samples.hmc <- Hmc(theta0 = c(1,1), epsilon = 0.2, leap.nsteps = 10, L = L, M = 5000)
#' ts.plot(samples.hmc[,1])
#' ts.plot(samples.hmc[,2])
#' # Plot a histogram of the second variable, with true density
#' hist(samples.hmc[,2],freq=FALSE,breaks=50)
#' x = seq(-5,5,len=100)
#' lines(x,0.5*dexp(abs(x)),col="red") #true density in red
#' # Autocorrelation function
#' acf(samples.hmc)
#'
#'
#' ### For HmcDual plot the MCMC chains and densities
#' samples.hmcdual <- HmcDual(theta0=c(1,1), delta = 0.65, lambda = 1.5, L = L, M = 5000, Madapt = 2500)
#' ts.plot(samples.hmcdual$samples[,1])
#' ts.plot(samples.hmcdual$samples[,2])
#' # Plot a histogram of the second variable, with true density
#' hist(samples.hmcdual$samples[,2],freq=FALSE,breaks=50)
#' x = seq(-5,5,len=100)
#' lines(x,0.5*dexp(abs(x)),col="red") #true density in red
#' # Autocorrelation function
#' acf(samples.hmcdual$samples)
#'
#' ### For NUTS plot the MCMC chains and densities
#' samples.nuts <- NutsDual(theta0 = c(1,1), delta = 0.6, L = L, M = 5000, Madapt = 2500)
#' ts.plot(samples.nuts$samples[,1])
#' ts.plot(samples.nuts$samples[,2])
#' # Plot a histogram of the second variable, with true density
#' hist(samples.nuts$samples[,2],freq=FALSE,breaks=50)
#' x = seq(-5,5,len=100)
#' lines(x,0.5*dexp(abs(x)),col="red") #true density in red
#' # Autocorrelation function
#' acf(samples.nuts$samples)
#'
#' ################################################################
#' #Comparing histograms
#'
#' hist(samples.hmc[,1],freq=FALSE,breaks=50,col=rgb(1,0,0,0.5),ylim=c(0,0.5), xlab = "Sampled values",main = "Histograms of the first variable" )
#' hist(samples.hmcdual$samples[,1],freq=FALSE,breaks=50,add=T,col=rgb(0,0,1,0.5))
#' hist(samples.nuts$samples[,1],freq=FALSE,breaks=50,add=T,col=rgb(0,1,0,0.5))
#' x = seq(-5,5,len=100)
#' lines(x,0.5*dexp(abs(x)),col="black") #true density in black
#' library(rmutil)
#' lines(x,dlaplace(x),col="orange") #true density in black
#' box()
#' legend("topright",col = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0,1,0,0.5)), pch = 18, legend = c("Hmc","HmcDual","NutsDual"))
#'
#' ################################################################
#' ################################################################
#' ################################################################
#'
#' #' Comparing the performances of this version of HMC with respect to HybridMC::hybridMC
#' #'
#' #' This script contains the example from the function hybridMC (in the package HybridMC), that calls C code.
#' #' At the end an overview of the performances of the two functions is presented.
#'
#'
#' microbenchmark(hybridMC(y.start = c(1,1), n.samp=2000, logDens=logDens, dLogDens=dlogDens, epsilon=.2, LFsteps=10),
#'                Hmc(theta0 = c(1,1), epsilon = 0.2, leap.nsteps = 10, L = L, M = 2000))
#'
#' ################################################################
#' ################################################################
#' ################################################################
#'
#' time.hmc.doubexp <- proc.time()
#' plotmatHMC.doubexp <- PerfPlotHmc(theta0 = c(1, 1), lambda=lseq(0.5)[c(1,3,5,7,9)], L = L, M = 2000, Madapt = 1000)
#' time.nuts.doubexp<-proc.time()
#' plotmatNUTS.doubexp<-PerfPlotNuts(theta0 = c(1, 1), L = L, M = 2000, Madapt = 1000)
#' finalplot.dualav <- as.data.frame(rbind(plotmatNUTS.doubexp$DualAvMatrix,plotmatHMC.doubexp$DualAvMatrix))
#' finalplot.ESS <- as.data.frame(rbind(plotmatNUTS.doubexp$EssMatrix,plotmatHMC.doubexp$EssMatrix))
#'
#' plotmatHMC.doubexpSMALL <- PerfPlotHmc(theta0 = c(1, 1), lambda=sort(0.5*40^-(0:(9-1)/9)[c(3,5,7)]), L = L, M = 2000, Madapt = 1000)
#'
#'
#'
#'
#'
#' plotmatHMCMOD.doubexp <- PerfPlotHmcDualMod(theta0 = c(1, 1), lambda=c(0.22,0.5,1.135,13.275), L = L, M = 2000, Madapt = 1000)
#'
#'
#'
#'
#'
#'
#' #plotmatHMC.doubexp$ESSess <- ess(samples.hmcdual)
#' #plotmatNUTS.doubexp$ESSess<- ess(samples.nuts)
#' finalplot.dualav <- as.data.frame(rbind(plotmatNUTS.doubexp$DualAvMatrix,plotmatHMC.doubexpSMALL$DualAvMatrix,plotmatHMC.doubexp$DualAvMatrix))[-c(151:310,551:630,711:790),]
#' finalplot.ESS <- as.data.frame(rbind(plotmatNUTS.doubexp$EssMatrix,plotmatHMC.doubexpSMALL$EssMatrix,plotmatHMC.doubexp$EssMatrix))[-c(151:310,551:710),]
#'
#'
#' ### Plot the discrepancies as in Figure 3 in Hoffman & Gelman (2014)
#'
#' library(ggplot2)
#' sp <- ggplot(finalplot.dualav, aes(x=Delta, y=Discrepancy)) + ylim(-0.1,0.2) + geom_point(shape=1)
#' sp <- sp + facet_grid( ~ Model)
#' sp <- sp + stat_summary(fun.y=mean, colour="red", geom="line")
#' sp <- sp + theme_bw()
#' sp <- sp + ggtitle("Double Exponential") + theme(plot.title = element_text(hjust = 0.5))
#' sp <- sp + labs(x = expression("Target acceptance rate statistics "*delta), y = expression("h - "*delta))
#' sp
#'
#' ### Plot the discrepancies (not exactly) like in Figure 6 in Hoffman & Gelman (2014)
#'
#' library(ggplot2)
#' spESS <- ggplot(finalplot.ESS, aes(x=Delta, y=ESS)) + ylim(0,2000) + geom_point(shape=1)
#' spESS <- spESS + facet_grid( ~ Model)
#' spESS <- spESS + stat_summary(fun.y=mean, colour="red", geom="line")
#' spESS <- spESS + theme_bw()
#' spESS <- spESS + ggtitle("Double Exponential") + theme(plot.title = element_text(hjust = 0.5))
#' spESS <- spESS + labs(x = expression("Target acceptance rate statistics "*delta), y = "ESS")
#' spESS
