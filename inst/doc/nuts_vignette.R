### R code from vignette source 'nuts_vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: chu (eval = FALSE)
###################################################
##  L <- function(x) {
##    return(list(- sum(abs(x)), - sign(x)))
##  }
##  samples.hmc <- Hmc(theta0 = c(1,1), epsilon = 0.2, leap.nsteps = 10, L = L, M = 5000)
##  samples.hmcdual <- HmcDual(theta0=c(1,1), delta = 0.65, lambda = 1.5, L = L, M = 5000, Madapt = 2500)
##  samples.nuts <- NutsDual(theta0 = c(1,1), delta = 0.6, L = L, M = 5000, Madapt = 2500)


###################################################
### code chunk number 2: doublePrint
###################################################
  load("/homes/palma/Documents/NUTS/data/SAMPLE.RData")
   #load("https://github.com/jodie0399/NUTS/DoubExp.RData")
  print(samples.nuts)


