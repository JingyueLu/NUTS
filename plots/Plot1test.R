###PLOTS: average acceptance rates

##Normal Case
L <- function(theta){
  cov = c(1,4)
  grad0 = - theta/cov
  log0 = 0.5*t(theta) %*%grad0
  return(list(log0, grad0))
}

theta0 <- rnorm(2)
r <- rnorm(2)
# Testing Find epsilon


M = 5000
Madapt = 1000

g(samples, acceprate, epsilonNor, height) %=% NutsDual(theta0, 0.6, L, M, Madapt)

deltaSeq <-seq(0.25, 0.95, by= 0.10 )
NumSeed = 3
DualAvMatrix = matrix(0, nrow =length(deltaSeq)*NumSeed, ncol=4)
EssMatrix = matrix(0, nrow =length(deltaSeq)*NumSeed, ncol=4)

for (i in 1:length(deltaSeq) ){
  for (j in 1:NumSeed){
    set.seed(j)
    g(samples, acceprate, epsilonNor, steps) %=% NutsDual(theta0, deltaSeq[i], L, M, Madapt)
    DualAvMatrix[(NumSeed*(i-1)+j),] <-  c(j, deltaSeq[i], (mean(acceprate)-deltaSeq[i]),9999)
    EssMatrix[(NumSeed*(i-1)+j),] <- c(j, deltaSeq[i], multiESS(samples),9999)
    print(NumSeed*(i-1)+j)
  }
}

essbar <- c()
for (i in 1:length(deltaSeq)){
  a <-EssMatrix[((i-1)*NumSeed+1):((i-1)*NumSeed+NumSeed), 3]
  essbar <- c(essbar, mean(a))
}
plot(deltaSeq, essbar)

#EssMatrix = matrix(0, nrow =150, ncol=4)
