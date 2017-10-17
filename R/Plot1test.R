###PLOTS: average acceptance rates
M = 2000
Madapt = 1000
DualAvMatrix = matrix(0, nrow =150, ncol=4)
for (i in seq(0.25, 0.95, by= 0.05)){
  for (j in 1:10){
    set.seed(j)
    g(samples, acceprate, epsilonNor, steps) %=% NutsDual(theta0, delta, L1, M, Madapt)
    DualAvMatrix[i,j] <  c(j, i, (mean(acceprate)-i),"NUTS")
  }
}


#EssMatrix = matrix(0, nrow =150, ncol=4)
