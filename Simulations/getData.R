rm(list=ls())

#: Define the setting:
n = 1200
nseg = 4

omega = list()
omega[[1]] = c(1/25, 1/10)
omega[[2]] = c(1/20)
omega[[3]] = c(1/30, 1/15)
omega[[4]] = c(1/20)

beta = list()
beta[[1]] = matrix(c(0.5, 1.5, 
                     2.0, 1.0),2,2)
beta[[2]] = matrix(c(2.0,
                     1.0),1,2)
beta[[3]] = matrix(c(1.5, 2.0, 
                     0.5, 1.0),2,2)
beta[[4]] = matrix(c(0.5, 
                     2.0),1,2)

#: Create data
Nsim = 20
ncp = nseg - 1
minLen = 100
seeds = sample(1:1e4,Nsim)
cp = matrix(NA,Nsim,ncp)

for (it in 1:Nsim) {
  set.seed(seeds[it])
  
  tau = rep(2,ncp)
  while (min(diff(c(1,tau,n))) < minLen) tau = sample(1:n,ncp)
  cp[it,] = tau
  
  knots = c(1,tau,n+1)
  y = NULL
  s = NULL
  
  for (i in 1:nseg) {
    x = seq(knots[i],knots[i+1]-1)
    z = 0
    for (k in 1:length(omega[[i]])) z = z +
      beta[[i]][k,1]*sin(2*pi*x*omega[[i]][k])+beta[[i]][k,2]*cos(2*pi*x*omega[[i]][k])
    s = c(s,z)
    y = c(y,z+rnorm(length(x)))
  }
  
  write.table(y,file=paste0("data_sim/setting_1/data",it,".txt"),row.names=F)
  write.table(s,file=paste0("data_sim/setting_1/signal",it,".txt"),row.names=F)
}

plot(y)

save(list=c("n","nseg","omega","beta","cp",
            "Nsim","ncp","minLen","seeds"),file="data_sim/setting_1.RData")

