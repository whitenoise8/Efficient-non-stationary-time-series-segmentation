rm(list=ls())
library(CPSuSiE)
library(TSA)
library(BayesSpec)
library(astsa)

######################
#: Simulation study :#
######################

#: Setting
load("data_sim/setting_1.RData")

#: Collect results
cp_mat = array(NA,dim=c(Nsim,length(cp[1,]),4))
s_mat = array(NA,dim=c(n,Nsim,4))
ctime = matrix(NA,Nsim,4)

#: Params
l = 4
n_om = 50
f = 1

for (it in 1:10) {
  
  y = read.table(paste0("./data_sim/setting_1/data",it,".txt"),header=T)$x
  s = read.table(paste0("./data_sim/setting_1/signal",it,".txt"),header=T)$x
  
  #: Set omega grid and define covariates
  p = periodogram(y,plot=F)
  p_sort = p$freq[sort.int(p$spec,decreasing=T,index.return=T)$ix]
  om_grid = p_sort[1:n_om]

  x = seq(1,n)
  Xout = getX(x,om_grid)
  Xsin = Xout$Xsin
  Xcos = Xout$Xcos
  
  #: Bin Segm
  # ctime[it,1] = system.time(out_bs <- BS_MCP(y,Xsin,Xcos,om_grid,L=l,f=f,cost=10))[3]

  #: Naive OS
  ctime[it,1] = system.time(out_nos <- nOS_MCP(y,Xsin,Xcos,L=l,f=f,cost=10,nu=1/3))[3]
  
  ft = fit_segments(y,Xsin,Xcos,L=l,f=f,cp=out_nos,pl=F)
  if (length(out_nos) == length(cp[it,])) cp_mat[it,,1] = out_nos
  s_mat[,it,1] = ft$signal_estimate
  
  print(it)
}

cbind(cp_mat[1:it,,1],cp[1:it,])
cbind(ctime[1:it,1])


#: Individual Checks
rm(list=ls())
library(CPSuSiE)
library(TSA)
library(BayesSpec)
library(astsa)

load("data_sim/setting_1.RData")

it = 1

y = read.table(paste0("./data_sim/setting_1/data",it,".txt"),header=T)$x
s = read.table(paste0("./data_sim/setting_1/signal",it,".txt"),header=T)$x

#: Set omega grid and define covariates
plot(y)

p = periodogram(y,plot=T)
abline(v=unlist(omega),col=2)

p_sort = p$freq[sort.int(p$spec,decreasing=T,index.return=T)$ix]
n_om = 50
om_grid = p_sort[1:n_om]

plot(om_grid,rep(0,n_om))
abline(v=unlist(omega),col=2)

x = seq(1,n)
Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

l = 4
f = 1

system.time(out_nos <- nOS_MCP(y,Xsin,Xcos,L=l,f=f,cost=10,nu=1/3))

ft = fit_segments(y,Xsin,Xcos,L=l,f=f,cp=out_nos,pl=TRUE)
abline(v=cp[it,],col=4,lwd=3,lty="dashed")
lines(s,lwd=2,col=4)

par(mfrow=c(2,2),mar=c(2,2,2,2))
for (sgm in 1:4) {
plot(om_grid,1-apply(ft$models[[sgm]]$alpha,1,function(x) prod(1-x)))
abline(v=omega[[sgm]],col=2,lwd=1,lty='dashed')
}
par(mfrow=c(1,1),mar=c(2,2,2,2))


# Prova adaptspec e autoparm ma boh
plot(y)
fit_ada = adaptspec(1000, 500, 4, y)
str(fit_ada)

fit_ap = autoParm(y)
fit_ap$breakpoints
