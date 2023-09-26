rm(list=ls())
library(CPSuSiE)
library(TSA)

######################
#: Simulation study :#
######################

#: Setting
load("data_sim/setting_1.RData")

#: Collect results
cp_mat = array(NA,dim=c(Nsim,length(cp[1,]),2))
s_mat = array(NA,dim=c(n,Nsim,2))
ctime = matrix(NA,Nsim,2)

#: Params
l = 4
n_om = 50

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
  # ctime[it,1] = system.time(out_bs <- BS_MCP(y,Xsin,Xcos,L=l))[3]

  #: Naive OS
  ctime[it,2] = system.time(out_nos <- nOS_MCP(y,Xsin,Xcos,L=l))[3]
  ft = fit_segments(y,Xsin,Xcos,L=l,cp=out_nos,pl=F)
  
  if (length(out_nos) == length(cp[it,])) cp_mat[it,,2] = out_nos
  s_mat[,it,2] = ft$s_hat
  
  print(it)
}

cbind(cp_mat[1:it,,2],cp[1:it,])
ctime[1:it,2]

#: Single case
it = 1

y = read.table(paste0("./data_sim/setting_1/data",it,".txt"),header=T)$x
s = read.table(paste0("./data_sim/setting_1/signal",it,".txt"),header=T)$x

plot(y)
abline(v=cp[it,],lwd=3,lty='dashed')

#: Set omega grid and define covariates
p = periodogram(y,plot=T)
p_sort = p$freq[sort.int(p$spec,decreasing=T,index.return=T)$ix]
om_grid = p_sort[1:n_om]

x = seq(1,n)
Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

#: Naive OS
l = 6
system.time(out_nos <- nOS_MCP(y,Xsin,Xcos,L=l))

ft = fit_segments(y,Xsin,Xcos,L=l,cp=out_nos,pl=T)
lines(s,col=4)

out_nos
cp[it,]

sum(dnorm(y,ft$s_hat,1,log=T))
