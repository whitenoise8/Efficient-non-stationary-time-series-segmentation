#: Individual Checks
rm(list=ls())
library(CPSuSiE)
library(TSA)
library(astsa)
library(BayesSpec)

n = 1000
ncp = 4

load(paste0("data_sim/setting_n",n,"_ncp",ncp,".RData"))

it = 1

y = read.table(paste0("./data_sim/setting_n",n,"_ncp",ncp,"/data",it,".txt"),header=T)$x
s = read.table(paste0("./data_sim/setting_n",n,"_ncp",ncp,"/signal",it,".txt"),header=T)$x

#: Set omega grid and define covariates
par(mfrow=c(2,1),mar=c(2,2,2,2))
plot(y,type='l')
abline(v=cp[it,],lty='dashed',lwd=2,col=2)
plot(s,type='l')
abline(v=cp[it,],lty='dashed',lwd=2,col=2)
par(mfrow=c(1,1),mar=c(2,2,2,2))

#: Grid using periodogram
p = periodogram(y,plot=T)

p_sort = p$freq[sort.int(p$spec,decreasing=T,index.return=T)$ix]
n_om = 100
om_grid = sort(p_sort[1:n_om])
plot(om_grid,rep(0,n_om))

x = seq(1,n)
Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

l = 1
cp_set <- getCandidateSet(y,Xsin,Xcos,L=l,nu=1/2,Trace=1)
cp_hat <- pruneSet(cp_set,y,Xsin,Xcos,L=l,IC="MDL",option='joint',Trace=1)

ft = fit_segments(y,Xsin,Xcos,L=4,cp=cp_hat,pl=TRUE)
abline(v=cp[it,],col=4,lwd=3,lty="dashed")
lines(s,lwd=2,col=4)

cp_set <- getCandidateSet(y,Xsin,Xcos,L=l,sigma2=9,nu=1/2,Trace=1)
cp_hat <- pruneSet(cp_set,y,Xsin,Xcos,L=l,sigma2=9,IC="MDL",option='joint',Trace=1)

ft = fit_segments(y,Xsin,Xcos,L=4,cp=cp_hat,sigma2=9,pl=TRUE)
abline(v=cp[it,],col=4,lwd=3,lty="dashed")
lines(s,lwd=2,col=4)

#: Grid using seq
n_om = 200
om_grid = seq(1e-3,0.5,length.out=n_om)
plot(om_grid,rep(0,n_om))

x = seq(1,n)
Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

l = 1
cp_set <- getCandidateSet(y,Xsin,Xcos,L=l,nu=1/2,Trace=1)
cp_hat <- pruneSet(cp_set,y,Xsin,Xcos,L=l,IC="MDL",option='joint',Trace=1)

ft = fit_segments(y,Xsin,Xcos,L=4,cp=cp_hat,pl=TRUE)
abline(v=cp[it,],col=4,lwd=3,lty="dashed")
lines(s,lwd=2,col=4)



#: AUTOPARM (Davis et al., 2006)
system.time(fit_ap <- autoParm(y))
cp_ap = fit_ap$breakpoints

s_ap = fit_autoParm(y,fit_ap)
abline(v=cp[it,],col=4,lwd=3,lty="dashed")
lines(s,lwd=2,col=4)


#: AdaptSpec (Rosen, Wood and Stoffer, 2012)
fit_ada <- adaptspec(2000, 500, 4, y)

n_seg = round(mean(fit_ada$nexp_curr))
n_cp = n_seg-1
if (n_cp > 0) cp_hat = round(rowMeans(fit_ada$xi[[n_seg]][1:n_cp,]))
if (n_cp == 0) cp_hat = NULL
cp_hat
