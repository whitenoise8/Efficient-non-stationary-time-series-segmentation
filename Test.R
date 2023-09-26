rm(list=ls())
library(CPSuSiE)
library(TSA)

#######################################################
#: One change point (first part of Beniamino's data) :#
#######################################################

y = read.table("AutoNOM-Beniamino/Data/data1.txt")$V1
y = y[1:600]

plot(y)
abline(v=300,col=1,lty='dashed',lwd=2)

#: Set omega grid and define covariates
n = length(y)
p = periodogram(y,plot=T)
om_grid = sort(p$freq[p$spec>10])
#om_grid = seq(1e-3,0.5,length.out=100)
x = seq(1,n)
Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

#: Bin Segm
system.time(out_bs <- BS(y,Xsin,Xcos,L=4))
plot(y)
abline(v=out_bs$cp_hat,col=2,lty='dashed',lwd=2)
out_bs$cp_hat
plot(out_bs$IC_seq,ylim=c(min(out_bs$IC_seq,out_bs$IC0,na.rm=T),
                          max(out_bs$IC_seq,out_bs$IC0,na.rm=T)))
abline(h=out_bs$IC0,lwd=3,lty='dashed')

#: Naive OS
system.time(out_nos <- nOS_MCP(y,Xsin,Xcos,L=4))
plot(y)
abline(v=out_nos,col=3,lty='dashed',lwd=2)
out_nos

# #: OS
# system.time(out_os <- aOS_MCP(y,Xsin,Xcos,L=4,thres=2))
# abline(v=out_os,col=4,lty='dashed',lwd=2)
# out_os

segment_fit = fit_segments(y,Xsin,Xcos,L=4,cp=out_nos,pl=TRUE)



##############################################
#: Multiple change point (Beniamino's data) :#
##############################################

y = read.table("AutoNOM-Beniamino/Data/data1.txt")$V1

plot(y)
abline(v=300,col=1,lty='dashed',lwd=2)
abline(v=650,col=1,lty='dashed',lwd=2)

#: Set omega grid and define covariates
n = length(y)
p = periodogram(y,plot=T)
om_grid = sort(p$freq[p$spec>10])
x = seq(1,n)
Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

#: Bin Segm
system.time(out_bs <- BS_MCP(y,Xsin,Xcos,L=4))
plot(y)
abline(v=out_bs,col=2,lty='dashed',lwd=2)
out_bs

#: Naive OS
system.time(out_nos <- nOS_MCP(y,Xsin,Xcos,L=4,nu=1/3))
plot(y)
abline(v=out_nos,col=3,lty='dashed',lwd=2)
out_nos

# #: OS
# system.time(out_os <- aOS_MCP(y,Xsin,Xcos,L=4))
# abline(v=out_os,col=4,lty='dashed',lwd=2)
# out_os

segment_fit = fit_segments(y,Xsin,Xcos,L=4,cp=out_nos,pl=TRUE)

