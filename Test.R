rm(list=ls())
library(CPSuSiE)
library(TSA)

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
om_grid = sort(p$freq[p$spec>100])
x = seq(1,n)
Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos


BS_MCP(y,Xsin,Xcos,om_grid,L=4,f=1,w=10,cost=1)

f = 1
cp_hat = nOS_MCP(y,Xsin,Xcos,om_grid,L=4,f=f,cost=1,nu=1/2,w=10)
fit = fit_segments(y,Xsin,Xcos,L=4,f=f,cp=cp_hat,om_grid=om_grid,pl=TRUE)
