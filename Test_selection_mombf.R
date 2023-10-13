rm(list=ls())
library(CPSuSiE)
library(TSA)
library(mombf)

#: Frequencies (same)
om = 0.032
om2 = 0.086

beta1 = 1.4
beta2 = 0.8

#: True n=100,200,500
n = 200
x = seq(1,n)

set.seed(1234)
y = beta1*(sin(om*2*pi*x)+cos(om*2*pi*x)) +
  beta2*(sin(om2*2*pi*x)+cos(om2*2*pi*x)) + rnorm(length(x))

plot(y)

#: True omegas not in the grid 
p = periodogram(y,plot=T)
p_sort = p$freq[sort.int(p$spec,decreasing=T,index.return=T)$ix]
n_om = 50
om_grid = p_sort[1:n_om]

Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos
X = Xsin+Xcos

R1 = cor(Xsin)
R2 = cor(Xcos)

system.time(fit <- SuSiE_group_refined(y,Xsin,Xcos,L=4,f=1,R1,R2))

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

priorCoef = momprior(tau=0.348)
priorDelta = modelbbprior(alpha.p=1,beta.p=1)
system.time(fit1 <- modelSelection(y=y, x=X, center=FALSE, scale=FALSE,
                       priorCoef=priorCoef, priorDelta=priorDelta))
idmombf = which(fit1$margpp > 0.5)
abline(v = om_grid[idmombf],col=4,lwd=2,lty="dashed")

#: Add the true omegas 
om_grid = sort(c(om,om2,p_sort[1:n_om]))

Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos
X = Xsin+Xcos

R1 = cor(Xsin)
R2 = cor(Xcos)

fit = SuSiE_group_refined(y,Xsin,Xcos,L=4,f=1,R1,R2)

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

priorCoef = momprior(tau=0.348)
priorDelta = modelbbprior(alpha.p=1,beta.p=1)
fit1 = modelSelection(y=y, x=X, center=FALSE, scale=FALSE,
                      priorCoef=priorCoef, priorDelta=priorDelta)
idmombf = which(fit1$margpp > 0.5)
abline(v = om_grid[idmombf],col=4,lwd=2,lty="dashed")

