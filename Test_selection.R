#: This code compares the selection of the frequencies within a segment using SuSiE and
#: mombf package of David Rossell for Bayesian variable selection with non-local priors

rm(list=ls())
library(CPSuSiE)
library(TSA)
library(mombf)

#: Frequencies
om = 0.032
om2 = 0.086

beta1 = 1.4
beta2 = 0.8

n = 500
x = seq(1,n)

set.seed(1234)
y = beta1*(sin(om*2*pi*x)+cos(om*2*pi*x)) +
  beta2*(sin(om2*2*pi*x)+cos(om2*2*pi*x)) + rnorm(length(x))

plot(y)

#: True omegas not in the grid 
om_grid = seq(1e-3,0.5,length.out=100)

Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos
X = Xsin+Xcos

system.time(fit <- SuSiE_group(y,Xsin,Xcos,L=4,f=1,sigma2=0))

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

priorCoef = momprior(tau=0.348)
priorDelta = modelbbprior(alpha.p=1,beta.p=1)
system.time(fit1 <- modelSelection(y=y, x=X, center=FALSE, scale=FALSE,
                       priorCoef=priorCoef, priorDelta=priorDelta))
idmombf = which(fit1$margpp > 0.5)
plot(om_grid,fit1$margpp)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

