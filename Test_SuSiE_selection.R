#: This code explores the selection of the frequencies within a segment using SuSiE

rm(list=ls())
library(CPSuSiE)
library(TSA)

#: Frequencies (same)
om = 0.032
om2 = 0.086

beta_sin = 1.4
beta2_sin = 0.8

beta_cos = 2.0
beta2_cos = 1.2

n = 200
x = seq(1,n)

set.seed(123)
y = beta_sin*sin(om*2*pi*x) +
  beta2_sin*sin(om2*2*pi*x) +
  beta_cos*cos(om*2*pi*x) +
  beta2_cos*cos(om2*2*pi*x) + rnorm(length(x))

plot(y)

#: 1) Periodogram
p = periodogram(y,plot=T)
p_sort = p$freq[sort.int(p$spec,decreasing=T,index.return=T)$ix]
n_om = 100
om_grid = sort(p_sort[1:n_om])

Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

#: f coefficient for the fractional posterior
#: if sigma2 > 0 then sigma2 is fixed in SuSiE
fit = SuSiE_group(y,Xsin,Xcos,L=4,f=1,sigma2=0)

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a,xlim=c(0.07,0.1))
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

fit$sigma2

plot(y)
lines(fit$yhat,col=2,lwd=2)


#: 1) Equi-spaced grid
om_grid = seq(1e-3,0.5,length.out=1000)

Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

fit = SuSiE_group(y,Xsin,Xcos,L=4,f=1,sigma2=0)

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a,xlim=c(0.07,0.1))
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

fit$sigma2

plot(y)
lines(fit$yhat,col=2,lwd=2)



