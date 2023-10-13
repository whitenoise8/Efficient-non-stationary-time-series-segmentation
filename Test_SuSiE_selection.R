rm(list=ls())
library(CPSuSiE)
library(TSA)

getX = function(x,om){
  
  n = length(x)
  Xsin = matrix(0,nrow=n,ncol=length(om))
  Xcos = matrix(0,nrow=n,ncol=length(om))
  
  for (j in 1:length(om)){
    Xsin[,j] = sin(om[j]*2*pi*x)
    Xcos[,j] = cos(om[j]*2*pi*x)
  }
  
  list(Xsin=Xsin, Xcos=Xcos)
}

getCS = function(fit,R1,R2,rho_CS,thres){
  CS_list = list()
  
  k_0 = apply(fit$alpha,2,function(x) min(which(cumsum(sort(x,decreasing=T)) > rho_CS)))
  
  for (l in 1:ncol(fit$alpha)) {
    CS = sort.int(fit$alpha[,l],index.return=T,decreasing=T)$ix[1:k_0[l]]
    
    if (min(abs(c(R1[CS,CS],R2[CS,CS]))) < thres) {
      fit$alpha[,l] = 0
      CS = NULL
    }
    
    CS_list[[l]] = CS
  }
  
  list(CS=CS_list,fit=fit)
}



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

#: 1) Equi-spaced grid
p = periodogram(y,plot=T)
p_sort = p$freq[sort.int(p$spec,decreasing=T,index.return=T)$ix]
n_om = 10
om_grid = p_sort[1:n_om]

om_grid = seq(1e-3,0.5,length.out=1000)

Xout = getX(x,om_grid)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

fit = SuSiE_group(y,Xsin,Xcos,L=4,f=1)

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a,xlim=c(0.07,0.1))
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

R1 = cor(Xsin)
R2 = cor(Xcos)

out = getCS(fit,R1,R2,0.95,0.5)
fit = out$fit

a = 1-apply(fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

L_hat = length(out$CS)
om_hat = rep(NA,L_hat)
for (l in 1:L_hat) {
  wei = a[out$CS[[l]]]/sum(a[out$CS[[l]]])
  om_hat[l] = sum(wei*om_grid[out$CS[[l]]])
}
om_hat

Xout = getX(x,om_hat)
Xsin = Xout$Xsin
Xcos = Xout$Xcos

fit = SuSiE_group(y,Xsin,Xcos,L=L_hat)

plot(c(om,om2,om_hat),c(beta_sin,beta2_sin,fit$B1),type='n')
points(c(om,om2),c(beta_sin,beta2_sin),pch=16,col=2)
points(om_hat,fit$B1,pch=16)

plot(c(om,om2,om_hat),c(beta_cos,beta2_cos,fit$B2),type='n')
points(c(om,om2),c(beta_cos,beta2_cos),pch=16,col=2)
points(om_hat,fit$B2,pch=16)

fit$sigma2

plot(y)
lines(fit$yhat,col=2)



R1 = cor(Xsin)
R2 = cor(Xcos)

fit1 = SuSiE_group(y,Xsin,Xcos,L=4,f=1)
fit2 = SuSiE_group_refined(y,Xsin,Xcos,L=4,f=1,R1,R2,om_grid)

a = 1-apply(fit1$alpha,1,function(x) prod(1-x))
plot(om_grid,a)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

a = 1-apply(fit2$outCS$fit$alpha,1,function(x) prod(1-x))
plot(om_grid,a)
abline(v=c(om,om2),col=2,lwd=2,lty='dashed')

plot(fit1$B1,fit2$B1)
abline(0,1)
plot(fit1$B2,fit2$B2)
abline(0,1)
c(fit1$sigma2,fit2$sigma2)
c(fit1$elbo,fit2$elbo)



