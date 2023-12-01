rm(list=ls())

#: Define the setting and the signal 
signal = function(x,beta1,beta2,omega) {
  n = length(x)
  M = length(omega)
  S = matrix(NA,n,M)
  for (j in 1:M) S[,j] = beta1[j]*sin(2*pi*x*omega[j])+beta2[j]*cos(2*pi*x*omega[j])
  rowSums(S)
}

setting = list()
setting[[1]] = matrix(c(2.0, 3.5, 1/30), 1, 3, byrow=TRUE)
setting[[2]] = matrix(c(4.0, 3.0, 1/12), 1, 3, byrow=TRUE)

setting[[3]] = matrix(c(1.5, 2.0, 1/35,
                        2.5, 3.0, 1/20), 2, 3, byrow=TRUE)
setting[[4]] = matrix(c(2.5, 4.0, 1/22,
                        4.0, 2.0, 1/15), 2, 3, byrow=TRUE)

setting[[5]] = matrix(c(3.0, 2.0, 1/40,
                        2.0, 4.0, 1/20,
                        1.0, 1.5, 1/10), 3, 3, byrow=TRUE)
setting[[6]] = matrix(c(2.0, 3.0, 1/24,
                        4.0, 5.0, 1/15,
                        1.0, 2.5, 1/7), 3, 3, byrow=TRUE)

colnames(setting[[1]]) = colnames(setting[[2]]) = colnames(setting[[3]]) = c("beta_1","beta_2","omega")
colnames(setting[[4]]) = colnames(setting[[5]]) = colnames(setting[[6]]) = c("beta_1","beta_2","omega")

par(mfrow=c(3,2),mar=c(2,2,2,2))
for (j in 1:6) {
  plot(signal(1:200,setting[[j]][,"beta_1"],setting[[j]][,"beta_2"],setting[[j]][,"omega"]),type='l',lwd=2)
}
par(mfrow=c(1,1),mar=c(2,2,2,2))

library(ggplot2)
dfpl = data.frame(x=rep(1:200,6),
                  signal=rep(1:6,each=200))
y = NULL
for (j in 1:6) y = c(y,signal(1:200,setting[[j]][,"beta_1"],setting[[j]][,"beta_2"],setting[[j]][,"omega"]))
dfpl$y = y

signal_names <- list(
  '1'=expression(paste(beta[sin] == 2.0,", ",beta[cos] == 3.5,", ",omega == 1/30)),
  '2'=expression(paste(beta[sin] == 4.0,", ",beta[cos] == 3.0,", ",omega == 1/12)),
  '3'=expression(paste(beta[sin],"=(1.5,2.5), ",beta[cos],"=(2.0,3.0), ",omega,"=(1/35,1/20)")),
  '4'=expression(paste(beta[sin],"=(2.5,4.0), ",beta[cos],"=(4.0,2.0), ",omega,"=(1/22,1/25)")),
  '5'=expression(paste(beta[sin],"=(3.0,2.0,1.0), ",beta[cos],"=(2.0,4.0,1.5), ",omega,"=(1/40,1/20,1/10)")),
  '6'=expression(paste(beta[sin],"=(2.0,4.0,1.0), ",beta[cos],"=(3.0,5.0,2.5), ",omega,"=(1/24,1/15,1/7)"))
  )
signal_labeller = function(variable,value){
  return(signal_names[value])
}

ggplot(data=dfpl,aes(x=x,y=y)) + theme_bw() +
    geom_line(size=1) + xlab("") + ylab(expression(f(t,beta,omega))) +
    theme(text=element_text(size=15)) + 
  facet_wrap(signal ~ ., labeller=signal_labeller,ncol=2,scales='free')

ggsave("Figures/sim_signal.pdf",height=7,width=10)


#: Create data

n = 100
ncp = 1
nseg = ncp + 1
minLen = 20
minLen*nseg

sigma = 3

Nsim = 50
seeds = sample(1:1e4,Nsim)

cp = matrix(NA,Nsim,ncp)
id = matrix(NA,Nsim,nseg)
sets = matrix(NA,Nsim,ncp+2)

for (it in 1:Nsim) {
  set.seed(seeds[it])
  
  id[it,1:nseg] = 1
  while (any(diff(id[it,1:nseg])==0)) id[it,1:nseg] = sample(1:6,nseg,replace=T)
  
  tau = rep(2,ncp)
  while (min(diff(c(1,tau,n))) < minLen) tau = sort(sample(1:n,ncp))
  cp[it,1:ncp] = tau
  sets[it,1:(nseg+1)] = c(0,tau,n)
  
  knots = c(0,tau,n)
  y = NULL
  s = NULL
  for (i in 1:nseg) {
    x = seq((knots[i]+1),knots[i+1])
    s_i = signal(x,setting[[id[it,i]]][,"beta_1"],setting[[id[it,i]]][,"beta_2"],setting[[id[it,i]]][,"omega"])
    s = c(s,s_i)
    y = c(y,s_i+rnorm(length(x),0,sigma))
  }
  
  write.table(y,file=paste0("data_sim/setting_n",n,"_ncp",ncp,"/data",it,".txt"),row.names=F)
  write.table(s,file=paste0("data_sim/setting_n",n,"_ncp",ncp,"/signal",it,".txt"),row.names=F)
}

save(list=c("setting","n","Nsim","minLen","seeds","ncp","nseg","id","cp","sets"),file=paste0("data_sim/setting_n",n,"_ncp",ncp,".RData"))

par(mfrow=c(2,1),mar=c(2,2,2,2))
plot(y,type='l')
abline(v=tau,lty='dashed',lwd=2,col=2)
plot(s,type='l')
abline(v=tau,lty='dashed',lwd=2,col=2)
par(mfrow=c(1,1),mar=c(2,2,2,2))
