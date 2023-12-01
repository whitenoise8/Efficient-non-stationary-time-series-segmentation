#: Individual Checks
rm(list=ls())
library(ggplot2)

n = 500
ncp = 2
Nsim = 50


load(paste0("results/results_CPSuSiE_setting_n",n,"_ncp",ncp,".RData"))
load(paste0("results/results_autoNom_setting_n",n,"_ncp",ncp,".RData"))

cov_mat = cbind(cov_mat,cov_autoNom)
bias_mat = cbind(bias_mat,bias_autoNom)
mse_mat = cbind(mse_mat,mse_autoNom)
time_mat = cbind(time_mat,time_autoNom)

colMeans(cov_mat)
colMeans(bias_mat)
apply(bias_mat,2,function(x) mean(x==0))
colMeans(mse_mat)
colMeans(time_mat)

dfpl = data.frame(x=rep(colnames(cov_mat),each=Nsim),
                  y=as(cov_mat,"vector"))
ggplot(data=dfpl,aes(x=x,y=y)) + theme_bw() +
  geom_boxplot()

dfpl = data.frame(x=rep(colnames(cov_mat),each=Nsim),
                  y=as((time_mat),"vector"))
ggplot(data=dfpl,aes(x=x,y=y,fill=factor(x))) + theme_bw() +
  geom_boxplot()

