#: Individual Checks
rm(list=ls())
library(CPSuSiE)

n = 500
ncp = 2

load(paste0("data_sim/setting_n",n,"_ncp",ncp,".RData"))

cov_autoNom = bias_autoNom = mse_autoNom = time_autoNom = matrix(NA,Nsim,3)
colnames(cov_autoNom) = colnames(bias_autoNom) = colnames(mse_autoNom) = colnames(time_autoNom) = c("AutoNOM_lam_0.1","AutoNOM_lam_1","AutoNOM_lam_10")

for (lam_id in 1:3) {
  
  cp_hat = as(read.table(paste0("./res_julia/cp_setting_n",n,"_ncp",ncp,"_lambda_setting",lam_id,".txt"),header=F,sep=','),"matrix")
  time_hat = read.table(paste0("./res_julia/time_setting_n",n,"_ncp",ncp,"_lambda_setting",lam_id,".txt"),header=F,sep=',')$V1
  s_hat = as(read.table(paste0("./res_julia/signal_setting_n",n,"_ncp",ncp,"_lambda_setting",lam_id,".txt"),header=F,sep=','),"matrix")
  
  for (it in 1:Nsim) {
    
    s = read.table(paste0("./data_sim/setting_n",n,"_ncp",ncp,"/signal",it,".txt"),header=T)$x
    cp_true = cp[it,]
    cp_hat_auto = round(cp_hat[it,cp_hat[it,]!=0])
    
    cov_autoNom[it,lam_id] = covering_metric(cp_true,cp_hat_auto,n)
    bias_autoNom[it,lam_id] = abs(ncp-length(cp_hat_auto))
    mse_autoNom[it,lam_id] = mean((s_hat[it,]-s)^2)
    time_autoNom[it,lam_id] = time_hat[it]
    
    cat("-----------------------------------\n")
    cat("Simulation n°:",it,"/",Nsim," done!\n")
    cat("-----------------------------------\n")
  }
}

save(list=c("cov_autoNom","bias_autoNom","mse_autoNom","time_autoNom"),
     file=paste0("results/results_autoNOM_setting_n",n,"_ncp",ncp,".RData"))
