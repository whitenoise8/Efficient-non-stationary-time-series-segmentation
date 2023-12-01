rm(list=ls())
library(CPSuSiE)
library(BayesSpec)
library(astsa)

######################
#: Simulation study :#
######################

#: Setting
n = 1000
ncp = 3

load(paste0("data_sim/setting_n",n,"_ncp",ncp,".RData"))

Ndraws = 5000
Nburn = 1000
sample_id = (Nburn+1):Ndraws

#: Collect results
cov_mat = matrix(NA,Nsim,3)
bias_mat = matrix(NA,Nsim,3)
mse_mat = matrix(NA,Nsim,3)
time_mat = matrix(NA,Nsim,3)
autoParm_list = ada_list7 = ada_list15 = list()

methods = c("autoParm","adaptspec-J7","adaptspec-J15")
colnames(cov_mat) = colnames(bias_mat) = colnames(mse_mat) = colnames(time_mat) = methods

#: Params

for (it in 1:Nsim) {
  
  y = read.table(paste0("./data_sim/setting_n",n,"_ncp",ncp,"/data",it,".txt"),header=T)$x
  s = read.table(paste0("./data_sim/setting_n",n,"_ncp",ncp,"/signal",it,".txt"),header=T)$x
  
  cat("----------------------- AUTO-PARM ---------------------------\n")
  
  time_mat[it,1] = system.time(fit_ap <- autoParm(y))[3]
  
  cp_hat = fit_ap$breakpoints
  cp_hat = cp_hat[cp_hat > 1]
  cp_hat = cp_hat[cp_hat < n]
  
  cp_true = cp[it,]
  cov_mat[it,1] = covering_metric(cp_true,cp_hat,n)
  bias_mat[it,1] = abs(ncp-length(cp_hat))
  mse_mat[it,1] = mean((fit_autoParm(y,fit_ap,pl=FALSE)-s)^2)
  
  autoParm_list[[it]] = list(model = fit_ap)
  
  cat("----------------------- ADAPTSPEC-nBasis=7 ---------------------------\n")
  
  time_mat[it,2] = system.time(fit_ada <- adaptspec(Ndraws, Nburn, 4, y, nbasis=7))[3]
  
  nseg_ada = round(mean(fit_ada$nexp_curr[sample_id]))
  ncp_ada = nseg_ada - 1
  if (ncp_ada > 0) {
    cp_hat = round(rowMeans(fit_ada$xi[[nseg_ada]][1:ncp_ada,sample_id]))
  } else {
    cp_hat = NULL
  }
  
  cp_true = cp[it,]
  cov_mat[it,2] = covering_metric(cp_true,cp_hat,n)
  bias_mat[it,2] = abs(ncp-length(cp_hat))
  mse_mat[it,2] = NA
  
  if (nseg_ada > 1) {
    lspect = apply(fit_ada$log_spec_hat[[nseg_ada]][,,sample_id],c(1,2),mean)
  } else {
    lspect = rowMeans(fit_ada$log_spec_hat[[nseg_ada]][,1,sample_id])
  }
  ada_list7[[it]] = list(cp_hat = cp_hat,
                         lspect = lspect)
  
  cat("----------------------- ADAPTSPEC-nBasis=15 ---------------------------\n")
  
  time_mat[it,3] = system.time(fit_ada <- adaptspec(Ndraws, Nburn, 4, y, nbasis=15))[3]
  
  nseg_ada = round(mean(fit_ada$nexp_curr[sample_id]))
  ncp_ada = nseg_ada - 1
  if (ncp_ada > 0) {
    cp_hat = round(rowMeans(fit_ada$xi[[nseg_ada]][1:ncp_ada,sample_id]))
  } else {
    cp_hat = NULL
  }
  
  cp_true = cp[it,]
  cov_mat[it,3] = covering_metric(cp_true,cp_hat,n)
  bias_mat[it,3] = abs(ncp-length(cp_hat))
  mse_mat[it,3] = NA
  
  if (nseg_ada > 1) {
    lspect = apply(fit_ada$log_spec_hat[[nseg_ada]][,,sample_id],c(1,2),mean)
  } else {
    lspect = rowMeans(fit_ada$log_spec_hat[[nseg_ada]][,1,sample_id])
  }
  ada_list15[[it]] = list(cp_hat = cp_hat,
                          lspect = lspect)
  
  cat("-----------------------------------\n")
  cat("Simulation n°:",it,"/",Nsim," done!\n")
  cat("-----------------------------------\n")
}

save(list=c("cov_mat","bias_mat","mse_mat","time_mat"),
     file=paste0("results/results_other_setting_n",n,"_ncp",ncp,".RData"))
save(list=c("autoParm_list","ada_list7","ada_list15"),
     file=paste0("results/results_otherModels_setting_n",n,"_ncp",ncp,".RData"))

