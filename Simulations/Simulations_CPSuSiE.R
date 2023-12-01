rm(list=ls())
library(CPSuSiE)
library(TSA)

######################
#: Simulation study :#
######################

#: Setting
n = 5000
ncp = 8

load(paste0("data_sim/setting_n",n,"_ncp",ncp,".RData"))


#: Collect results
cov_mat = matrix(NA,Nsim,16)
bias_mat = matrix(NA,Nsim,16)
mse_mat = matrix(NA,Nsim,16)
time_mat = matrix(NA,Nsim,16)

methods = expand.grid(c("joint","conditional"),c("MDL","mBIC"),c("L1","L2"),c("grid-p","grid-eq"))
colnames(cov_mat) = colnames(bias_mat) = colnames(mse_mat) = colnames(time_mat) = apply(methods,1,function(x) paste(x,collapse='-'))

#: Params
n_om = 100
l_grid = 200

for (it in 1:Nsim) {
  
  y = read.table(paste0("./data_sim/setting_n",n,"_ncp",ncp,"/data",it,".txt"),header=T)$x
  s = read.table(paste0("./data_sim/setting_n",n,"_ncp",ncp,"/signal",it,".txt"),header=T)$x
  sd_y = sd(y)
  y_std = y

  #: Set omega grid and define covariates
  p = periodogram(y_std,plot=F)
  p_sort = p$freq[sort.int(p$spec,decreasing=T,index.return=T)$ix]
  om_grid = sort(p_sort[1:n_om])

  x = seq(1,n)
  Xout = getX(x,om_grid)
  Xsin = Xout$Xsin
  Xcos = Xout$Xcos
  
  time_mat[it,1:4] = system.time(cp_set <- getCandidateSet(y_std,Xsin,Xcos,L=1,f=1,nu=1/2))[3]
  
  for (m in 1:4) {
    time_mat[it,m] = time_mat[it,m] + system.time(cp_hat <- pruneSet(cp_set,y_std,Xsin,Xcos,L=1,f=1,IC=methods[m,2],option=methods[m,1]))[3]
    
    ft = fit_segments(y,Xsin,Xcos,L=4,f=1,cp=cp_hat,pl=F)

    cp_true = cp[it,]
    cov_mat[it,m] = covering_metric(cp_true,cp_hat,n)
    bias_mat[it,m] = abs(ncp-length(cp_hat))
    mse_mat[it,m] = mean((ft$signal_estimate-s)^2)
  }
  
  time_mat[it,5:8] = system.time(cp_set <- getCandidateSet(y_std,Xsin,Xcos,L=2,f=1,nu=1/2))[3]
  
  for (m in 5:8) {
    time_mat[it,m] = time_mat[it,m] + system.time(cp_hat <- pruneSet(cp_set,y_std,Xsin,Xcos,L=2,f=1,IC=methods[m,2],option=methods[m,1]))[3]
    
    ft = fit_segments(y,Xsin,Xcos,L=4,f=1,cp=cp_hat,pl=F)
    
    cp_true = cp[it,]
    cov_mat[it,m] = covering_metric(cp_true,cp_hat,n)
    bias_mat[it,m] = abs(ncp-length(cp_hat))
    mse_mat[it,m] = mean((ft$signal_estimate-s)^2)
  }
  
  #: Set omega grid equally spaced
  om_grid = seq(1e-3,0.5,length.out=l_grid)
  
  x = seq(1,n)
  Xout = getX(x,om_grid)
  Xsin = Xout$Xsin
  Xcos = Xout$Xcos
  
  time_mat[it,9:12] = system.time(cp_set <- getCandidateSet(y_std,Xsin,Xcos,L=1,f=1,nu=1/2))[3]
  
  for (m in 9:12) {
    time_mat[it,m] = time_mat[it,m] + system.time(cp_hat <- pruneSet(cp_set,y_std,Xsin,Xcos,L=1,f=1,IC=methods[m,2],option=methods[m,1]))[3]
    
    ft = fit_segments(y,Xsin,Xcos,L=4,f=1,cp=cp_hat,pl=F)
    
    cp_true = cp[it,]
    cov_mat[it,m] = covering_metric(cp_true,cp_hat,n)
    bias_mat[it,m] = abs(ncp-length(cp_hat))
    mse_mat[it,m] = mean((ft$signal_estimate-s)^2)
  }
  
  time_mat[it,13:16] = system.time(cp_set <- getCandidateSet(y_std,Xsin,Xcos,L=2,f=1,nu=1/2))[3]
  
  for (m in 13:16) {
    time_mat[it,m] = time_mat[it,m] + system.time(cp_hat <- pruneSet(cp_set,y_std,Xsin,Xcos,L=2,f=1,IC=methods[m,2],option=methods[m,1]))[3]
    
    ft = fit_segments(y,Xsin,Xcos,L=4,f=1,cp=cp_hat,pl=F)
    
    cp_true = cp[it,]
    cov_mat[it,m] = covering_metric(cp_true,cp_hat,n)
    bias_mat[it,m] = abs(ncp-length(cp_hat))
    mse_mat[it,m] = mean((ft$signal_estimate-s)^2)
  }
  
  cat("-----------------------------------\n")
  cat("Simulation n°:",it,"/",Nsim," done!\n")
  cat("-----------------------------------\n")
}

save(list=c("cov_mat","bias_mat","mse_mat","time_mat"),
     file=paste0("results/results_CPSuSiE_setting_n",n,"_ncp",ncp,".RData"))

