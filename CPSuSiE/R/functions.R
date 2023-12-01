###########################
#-- Auxiliary functions --#
###########################

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
    
    if (min(abs(c(R1[CS,CS],R2[CS,CS]))) < thres) CS = NULL
    
    CS_list[[l]] = CS
  }
  
  CS_list
}


split_elbo = function(y,Xsin,Xcos,L,cp,f=1,sigma2=0){
  
  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  elbo = rep(NA,n_cp+1)
  
  for (i in 1:(n_cp+1)) {
    segm = (pts[i]+1):pts[i+1]
    yW = y[segm]
    XsinW = matrix(Xsin[segm,],nrow=length(yW))
    XcosW = matrix(Xcos[segm,],nrow=length(yW))
    
    fitW = SuSiE_group(yW,XsinW,XcosW,L,f=f,sigma2=sigma2) 
    elbo[i] = fitW$elbo
  }
  
  sum(elbo)
}

split_IC = function(y,Xsin,Xcos,L,cp,IC_name,f=1,sigma2=0){
  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  if (n_cp > 0) {
    pi = llik = rep(NA,n_cp+1)
    
    for (i in 1:(n_cp+1)) {
      segm = (pts[i]+1):pts[i+1]
      yW = y[segm]
      XsinW = matrix(Xsin[segm,],nrow=length(yW))
      XcosW = matrix(Xcos[segm,],nrow=length(yW))
      
      fitW = SuSiE_group(yW,XsinW,XcosW,L,f=f,sigma2=sigma2) 
      if (fitW$sigma2 <= 0) fitW$sigma2 = var(y)
      
      pi[i] =  2*sum(1-apply(fitW$alpha,1,function(x) prod(1-x))) + 1
      llik[i] = sum(dnorm(yW,fitW$yhat,sqrt(fitW$sigma2),log=T))
    }
    p = sum(pi)
    
  } else {
    segm = 1:n
    yW = y[segm]
    XsinW = matrix(Xsin[segm,],nrow=length(yW))
    XcosW = matrix(Xcos[segm,],nrow=length(yW))
    
    fitW = SuSiE_group(yW,XsinW,XcosW,L,f=f,sigma2=sigma2) 
    if (fitW$sigma2 <= 0) fitW$sigma2 = var(y)
    
    pi = 2*sum(1-apply(fitW$alpha,1,function(x) prod(1-x))) + 1
    p = pi
    llik = sum(dnorm(yW,fitW$yhat,sqrt(fitW$sigma2),log=T))
  }
  
  if (IC_name == 'MDL') IC = -sum(llik) + log(n_cp) + (n_cp+1)*log(n) + sum(log(pi)) + sum((pi/2+1)*log(diff(pts)))
  if (IC_name == 'BIC') IC = -sum(llik) + p/2*log(n) 
  if (IC_name == 'mBIC') IC = -sum(llik) + p/2*log(n) + 1/2*sum(pi*log(diff(pts/n)))
  
  IC
}


fit_segments = function(y,Xsin,Xcos,L,cp,f=1,sigma2=0,pl=TRUE){
  
  out = list()
  yhat = NULL

  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  for (i in 1:(n_cp+1)) {
    segm = (pts[i]+1):pts[i+1]
    
    yW = y[segm]
    XsinW = Xsin[segm,] 
    XcosW = Xcos[segm,]
    
    fitW = SuSiE_group(yW,XsinW,XcosW,L,f=f,sigma2=sigma2)

    out[[i]] = fitW
    yhat = c(yhat,fitW$yhat)
  }
  
  if (pl) {
    plot(y)
    abline(v=cp,lty='dashed',lwd=2)
    lines(yhat,col=2,lwd=2)
  }
  
  list(models = out,
       signal_estimate = yhat)
}

fit_autoParm = function(y,fit,pl=TRUE){
  
  yhat = NULL
  
  n = length(y)
  cp = fit$breakpoints[-c(1,length(fit$breakpoints))]
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  for (i in 1:(n_cp+1)) {
    segm = (pts[i]+1):pts[i+1]
    
    fitW = arima(y[segm],order=c(fit$segment_AR_orders[i],0,0))
    yhat = c(yhat,y[segm]-fitW$residuals)
  }
  
  if (pl) {
    plot(y)
    abline(v=cp,lty='dashed',lwd=2)
    lines(yhat,col=2,lwd=2)
  }
  
  yhat
}

covering_metric = function(cp_true,cp_hat,n) {
  ncp = length(cp_true)
  nseg = ncp+1
  sets = c(0,cp_true,n)
  
  ncp_hat = length(cp_hat)
  nseg_hat = ncp_hat+1
  sets_hat = c(0,cp_hat,n)
  
  C = 0
  for (i in 1:nseg) {
    A = (sets[i]+1):sets[i+1]
    
    J = rep(NA,nseg_hat)
    for (j in 1:nseg_hat) {
      B = (sets_hat[j]+1):sets_hat[j+1]
      J[j] = length(intersect(A,B))/length(union(A,B))
    }
    C = C + 1/n*length(A)*max(J)
  }
  
  C
}

#########################
#-- Optimistic search --#
#########################

nOS_internal = function(y,Xsin,Xcos,L,lt,t,rt,l,r,nu=1/2,f=1,sigma2=0){
  
  n = nrow(Xsin)
  
  if ((rt-lt) <= 5){
    elbo = c()
    if (lt == l) lt = l+1
    if (rt == r) rt = r-1
    for (i in lt:rt){
      elbo = c(elbo,split_elbo(y,Xsin,Xcos,L,cp=i,f=f,sigma2=sigma2))
    }
    cp = which.max(elbo)+lt-1
    out = list(cp=cp,
               elbo=max(elbo,na.rm=T))
    return(out)
  } else {
    
    if ((rt-t) > (t-lt)) {
      w = floor(rt - (rt-t)*nu)
      elbo_w = split_elbo(y,Xsin,Xcos,L,cp=w,f=f,sigma2=sigma2)
      elbo_t = split_elbo(y,Xsin,Xcos,L,cp=t,f=f,sigma2=sigma2)
      if (elbo_w > elbo_t) { nOS_internal(y,Xsin,Xcos,L,t,w,rt,l,r,nu,f=f,sigma2=sigma2) } else { nOS_internal(y,Xsin,Xcos,L,lt,t,w,l,r,nu,f=f,sigma2=sigma2) }
    } else {
      w = floor(lt + (t-lt)*nu)
      elbo_w = split_elbo(y,Xsin,Xcos,L,cp=w,f=f,sigma2=sigma2)
      elbo_t = split_elbo(y,Xsin,Xcos,L,cp=t,f=f,sigma2=sigma2)
      if (elbo_w > elbo_t) { nOS_internal(y,Xsin,Xcos,L,lt,w,t,l,r,nu,f=f,sigma2=sigma2) } else { nOS_internal(y,Xsin,Xcos,L,w,t,rt,l,r,nu,f=f,sigma2=sigma2) }
    }
    
  }
}

nOS = function(y,Xsin,Xcos,L,l,r,nu=1/2,f=1,sigma2=0){
  
  n = length(y)
  t = floor((l+nu*r/(1+nu)))
  out_nOS = nOS_internal(y,Xsin,Xcos,L=L,
                         lt=l,t=t,rt=r,l=l,r=r,nu=nu,f=f,sigma2=sigma2)
  
  cp_hat = NA
  
  fit0 = SuSiE_group(y,Xsin,Xcos,L,f=f,sigma2=sigma2)
  
  if (out_nOS$elbo > fit0$elbo)  cp_hat = out_nOS$cp
  
  c(cp_hat,out_nOS$elbo-fit0$elbo)
}

getCandidateSet = function(y,Xsin,Xcos,L,f=1,sigma2=0,nu=1/2,w=NULL,Trace=0) {
  n = length(y)
  if (is.null(w)) w = sqrt(n)
  
  it = 1
  if (Trace == 1) cat(paste0("* Searching candidate CP...\n"))
  eval_sgm = matrix(c(1,n),1,2)
  
  os = nOS(y,Xsin,Xcos,L,l=1,r=n,nu,f=f,sigma2=sigma2)
  cp_hat = os[1]
  elbo_hat = os[2]
  
  cp_set = NULL
  cp_list = NULL
  if (!is.na(cp_hat)) {
    cp_set = cp_hat
    cp_list = list(cp_hat)
    s_set = 1
  } else {
    s_set = 0
  }
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      it = it + 1
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      elbs = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        
        if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
          eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
          
          os = nOS(y[x],Xsin[x,],Xcos[x,],L,1,length(x),nu,f=f,sigma2=sigma2)
          cp_hat = os[1] + pts[k]
          elbo_hat = os[2]
          
          if (!is.na(cp_hat)) if(min(abs(cp_set-cp_hat)) > w) if (cp_hat-w > 1) if (cp_hat+w < n) {
            cps[k] = cp_hat 
            elbs[k] = elbo_hat
          } 
        }
      }
      
      elbs = elbs[!is.na(cps)]
      cps = cps[!is.na(cps)]
      cp_set = unique(sort(c(cp_set,cps)))
      if (length(cps) > 0) cp_list[[it]] = cps[sort.int(elbs,decreasing=T,index.return=T)$ix]
      
      if (length(cp_set) > s_set) {
        s_set = length(cp_set)
      } else {
        conv = 1
      }
      
    }
  }
  
  if (Trace == 1) cat(paste0("Candidate CP set: {",paste(cp_set,collapse=','),"}\n"))
  cp_list
}


pruneSet = function(cp_set,y,Xsin,Xcos,L,f=1,sigma2=0,IC="MDL",option="joint",Trace=0) {
  
  if (length(cp_set) == 0) out = NULL
  
  if (length(cp_set) > 0) {
    if (option == "joint") {
      cp_set = unlist(cp_set)
      if (length(cp_set) > 8) cp_set = cp_set[1:8]
      
      n = length(y)
      
      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      if (length(cp_set) > 1) {
        cp_mat = as(expand.grid(rep(list(0:1),length(cp_set))),"matrix")
        cp_mat = cp_mat[-1,]
        IC_star = rep(NA,nrow(cp_mat))
        for (i in 1:nrow(cp_mat)) {
          IC_star[i] = split_IC(y,Xsin,Xcos,L,cp=sort(cp_set[cp_mat[i,]==1]),IC,f=f,sigma2=sigma2)
        }
        
        out = sort(cp_set[cp_mat[which.min(IC_star),]==1])
      } else {
        out = sort(cp_set)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
    if (option == 'conditional') {
      n = length(y)
      
      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      lev = length(cp_set)
      cp_hat = cp_set[[1]]
      IC_hat = split_IC(y,Xsin,Xcos,L,cp=cp_hat,IC,f=f,sigma2=sigma2)
      
      if (lev > 1) {
        conv = 0
        lv = 1
        while (conv == 0) {
          lv = lv + 1
          cps = cp_set[[lv]]
          for (i in 1:length(cps)) {
            cp_star = unique(sort(c(cp_hat,cps[1:i])))
            IC_star = split_IC(y,Xsin,Xcos,L,cp=cp_star,IC,f=f,sigma2=sigma2)
            if (IC_star < IC_hat) {
              cp_hat = cp_star
              IC_hat = IC_star
            } else {
              conv = 1
              break
            }
          }
          if (lv == lev) conv = 1
        }
        
        out = sort(cp_hat)
      } else {
        out = sort(cp_hat)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
  }
  
  out
}

