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

SuSiE_group_refined = function(y,Xsin,Xcos,L,f,R1,R2) {
  fit = SuSiE_group(y,Xsin,Xcos,L,f)

  outCS = getCS(fit,R1,R2,0.95,0.5)
  
  L_hat = length(outCS)
  fit = SuSiE_group(y,Xsin,Xcos,L=L_hat)
  fit$L_hat = L_hat
  fit
}


split_elbo = function(y,Xsin,Xcos,L,f,R1,R2,cp){
  
  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  elbo = rep(NA,n_cp+1)
  
  for (i in 1:(n_cp+1)) {
    segm = pts[i]:pts[i+1]
    yW = y[segm]
    XsinW = matrix(Xsin[segm,],nrow=length(yW))
    XcosW = matrix(Xcos[segm,],nrow=length(yW))
    
    fitW = SuSiE_group_refined(yW,XsinW,XcosW,L,f,R1,R2) 
    elbo[i] = fitW$elbo
  }
  
  sum(elbo)
}

multiple_split_IC = function(y,Xsin,Xcos,L,f,cp,cost,R1,R2){
  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  if (n_cp > 0) {
    p = llik = rep(NA,n_cp+1)
    
    for (i in 1:(n_cp+1)) {
      segm = (pts[i]+1):pts[i+1]
      yW = y[segm]
      XsinW = matrix(Xsin[segm,],nrow=length(yW))
      XcosW = matrix(Xcos[segm,],nrow=length(yW))
      
      fitW = SuSiE_group_refined(yW,XsinW,XcosW,L,f,R1,R2) 
      if (fitW$sigma2 <= 0) fitW$sigma2 = var(y)
      
      p[i] =  2*fitW$L_hat + 1
      llik[i] = sum(dnorm(yW,fitW$yhat,sqrt(fitW$sigma2),log=T))
    }
    p = sum(p)
    
  } else {
    segm = 1:n
    yW = y[segm]
    XsinW = matrix(Xsin[segm,],nrow=length(yW))
    XcosW = matrix(Xcos[segm,],nrow=length(yW))
    
    fitW = SuSiE_group_refined(yW,XsinW,XcosW,L,f,R1,R2) 
    if (fitW$sigma2 <= 0) fitW$sigma2 = var(y)
    
    p = 2*fitW$L_hat + 1
    llik = sum(dnorm(yW,fitW$yhat,sqrt(fitW$sigma2),log=T))
  }
  
  IC = -2*sum(llik) + (n_cp+1)*p*log(n) + cost*sum((diff(pts)/n-1/(n_cp+1))^2)*log(n)
  IC
}


fit_segments = function(y,Xsin,Xcos,L,f,cp,pl=TRUE){
  
  out = list()
  yhat = NULL
  
  R1 = cor(Xsin)
  R2 = cor(Xcos)

  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  for (i in 1:(n_cp+1)) {
    segm = (pts[i]+1):pts[i+1]
    
    yW = y[segm]
    XsinW = Xsin[segm,] 
    XcosW = Xcos[segm,]
    
    fitW = SuSiE_group_refined(yW,XsinW,XcosW,L,f,R1,R2)
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

###########################
#-- Binary segmentation --#
###########################

BS = function(y,Xsin,Xcos,L,f,w,R1,R2) {
  
  n = length(y)
  cp_hat = NA
  elbo = rep(NA,n+1)
  
  fit0 = SuSiE_group_refined(y,Xsin,Xcos,L,f,R1,R2)
  elbo[1] = fit0$elbo
  
  for (i in w:(n-w)){
    yA = y[1:i]
    yB = y[i:n]
    
    XsinA = Xsin[1:i,]
    XcosA = Xcos[1:i,]
    XsinB = Xsin[i:n,]
    XcosB = Xcos[i:n,]
    
    elbo[i+1] = split_elbo(y,Xsin,Xcos,L,f,R1,R2,cp=i)
  }
  
  maxelbo = max(elbo,na.rm=TRUE)
  cp_star = which.max(elbo)-1
  
  if (cp_star > 0) cp_hat = cp_star
  
  c(cp_hat,maxelbo-elbo[1])
}


BS_MCP = function(y,Xsin,Xcos,L,f=1,cost=1,w=10) {
  n = length(y)
  
  R1 = cor(Xsin)
  R2 = cor(Xcos)
  
  it = 1
  cat(paste0("* Searching CP - Split stage ",it,":\n"))
  eval_sgm = matrix(c(1,n),1,2)
  
  bs = BS(y,Xsin,Xcos,L,f,w,R1,R2)
  cp_hat = bs[1]
  elbo_hat = bs[2]
  
  if (!is.na(cp_hat)) {
    cat(paste0("Segment: ",1," - ",n,"...candidate CP found: ",cp_hat,"\n"))
    IC0 = multiple_split_IC(y,Xsin,Xcos,L,f,cp=NULL,cost=cost,R1,R2)
    IC_hat = multiple_split_IC(y,Xsin,Xcos,L,f,cp=cp_hat,cost=cost,R1,R2)
    if (IC_hat < IC0) {
      cp_set = cp_hat
      s_set = 1
      cat(paste0("CP found via IC: ",cp_hat,".\n Current CP set: {",paste(cp_set,collapse=','),"}\n"))
    } else {
      cat(paste0("CP found via IC: None.\n Current CP set: {}\n"))
    }
  } else {
    cat(paste0("Segment: ",1," - ",n,"...candidate CP not found!\n"))
    s_set = 0
  }
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      it = it + 1
      cat(paste0("* Searching CP - Split stage ",it,":\n"))
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      elbs = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        
        if ((pts[k+1]-pts[k]) > (2*w+1)) {
          if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
            eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
            
            bs = BS(y[x],Xsin[x,],Xcos[x,],L,f,w,R1,R2)
            cp_hat = bs[1]
            elbo_hat = bs[2]
            
            if (!is.na(cp_hat)) {
              cat(paste0("Segment: ",pts[k], " - ",pts[k+1],"...candidate CP found: ",cp_hat + pts[k],"\n"))
              cps[k] = cp_hat + pts[k]
              elbs[k] = elbo_hat
            } else {
              cat(paste0("Segment: ",pts[k], " - ",pts[k+1],"...candidate CP not found!\n"))
            }
            
          }
        }
      }
      
      elbs = elbs[!is.na(cps)]
      cps = cps[!is.na(cps)]
      if (length(elbs) > 0) {
        idx = sort.int(elbs,decreasing=TRUE,index.return=TRUE)$ix
        
        IC_star = rep(NA,length(idx))
        cps_eval = NULL
        for (i in idx) {
          cps_eval = c(cps_eval,cps[i])
          cp_star = unique(sort(c(cp_set,cps_eval)))
          IC_star[i] = multiple_split_IC(y,Xsin,Xcos,L,f,cp=cp_star,cost=cost,R1,R2)
        }
        cp_hat = cps[IC_star < IC_hat]
        IC_hat = min(IC_star)
        
        cp_set = unique(sort(c(cp_set,cp_hat)))
        if (length(cp_hat) > 0) {
          cat(paste0("CP found via IC: ",cp_hat,".\n Current CP set: {",paste(cp_set,collapse=','),"}\n"))
        } else {
          cat(paste0("CP found via IC: None.\n Current CP set: {",paste(cp_set,collapse=','),"}\n"))
        }
      }
      
      if (length(cp_set) > s_set) {
        s_set = length(cp_set)
      } else {
        conv = 1
      }
      
    }
  }
  
  cp_set
}

#########################
#-- Optimistic search --#
#########################

nOS_internal = function(y,Xsin,Xcos,L,f,R1,R2,lt,t,rt,l,r,nu=1/2){
  
  n = nrow(Xsin)
  
  if ((rt-lt) <= 5){
    elbo = c()
    if (lt == l) lt = l+1
    if (rt == r) rt = r-1
    for (i in lt:rt){
      elbo = c(elbo,split_elbo(y,Xsin,Xcos,L,f,R1,R2,cp=i))
    }
    cp = which.max(elbo)+lt-1
    out = list(cp=cp,
               elbo=max(elbo))
    return(out)
  } else {
    
    if ((rt-t) > (t-lt)) {
      w = floor(rt - (rt-t)*nu)
      elbo_w = split_elbo(y,Xsin,Xcos,L,f,R1,R2,cp=w)
      elbo_t = split_elbo(y,Xsin,Xcos,L,f,R1,R2,cp=t)
      if (elbo_w > elbo_t) { nOS_internal(y,Xsin,Xcos,L,f,R1,R2,t,w,rt,l,r,nu) } else { nOS_internal(y,Xsin,Xcos,L,f,R1,R2,lt,t,w,l,r,nu) }
    } else {
      w = floor(lt + (t-lt)*nu)
      elbo_w = split_elbo(y,Xsin,Xcos,L,f,R1,R2,cp=w)
      elbo_t = split_elbo(y,Xsin,Xcos,L,f,R1,R2,cp=t)
      if (elbo_w > elbo_t) { nOS_internal(y,Xsin,Xcos,L,f,R1,R2,lt,w,t,l,r,nu) } else { nOS_internal(y,Xsin,Xcos,L,f,R1,R2,w,t,rt,l,r,nu) }
    }
    
  }
}

nOS = function(y,Xsin,Xcos,L,f,R1,R2,l,r,nu=1/2){
  
  n = length(y)
  t = floor((l+nu*r/(1+nu)))
  out_nOS = nOS_internal(y,Xsin,Xcos,L=L,f=f,R1=R1,R2=R2,
                         lt=l,t=t,rt=r,l=l,r=r,nu=nu)
  
  cp_hat = NA
  
  fit0 = SuSiE_group_refined(y,Xsin,Xcos,L,f,R1,R2)
  
  if (out_nOS$elbo > fit0$elbo)  cp_hat = out_nOS$cp
  
  c(cp_hat,out_nOS$elbo-fit0$elbo)
}


nOS_MCP = function(y,Xsin,Xcos,L,f=1,cost=1,nu=1/2,w=10) {
  n = length(y)
  
  R1 = cor(Xsin)
  R2 = cor(Xcos)
  
  it = 1
  cat(paste0("* Searching CP - Split stage ",it,":\n"))
  eval_sgm = matrix(c(1,n),1,2)
  
  os = nOS(y,Xsin,Xcos,L,f,R1,R2,l=1,r=n,nu)
  cp_hat = os[1]
  elbo_hat = os[2]
  
  if (!is.na(cp_hat)) {
    cat(paste0("Segment: ",1," - ",n,"...candidate CP found: ",cp_hat,"\n"))
    IC0 = multiple_split_IC(y,Xsin,Xcos,L,f,cp=NULL,cost=cost,R1,R2)
    IC_hat = multiple_split_IC(y,Xsin,Xcos,L,f,cp=cp_hat,cost=cost,R1,R2)
    if (IC_hat < IC0) {
      cp_set = cp_hat
      s_set = 1
      cat(paste0("CP found via IC: ",cp_hat,".\n Current CP set: {",paste(cp_set,collapse=','),"}\n"))
    } else {
      cat(paste0("CP found via IC: None.\n Current CP set: {}\n"))
    }
  } else {
    cat(paste0("Segment: ",1," - ",n,"...candidate CP not found!\n"))
    s_set = 0
  }
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      it = it + 1
      cat(paste0("* Searching CP - Split stage ",it,":\n"))
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      elbs = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        
        if ((pts[k+1]-pts[k]) > (2*w+1)) {
          if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
            eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
            
            os = nOS(y[x],Xsin[x,],Xcos[x,],L,f,R1,R2,1,length(x),nu)
            cp_hat = os[1]
            elbo_hat = os[2]
            
            if (!is.na(cp_hat)) {
              cat(paste0("Segment: ",pts[k], " - ",pts[k+1],"...candidate CP found: ",cp_hat + pts[k],"\n"))
              cps[k] = cp_hat + pts[k]
              elbs[k] = elbo_hat
            } else {
              cat(paste0("Segment: ",pts[k], " - ",pts[k+1],"...candidate CP not found!\n"))
            }
            
          }
        }
      }
      
      elbs = elbs[!is.na(cps)]
      cps = cps[!is.na(cps)]
      if (length(elbs) > 0) {
        idx = sort.int(elbs,decreasing=TRUE,index.return=TRUE)$ix
        
        IC_star = rep(NA,length(idx))
        cps_eval = NULL
        for (i in idx) {
          cps_eval = c(cps_eval,cps[i])
          cp_star = unique(sort(c(cp_set,cps_eval)))
          IC_star[i] = multiple_split_IC(y,Xsin,Xcos,L,f,cp=cp_star,cost=cost,R1,R2)
        }
        cp_hat = cps[IC_star < IC_hat]
        IC_hat = min(IC_star)
        
        cp_set = unique(sort(c(cp_set,cp_hat)))
        if (length(cp_hat) > 0) {
          cat(paste0("CP found via IC: ",cp_hat,".\n Current CP set: {",paste(cp_set,collapse=','),"}\n"))
        } else {
          cat(paste0("CP found via IC: None.\n Current CP set: {",paste(cp_set,collapse=','),"}\n"))
        }
      }
      
      if (length(cp_set) > s_set) {
        s_set = length(cp_set)
      } else {
        conv = 1
      }
      
    }
  }
  
  cp_set
}


#########################
#: Dynamic programming :#
#########################

PELT = function(y,Xsin,Xcos,om_grid,L,fr=1) {
  n = length(y)
  
  R1 = cor(Xsin)
  R2 = cor(Xcos)
  
  beta = (L+1)*log(n)
  
  f = rep(NA,n+1)
  f[1] = -beta
  
  cp = list()
  cp[[1]] = 0
  
  R = list()
  R[[1]] = 0
  
  
  for(t_star in 1:n) {

    f_hat = rep(NA,length(R[[t_star]]))
    it = 0
    for (t in R[[t_star]]) {
      it = it + 1
      
      segm = (t+1):t_star
      yW = y[segm]
      XsinW = matrix(Xsin[segm,],nrow=length(yW))
      XcosW = matrix(Xcos[segm,],nrow=length(yW))
      
      fitW = SuSiE_group_refined(yW,XsinW,XcosW,L,fr,R1,R2,om_grid) 
      if (fitW$sigma2 <= 0) fitW$sigma2 = var(y)
      f_hat[it] = f[t+1] - sum(dnorm(yW,fitW$yhat,sqrt(fitW$sigma2),log=T)) + beta
    }
    
    f[t_star+1] = min(f_hat)
    t_prime = R[[t_star]][which.min(f_hat)]
    cp[[t_star+1]] = unique(c(cp[[t_prime+1]],t_prime))
    R[[t_star+1]] = unique(c(t_star, R[[t_star]][f_hat - beta < f[t_star+1]]))
    
  }
  
  cp[[n]][-1]+1
}



