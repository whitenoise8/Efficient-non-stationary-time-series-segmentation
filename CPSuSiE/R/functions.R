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

single_split_IC = function(y,Xsin,Xcos,L,cp,cost){
  
  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  p = matrix(NA,ncol(Xsin),2)
  llik = rep(NA,2)
  
  for (i in 1:(n_cp+1)) {
    segm = (pts[i]+1):pts[i+1]
    yW = y[segm]
    XsinW = matrix(Xsin[segm,],nrow=length(yW))
    XcosW = matrix(Xcos[segm,],nrow=length(yW))
    
    fitW = SuSiE_group(yW,XsinW,XcosW,L) 
    
    p[,i] = 1-apply(fitW$alpha,1,function(x) prod(1-x))
    llik[i] = sum(dnorm(yW,fitW$yhat,sqrt(fitW$sigma2),log=T))
  }
  
  p = 2*(sum(p > 0.5) + 1)
  IC = -2*sum(llik) + 2*p*log(n) + cost*sum((diff(pts)/n-0.5)^2)*log(n)
  
  IC
}

fit_segments = function(y,Xsin,Xcos,L,cp,pl=FALSE){
  
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
    
    fitW = SuSiE_group(yW,XsinW,XcosW,L) 
    out[[i]] = fitW
    yhat = c(yhat,fitW$yhat)
  }
  
  if (pl) {
    plot(y)
    abline(v=cp,lty='dashed',lwd=2)
    lines(yhat,col=2,lwd=2)
  }
  
  list(fits=out,
       s_hat=yhat)
}

###########################
#-- Binary segmentation --#
###########################

BS = function(y,Xsin,Xcos,L,w=10,cost=10) {
  
  n = length(y)
  cp_hat = NULL
  
  fit0 = SuSiE_group(y,Xsin,Xcos,L)
  p0 = sum(1-apply(fit0$alpha,1,function(x) prod(1-x)) > 0.5)
  IC0 = -2*sum(dnorm(y,fit0$yhat,sqrt(fit0$sigma2),log=T)) + (2*p0+1)*log(n)
  
  IC = rep(NA,n)
  
  for (i in w:(n-w)){
    yA = y[1:i]
    yB = y[(i+1):n]
    
    XsinA = Xsin[1:i,]
    XcosA = Xcos[1:i,]
    XsinB = Xsin[(i+1):n,]
    XcosB = Xcos[(i+1):n,]
    
    IC[i] = single_split_IC(y,Xsin,Xcos,L,cp=i,cost=cost)
  }
  
  IC_min = min(IC,na.rm=T)
  cp_star = which.min(IC)
  
  if (IC_min < IC0) cp_hat = cp_star
  
  list(cp_hat = cp_hat,
       IC_seq = IC,
       IC0 = IC0)
}

BS_MCP = function(y,Xsin,Xcos,L,w=10,cost=10) {
  n = length(y)
  
  cp_set = NULL
  eval_sgm = matrix(c(1,n),1,2)
  
  print(paste0("Segment: ",1," - ",n))
  out = BS(y,Xsin,Xcos,L,w,cost)
  if (!is.na(out$cp)) print(paste0("CP found: ",out$cp))
  print("---------------------------------------")
  
  cp_set = sort(c(cp_set,out$cp_hat))
  s_set = length(cp_set)
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
          print(paste0("Segment: ",pts[k], " - ",pts[k+1]))
          eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
          
          out_ = BS(y[x],Xsin[x,],Xcos[x,],L,w,cost)
          if (!is.null(out_$cp_hat)) {
            cps[k] = out_$cp_hat + pts[k]
            print(paste0("CP found: ",cps[k]))
          }
          print("---------------------------------------")
        }
      }
      
      cp_set = sort(c(cp_set,cps[!is.na(cps)]))
      
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

nOS_internal = function(y,Xsin,Xcos,L,lt,t,rt,l,r,nu=1/2,cost=10){
  
  n = nrow(Xsin)
  
  if ((rt-lt) <= 5){
    IC = c()
    if (lt == l) lt = l+1
    if (rt == r) rt = r-1
    for (i in lt:rt){
      IC = c(IC,single_split_IC(y,Xsin,Xcos,L,i,cost))
    }
    cp = which.min(IC)+lt-1
    out = list(cp=cp,
               IC=min(IC))
    return(out)
  } else {
    
    if ((rt-t) > (t-lt)) {
      w = floor(rt - (rt-t)*nu)
      IC_w = single_split_IC(y,Xsin,Xcos,L,w,cost)
      IC_t = single_split_IC(y,Xsin,Xcos,L,t,cost)
      if (IC_w < IC_t) { nOS_internal(y,Xsin,Xcos,L,t,w,rt,l,r,nu=1/2,cost=cost) } else { nOS_internal(y,Xsin,Xcos,L,lt,t,w,l,r,nu=1/2,cost=cost) }
    } else {
      w = floor(lt + (t-lt)*nu)
      IC_w = single_split_IC(y,Xsin,Xcos,L,w,cost)
      IC_t = single_split_IC(y,Xsin,Xcos,L,t,cost)
      if (IC_w < IC_t) { nOS_internal(y,Xsin,Xcos,L,lt,w,t,l,r,nu=1/2,cost=cost) } else { nOS_internal(y,Xsin,Xcos,L,w,t,rt,l,r,nu=1/2,cost=cost) }
    }
    
  }
}

nOS = function(y,Xsin,Xcos,L,l,r,nu=1/2,cost=10){
  
  t = floor((l+nu*r/(1+nu)))
  out_nOS = nOS_internal(y,Xsin,Xcos,L=L,lt=l,t=t,rt=r,l=l,r=r,nu=nu,cost=cost)
  
  cp_hat = NA
  
  fit0 = SuSiE_group(y,Xsin,Xcos,L)
  p0 = sum(1-apply(fit0$alpha,1,function(x) prod(1-x)) > 0.5)
  IC0 = -2*sum(dnorm(y,fit0$yhat,sqrt(fit0$sigma2),log=T)) + (2*p0+1)*log(n)
  
  if (out_nOS$IC < IC0)  cp_hat = out_nOS$cp

  list(cp_hat = cp_hat,
       IC = out_nOS$IC,
       IC0 = IC0)
}

nOS_MCP = function(y,Xsin,Xcos,L,cost=10,nu=1/2) {
  n = length(y)
  
  cp_set = NULL
  eval_sgm = matrix(c(1,n),1,2)
  
  print(paste0("Segment: ",1," - ",n))
  out = nOS(y,Xsin,Xcos,L,1,n,nu,cost)
  if (!is.na(out$cp)) print(paste0("CP found: ",out$cp))
  print("---------------------------------------")
  
  cp_set = sort(c(cp_set,out$cp))
  s_set = length(cp_set)
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        if (pts[k+1]-pts[k] > 100) {
          if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
            print(paste0("Segment: ",pts[k], " - ",pts[k+1]))
            eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
            
            out_ = nOS(y[x],Xsin[x,],Xcos[x,],L,1,length(x),nu,cost)
            if (!is.na(out_$cp)) {
              cps[k] = out_$cp + pts[k] - 1
              print(paste0("CP found: ",cps[k]))
            }
            print("---------------------------------------")
          }
        }
      }
      
      cp_set = sort(c(cp_set,cps[!is.na(cps)]))
      
      if (length(cp_set) > s_set) {
        s_set = length(cp_set)
      } else {
        conv = 1
      }
    }
  }
  
  cp_set
}



#: ALTRO

BS2 = function(y,Xsin,Xcos,L,w=10,thres=1) {
  
  n = length(y)
  cp_hat = NULL
  fit0 = SuSiE_group(y,Xsin,Xcos,L)
  
  marglik = rep(NA,n)
  
  for (i in w:(n-w)){
    yA = y[1:i]
    yB = y[(i+1):n]
    
    XsinA = Xsin[1:i,]
    XcosA = Xcos[1:i,]
    XsinB = Xsin[(i+1):n,]
    XcosB = Xcos[(i+1):n,]
    
    fitA = SuSiE_group(yA,XsinA,XcosA,L)
    fitB = SuSiE_group(yB,XsinB,XcosB,L)
    
    marglik[i] = fitA$elbo + fitB$elbo
  }
  
  elbo_max = max(marglik,na.rm=T)
  cp_star = which.max(marglik)
  
  if (elbo_max > fit0$elbo) if ((elbo_max-fit0$elbo)/abs(fit0$elbo) > thres) cp_hat = cp_star
  
  list(cp_hat = cp_hat,
       elbo_seq = marglik,
       elbo_0 = fit0$elbo)
}

BS_MCP2 = function(y,Xsin,Xcos,L,w=10,thres=1) {
  n = length(y)
  
  cp_set = NULL
  eval_sgm = matrix(c(1,n),1,2)
  
  print(paste0("Segment: ",1," - ",n))
  out = BS(y,Xsin,Xcos,L,w,thres)
  if (!is.na(out$cp)) print(paste0("CP found: ",out$cp))
  print(paste0("Segment: ",1," - ",n,"; elbo-variation: ",round((max(out$elbo_seq,na.rm=T)-out$elbo_0)/abs(out$elbo_0),2)))
  print("---------------------------------------")
  
  cp_set = sort(c(cp_set,out$cp_hat))
  s_set = length(cp_set)
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
          print(paste0("Segment: ",pts[k], " - ",pts[k+1]))
          eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
          
          out_ = BS(y[x],Xsin[x,],Xcos[x,],L,w,thres)
          if (!is.null(out_$cp_hat)) {
            cps[k] = out_$cp_hat + pts[k]
            print(paste0("CP found: ",cps[k]))
          }
          print(paste0("Segment: ",pts[k], " - ",pts[k+1],"; elbo-variation: ",round((max(out_$elbo_seq,na.rm=T)-out_$elbo_0)/abs(out_$elbo_0),2)))
          print("---------------------------------------")
        }
      }
      
      cp_set = sort(c(cp_set,cps[!is.na(cps)]))
      
      if (length(cp_set) > s_set) {
        s_set = length(cp_set)
      } else {
        conv = 1
      }
      
    }
  }
  
  cp_set
}

nOS_internal2 = function(y,Xsin,Xcos,L,lt,t,rt,l,r,nu=1/2){
  
  n = nrow(Xsin)
  
  if ((rt-lt) <= 5){
    marglik = c()
    for (i in lt:rt){
      marglik_i = split_SuSiE_group(y,Xsin,Xcos,L,i)
      marglik = c(marglik,marglik_i)
    }
    cp = which.max(marglik)+lt-1
    out = list(cp=cp,
               score=max(marglik))
    return(out)
  } else {
    
    if ((rt-t) > (t-lt)) {
      w = floor(rt - (rt-t)*nu)
      elbo_w = split_SuSiE_group(y,Xsin,Xcos,L,w)
      elbo_t = split_SuSiE_group(y,Xsin,Xcos,L,t)
      if (elbo_w > elbo_t) { nOS_internal(y,Xsin,Xcos,L,t,w,rt,l,r,nu=1/2) } else { nOS_internal(y,Xsin,Xcos,L,lt,t,w,l,r,nu=1/2) }
    } else {
      w = floor(lt + (t-lt)*nu)
      elbo_w = split_SuSiE_group(y,Xsin,Xcos,L,w)
      elbo_t = split_SuSiE_group(y,Xsin,Xcos,L,t)
      if (elbo_w > elbo_t) { nOS_internal(y,Xsin,Xcos,L,lt,w,t,l,r,nu=1/2) } else { nOS_internal(y,Xsin,Xcos,L,w,t,rt,l,r,nu=1/2) }
    }
    
  }
}

nOS2 = function(y,Xsin,Xcos,L,l,r,nu=1/2,thres=1){
  t = floor((l+nu*r/(1+nu)))
  out_nOS = nOS_internal(y,Xsin,Xcos,L=L,lt=l,t=t,rt=r,l=l,r=r,nu=nu)
  
  cp_hat = NA
  delta_elbo = 0
  elbo_0 = SuSiE_group(y,Xsin,Xcos,L)$elbo
  
  if (out_nOS$score > elbo_0) {
    delta_elbo = (out_nOS$score-elbo_0)/abs(elbo_0)
    if ((out_nOS$score-elbo_0)/abs(elbo_0) > thres) cp_hat = out_nOS$cp
  }
  
  list(cp_hat = cp_hat,
       elbo_OS = out_nOS$score,
       elbo_0 = elbo_0,
       delta_elbo = delta_elbo)
}

nOS_MCP2 = function(y,Xsin,Xcos,L,thres=1,nu=1/2) {
  n = length(y)
  
  cp_set = NULL
  eval_sgm = matrix(c(1,n),1,2)
  
  print(paste0("Segment: ",1," - ",n))
  out = nOS(y,Xsin,Xcos,L,1,n,nu,thres)
  if (!is.na(out$cp)) print(paste0("CP found: ",out$cp))
  print(paste0("Segment: 1 - ",n,"; elbo-variation: ",round(out$delta_elbo,2)))
  print("---------------------------------------")
  
  cp_set = sort(c(cp_set,out$cp))
  s_set = length(cp_set)
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        if (pts[k+1]-pts[k] > 100) {
          if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
            print(paste0("Segment: ",pts[k], " - ",pts[k+1]))
            eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
            
            out_ = nOS(y[x],Xsin[x,],Xcos[x,],L,1,length(x),nu,thres)
            if (!is.na(out_$cp)) {
              cps[k] = out_$cp + pts[k] - 1
              print(paste0("CP found: ",cps[k]))
            }
            print(paste0("Segment: ",pts[k], " - ",pts[k+1],"; elbo-variation: ",round(out_$delta_elbo,2)))
            print("---------------------------------------")
          }
        }
      }
      
      cp_set = sort(c(cp_set,cps[!is.na(cps)]))
      
      if (length(cp_set) > s_set) {
        s_set = length(cp_set)
      } else {
        conv = 1
      }
      #if (length(cp_set) >= 2) conv = 1
    }
  }
  
  cp_set
}

aOS2 = function(y,Xsin,Xcos,L,l,r,nu=1/2,thres=1){
  k = floor(log((r-l)/2,2))
  set_D = sort(unique(floor(c(l+(0.5)^(1:k)*(r-l),r-(0.5)^(1:k)*(r-l)))))
  
  marglik = c()
  for (i in set_D){
    marglik_i = split_SuSiE_group(y,Xsin,Xcos,L,i)
    marglik = c(marglik,marglik_i)
  }
  t_star = set_D[which.max(marglik)]
  
  if (t_star <= (r+l)/2) {
    lt = floor(t_star - (t_star-l)/2)
    rt = floor(t_star + (t_star-l))
  } else {
    lt = floor(t_star - (r-t_star))
    rt = floor(t_star + (r-t_star)/2)
  }
  
  out_aOS = nOS_internal(y,Xsin,Xcos,L=L,lt=lt,t=t_star,rt=rt,l=l,r=r,nu=nu)
  
  cp_hat = NA
  delta_elbo = 0
  elbo_0 = SuSiE_group(y,Xsin,Xcos,L)$elbo
  
  if (out_aOS$score > elbo_0) {
    delta_elbo = (out_aOS$score-elbo_0)/abs(elbo_0)
    if ((out_aOS$score-elbo_0)/abs(elbo_0) > thres) cp_hat = out_aOS$cp
  }
  
  list(cp_hat = cp_hat,
       elbo_OS = out_aOS$score,
       elbo_0 = elbo_0,
       delta_elbo = delta_elbo)
}

aOS_MCP2 = function(y,Xsin,Xcos,L,thres=0.5,nu=1/2) {
  n = length(y)
  
  cp_set = NULL
  eval_sgm = matrix(c(1,n),1,2)
  
  print(paste0("Segment: ",1," - ",n))
  out = aOS(y,Xsin,Xcos,L,1,n,nu,thres)
  if (!is.na(out$cp)) print(paste0("CP found: ",out$cp))
  print(paste0("Segment: 1 - ",n,"; elbo-variation: ",round(out$delta_elbo,2)))
  print("---------------------------------------")
  
  cp_set = sort(c(cp_set,out$cp))
  s_set = length(cp_set)
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
          print(paste0("Segment: ",pts[k], " - ",pts[k+1]))
          eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
          
          out_ = aOS(y[x],Xsin[x,],Xcos[x,],L,1,length(x),nu,thres)
          if (!is.na(out_$cp)) {
            cps[k] = out_$cp + pts[k] - 1
            print(paste0("CP found: ",cps[k]))
          }
          print(paste0("Segment: ",pts[k], " - ",pts[k+1],"; elbo-variation: ",round(out_$delta_elbo,2)))
          print("---------------------------------------")
        }
      }
      
      cp_set = sort(c(cp_set,cps[!is.na(cps)]))
      
      if (length(cp_set) > s_set) {
        s_set = length(cp_set)
      } else {
        conv = 1
      }
      
    }
  }
  
  cp_set
}
