solqp=function(H=NULL,f,A,B,X,neqtl,nctl,numVar){
  f=as.matrix(f)
  B=as.matrix(B)
  numVar=length(X)
  flag = 1
  iter = 0
  iterflag = 0
  if (is.null(H)){
    qpflag=0
  } else if (norm(H,'I')==0){
    qpflag=0
  }else{
    qpflag=1
  }
  
  imax = 10*numVar
  tolcon  = 10e-5
  Astan = matrix(sqrt(apply(A^2,1,sum)),ncol=1)
  B=B/Astan
  A=A/matrix(rep(Astan,ncol(A)),ncol=ncol(A))
  tolerr = 0.01*sqrt(.Machine$double.eps)
  tolf = 100*numVar*.Machine$double.eps
  
  orignctl = nctl;
  lambda = matrix(0,orignctl,1)
  if(neqtl>0){
    indxact=ieq = 1:neqtl
  }else{
    indxact=ieq = NULL
  }
  
  xieq = matrix(0,nctl,1)
  xieq[indxact] = 1
  iact = length(indxact)
  Act = A[indxact,,drop=F]
  indall = 1:nctl
  idel = NULL
  if (iact > 0 ){
    tolcons = 1e-10
    Z=NULL
    s=qr(A[ieq,])
    Qa=qr.Q(s)
    Ra=qr.R(s)
    idxdep =which(abs(diag(Ra))<tolf)
    if (neqtl > numVar){
      idxdep = rbind(idxdep, matrix((numVar+1):neqtl,ncol=1))
    }
    if (length(idxdep)>0){
      flag = 1
      bidxdep =  abs(t(Qa[,idxdep])%*%B[ieq]) >= tolf ;
      if (any( bidxdep )){  
        flag = -2;
      } else {
        st=qr(t(A[ieq,]))
        Qat=qr.Q(st)
        Rat=qr.R(st)
        idel=st$pivot[idxdep]
        numDepend = length(idel)
        A=A[-idel,]
        B=B[-idel,,drop=F]
        neqtl = neqtl - numDepend;
        nctl = nctl - numDepend;
        ieq = 1:neqtl;
        xieq=xieq[-idel,] 
        indxact=indxact[-(1:numDepend) ]
        indxact = indxact - numDepend;      
        Act = A[indxact,,drop=F]
        iact = iact - numDepend;       
      }    
    }
    
    if (iact >= numVar){
      iact = max(neqtl, numVar-1);
      indxact = indxact[1:iact];
      Act = A[indxact,,drop=F];
      xieq = matrix(0,nctl,1);
      xieq[indxact] = 1;
    }
    
    if (iact > neqtl){
      sa=qr(t(A[ieq,]))
      Qat=qr.Q(sa)
      Rat=qr.R(sa)
      idxdep =which(abs(diag(Ra))<tolf)
      if (length(idxdep)>0){       
        idel2=st$pivot[idxdep]
        idelEq   = idel2[which(idel2 <= neqtl)]
        idelIneq = idel2[which(idel2 > neqtl)]
        if (length(idelEq)>0){
          if(neqtl){
            indxact = 1:neqtl
          }else{indxact=NULL}
        } else{
          indxact=indxact[-idelIneq]
        }
        xieq = matrix(0,nctl,1);
        xieq[indxact] = 1;
        Act = A[indxact,,drop=F]
        iact = length(indxact); 
      }
    }
    
    s=qr(t(Act))
    Q=qr.Q(s,complete=T)
    R=qr.R(s,complete=T)
    Z = Q[,(iact+1):numVar]
    
    if(flag != -2 & iact > 0){       
      deltax = Q[,1:iact]%*%(solve(t(R[1:iact,1:iact]), (B[indxact] - Act%*%X)))
      X = X + deltax
      err = A%*%X - B;
      err[ieq] = abs(err[ieq])
      if (any(err > .Machine$double.eps)){
        Xb = ginv(Act)%*%B[indxact,,drop=F];
        errb = A%*%Xb - B;
        errb[ieq] = abs(errb[ieq])
        if (max(errb) < max(err)){
          X = Xb
        }
      }
    }
    
    if (length(idel)>0){
      indall=indall[-idel]
      Astan = Astan[indall,,drop=F]
    }
    
    if (flag == -2){   
      indxact = indall[indxact]
      flag = -2
      return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
      
    }
    
    err = 0
    if (neqtl >= numVar){
      err = max(abs(A[ieq,]%*%X-B[ieq]))
      if(err > 1e-8){
        flag = -2
        indxact = indall[indxact]
        return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
        
        
      } else {   
        if(max(A%*%X-B) > 1e-8){
          flag = -2
        }     
      }
      
      if(qpflag){
        lambdact = -ginv(R)%*%(t(Q)%*%(H%*%X+f))
      } else{
        lambdact = -ginv(R)%*%(t(Q)%*%f)
      }
      
      lambda[indall[indxact],,drop=F] = (lambdact/Astan[indxact]);
      indxact = indall[indxact];
      return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))    
    }      
  } else {
    if(iact == 0){
      Q = diag(1,numVar)
      R = NULL
      Z = 1
    }else{
      s=qr(t(Act))
      Q=qr.Q(s,complete=T)
      R=qr.R(s,complete=T)
      Z = Q[,(iact+1):numVar,drop=F]
    }   
  }
  
  const = A%*%X-B;
  constdf = 0
  tolconnorm = .Machine$double.eps
  
  if(nctl > neqtl){
    constdf=max(const[(neqtl+1):nctl]) 
    id=which.max(const[(neqtl+1):nctl])
    tolconnorm = tolcon/Astan[neqtl+id]
  }
  
  if (constdf > tolconnorm){   
    if(neqtl){
      indxact2 = 1:neqtl
    }else{indxact2=NULL}
    A2=cbind(rbind(A,matrix(0,1,numVar)),rbind(matrix(0,neqtl,1),-matrix(1,nctl+1-neqtl,1)))
    st=solqp(NULL,matrix(c(rep(0,numVar),1),ncol=1),A2,rbind(B,1e-5),rbind(X,constdf+1),
             neqtl,nrow(A2),numVar+1)
    Xp=st$X;lambdaS=st$lambda;
    slack=Xp[numVar+1]
    X = Xp[1:numVar,,drop=F]
    const = A%*%X-B
    constdf=max(const[(neqtl+1):nctl]) 
    id=which.max(const[(neqtl+1):nctl])
    tolconnorm = .Machine$double.eps
    if(slack > tolconnorm){
      if(slack>1e-8){
        flag=-2
      } else{
        flag=-4
      }
      lambda[indall,] <- lambdaS[(1:nctl),]/Astan; 
      if(neqtl){
        indxact = 1:neqtl
      }else{indxact=NULL}
      indxact = indall[indxact]
      return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))     
    } else {
      if(neqtl){
        indxact = 1:neqtl
      }else{indxact=NULL}
      Act = A[indxact,,drop=F]
      iact = length(indxact)
      xieq = matrix(0,nctl,1)
      xieq[indxact] = 1
      if(iact == 0){
        Q = matrix(0,numVar,numVar)
        R=NULL
        Z=1
      } else{
        s=qr(t(Act))
        Q=qr.Q(s,complete=T)
        R=qr.R(s,complete=T)
        Z = Q[,(iact+1):numVar]
      }   
    }
  } 
  
  if (iact >= numVar - 1) {
    iterflag = 1
  } 
  
  if(qpflag){
    ggf=H%*%X+f 
    ss=getDir(Z,H,ggf,numVar,f)
    Dir=ss$Dir
    searchDir=getElement(ss,"searchDir")
  }else{
    ggf=f
    if(length(Z)>1){
      Dir=-Z%*%t(Z)%*%ggf
    } else{
      Dir=-Z*Z*ggf  
    }
    searchDir='SteepDescent'
    if(norm(Dir)<1e-10 & neqtl){
      lambdact = -ginv(R)%*%(t(Q)%*%ggf)
      lambda[indall[indxact],,drop=F] = (lambdact/Astan[indxact]);
      indxact = indall[indxact];
      return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))    
    }
  }
  
  
  ind_old = 0
  while (iter < imax){
    iter = iter + 1
    GDir=A%*%Dir
    indf = which((GDir > tolerr * norm(Dir,"2"))  &  !as.matrix(xieq))
    if(length(indf)==0){
      minalpha = 1e16
      ind=NULL
    }else{
      dist = abs(const[indf,,drop=F])/GDir[indf,,drop=F];
      minalpha = min(dist);
      ind2 = which(dist == minalpha);
      ind = indf[min(ind2)]
    }
    
    flagremove = 0
    if(length(indf)>0 & is.finite(minalpha)){
      if(identical(searchDir,'Newton')){
        if(minalpha>1){
          minalpha=1
          flagremove = 1
        }
        X = X+minalpha*Dir
      }else{ 
        X = X+minalpha*Dir        
      } 
    } else { 
      if(identical(searchDir,'Newton')){
        minalpha = 1
        X = X + Dir
        flagremove = 1
        
      } else{
        if(!qpflag | identical(searchDir,'NegCurv')){
          if(norm(Dir,"2")>tolerr){
            minalpha = 1e16
            X = X + minalpha*Dir
            flag = -3                   
          } else { 
            flag = -7 
          }
          indxact = indall(indxact)
          return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
          
        }else{
          if(qpflag){
            ZHZ = t(Z)%*%H%*%Z; 
            Zggf = t(Z)%*%ggf;
            Dirproj = ginv(ZHZ)%*%(-Zggf)
          }
          
          if (qpflag & (norm(ZHZ%*%Dirproj+Zggf,"2") > 10*.Machine$double.eps*(norm(ZHZ,"2") + norm(Zggf,"2")))){
            if(norm(Dir,"2") > tolerr){
              minalpha = 1e16
              X = X + minalpha*Dir
              flag = -3
            } else {
              flag = -7
            }
            indxact = indall[indxact]
            return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
            
          } else {
            Dir = Z%*%Dirproj   
            if(t(ggf)%*%Dir > 0){
              Dir = -Dir
            }
            searchDir = 'singular'
            GDir=A%*%Dir
            indf = which((GDir > tolerr * norm(Dir,"2"))  &  !as.matrix(xieq))
            if(length(indf)==0){
              minalpha = 1e16
              ind=NULL
            }else{
              dist = abs(const[indf,,drop=F])/GDir[indf,,drop=F];
              minalpha = min(dist);
              ind2 = which(dist == minalpha);
              ind = indf[min(ind2)]
            }
            if (minalpha > 1){
              minalpha = 1;
              flagremove = 1
            }
            X = X + minalpha*Dir           
          }     
        }      
      } 
    } 
    
    if(qpflag){
      ggf=H%*%X+f
    }
    
    const = A%*%X-B;
    const[ieq] = abs(const[ieq,,drop=F]);
    if(neqtl<nctl){
      constdf = max(const[(neqtl+1):nctl])
    } else{constdf=0}
    
    if (flagremove){
      if (iact>0){
        rlambda = -ginv(R)%*%(t(Q)%*%ggf)
        lambdact = rlambda;
        lambdact[ieq] = abs(rlambda[ieq,,drop=F]);
        indlam = which(lambdact < 0)
        if(!length(indlam)){
          lambda[indall[indxact]] = (rlambda/Astan[indxact]);
          indxact = indall[indxact]
          return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
          
        }
        
        minindact = which(indxact == min(indxact[indlam]))[1]
        if(!is.na(minindact)){
          Act=Act[-minindact,,drop=F] 
          xieq[indxact[minindact]] = 0;
          
          st=qprm(Q,R,minindact);
          Q=st$Q
          R=st$R
          indxact=indxact[-minindact]
          iact = length(indxact)
          iterflag = 0;
          ind = 0;
        }
      }else {
        return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))   
      } 
      flagremove = 0
    }
    
    if (constdf > 1e5*tolerr){
      flag = -1
    }else{ 
      flag = 1;
    } 
    
    if(!is.null(ind) ){
      xieq[ind,]=1;
      indplus = length(indxact) + 1
      Act=rbind(Act,A[ind,])
      indxact=c(indxact,ind)
      m=nrow(Act)
      n=ncol(Act)
      
      st = qrin(Q,R,indplus,t(A[ind,,drop=F]))
      Q=st$Q
      R=st$R
      iact = length(indxact)
    }
    
    if(!iterflag){
      m=nrow(Act)
      n=ncol(Act)
      Z=Q[,(m+1):n,drop=F]
      if(iact == numVar - 1){
        iterflag = 1
      }
      ind_old = 0
    } else {
      rlambda = -ginv(R)%*%(t(Q)%*%ggf)
      if (is.infinite(rlambda[1]) && rlambda[1] < 0){
        rlambda = -ginv(t(Act + sqrt(.Machine$double.eps)*matrix(rnorm(m*n),nrow=m)))%*%ggf  
      }
      lambdact = rlambda;
      lambdact[ieq]=abs(lambdact[ieq,,drop=F]);
      indlam =which(lambdact<0)
      if(length(indlam)){ 
        if (minalpha > tolerr){
          minl=min(lambdact)
          minindact=which(lambdact == min(lambdact))[1]
        }else{
          minindact = which(indxact == min(indxact[indlam,,drop=F]))[1]
        }
        Act=Act[-minindact,,drop=F]
        xieq[indxact[minindact]] = 0
        st=qprm(Q,R,minindact);
        Q=st$Q
        R=st$R
        Z = Q[,numVar,drop=F]
        ind_old = indxact[minindact]
        indxact=indxact[-minindact]
        iact = length(indxact)
        
      }else{ 
        lambda[indall[indxact]]= (rlambda/Astan[indxact,,drop=F])
        indxact = indall[indxact]
        return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
      }
    }
    
    if(qpflag){
      Zggf = t(Z)%*%ggf
      if (length(Zggf)>0 & (norm(Zggf,"2") < 1e-15)){
        Dir = matrix(0,numVar,1)
        searchDir = 'ZeroStep'
      }else{
        ss=getDir(Z,H,ggf,numVar,f)
        Dir=ss$Dir
        searchDir=getElement(ss,"searchDir")
      }  
    }else{
      if(!iterflag){
        Dir=-Z%*%(t(Z)%*%ggf)
        gradDir = norm(Dir,'2')
      }else{
        gradDir = t(Z)%*%ggf
        if(gradDir>0){
          Dir=-Z
        }else{
          Dir = Z
        }
      } 
      
      if(abs(gradDir) < 1e-10){
        if(!ind_old){
          rlambda = -ginv(R)%*%(t(Q)%*%ggf)
          indxacttmp = indxact; Qtmp = Q; Rtmp = R;
        }else{ 
          indxacttmp = indxact
          indxacttmp[(minindact+1):(iact+1)] = indxact[minindact:iact];
          indxacttmp[minindact] = ind_old;
          st = qrin(Q,R,minindact,t(A[ind_old,]))
          Qtmp=st$Q
          Rtmp=st$R
        } 
        lambdact = rlambda;
        if(neqtl){
          lambdact[1:neqtl] = abs(lambdact[1:neqtl])
        }
        indlam = which(lambdact < tolerr);
        lambda[indall[indxacttmp]] = (rlambda/Astan[indxacttmp,,drop=F]);
        if(!length(indlam)){
          indxact = indall[indxact]
          return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
        }
        
        
        indplusmax = length(indlam);
        indpluscnt = 0;
        m = length(indxacttmp)
        while ((abs(gradDir) < 1e-10) && (indpluscnt < indplusmax)){
          indpluscnt = indpluscnt + 1;
          minindact = indlam[indpluscnt]
          st=qprm(Qtmp,Rtmp,minindact);
          Q=st$Q
          R=st$R
          Z = Q[,m:numVar]
          if(m!=numVar){
            if(length(Z)>1){
              Dir=-Z%*%t(Z)%*%ggf
            } else{
              Dir=-Z*Z*ggf  
            }
            gradDir = norm(Dir,'2')
          }else{
            gradDir = t(Z)%*%ggf
            if(gradDir > 0){
              Dir = -Z
            }else{
              Dir = Z
            }
          }     
        }
        if(abs(gradDir) < 1e-10){
          indxact = indall[indxact]
          return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
        }else{ 
          indxact = indxacttmp
          indxact=indxact[-minindact]
          xieq = matrix(0,nctl,1)
          iact = length(indxact)
          Act = A(indxact,,drop=F)
        } 
        lambda = matrix(0,orignctl,1)
      }
    } 
  } 
  if(iter >= imax){
    indxact = indall[indxact]
    flag = 0
  }
  return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
  
}

getDir=function(Z,H,ggf,numVar,f){
  if(length(Z)==1){flag=1}else{flag=0}
  if(flag==1){
    ZHZ = Z*H*Z
  }else{
    ZHZ = t(Z)%*%H%*%Z
  }
  
  if(any(eigen(ZHZ)$values>0)){
    tep=chol(ZHZ)
    if(flag==1){
      Dir = - Z*(ginv(tep)%*%( ginv(t(tep))%*%(Z*ggf)))
    }else{
      Dir = - Z%*%(ginv(tep)%*%( ginv(t(tep))%*%(t(Z)%*%ggf)))
    }
    searchDir = "Newton"
    
  } else {
    ss=cholfac(ZHZ)
    sneg=ss$sneg
    if((!is.null(sneg)) & t(sneg)%*%ZHZ%*%sneg < -sqrt(.Machine$double.eps)) {#negative enough
      if(flag==1){
        Dir=Z*sneg
      }else{
        Dir=Z%*%sneg
      }
      searchDir = 'NegCurv'
    }else{
      if(flag==1){
        Dir = - Z*(Z*ggf)
      }else{
        Dir = - Z%*%(t(Z)%*%ggf)
      }
      searchDir = 'SteepDescent'
    }
  }
  return(list(Dir=Dir,searchDir=searchDir))
} 

cholfac=function(ZHZ) {
  sneg=NULL
  n=nrow(ZHZ)
  L=diag(n)
  tol=0
  for(k in 1:(n-1)){
    if(ZHZ[k,k]<=tol){
      elem = matrix(0,length(ZHZ),1)
      elem[k,1] = 1
      sneg = t(L) %*% ginv(elem)
      return(list(L=L,sneg=sneg))
    }else{
      L[k,k] = sqrt(ZHZ[k,k])
      s = (k+1):n;
      L[s,k] = ZHZ[s,k]/L[k,k]
      ZHZ[(k+1):n,(k+1):n]=ZHZ[(k+1):n,(k+1):n]-lower.tri(L[(k+1):n,k]%*%t(L[(k+1):n,k]))
    }
  }
  if(ZHZ[n,n] <= tol){
    elem = matrix(0,length(ZHZ),1)
    elem[n,1] = 1
    sneg = t(L) %*% ginv(elem)
  }else{
    L[n,n] = sqrt(ZHZ[n,n])
  }
  return(list(L=L,sneg=sneg))
} 


qprm=function(Q,R,j){ 
  R=R[,-j,drop=F]
  m=nrow(R)
  n=ncol(R)
  
  for (k in j:min(n,m-1)){
    p=k:(k+1)
    teppl=rotateplan(R[p,k])
    R[p,k]=teppl$x
    G=teppl$G
    if(k<n){
      R[p,((k+1):n)]=G%*%R[p,((k+1):n)]
    }
    Q[,p]=Q[,p]%*%t(G)  
  }
  return(list(Q=Q,R=R)) 
}

rotateplan=function(x){
  if (x[2]!=0){
    r=norm(x,"2")
    G=rbind(x,cbind(-x[2],x[1]))/r
    x=rbind(r,0)
  }else{
    G=diag(length(x))
  }
  return(list(G=G,x=x))
}

qrin=function(Q,R,j,x){
  x=as.matrix(x)
  if(length(R)>0){
    m=nrow(R)
    n=ncol(R)
  }else{
    m=0;n=0;
  }
  
  if(n==0){
    ss=qr(x)
    Q=qr.Q(ss,complete=T)
    R=qr.R(ss,complete=T)
    return(list(Q=Q,R=R))
  }
  R=cbind(R,matrix(0,nrow(R),1))
  if(j<n+1){
    R[,(j+1):(n+1)] = R[,j:n]
    R[,j] = t(Q)%*%x
  }else if(j==n+1){
    R[,j] = t(Q)%*%x
  }
  n=n+1
  
  if(m>j){
    for (k in (m-1):j){
      p=k:(k+1)
      teppl=rotateplan(R[p,j])
      R[p,j]=teppl$x
      G=teppl$G
      if(k<n){
        R[p,((k+1):n)]=G%*%R[p,((k+1):n)]
      }
      Q[,p]=Q[,p]%*%t(G) 
    }
  }
  return(list(Q=Q,R=R)) 
}


estgradient=function(Xin,objfun=NULL,confun=NULL,lb,ub,fin,cIneq,cEq,variables, 
                     medx,dims,gradf,estgIn=NULL,estgeq=NULL){
  fin= as.vector(fin)
  lbflag = length(lb)>0;
  ubflag = length(ub)>0;
  Xin = as.vector(Xin)
  deltaX = sqrt(.Machine$double.eps)*ifelse(Xin>=0,1,-1)*max(abs(Xin),abs(medx))
  
  for (i in variables){
    if ((lbflag & is.finite(lb[i])) || (ubflag & is.finite(ub[i]))){
      deltaX[i] = boundflag(Xin[i],lb[i],ub[i],deltaX[i],i);
    }
    xfix = Xin[i];
    Xin[i] = xfix + deltaX[i];
    if (!is.null(objfun) ){
      fplus = objfun(matrix(Xin,dims$nxrows,dims$nxcols), N = N, T = T, f_hat_vec = f_bar, u_mat = u_mat, asset_group_vec = asset_group_vec, 
                     n_vec = n_vec, ind_t_dist = 1, ind_Rt = 0, ind_sim = 0)
      gradf[i] = (fplus - fin)/deltaX[i];
    }
    if (!is.null(confun) ){
      cons=confun(matrix(Xin,dims$nxrows,dims$nxcols))
      cIneqPlus=cons$c
      cEqPlus=cons$ceq   
      cIneqPlus = as.vector(cIneqPlus); cEqPlus = as.vector(cEqPlus);
      estgIn[i,] = t(cIneqPlus - cIneq) / deltaX[i];
      estgeq[i,] = (cEqPlus - cEq) / deltaX[i];            
    }
    Xin[i] = xfix;
  }
  nvals = length(variables)
  return(list(gradf=gradf,estgIn=estgIn,estgeq=estgeq,nvals=nvals)) 
}


boundflag=function(x,lb,ub, delta,i){
  if (lb != ub & x >= lb & x <= ub){
    if  ((x + delta > ub) || (x + delta < lb)){
      delta = -delta
      if((x + delta) > ub || (x + delta) < lb) {
        delta=ifelse(x-lb>ub-x, -(x-lb),ub-x)
      }
    }
  }  
  return(delta)
}


function (X = NULL, objfun = NULL, confun = NULL, A = NULL, B = NULL, 
          Aeq = NULL, Beq = NULL, lb = NULL, ub = NULL, tolX = 1e-05, 
          tolFun = 1e-06, tolCon = 1e-06, maxnFun = 1e+07, maxIter = 4000) 
{
  
  X =  para_start_Mf_LT_2
  objfun = LogLik_Copula_LT_factor_given_omega_Optim
  N = N
  T = T
  f_hat_vec = f_bar_sim
  u_mat = x_mat_sim
  asset_group_vec = asset_group_vec
  
  n_vec = n_vec
  ind_t_dist = 1
  ind_Rt = 0
  ind_sim = 1
  confun = NULL
  A = NULL
  B = NULL
  Aeq = NULL
  Beq = NULL
  lb = Lb_Mf_LT_2
  ub = Rb_Mf_LT_2
  tolX = 1e-05
  tolFun = 1e-06
  tolCon = 1e-06
  maxnFun = 1e+07
  maxIter = 4000
  
  if (is.null(X)) 
    stop("please input initial value")
  if (is.null(objfun)) 
    stop("please write objective function")
  X = as.matrix(X)
  Xtarget = as.vector(X)
  dims = NULL
  dims$nxrows = nrow(X)
  dims$nxcols = ncol(X)
  dims$nvar = length(Xtarget)
  B = as.vector(B)
  Beq = as.vector(Beq)
  if (is.null(Aeq)) {
    Aeq = matrix(NA, nrow = 0, ncol = dims$nvar)
  }
  if (is.null(A)) {
    A = matrix(NA, 0, dims$nvar)
  }
  nLineareq = nrow(Aeq)
  ncolAeq = ncol(Aeq)
  nLinearIneq = nrow(A)
  Ainput = A
  lenghlb = length(lb)
  lenghub = length(ub)
  if (lenghlb > dims$nvar) {
    print("invalid bounds")
    lb = lb(1:dims$nvar)
    lenghlb = dims$nvar
  }
  else if (lenghlb < dims$nvar) {
    if (lenghlb > 0) {
      print("invalid bounds")
    }
    lb = rbind(lb, matrix(rep(-Inf, dims$nvar - lenghlb), 
                          ncol = 1))
    lenghlb = dims$nvar
  }
  if (lenghub > dims$nvar) {
    print("invalid bounds")
    ub = ub(1:dims$nvar)
    lenghub = dims$nvar
  }
  else if (lenghub < dims$nvar) {
    if (lenghub > 0) {
      print("invalid bounds")
    }
    ub = rbind(ub, matrix(rep(Inf, dims$nvar - lenghub), 
                          ncol = 1))
    lenghub = dims$nvar
  }
  len = min(lenghlb, lenghub)
  if (sum(lb > ub) > 0) {
    return("please check bounds")
  }
  start = NULL
  start$xform = X
  medx = matrix(rep(1, dims$nvar), ncol = 1)
  Xtarget[Xtarget < lb] = lb[Xtarget < lb]
  Xtarget[Xtarget > ub] = ub[Xtarget > ub]
  X = matrix(Xtarget, dim(start$xform))
  start$g = matrix(0, dims$nvar, 1)
  start$f = objfun(params = X, N = N, T = T, f_hat_vec = f_bar, u_mat = u_mat, asset_group_vec = asset_group_vec, 
                   n_vec = n_vec, ind_t_dist = 1, ind_Rt = 0, ind_sim = 0)
  if (!is.null(confun)) {
    conf = confun(X)
    ctmp = conf$c
    ceqtmp = conf$ceq
    start$ncineq = as.vector(ctmp)
    start$nceq = as.vector(ceqtmp)
    start$gnc = matrix(0, dims$nvar, length(start$ncineq))
    start$gnceq = matrix(0, dims$nvar, length(start$nceq))
  }
  else {
    conf = ctmp = ceqtmp = start$ncineq = start$nceq = start$gnc = start$gnceq = NULL
  }
  if (is.null(start$ncineq)) 
    start$ncineq = matrix(0, 0, 1)
  if (is.null(start$nceq)) 
    start$nceq = matrix(0, 0, 1)
  if (is.null(start$gnc)) 
    start$gnc = matrix(0, dims$nvar, 0)
  if (is.null(start$gnceq)) 
    start.gnceq = matrix(0, dims$nvar, 0)
  fval = NULL
  lambda_out = NULL
  lambda_nc = NULL
  GRADIENT = NULL
  xform = start$xform
  iter = 0
  Xtarget = as.vector(X)
  numVar = length(Xtarget)
  DIR = matrix(1, numVar, 1)
  finalf = Inf
  steplength = 1
  HESS = diag(numVar)
  finishflag = F
  lbflag = is.finite(lb)
  ubflag = is.finite(ub)
  boundM = diag(max(lenghub, lenghlb))
  if (sum(lbflag) > 0) {
    lbM = -boundM[lbflag, 1:numVar]
    lbright = -as.matrix(lb)[lbflag, , drop = F]
  }else {
    lbM = NULL
    lbright = NULL
  }
  if (sum(ubflag) > 0) {
    ubM = boundM[ubflag, 1:numVar]
    ubright = as.matrix(ub)[ubflag, , drop = F]
  }else {
    ubM = NULL
    ubright = NULL
  }
  A = rbind(lbM, ubM, A)
  B = as.vector(c(lbright, ubright, B))
  if (length(A) == 0) {
    A = matrix(0, 0, numVar)
    B = matrix(0, 0, 1)
  }
  if (length(Aeq) == 0) {
    Aeq = matrix(0, 0, numVar)
    Beq = matrix(0, 0, 1)
  }
  LAMBDA_new = NULL
  LAMBDA = NULL
  LAMBDA_old = NULL
  X = matrix(Xtarget, dim(start$xform))
  f = start$f
  nceq = start$nceq
  ncineq = start$ncineq
  nctmp = ncineq
  nc = rbind(as.matrix(nceq), as.matrix(ncineq))
  c = rbind(rbind(rbind(Aeq %*% Xtarget - Beq, as.matrix(nceq)), 
                  as.matrix(A %*% Xtarget) - matrix(B, ncol = 1)), as.matrix(ncineq))
  nonlin_eq = length(nceq)
  nonlin_Ineq = length(ncineq)
  nLineareq = nrow(Aeq)
  nLinearIneq = nrow(A)
  eq = nonlin_eq + nLineareq
  ineq = nonlin_Ineq + nLinearIneq
  nctl = ineq + eq
  if (nonlin_eq > 0) {
    nleq_i = (1:nonlin_eq)
  }else {
    nleq_i = NULL
  }
  if (nonlin_Ineq > 0) {
    nlineq_i = ((nonlin_eq + 1):(nonlin_eq + nonlin_Ineq))
  }else {
    nlineq_i = NULL
  }
  if (eq > 0 & ineq > 0) {
    ga = rbind(abs(c[1:eq, , drop = F]), c[(eq + 1):nctl, 
                                           , drop = F])
  }else if (eq > 0) {
    ga = abs(c[1:eq, , drop = F])
  }else if (ineq > 0) {
    ga = c[(eq + 1):nctl, , drop = F]
  }else {
    ga = NULL
  }
  if (length(c) > 0) {
    maxga = max(ga)
  }else maxga = 0
  x_old = Xtarget
  c_old = c
  ggf_old = matrix(0, numVar, 1)
  ggf = start$g
  gnc = cbind(start$gnceq, start$gnc)
  tgc_old = matrix(0, nctl, numVar)
  LAMBDA = matrix(0, nctl, 1)
  lambda_nc = matrix(0, nctl, 1)
  nfval = 1
  ngval = 1
  errfloat = NULL
  while (!finishflag) {
    len_nc = length(nc)
    nctl = nLineareq + nLinearIneq + len_nc
    resfd = estgradient(Xtarget, objfun, confun, lb, ub, 
                        f, nc[nlineq_i], nc[nleq_i], 1:numVar, medx, dims, 
                        ggf, gnc[, nlineq_i, drop = F], gnc[, nleq_i, drop = F])
    ggf = resfd$gradf
    gnc[, nlineq_i] = resfd$estgIn
    gnc[, nleq_i] = resfd$estgeq
    nvals = resfd$nvals
    nfval = nfval + nvals
    if (length(gnc) > 0) {
      gc = cbind(cbind(cbind(t(Aeq), gnc[, nleq_i]), t(A)), 
                 gnc[, nlineq_i])
    }else if (length(Aeq) > 0 || length(A) > 0) {
      gc = cbind(t(Aeq), t(A))
    }else {
      gc = matrix(0, numVar, 0)
    }
    tgc = t(gc)
    if (eq > 0) {
      for (i in 1:eq) {
        iopp = tgc[i, ] %*% ggf
        if (iopp > 0) {
          tgc[i, ] = -tgc[i, ]
          c[i] = -c[i]
        }
      }
    }
    if (iter > 0) {
      maxgrad = norm(as.matrix(ggf + t(tgc) %*% lambda_nc), 
                     "I")
      if (nctl > eq) {
        maxc = norm(as.matrix(lambda_nc[(eq + 1):nctl] * 
                                c[(eq + 1):nctl]), "I")
      }else {
        maxc = 0
      }
      if (is.finite(maxgrad) & is.finite(maxc)) {
        errfloat = max(maxgrad, maxc)
      }else {
        errfloat = Inf
      }
      errga = maxga
      if (errfloat < tolFun & errga < tolCon) {
        finishflag = TRUE
      }else {
        if (nfval > maxnFun) {
          Xtarget = xtrial
          f = f_old
          ggf = ggf_old
          finishflag = TRUE
        }
        if (iter >= maxIter) {
          finishflag = TRUE
        }
      }
    }
    if (!finishflag) {
      iter = iter + 1
      if (ngval > 1) {
        LAMBDA_new = LAMBDA
        g_new = ggf + t(tgc) %*% LAMBDA_new
        g_old = ggf_old + t(tgc_old) %*% LAMBDA
        gdif = g_new - g_old
        xdif = Xtarget - x_old
        if (t(gdif) %*% xdif < steplength^2 * 0.001) {
          while (t(gdif) %*% xdif < -1e-05) {
            gdif[order(gdif * xdif)[1], ] = gdif[order(gdif * 
                                                         xdif)[1], ]/2
          }
          if (t(gdif) %*% xdif < (.Machine$double.eps * 
                                  norm(HESS, "F"))) {
            tgccdif = t(tgc) %*% c - t(tgc_old) %*% c_old
            tgccdif = tgccdif * (xdif * tgccdif > 0) * 
              (gdif * xdif <= .Machine$double.eps)
            weight = 0.01
            if (max(abs(tgccdif)) == 0) {
              tgccdif = 1e-05 * sign(xdif)
            }
            while (t(gdif) %*% xdif < (.Machine$double.eps * 
                                       norm(HESS, "F")) & (weight < 1/.Machine$double.eps)) {
              gdif = gdif + weight * tgccdif
              weight = weight * 2
            }
          }
        }
        if (t(gdif) %*% xdif > .Machine$double.eps) {
          HESS = HESS + (gdif %*% t(gdif))/c(t(gdif) %*% 
                                               xdif) - ((HESS %*% xdif) %*% (t(xdif) %*% 
                                                                               t(HESS)))/c(t(xdif) %*% HESS %*% xdif)
        }
      }else {
        LAMBDA_old = matrix(.Machine$double.eps + t(ggf) %*% 
                              ggf, nctl, 1)/(apply(t(tgc) * t(tgc), 2, sum) + 
                                               .Machine$double.eps)
        iact = 1:eq
      }
      ngval = ngval + 1
      L_old = LAMBDA
      tgc_old = tgc
      ggf_old = ggf
      c_old = c
      f_old = f
      x_old = Xtarget
      xint = matrix(0, numVar, 1)
      HESS = (HESS + t(HESS)) * 0.5
      tryCatch({
        resqp = solve.QP(HESS, -ggf, -t(tgc), c, meq = eq)
        DIR = resqp$solution
        lambda = resqp$Lagrangian
        iact = resqp$iact
      }, error = function(e) {
        resqp = solqp(HESS, ggf, tgc, -c, xint, eq, nrow(tgc), 
                      numVar)
        DIR = resqp$X
        lambda = resqp$lambda
        iact = resqp$indxact
      })
      lambda_nc[, 1] = 0
      lambda_nc[iact, ] = lambda[iact]
      lambda[(1:eq)] = abs(lambda[(1:eq)])
      if (eq > 0 & ineq > 0) {
        ga = rbind(abs(c[1:eq, , drop = F]), c[(eq + 
                                                  1):nctl, , drop = F])
      }else if (eq > 0) {
        ga = abs(c[1:eq, , drop = F])
      }else if (ineq > 0) {
        ga = c[(eq + 1):nctl, , drop = F]
      }
      if (length(c) > 0) {
        maxga = max(ga)
      }else maxga = 0
      LAMBDA = lambda[(1:nctl)]
      LAMBDA_old = apply(cbind(LAMBDA, 0.5 * (LAMBDA + 
                                                LAMBDA_old)), 1, max)
      ggfDIR = t(ggf) %*% DIR
      xtrial = Xtarget
      commerit = f + sum(LAMBDA_old * (ga > 0) * ga) + 
        1e-30
      if (maxga > 0) {
        commerit2 = maxga
      }else if (f >= 0) {
        commerit2 = -1/(f + 1)
      }else {
        commerit2 = 0
      }
      if (f < 0) {
        commerit2 = commerit2 + f - 1
      }
      if ((maxga < .Machine$double.eps) & (f < finalf)) {
        finalf = f
        finalx = Xtarget
        finalHess = HESS
        finalgrad = ggf
        finallambda = lambda
        finalmaxga = maxga
        finalerrfloat = errfloat
      }
      search = T
      alpha = 2
      while (search & nfval < maxnFun) {
        alpha = alpha/2
        if (alpha < 1e-04) {
          alpha = -alpha
        }
        if ((norm(as.matrix(DIR), "I") < 2 * tolX || 
             abs(alpha * ggfDIR) < tolFun) & (maxga < tolCon)) {
          finishflag = T
        }
        Xtarget = xtrial + alpha * DIR
        X = matrix(Xtarget, dim(start$xform))
        f = objfun(X, N = N, T = T, f_hat_vec = f_bar, u_mat = u_mat, asset_group_vec = asset_group_vec, 
                   n_vec = n_vec, ind_t_dist = 1, ind_Rt = 0, ind_sim = 0)
        if (!is.null(confun)) {
          conf = confun(X)
          nctmp = conf$c
          nceqtmp = conf$ceq
        }else {
          conf = nctmp = nceqtmp = NULL
        }
        nctmp = as.vector(nctmp)
        nceqtmp = as.vector(nceqtmp)
        nfval = nfval + 1
        nc = c(nceqtmp, nctmp)
        c = matrix(c(Aeq %*% Xtarget - Beq, nceqtmp, 
                     A %*% Xtarget - matrix(B, ncol = 1), nctmp), 
                   ncol = 1)
        if (eq > 0 & ineq > 0) {
          ga = rbind(abs(c[1:eq, , drop = F]), c[(eq + 
                                                    1):nctl, , drop = F])
        }else if (eq > 0) {
          ga = abs(c[1:eq, , drop = F])
        }else if (ineq > 0) {
          ga = c[(eq + 1):nctl, , drop = F]
        }
        if (length(c) > 0) {
          maxga = max(ga)
        }else maxga = 0
        merit = f + sum(LAMBDA_old * (ga > 0) * ga)
        if (maxga > 0) {
          merit2 = maxga
        }else if (f >= 0) {
          merit2 = -1/(f + 1)
        }else {
          merit2 = 0
        }
        if (f < 0) {
          merit2 = merit2 + f - 1
        }
        search = (merit2 > commerit2) & (merit > commerit)
      }
      steplength = alpha
      if (!finishflag) {
        mf = abs(steplength)
        LAMBDA = mf * LAMBDA + (1 - mf) * L_old
        X = matrix(Xtarget, dim(start$xform))
      }
    }
  }
  GRADIENT = ggf
  if (f > finalf) {
    Xtarget = finalx
    f = finalf
    HESS = finalHess
    GRADIENT = finalgrad
    lambda = finallambda
    maxga = finalmaxga
    ggf = finalgrad
    errfloat = finalerrfloat
  }
  fval = f
  X = matrix(Xtarget, dim(start$xform))
  nLinearIneq = nrow(Ainput)
  lambda_out = NULL
  lambda_out$lower = matrix(0, lenghlb, 1)
  lambda_out$upper = matrix(0, lenghub, 1)
  if (nLineareq > 0) {
    lambda_out$eqlin = lambda_nc[1:nLineareq]
  }
  ii = nLineareq
  if (nonlin_eq > 0) {
    lambda_out$eqnonlin = lambda_nc[(ii + 1):(ii + nonlin_eq)]
  }
  ii = ii + nonlin_eq
  if (sum(lbflag != 0) > 0) {
    lambda_out$lower[lbflag] = lambda_nc[(ii + 1):(ii + sum(lbflag != 
                                                              0))]
  }
  ii = ii + sum(lbflag != 0)
  if (sum(ubflag != 0) > 0) {
    lambda_out$upper[ubflag] = lambda_nc[(ii + 1):(ii + sum(ubflag != 
                                                              0))]
  }
  ii = ii + sum(ubflag != 0)
  if (nLinearIneq > 0) {
    lambda_out$ineqlin = lambda_nc[(ii + 1):(ii + nLinearIneq)]
  }
  ii = ii + nLinearIneq
  if (nonlin_Ineq > 0) {
    lambda_out$ineqnonlin = lambda_nc[(ii + 1):length(lambda_nc)]
  }
  return(list(par = X, fn = fval, counts = cbind(nfval, ngval), 
              lambda = lambda_out, grad = GRADIENT, hessian = HESS))
}