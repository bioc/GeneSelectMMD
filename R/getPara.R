"getPara" <-
function(X, memSubjects, memGenes, eps=1.0e-6)
{
  Xupp<-X[memGenes==1,, drop=FALSE] 
  Xlow<-X[memGenes==3,, drop=FALSE] 
  Xnon<-X[memGenes==2,, drop=FALSE]


  n.upp<-nrow(Xupp)
  n.low<-nrow(Xlow)
  n.non<-nrow(Xnon)

  n<-nrow(X)
  memGenes2<-rep(1, n)
  memGenes2[memGenes==2]<-0

  pi.1<-n.upp/n
  pi.3<-n.low/n
  pi.2<-1-pi.1-pi.3

  Xuppc<-Xupp[,memSubjects==1, drop=FALSE]
  Xuppn<-Xupp[,memSubjects==0, drop=FALSE]

  Xlowc<-Xlow[,memSubjects==1, drop=FALSE]
  Xlown<-Xlow[,memSubjects==0, drop=FALSE]


  # within gene variation
  v.uppc.w<-mean(apply(Xuppc, 1, var))
  v.uppn.w<-mean(apply(Xuppn, 1, var))

  v.non.w<-mean(apply(Xnon, 1, var))

  v.lowc.w<-mean(apply(Xlowc, 1, var))
  v.lown.w<-mean(apply(Xlown, 1, var))

  # between gene variation
  v.uppc.b<-var(apply(Xuppc, 1, mean))
  v.uppn.b<-var(apply(Xuppn, 1, mean))

  v.non.b<-var(apply(Xnon, 1, mean))

  v.lowc.b<-var(apply(Xlowc, 1, mean))
  v.lown.b<-var(apply(Xlown, 1, mean))

  # marginal correlation
  rho.c1<-v.uppc.b/(v.uppc.b+v.uppc.w)
  rho.n1<-v.uppn.b/(v.uppn.b+v.uppn.w)
  rho.2<-v.non.b/(v.non.b+v.non.w)
  rho.c3<-v.lowc.b/(v.lowc.b+v.lowc.w)
  rho.n3<-v.lown.b/(v.lown.b+v.lown.w)

  mu.c1<-mean(apply(Xuppc, 1, mean))
  mu.n1<-mean(apply(Xuppn, 1, mean))
  sigma2.c1<-v.uppc.w
  sigma2.n1<-v.uppn.w

  mu.2<-mean(apply(Xnon, 1, mean))
  sigma2.2<-v.non.w

  mu.c3<-mean(apply(Xlowc, 1, mean))
  mu.n3<-mean(apply(Xlown, 1, mean))
  sigma2.c3<-v.lowc.w
  sigma2.n3<-v.lown.w

  paraIni<-c(pi.1, pi.2, pi.3, 
             mu.c1, sigma2.c1, rho.c1, 
             mu.n1, sigma2.n1, rho.n1, 
             mu.2, sigma2.2, rho.2, 
             mu.c3, sigma2.c3, rho.c3, 
             mu.n3, sigma2.n3, rho.n3)
  names(paraIni)<-paraNames

  mat<-t(apply(X, 1, wiFun, Psi.m=paraIni, memSubjects=memSubjects, eps=eps))

  if(mu.c1<mu.n1 || mu.c3>mu.n3 )
  { llkh<- -Inf }
  else { 
    llkh<-QFunc(paraIni, X, mat[,1:3], memSubjects, memGenes, eps=eps)
  }

  return(list(para=paraIni, llkh=llkh))
}

