"getPara" <-
function(X, memSubjects, memGenes, eps=1.0e-6)
{
  Xupp<-X[memGenes==1,, drop=FALSE] 
  Xlow<-X[memGenes==3,, drop=FALSE] 
  Xnon<-X[memGenes==2,, drop=FALSE]

  nCases<-sum(memSubjects==1, na.rm=TRUE) # number of cases
  nControls<-sum(memSubjects==0, na.rm=TRUE) # number of controls
  nc<-nCases
  nn<-nControls
  nSubj<-nCases+nControls
  n<-nSubj
  
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
  v.uppc.w<-mean(apply(Xuppc, 1, var, na.rm=TRUE), na.rm=TRUE)
  v.uppn.w<-mean(apply(Xuppn, 1, var, na.rm=TRUE), na.rm=TRUE)

  v.non.w<-mean(apply(Xnon, 1, var, na.rm=TRUE), na.rm=TRUE)

  v.lowc.w<-mean(apply(Xlowc, 1, var, na.rm=TRUE), na.rm=TRUE)
  v.lown.w<-mean(apply(Xlown, 1, var, na.rm=TRUE), na.rm=TRUE)

  # between gene variation
  v.uppc.b<-var(apply(Xuppc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  v.uppn.b<-var(apply(Xuppn, 1, mean, na.rm=TRUE), na.rm=TRUE)

  v.non.b<-var(apply(Xnon, 1, mean, na.rm=TRUE), na.rm=TRUE)

  v.lowc.b<-var(apply(Xlowc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  v.lown.b<-var(apply(Xlown, 1, mean, na.rm=TRUE), na.rm=TRUE)

  # marginal correlation
  rho.c1<-v.uppc.b/(v.uppc.b+v.uppc.w)
  rho.n1<-v.uppn.b/(v.uppn.b+v.uppn.w)
  rho.2<-v.non.b/(v.non.b+v.non.w)
  rho.c3<-v.lowc.b/(v.lowc.b+v.lowc.w)
  rho.n3<-v.lown.b/(v.lown.b+v.lown.w)

  r.c1<-log((1+(nCases-1)*rho.c1)/((1-rho.c1)*(nCases-1)))
  r.n1<-log((1+(nControls-1)*rho.n1)/((1-rho.n1)*(nControls-1)))
  r.2<-log((1+(nSubj-1)*rho.2)/((1-rho.2)*(nSubj-1)))
  r.c3<-log((1+(nCases-1)*rho.c3)/((1-rho.c3)*(nCases-1)))
  r.n3<-log((1+(nControls-1)*rho.n3)/((1-rho.n3)*(nControls-1)))

  mu.c1<-mean(apply(Xuppc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  mu.n1<-mean(apply(Xuppn, 1, mean, na.rm=TRUE), na.rm=TRUE)
  sigma2.c1<-v.uppc.w
  sigma2.n1<-v.uppn.w

  mu.2<-mean(apply(Xnon, 1, mean, na.rm=TRUE), na.rm=TRUE)
  sigma2.2<-v.non.w

  mu.c3<-mean(apply(Xlowc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  mu.n3<-mean(apply(Xlown, 1, mean, na.rm=TRUE), na.rm=TRUE)
  sigma2.c3<-v.lowc.w
  sigma2.n3<-v.lown.w

  if(mu.c1<=mu.n1 || mu.c3>=mu.n3 )
  { llkh<- -Inf }
  else { 
    delta.n1<-log(mu.c1-mu.n1)
    delta.n3<-log(mu.n3-mu.c3)
    paraIni<-c(pi.1, pi.2, 
               mu.c1, log(sigma2.c1), r.c1, 
               delta.n1, log(sigma2.n1), r.n1, 
               mu.2, log(sigma2.2), r.2, 
               mu.c3, log(sigma2.c3), r.c3, 
               delta.n3, log(sigma2.n3), r.n3)
    names(paraIni)<-paraNamesRP
 

#   call the wiFun Fortran routine below 
#    mat<-t(apply(X, 1, wiFun, Psi.m=paraIni, memSubjects=memSubjects, eps=eps))

# initialize the wi arrays
    wi1<-rep(0,nrow(X))
    wi2<-rep(0,nrow(X))
    wi3<-rep(0,nrow(X))

#   checkpara was called within wiFun.R, however this will be handled
#   in Fortran, so checking the dimension of Psi.m may not be included
#   in the future

# replace NAs with zeros and pass this matrix to wifun
    Xnna<-X
    Xnna[is.na(Xnna)] <- 0

    res <- .Fortran("wifun", wi1=as.double(wi1), wi2=as.double(wi2), wi3=as.double(wi3), 
            as.double(Xnna), as.double(paraIni),
            as.integer(memSubjects), as.double(eps),
            as.integer(nrow(X)), as.integer(ncol(X)),
            as.integer(sum(memSubjects==1, na.rm=TRUE)),
            as.integer(sum(memSubjects==0, na.rm=TRUE))     )

    mat<-cbind(res$wi1,res$wi2,res$wi3)


    w1<-as.numeric(mat[,1])
    w2<-as.numeric(mat[,2])
    w3<-as.numeric(mat[,3])
    xcMat<-X[,memSubjects==1,drop=FALSE]
    xnMat<-X[,memSubjects==0,drop=FALSE]
 
    xcTxc<-apply(xcMat, 1, function(x) {sum(x^2, na.rm=TRUE)})
    xcT1<-apply(xcMat, 1, function(x) {sum(x, na.rm=TRUE)})
 
    xnTxn<-apply(xnMat, 1, function(x) {sum(x^2, na.rm=TRUE)})
    xnT1<-apply(xnMat, 1, function(x) {sum(x, na.rm=TRUE)})
 
    xTx<-apply(X, 1, function(x) {sum(x^2, na.rm=TRUE)})
    xT1<-apply(X, 1, function(x) {sum(x, na.rm=TRUE)})

    sumw1<-mean(mat[,1], na.rm=TRUE)
    sumw2<-mean(mat[,2], na.rm=TRUE)
    sumw3<-mean(mat[,3], na.rm=TRUE)
 
    sumw1xcTxc<-sum(w1*xcTxc, na.rm=TRUE)
    sumw1xcT1<-sum(w1*xcT1, na.rm=TRUE)
    sumw1xcT1.sq<-sum(w1*xcT1^2, na.rm=TRUE)
  
    sumw1xnTxn<-sum(w1*xnTxn, na.rm=TRUE)
    sumw1xnT1<-sum(w1*xnT1, na.rm=TRUE)
    sumw1xnT1.sq<-sum(w1*xnT1^2, na.rm=TRUE)
  
    sumw2xTx<-sum(w2*xTx, na.rm=TRUE)
    sumw2xT1<-sum(w2*xT1, na.rm=TRUE)
    sumw2xT1.sq<-sum(w2*xT1^2, na.rm=TRUE)
  
    sumw3xcTxc<-sum(w3*xcTxc, na.rm=TRUE)
    sumw3xcT1<-sum(w3*xcT1, na.rm=TRUE)
    sumw3xcT1.sq<-sum(w3*xcT1^2, na.rm=TRUE)
  
    sumw3xnTxn<-sum(w3*xnTxn, na.rm=TRUE)
    sumw3xnT1<-sum(w3*xnT1, na.rm=TRUE)
    sumw3xnT1.sq<-sum(w3*xnT1^2, na.rm=TRUE)
  

# call the Fortran routine negqfunc (converted from negQFunc)
  llkh <- 0

  res <- .Fortran("negqfunc", llkh=as.double(llkh), as.double(paraIni),
            as.double(nc), as.double(nn), as.double(n),
            as.double(sumw1), as.double(sumw2), as.double(sumw3),
            as.double(sumw1xcTxc), as.double(sumw1xcT1), as.double(sumw1xcT1.sq),
            as.double(sumw1xnTxn), as.double(sumw1xnT1), as.double(sumw1xnT1.sq),
            as.double(sumw2xTx), as.double(sumw2xT1), as.double(sumw2xT1.sq),
            as.double(sumw3xcTxc), as.double(sumw3xcT1), as.double(sumw3xcT1.sq),
            as.double(sumw3xnTxn), as.double(sumw3xnT1), as.double(sumw3xnT1.sq)       )

  llkh <- -res$llkh
  }

  return(list(para=paraIni, llkh=llkh))
}

