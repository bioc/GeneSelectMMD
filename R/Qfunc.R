# sum_{i=1}^{nGenes} log10(pi.1*f1(xi)+pi.2*f2(xi)+pi.3*f3(xi))

# M-step's objective function
# Q(Psi)=sum_{i=1}^{nGenes} wi*V(pi)+\sum_{i=1}^{nGenes) wi*U(theta)
"QFunc" <-
function(para, dat, wiMat, memSubjects, memGenes, eps=1.0e-6)
{
  # check if parameters are in appropriate ranges
  checkPara(para)

  pi.1<-para[1]; pi.2<-para[2]; pi.3<-para[3];

  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-para[4]; 

  # mean expression level for cluster 1 for normal subjects
  mu.n1<-para[7]; 

  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-para[13]; 

  # mean expression level for cluster 3 for normal subjects
  mu.n3<-para[16]; 

  if(mu.c1<mu.n1 || mu.c3>mu.n3)
  { 
    return(-Inf)
  }

  nGenes<-nrow(dat)

  tmp<-tapply(X=1:nGenes, INDEX=1:nGenes, FUN=QFunc.default, dat=dat, 
              para=para, wiMat=wiMat, memSubjects=memSubjects, memGenes=memGenes, eps=eps)
  res<-sum(tmp)

  return(res)
}


# log(pi.1*f1(xi)+pi.2*f2(xi)+pi.3*f3(xi))
"QFunc.default" <-
function(i, dat, para, wiMat=wiMat, memSubjects, memGenes, eps=1.0e-6)
{
  Xi<-dat[i,]

  # mixture proportions
  pi.1<-para[1]; pi.2<-para[2]; pi.3<-para[3];

  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-para[4]; 
  # variance of expression levels for cluster 1 for diseased subjects
  sigma2.c1<-para[5]; 
  # correlation among expression levels for cluster 1 for diseased subjects
  rho.c1<-para[6]; 

  # mean expression level for cluster 1 for normal subjects
  mu.n1<-para[7]; 
  # variance of expression levels for cluster 1 for normal subjects
  sigma2.n1<-para[8]; 
  # correlation among expression levels for cluster 1 for normal subjects
  rho.n1<-para[9]; 

  # mean expression level for cluster 2
  mu.2<-para[10]; 
  # variance of expression levels for cluster 2
  sigma2.2<-para[11]; 
  # correlation among expression levels for cluster 2
  rho.2<-para[12]; 

  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-para[13]; 
  # variance of expression levels for cluster 3 for diseased subjects
  sigma2.c3<-para[14]; 
  # correlation among expression levels for cluster 3 for diseased subjects
  rho.c3<-para[15]; 

  # mean expression level for cluster 3 for normal subjects
  mu.n3<-para[16]; 
  # variance of expression levels for cluster 3 for normal subjects
  sigma2.n3<-para[17]; 
  # correlation among expression levels for cluster 3 for normal subjects
  rho.n3<-para[18]; 

  ##############################
  n<-length(Xi) # number of patients/subjects
  nc<-sum(memSubjects==1)  # number of cases
  nn<-sum(memSubjects==0)  # number of controls
  if(nc>=n){
    stop("Number of cases >= Total number of patients!\n")
  }
  if(nn>=n){
    stop("Number of controls >= Total number of patients!\n")
  }

  # expression levels of the gene for diseased subjects
  Xci<-Xi[memSubjects==1]
  # expression levels of the gene for non-diseased subjects
  Xni<-Xi[memSubjects==0]

 #########################
  ###
  # density for genes in cluster 1 (over-expressed)
  ###
  XiTXi.c<-sum(Xci^2)
  XiT1.c<-sum(Xci)
  aiTai.c1<-XiTXi.c-2*mu.c1*XiT1.c+nc*mu.c1^2
  aiT12.c1<-XiT1.c^2-2*nc*mu.c1*XiT1.c+nc^2*mu.c1^2

  XiTXi.n<-sum(Xni^2)
  XiT1.n<-sum(Xni)
  aiTai.n1<-XiTXi.n-2*mu.n1*XiT1.n+nn*mu.n1^2
  aiT12.n1<-XiT1.n^2-2*nn*mu.n1*XiT1.n+nn^2*mu.n1^2

  part.c1<-( aiTai.c1  - rho.c1* aiT12.c1 / (1+(nc-1)*rho.c1) )/ (sigma2.c1 * (1-rho.c1)) 
  part.n1<-( aiTai.n1  - rho.n1* aiT12.n1 / (1+(nn-1)*rho.n1) )/ (sigma2.n1 * (1-rho.n1)) 
  delta1<- (part.c1+part.n1)

  log.detSigma1<-nc*log(sigma2.c1)+(nc-1)*log(1-rho.c1)+log(1+(nc-1)*rho.c1)+
                 nn*log(sigma2.n1)+(nn-1)*log(1-rho.n1)+log(1+(nn-1)*rho.n1)

  ###
  # density for genes in cluster 2 (non-expressed)
  ###
  XiTXi<-sum(Xi^2)
  XiT1<-sum(Xi)
  aiTai.2<-XiTXi-2*mu.2*XiT1+n*mu.2^2
  aiT12.2<-XiT1^2-2*n*mu.2*XiT1+n^2*mu.2^2

  n<-length(Xi)
  one<-rep(1,n)

  Sigma.2<-sigma2.2*(1-rho.2)*(diag(n)+rho.2/(1-rho.2)*one%*%t(one))
  egvalues.2<-eigen(Sigma.2)$values
  tmp.2<-(try( solve(Sigma.2), silent=TRUE))
  if(class(tmp.2) == "try-error")
  { 
    tmppos.2<-which(egvalues.2==min(egvalues.2))
    egvalues.2<-egvalues.2[-tmppos.2]
    det.Sigma.2<-prod(egvalues.2)
    iSigma.2<-ginv(Sigma.2)
    ai.2<-Xi-mu.2*one
    delta2<-t(ai.2)%*%iSigma.2%*%ai.2
    log.detSigma2<-log(det.Sigma.2)
  } else {
    delta2<-( aiTai.2  - rho.2* aiT12.2 / (1+(n-1)*rho.2) )/ (sigma2.2 * (1-rho.2)) 
    log.detSigma2<-n*log(sigma2.2)+(n-1)*log(1-rho.2)+log(1+(n-1)*rho.2)
  }
 
  ###
  # density for genes in cluster 3 (under-expressed)
  ###
  aiTai.c3<-XiTXi.c-2*mu.c3*XiT1.c+nc*mu.c3^2
  aiT12.c3<-XiT1.c^2-2*nc*mu.c3*XiT1.c+nc^2*mu.c3^2

  aiTai.n3<-XiTXi.n-2*mu.n3*XiT1.n+nn*mu.n3^2
  aiT12.n3<-XiT1.n^2-2*nn*mu.n3*XiT1.n+nn^2*mu.n3^2

  part.c3<-( aiTai.c3  - rho.c3* aiT12.c3 / (1+(nc-1)*rho.c3) )/ (sigma2.c3 * (1-rho.c3)) 
  part.n3<-( aiTai.n3  - rho.n3* aiT12.n3 / (1+(nn-1)*rho.n3) )/ (sigma2.n3 * (1-rho.n3)) 
  delta3<- (part.c3+part.n3)

  log.detSigma3<-nc*log(sigma2.c3)+(nc-1)*log(1-rho.c3)+log(1+(nc-1)*rho.c3)+
                 nn*log(sigma2.n3)+(nn-1)*log(1-rho.n3)+log(1+(nn-1)*rho.n3)

 #########################
  const<- -n/2*log(2*PI)
  logf1<- const-(log.detSigma1+delta1)/2
  logf2<- const-(log.detSigma2+delta2)/2
  logf3<- const-(log.detSigma3+delta3)/2


  wi.vec<-wiMat[i,]
  v.vec<-c(log(pi.1), log(pi.2), log(pi.3))
  u.vec<-c(logf1, logf2, logf3)

  part1<-sum(wi.vec*v.vec)
  part2<-sum(wi.vec*u.vec)

  res<-part1+part2

  return(res)
}

