# get wij=pi.j*f.j(xi)/[pi.1*f.1(xi)+pi.2*f.2(xi)+pi.3*f.3(xi)]
# j=1,2,3
# Psi.m=(pi.1, pi.2, pi.3, 
#        mu.c1, sigma2.c1, rho.c1, 
#        mu.n1, sigma2.n1, rho.n1, 
#        mu.2, sigma2.2, rho.2, 
#        mu.c3, sigma2.c3, rho.c3, 
#        mu.n3, sigma2.n3, rho.n3) 

"wiFun" <-
function(Xi, Psi.m, memSubjects, eps=1.0e-6)
{
  # check if parameters are in appropriate ranges
  checkPara(Psi.m)

  # mixture proportions
  pi.1<-Psi.m[1]; pi.2<-Psi.m[2]; pi.3<-Psi.m[3];

  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-Psi.m[4]; 
  # variance of expression levels for cluster 1 for diseased subjects
  sigma2.c1<-Psi.m[5]; 
  # correlation among expression levels for cluster 1 for diseased subjects
  rho.c1<-Psi.m[6]; 

  # mean expression level for cluster 1 for normal subjects
  mu.n1<-Psi.m[7]; 
  # variance of expression levels for cluster 1 for normal subjects
  sigma2.n1<-Psi.m[8]; 
  # correlation among expression levels for cluster 1 for normal subjects
  rho.n1<-Psi.m[9]; 

  # mean expression level for cluster 2
  mu.2<-Psi.m[10]; 
  # variance of expression levels for cluster 2
  sigma2.2<-Psi.m[11]; 
  # correlation among expression levels for cluster 2
  rho.2<-Psi.m[12]; 

  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-Psi.m[13]; 
  # variance of expression levels for cluster 3 for diseased subjects
  sigma2.c3<-Psi.m[14]; 
  # correlation among expression levels for cluster 3 for diseased subjects
  rho.c3<-Psi.m[15]; 

  # mean expression level for cluster 3 for normal subjects
  mu.n3<-Psi.m[16]; 
  # variance of expression levels for cluster 3 for normal subjects
  sigma2.n3<-Psi.m[17]; 
  # correlation among expression levels for cluster 3 for normal subjects
  rho.n3<-Psi.m[18]; 

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


 
  tt<-c(delta1, delta2, delta3)
  names(tt)<-c("delta1", "delta2", "delta3")
  pos<-which(tt==min(tt))
  pos<-pos[1]
  ################


  if(pos==1) {
    # wi1=pi.1*f1/(pi.1*f1+pi.2*f2+pi.3*f3)
    #    = pi.1/(pi.1+pi.2*f2/f1 + pi.3 * f3/f1)
    # wi2=pi.2*f2/f1/[pi.1+pi.2*f2/f1 + pi.3 * f3/f1]
    # wi3=pi.3*f3/f1/[pi.1+pi.2*f2/f1 + pi.3 * f3/f1]
    # f2/f1
    # log(f2/f1)=[log(|Sigma1|)-log(|Sigma2|)+delta1-delta2]/2
    logf2df1<-(log.detSigma1-log.detSigma2+delta1-delta2)/2
    f2df1<-exp(logf2df1)
    
    # f3/f1
    # log(f3/f1)=[log(|Sigma1|)-log(|Sigma3|)+delta1-delta3]/2
    logf3df1<-(log.detSigma1-log.detSigma3+delta1-delta3)/2
    f3df1<-exp(logf3df1)

    denom.1<-(pi.1+pi.2*f2df1+pi.3*f3df1)
  
    wi1<-pi.1/denom.1
    wi2<-pi.2*f2df1/denom.1
    wi3<-pi.3*f3df1/denom.1
    
  } else if(pos==2){

    #wi2=pi.2*f2/(pi.1*f1+pi.2*f2+pi.3*f3)
    #   = pi.2/(pi.1*f1/f2+pi.2 + pi.3 * f3/f2)
    #wi1=pi.1 * f1/f2 / [pi.1*f1/f2+pi.2 + pi.3 * f3/f2]
    #wi3=pi.3 * f3/f2 / [pi.1*f1/f2+pi.2 + pi.3 * f3/f2]
    # f1/f2
    # log(f1/f2)=(log|Sigma2|-log|Sigma1|+detla2-delta1)/2
    logf1df2<-(log.detSigma2-log.detSigma1+delta2-delta1)/2
    f1df2<-exp(logf1df2)
    
    # f3/f2
    # log(f3/f2)=(log|Sigma2|-log|Sigma3|+delta2-detla3)/2
    logf3df2<-(log.detSigma2-log.detSigma3+delta2-delta3)/2
    f3df2<-exp(logf3df2)

    denom.2<-(pi.1*f1df2+pi.2+pi.3*f3df2)
  
    wi1<-pi.1*f1df2/denom.2
    wi2<-pi.2/denom.2
    wi3<-pi.3*f3df2/denom.2
  } else {
    #wi3=pi.3*f3/(pi.1*f1+pi.2*f2+pi.3*f3)
    #   = pi.3/(pi.1*f1/f3+pi.2*f2/f3 + pi.3)
    #wi1=pi.1*f1/f3 / [pi.1*f1/f3+pi.2*f2/f3 + pi.3]
    #wi2=pi.2*f2/f3 / [pi.1*f1/f3+pi.2*f2/f3 + pi.3]
    # f1/f3
    # log(f1/f3)=(log|Sigma3|-log|Sigma1|+delta3-delta1)/2
    logf1df3<-(log.detSigma3-log.detSigma1+delta3-delta1)/2
    f1df3<-exp(logf1df3)
 
    # f2/f3
    # log(f2/f3)=(log|Sigma3|-log|Sigma2|+delta3-delta2)/2
    logf2df3<-(log.detSigma3-log.detSigma2+delta3-delta2)/2
    f2df3<-exp(logf2df3)
  
    denom.3<-(pi.1*f1df3+pi.2*f2df3+pi.3)


    wi1<-pi.1*f1df3/denom.3
    wi2<-pi.2*f2df3/denom.3
    wi3<-pi.3/denom.3
  }


  res<-c(wi1, wi2, wi3)
  names(res)<-c("wi1", "wi2", "wi3")

  return(res)
}

checkPara<-function(para, eps=1.0e-6)
{
  if(length(para) !=TNumPara)
  { stop("Number of parameters is not correct!\n") }

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

  if(pi.1>=1 || pi.1<=0)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.1 should be in (0, 1)!\n") }
  if(pi.2>=1 || pi.2<=0)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.2 should be in (0, 1)!\n") }
  if(pi.3>=1 || pi.3<=0)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.3 should be in (0, 1)!\n") }
  if(abs((pi.1+pi.2+pi.3)-1)>eps)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.1+pi.2+pi.3 should be equal to 1!\n") }

  if(sigma2.c1<=0)
  { stop("sigma2.c1 should be positive!\n") }
  if(sigma2.n1<=0)
  { stop("sigma2.n1 should be positive!\n") }
  if(sigma2.2<=0)
  { stop("sigma2 should be positive!\n") }
  if(sigma2.c3<=0)
  { stop("sigma2.c3 should be positive!\n") }
  if(sigma2.n3<=0)
  { stop("sigma2.n3 should be positive!\n") }

  if(rho.c1< -1 || rho.c1>1)
  { stop("rho.c1 should be in (-1, 1)!\n") }
  if(rho.n1< -1 || rho.n1>1)
  { stop("rho.n1 should be in (-1, 1)!\n") }
  if(rho.2< -1 || rho.2>1)
  { stop("rho.2 should be in (-1, 1)!\n") }
  if(rho.c3< -1 || rho.c3>1)
  { stop("rho.c3 should be in (-1, 1)!\n") }
  if(rho.n3< -1 || rho.n3>1)
  { stop("rho.n3 should be in (-1, 1)!\n") }

}


