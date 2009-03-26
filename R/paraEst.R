"paraEst" <-
function(X, paraIni, memSubjects, maxFlag=TRUE, thrshPostProb=0.50, 
         geneNames=NULL, ITMAX=100, eps=1.0e-03, quiet=TRUE)
{
  # check if parameters are in appropriate ranges
  checkPara(paraIni)

  pi.1<-paraIni[1]; pi.2<-paraIni[2]; pi.3<-paraIni[3];

  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-paraIni[4]; 
  # variance of expression levels for cluster 1 for diseased subjects
  sigma2.c1<-paraIni[5]; 
  # correlation among expression levels for cluster 1 for diseased subjects
  rho.c1<-paraIni[6]; 

  # mean expression level for cluster 1 for normal subjects
  mu.n1<-paraIni[7]; 
  # variance of expression levels for cluster 1 for normal subjects
  sigma2.n1<-paraIni[8]; 
  # correlation among expression levels for cluster 1 for normal subjects
  rho.n1<-paraIni[9]; 

  # mean expression level for cluster 2
  mu.2<-paraIni[10]; 
  # variance of expression levels for cluster 2
  sigma2.2<-paraIni[11]; 
  # correlation among expression levels for cluster 2
  rho.2<-paraIni[12]; 

  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-paraIni[13]; 
  # variance of expression levels for cluster 3 for diseased subjects
  sigma2.c3<-paraIni[14]; 
  # correlation among expression levels for cluster 3 for diseased subjects
  rho.c3<-paraIni[15]; 

  # mean expression level for cluster 3 for normal subjects
  mu.n3<-paraIni[16]; 
  # variance of expression levels for cluster 3 for normal subjects
  sigma2.n3<-paraIni[17]; 
  # correlation among expression levels for cluster 3 for normal subjects
  rho.n3<-paraIni[18]; 

  n<-ncol(X) # number of patients/subjects
  nc<-sum(memSubjects==1)  # number of cases
  nn<-sum(memSubjects==0)  # number of controls
  if(nc>=n){
    stop("Number of cases >= Total number of patients!\n")
  }
  if(nn>=n){
    stop("Number of controls >= Total number of patients!\n")
  }


  # start em algorithm
  Psi.m<-paraIni
  nGenes<-nrow(X)

  loop<-0
  while(1)
  { 
    if(!quiet)
    {
      cat("************\n")
      cat("loop=", loop, "\n")
    }

    # E-step
    # obtain weights w_{ij}(\Psi^{(m)}), i.e. E(Z_{ij} | x_i, Psi^{(m)})
    mat<-t(apply(X, 1, wiFun, Psi.m=Psi.m, memSubjects=memSubjects, eps=eps))

    wi1<-mat[,1]
    wi2<-mat[,2]
    wi3<-mat[,3]
    if(!quiet)
    { cat("sumwi1=", sum(wi1), " sumwi2=", sum(wi2), " sumwi3=", sum(wi3), "\n") }

    # M-step

    # cluster 1
    # abnormal tissue samples
    res.c1<-updatePara3(dat=X[,memSubjects==1], wi=wi1)
    mu.c1<-res.c1$mu.hat
    sigma2.c1<-res.c1$sigma2.hat
    rho.c1<-res.c1$rho.hat

    # cluster 1
    # normal tissue samples
    res.n1<-updatePara3(dat=X[,memSubjects==0], wi=wi1)
    mu.n1<-res.n1$mu.hat
    sigma2.n1<-res.n1$sigma2.hat
    rho.n1<-res.n1$rho.hat

    # cluster 2
    res.2<-updatePara3(dat=X, wi=wi2)
    mu.2<-res.2$mu.hat
    sigma2.2<-res.2$sigma2.hat
    rho.2<-res.2$rho.hat

    # cluster 3
    # abnormal tissue samples
    res.c3<-updatePara3(dat=X[,memSubjects==1], wi=wi3)
    mu.c3<-res.c3$mu.hat
    sigma2.c3<-res.c3$sigma2.hat
    rho.c3<-res.c3$rho.hat

    # cluster 3
    # normal tissue samples
    res.n3<-updatePara3(dat=X[,memSubjects==0], wi=wi3)
    mu.n3<-res.n3$mu.hat
    sigma2.n3<-res.n3$sigma2.hat
    rho.n3<-res.n3$rho.hat

    # mixture proportions
    pi.1<-mean(wi1)
    pi.3<-mean(wi3)
    pi.2<-1-pi.1-pi.3
    pi.vec<-c(pi.1, pi.2, pi.3)

    if(mu.c1<mu.n1 || mu.c3>mu.n3)
    { 
      break 
    }

    # update parameters
    Psi.m.new<-c(pi.vec, 
                 mu.c1, sigma2.c1, rho.c1,
                 mu.n1, sigma2.n1, rho.n1,
                 mu.2, sigma2.2, rho.2,
                 mu.c3, sigma2.c3, rho.c3,
                 mu.n3, sigma2.n3, rho.n3)
   names(Psi.m.new)<-paraNames
   if(!quiet)
   { cat("Psi.m.new>>\n"); print(round(Psi.m.new,3)); cat("\n") }
      
    err.len<-sum(abs(Psi.m.new-Psi.m)<eps)
    if(sum((abs(Psi.m.new-Psi.m))<eps)==length(Psi.m))
      break

    Psi.m<-Psi.m.new
    loop<-loop+1

    if(loop>ITMAX)
    { cat("*********************\n") 
      cat("Warning! Number of looping exceeds ITMAX!\n")
      cat("EM algorithm did not converge!\n")
      cat("*********************\n") 
      break
    }
  }

  
  if(!quiet)
  { cat("Total iterations for EM algorithm=", loop, "\n") }
  # update gene cluster membership
  memGenes<-apply(mat[,1:3], 1, maxPosFun, maxFlag=maxFlag, 
                  thrshPostProb=thrshPostProb)
  llkh<-QFunc(Psi.m, X, mat[,1:3], memSubjects, memGenes, eps=eps)

  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0

  if(sum(is.null(geneNames)))
  {
    geneNames<-paste("gene", 1:nGenes, sep="")
  } 

  wiMat<-mat[,1:3]
  rownames(wiMat)<-geneNames
  colnames(wiMat)<-paste("cluster", 1:3, sep="")

  names(Psi.m)<-paraNames
  names(memGenes)<-geneNames
  names(memGenes2)<-geneNames

  invisible(list(para=Psi.m, llkh=llkh, memGenes=memGenes, 
                 memGenes2=memGenes2, wiMat=wiMat, loop=loop))
}

# find the position of the maximum element of 'x'
maxPosFun<-function(x, maxFlag=TRUE, thrshPostProb=0.50)
{
  if(maxFlag)
  { pos<-which(x==max(x))
    pos<-pos[1]
    return(pos)
  }
    
  pos<-which(x>thrshPostProb)
  len<-length(pos)
  if(len==0)
  { 
    return(2) 
  } # non-differentially expressed

  pos<-which(x==max(x))

  return(pos)
}

# update rho
# a=sum_{i=1}^{p} w_{i1}, say
# b=sum_{i=1}^{p} w_{i1}*(a_i^T*1)^2
# myc=sum_{i=1}^{p} w_{i1}*(a_i^T*a_i)
get.rho<-function(a, b, myc, n, sigma)
{
percent.1<- -(-24*b^3*myc-12*b^2*myc^2-3*n^6*a^4*sigma^8+3*n^7*a^4*sigma^8-
              36*b^2*n^2*myc^2+36*b^3*n*myc+48*b^2*n*myc^2+12*b*n^3*myc^3-
              24*b*n^2*myc^3+12*b*n*myc^3-12*b^4+24*n^2*a*sigma^2*myc^2*b-
              24*n*a*sigma^2*myc*b^2+12*n^3*a^2*sigma^4*myc*b-
              3*n^4*a^2*sigma^4*myc^2-6*n^5*a^3*sigma^6*myc-
              12*n^4*a^2*sigma^4*b*myc+12*n^3*a^2*sigma^4*b^2+
              3*n^5*a^2*sigma^4*myc^2+24*n^5*a^3*sigma^6*b-
              6*n^6*a^3*sigma^6*myc-12*n^2*a^2*sigma^4*b^2+
              72*b^2*n^3*myc*a*sigma^2-36*b*n^4*myc^2*a*sigma^2-
              84*b^2*n^2*myc*a*sigma^2+36*b*n^5*myc*a^2*sigma^4-
              36*b^3*n^2*a*sigma^2+72*b^3*n*a*sigma^2-
              36*b^2*n^4*a^2*sigma^4+12*b*myc^2*n^3*a*sigma^2-
              12*n^6*a^3*sigma^6*b)/(n-1)


percent.2<- (8*n^6*a^3*sigma^6-12*n^4*a^3*sigma^6+8*n^3*a^3*sigma^6- 
             12*n^5*a^3*sigma^6+8*myc^3+8*b^3-48*n^3*a*sigma^2*myc*b+ 
             72*n^2*a*sigma^2*b*myc-24*n*a*sigma^2*b*myc-
             24*n^3*a*sigma^2*myc^2- 24*n^2*a*sigma^2*myc^2+
             24*n*a*sigma^2*myc^2+ 48*n^3*a^2*sigma^4*b-
             12*n^3*a^2*sigma^4*myc+24*n^4*a^2*sigma^4*b+ 
             12*n^4*a^2*sigma^4*myc-48*n^2*a^2*sigma^4*b+ 
             24*n^2*a^2*sigma^4*myc+24*n^4*a*sigma^2*myc^2-
             24*n^5*a^2*sigma^4*myc+ 24*myc^2*n^2*b-24*myc*n*b^2-
             48*myc^2*n*b-8*myc^3*n^3+24*myc^3*n^2- 24*myc^3*n+ 
             24*myc*b^2+24*b*myc^2+24*n^2*a*sigma^2*b^2-48*n*a*sigma^2*b^2- 
             12*(percent.1)^(1/2)*sigma^2*a*n+
             12*(percent.1)^(1/2)*sigma^2*a*n^2) * (-2*n+1+n^2)^2*(n-1)^2 

percent.3<- (2*n*a*sigma^2*myc-n^3*a^2*sigma^4+n^2*a^2*sigma^4-
             2*b*n*myc+b^2+ 2*b*myc+myc^2*n^2- 2*myc^2*n+myc^2+
             2*n^2*a*sigma^2*b-2*n^3*a*sigma^2*myc-4*n*a*sigma^2*b+ 
             n^4*a^2*sigma^4) *(1+n^2-2*n) / 
             ((n-1)*n*a*sigma^2*(percent.2)^(1/3)) 

percent.4<- (percent.2)^(1/3)/ ((1+n^2-2*n)*sigma^2*a*n*(n-1))

percent.5<- (n^2*a*sigma^2-2*n*a*sigma^2-myc*n+b+myc)/((n-1)*n*a*sigma^2)

rho<- (1/6)*percent.4 + (2/3)*percent.3+(1/3)*percent.5

  return(rho)

}

# get initial parameter estimates
getIniParaAll<-function(dat, mat, memSubjects, maxFlag=TRUE, thrshPostProb=0.50)
{
  wi1<-mat[,1]
  wi2<-mat[,2]
  wi3<-mat[,3]

  pi.1<-mean(wi1)
  pi.3<-mean(wi3)
  pi.2<-1-pi.1-pi.3

  memGenes<-apply(mat[,1:3], 1, maxPosFun, maxFlag=maxFlag, 
                  thrshPostProb=thrshPostProb)

  datc1<-dat[memGenes==1, memSubjects==1]
  datn1<-dat[memGenes==1, memSubjects==0]
  dat2<-dat[memGenes==2, ]
  datc3<-dat[memGenes==3, memSubjects==1]
  datn3<-dat[memGenes==3, memSubjects==0]

  para.c1<-getIniPara(datc1)
  para.n1<-getIniPara(datn1)

  para.2<-getIniPara(dat2)

  para.c3<-getIniPara(datc3)
  para.n3<-getIniPara(datn3)

  para.est<-c(pi.1, pi.2, pi.3, para.c1, para.n1, para.2, para.c3, 
                    para.n3)

  return(para.est)
}

# Lagrange method to get optimal parameters
hFunc<-function(varCor, mat, pi.vec, mean.vec, memSubjects)
{
  pi.1<-pi.vec[1]; pi.2<-pi.vec[2]; pi.3<-pi.vec[3];
  mu.c1<-mean.vec[1]; mu.n1<-mean.vec[2];
  mu.2<-mean.vec[3];
  mu.c3<-mean.vec[4]; mu.n3<-mean.vec[5];

  n<-ncol(mat) # number of patients/subjects
  nc<-sum(memSubjects==1)  # number of cases
  nn<-sum(memSubjects==0)  # number of controls

  wi1<-mat[,1]
  wi2<-mat[,2]
  wi3<-mat[,3]

  xci.bar<-mat[,4]
  xni.bar<-mat[,5]
  x2.bar<-mat[,6]
  sumsq.c1<-mat[,7]
  sumsq.n1<-mat[,8]
  sumsq.2<-mat[,9]
  sumsq.c3<-mat[,10]
  sumsq.n3<-mat[,11]
  sum.c1<-mat[,12]
  sum.n1<-mat[,13]
  sum.2<-mat[,14]
  sum.c3<-mat[,15]
  sum.n3<-mat[,16]


  # variance of expression levels for cluster 1 for diseased subjects
  sigma2.c1<-varCor[1]; 
  # correlation among expression levels for cluster 1 for diseased subjects
  rho.c1<-varCor[2]; 

  # variance of expression levels for cluster 1 for normal subjects
  sigma2.n1<-varCor[3]; 
  # correlation among expression levels for cluster 1 for normal subjects
  rho.n1<-varCor[4]; 

  # variance of expression levels for cluster 2
  sigma2.2<-varCor[5]; 
  # correlation among expression levels for cluster 2
  rho.2<-varCor[6]; 

  # variance of expression levels for cluster 3 for diseased subjects
  sigma2.c3<-varCor[7]; 
  # correlation among expression levels for cluster 3 for diseased subjects
  rho.c3<-varCor[8]; 

  # variance of expression levels for cluster 3 for normal subjects
  sigma2.n3<-varCor[9]; 
  # correlation among expression levels for cluster 3 for normal subjects
  rho.n3<-varCor[10]; 

  # get part1=log(pi.1)*sum(wi1)+log(pi.2)*sum(wi2)+log(pi.3)*sum(wi3)
  sumwi1<-sum(wi1)
  sumwi2<-sum(wi2)
  sumwi3<-sum(wi3)

  part1<-log(pi.1)*sumwi1+log(1-pi.1-pi.3)*sumwi2+log(pi.3)*sumwi3

  const<- -n*log(2*PI)/2
  log.detSigma1<-nc*log(sigma2.c1)+(nc-1)*log(1-rho.c1)+log(1+(nc-1)*rho.c1)+
                 nn*log(sigma2.n1)+(nn-1)*log(1-rho.n1)+log(1+(nn-1)*rho.n1)
  part.c1<-( sum(wi1*sumsq.c1)  - rho.c1* sum(wi1*sum.c1^2) / (1+(nc-1)*rho.c1) )/ (sigma2.c1 * (1-rho.c1)) 
  part.n1<-( sum(wi1*sumsq.n1)  - rho.n1* sum(wi1*sum.n1^2) / (1+(nn-1)*rho.n1) )/ (sigma2.n1 * (1-rho.n1)) 
  delta1<- (part.c1+part.n1)

  wilogf1<-const*sumwi1-log.detSigma1*sumwi1/2-delta1/2

  delta2<-( sum(wi2*sumsq.2)  - rho.2* ( sum(wi2*sum.2^2) ) / (1+(n-1)*rho.2) )/ (sigma2.2 * (1-rho.2)) 

  log.detSigma2<-n*log(sigma2.2)+(n-1)*log(1-rho.2)+log(1+(n-1)*rho.2)

  wilogf2<-const*sumwi2-log.detSigma2*sumwi2/2-delta2/2

  log.detSigma3<-nc*log(sigma2.c3)+(nc-1)*log(1-rho.c3)+log(1+(nc-1)*rho.c3)+
                 nn*log(sigma2.n3)+(nn-1)*log(1-rho.n3)+log(1+(nn-1)*rho.n3)
  part.c3<-( sum(wi3*sumsq.c3)  - rho.c3* sum(wi3*sum.c3^2) / (1+(nc-1)*rho.c3) )/ (sigma2.c3 * (1-rho.c3)) 
  part.n3<-( sum(wi3*sumsq.n3)  - rho.n3* sum(wi3*sum.n3^2) / (1+(nn-1)*rho.n3) )/ (sigma2.n3 * (1-rho.n3)) 
  delta3<- (part.c3+part.n3)

  wilogf3<-const*sumwi3-log.detSigma3*sumwi3/2-delta3/2

  res<-part1+wilogf1+wilogf2+wilogf3

  return( -res)
}

getIniPara<-function(dat)
{
  nVariables<-ncol(dat)
  nSubjects<-nrow(dat)
  m.vec<-apply(dat, 2, sum)
  one<-rep(1, nVariables)
  mu.est<-sum(m.vec*one)/(nVariables*nSubjects)

  # within gene variation
  v2<-mean(apply(dat, 1, var))
  # between gene variation
  tau2<-var(apply(dat, 1, mean))

  sigma2.est<-v2+tau2

  rho.est<-tau2/sigma2.est
  
  para.est<-c(mu.est, sigma2.est, rho.est)
  names(para.est)<-c("mu.est", "sigma2.est", "rho.est")

  return(para.est)
}


updatePara<-function(wi, xiTxi, xiT1, n)
{
  sumwi<-sum(wi)
  sumwixiTxi<-sum(wi*xiTxi)
  sumwixiT1<-sum(wi*xiT1)
  sumwixiT12<-sum(wi*(xiT1)^2)

  delta<-sumwixiTxi*sumwi*n
  mu.hat<-sumwixiT1/(n*sumwi)

  numer.sigma2.hat<-sumwixiTxi*n^2*sumwi-delta+(sumwixiT1)^2-sumwixiT12*sumwi*n
  numer.sigma2.hat<-numer.sigma2.hat^2
  denom.sigma2.hat<-(sumwi)^2*n^2*(n-2)*(delta+(sumwixiT1)^2-2*sumwixiT12*sumwi)
  sigma2.hat<-numer.sigma2.hat/denom.sigma2.hat

  numer.rho.hat<- -delta+sumwixiT12*sumwi*n-(sumwixiT1)^2*n+(sumwixiT1)^2
  denom.rho.hat<- sumwixiTxi*n^3*sumwi-sumwixiT12*n^2*sumwi-2*sumwixiTxi*n^2*sumwi+delta+(sumwixiT1)^2*n+sumwixiT12*sumwi*n-(sumwixiT1)^2

  rho.hat<- - numer.rho.hat/denom.rho.hat

  res<-list(mu=mu.hat, sigma2=sigma2.hat, rho=rho.hat)
  return(res)
}

getIntermediate<-function(Xi)
{ 
  n<-length(Xi)
  # x_{c,i}^T * x_{c,i}
  xiTxi<-sum(Xi*Xi)

  # x_{c,i}^T * 1_{n_c}
  one<-rep(1, n)

  xiT1<-sum(Xi*one)

  return(c(xiTxi, xiT1))
}

updatePara3<-function(dat, wi)
{
  n<-ncol(dat)
  res.inter<-t(apply(dat, 1, getIntermediate))
  xiTxi<-res.inter[,1]
  xiT1<-res.inter[,2]

  sumwi<-sum(wi)
  sumwixiT1<-sum(wi*xiT1)
  mu.hat<-sumwixiT1/(n*sumwi)

  aiTai<-xiTxi-2*xiT1*mu.hat+n*mu.hat^2
  aiT12<-xiT1^2+n^2*mu.hat^2-2*n*mu.hat*xiT1

  sumwiaiTai<-sum(wi*aiTai)
  sumwiaiT12<-sum(wi*aiT12)

  rho.hat<-(sumwiaiT12-sumwiaiTai)/((n-1)*(sumwiaiTai))
  sigma2.hat<-sumwiaiTai/(n*sumwi)

  return(list(mu.hat=mu.hat, sigma2.hat=sigma2.hat, rho.hat=rho.hat))
}

