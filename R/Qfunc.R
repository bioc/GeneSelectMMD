# sum_{i=1}^{nGenes} log10(pi.1*f1(xi)+pi.2*f2(xi)+pi.3*f3(xi))

# negative of the M-step's objective function - the part relating to pi's
# Q(Psi)=sum_{i=1}^{nGenes} wi*V(pi)+\sum_{i=1}^{nGenes) wi*U(theta)
"negQFunc" <-
function(theta, dat, nc, nn, n, 
  xcMat, xnMat, xcTxc, xnTxn, xTx, xT1,
  sumw1, sumw2, sumw3, 
  sumw1xcTxc, sumw1xcT1, sumw1xcT1.sq,
  sumw1xnTxn, sumw1xnT1, sumw1xnT1.sq,
  sumw2xTx, sumw2xT1, sumw2xT1.sq,
  sumw3xcTxc, sumw3xcT1, sumw3xcT1.sq,
  sumw3xnTxn, sumw3xnT1, sumw3xnT1.sq
  )
{

  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-theta[1]; 
  # log(variance) of expression levels for cluster 1 for diseased subjects
  tau.c1<-theta[2]; 
  # modified logit of correlation among expression levels for cluster 1 for diseased subjects
  r.c1<-theta[3]; 
  rho.c1<-(exp(r.c1)-1/(nc-1))/(1+exp(r.c1))

  # mean expression level for cluster 1 for normal subjects
  delta.n1<-theta[4]; 
  mu.n1<-mu.c1-exp(delta.n1)
  # log(variance) of expression levels for cluster 1 for normal subjects
  tau.n1<-theta[5]; 
  # modified logit of correlation among expression levels for cluster 1 for normal subjects
  r.n1<-theta[6]; 
  rho.n1<-(exp(r.n1)-1/(nn-1))/(1+exp(r.n1))

  # mean expression level for cluster 2
  mu.2<-theta[7]; 
  # log(variance) of expression levels for cluster 2
  tau.2<-theta[8]; 
  # modified logit of correlation among expression levels for cluster 2
  r.2<-theta[9]; 
  rho.2<-(exp(r.2)-1/(n-1))/(1+exp(r.2))

  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-theta[10]; 
  # log(variance) of expression levels for cluster 3 for diseased subjects
  tau.c3<-theta[11]; 
  # modified logit of correlation among expression levels for cluster 3 for diseased subjects
  r.c3<-theta[12]; 
  rho.c3<-(exp(r.c3)-1/(nc-1))/(1+exp(r.c3))

  # mean expression level for cluster 3 for normal subjects
  delta.n3<-theta[13]; 
  mu.n3<-mu.c3+exp(delta.n3)
  # log(variance) of expression levels for cluster 3 for normal subjects
  tau.n3<-theta[14]; 
  # modified logit of correlation among expression levels for cluster 3 for normal subjects
  r.n3<-theta[15]; 
  rho.n3<-(exp(r.n3)-1/(nn-1))/(1+exp(r.n3))

  ###########
  sumw1acTac<-sumw1xcTxc-2*mu.c1*sumw1xcT1+nc*mu.c1^2*sumw1
  sumw1acT1.sq<-sumw1xcT1.sq-2*nc*mu.c1*sumw1xcT1+nc^2*mu.c1^2*sumw1
 
  sumw1anTan<-sumw1xnTxn-2*mu.n1*sumw1xnT1+nn*mu.n1^2*sumw1
  sumw1anT1.sq<-sumw1xnT1.sq-2*nn*mu.n1*sumw1xnT1+nn^2*mu.n1^2*sumw1
 
  sumw2aTa<-sumw2xTx-2*mu.2*sumw2xT1+n*mu.2^2*sumw2
  sumw2aT1.sq<-sumw2xT1.sq-2*n*mu.2*sumw2xT1+n^2*mu.2^2*sumw2
 
  sumw3acTac<-sumw3xcTxc-2*mu.c3*sumw3xcT1+nc*mu.c3^2*sumw3
  sumw3acT1.sq<-sumw3xcT1.sq-2*nc*mu.c3*sumw3xcT1+nc^2*mu.c3^2*sumw3
 
  sumw3anTan<-sumw3xnTxn-2*mu.n3*sumw3xnT1+nn*mu.n3^2*sumw3
  sumw3anT1.sq<-sumw3xnT1.sq-2*nn*mu.n3*sumw3xnT1+nn^2*mu.n3^2*sumw3

  ##########
  part1.1<- -nc*log(2*pi)/2-nc*tau.c1/2-nc*log(nc)/2+(nc-1)*log(nc-1)/2
  tt<-exp(r.c1)
  if(tt!=Inf)
  {
    part1.1<- part1.1+nc*log(1+exp(r.c1))/2-r.c1/2
  } else {
    part1.1<- part1.1+nc*r.c1/2-r.c1/2
  }

  part1.1<- part1.1 * sumw1

  part1.2<- -(nc-1)*(exp(-tau.c1)+exp(r.c1-tau.c1))/(2*nc)*sumw1acTac
  part1.3<- ((nc-2)*exp(-tau.c1)+exp(r.c1-tau.c1)*(nc-1)-exp(-r.c1-tau.c1))/(2*nc^2)*sumw1acT1.sq

  part1.4<- -nn*log(2*pi)/2-nn*tau.n1/2-nn*log(nn)/2+(nn-1)*log(nn-1)/2
  tt<-exp(r.n1)
  if(tt!=Inf)
  {
    part1.4<- part1.4+nn*log(1+exp(r.n1))/2-r.n1/2
  } else {
    part1.4<- part1.4+nn*r.n1/2-r.n1/2
  }
  part1.4<- part1.4 * sumw1

  tt<-exp(r.n1-tau.n1)
  if(tt!=Inf)
  { 
    part1.5<- -(nn-1)*(exp(-tau.n1)+exp(r.n1-tau.n1))/(2*nn)*sumw1anTan
    part1.6<- ((nn-2)*exp(-tau.n1)+exp(r.n1-tau.n1)*(nn-1)-exp(-r.n1-tau.n1))/(2*nn^2)*sumw1anT1.sq
    part1<-part1.1+part1.2+part1.3+part1.4+part1.5+part1.6
  } else {
    part1<-part1.1+part1.2+part1.3+part1.4
  }
#
#  part1<-part1.1+part1.2+part1.3+part1.4+part1.5+part1.6
#
  ##########
  part2.1<- -n*log(2*pi)/2-n*tau.2/2-n*log(n)/2+(n-1)*log(n-1)/2
  tt<-exp(r.2)
  if(tt!=Inf)
  {  
    part2.1<- part2.1+n*log(1+exp(r.2))/2-r.2/2
  } else {
    part2.1<- part2.1+n*r.2/2-r.2/2
  }
  part2.1<- part2.1 * sumw2

  part2.2<- -(n-1)*(exp(-tau.2)+exp(r.2-tau.2))/(2*n)*sumw2aTa
  part2.3<- ((n-2)*exp(-tau.2)+exp(r.2-tau.2)*(n-1)-exp(-r.2-tau.2))/(2*n^2)*sumw2aT1.sq

  if(!is.na(part2.2))
  { 
    part2<-part2.1+part2.2+part2.3
  } else {
    part2<-part2.1+part2.3
  }

  ##########
  part3.1<- -nc*log(2*pi)/2-nc*tau.c3/2-nc*log(nc)/2+(nc-1)*log(nc-1)/2
  tt<-exp(r.c3)
  if(tt!=Inf)
  {
    part3.1<- part3.1+nc*log(1+exp(r.c3))/2-r.c3/2
  } else {
    part3.1<- part3.1+nc*r.c3/2-r.c3/2
  }
  part3.1<- part3.1 * sumw3

  part3.2<- -(nc-1)*(exp(-tau.c3)+exp(r.c3-tau.c3))/(2*nc)*sumw3acTac
  part3.3<- ((nc-2)*exp(-tau.c3)+exp(r.c3-tau.c3)*(nc-1)-exp(-r.c3-tau.c3))/(2*nc^2)*sumw3acT1.sq

  part3.4<- -nn*log(2*pi)/2-nn*tau.n3/2-nn*log(nn)/2+(nn-1)*log(nn-1)/2
  tt<-exp(r.n3)
  if(tt!=Inf)
  {
    part3.4<- part3.4+nn*log(1+exp(r.n3))/2-r.n3/2
  } else {
    part3.4<- part3.4+nn*r.n3/2-r.n3/2
  }
  part3.4<- part3.4 * sumw3

  part3.5<- -(nn-1)*(exp(-tau.n3)+exp(r.n3-tau.n3))/(2*nn)*sumw3anTan
  part3.6<- ((nn-2)*exp(-tau.n3)+exp(r.n3-tau.n3)*(nn-1)-exp(-r.n3-tau.n3))/(2*nn^2)*sumw3anT1.sq

  part3<-part3.1+part3.2+part3.3+part3.4+part3.5+part3.6

  #Qfunc<-part0+part1+part2+part3
  Qfunc<-part1+part2+part3

  if(!is.na(Qfunc))
  { 
    return(-Qfunc)
  } else {
    return(1.0e+308)
  }
}

# first derivatives of the negative of the objective function
# Q(Psi)=sum_{i=1}^{nGenes} wi*V(pi)+\sum_{i=1}^{nGenes) wi*U(theta)
"dnegQFunc" <-
function(theta, dat, nc, nn, n,  
  xcMat, xnMat, xcTxc, xnTxn, xTx, xT1,
  sumw1, sumw2, sumw3, 
  sumw1xcTxc, sumw1xcT1, sumw1xcT1.sq,
  sumw1xnTxn, sumw1xnT1, sumw1xnT1.sq,
  sumw2xTx, sumw2xT1, sumw2xT1.sq,
  sumw3xcTxc, sumw3xcT1, sumw3xcT1.sq,
  sumw3xnTxn, sumw3xnT1, sumw3xnT1.sq
)
{

  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-theta[1]; 
  # log(variance) of expression levels for cluster 1 for diseased subjects
  tau.c1<-theta[2]; 
  # modified logit of correlation among expression levels for cluster 1 for diseased subjects
  r.c1<-theta[3]; 
  rho.c1<-(exp(r.c1)-1/(nc-1))/(1+exp(r.c1))

  # mean expression level for cluster 1 for normal subjects
  delta.n1<-theta[4]; 
  mu.n1<-mu.c1-exp(delta.n1)
  # log(variance) of expression levels for cluster 1 for normal subjects
  tau.n1<-theta[5]; 
  # modified logit of correlation among expression levels for cluster 1 for normal subjects
  r.n1<-theta[6]; 
  rho.n1<-(exp(r.n1)-1/(nn-1))/(1+exp(r.n1))

  # mean expression level for cluster 2
  mu.2<-theta[7]; 
  # log(variance) of expression levels for cluster 2
  tau.2<-theta[8]; 
  # modified logit of correlation among expression levels for cluster 2
  r.2<-theta[9]; 
  rho.2<-(exp(r.2)-1/(n-1))/(1+exp(r.2))

  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-theta[10]; 
  # log(variance) of expression levels for cluster 3 for diseased subjects
  tau.c3<-theta[11]; 
  # modified logit of correlation among expression levels for cluster 3 for diseased subjects
  r.c3<-theta[12]; 
  rho.c3<-(exp(r.c3)-1/(nc-1))/(1+exp(r.c3))

  # mean expression level for cluster 3 for normal subjects
  delta.n3<-theta[13]; 
  mu.n3<-mu.c3+exp(delta.n3)
  # log(variance) of expression levels for cluster 3 for normal subjects
  tau.n3<-theta[14]; 
  # modified logit of correlation among expression levels for cluster 3 for normal subjects
  r.n3<-theta[15]; 
  rho.n3<-(exp(r.n3)-1/(nn-1))/(1+exp(r.n3))

  ###########
  sumw1acTac<-sumw1xcTxc-2*mu.c1*sumw1xcT1+nc*mu.c1^2*sumw1
  sumw1acT1.sq<-sumw1xcT1.sq-2*nc*mu.c1*sumw1xcT1+nc^2*mu.c1^2*sumw1
 
  sumw1anTan<-sumw1xnTxn-2*mu.n1*sumw1xnT1+nn*mu.n1^2*sumw1
  sumw1anT1.sq<-sumw1xnT1.sq-2*nn*mu.n1*sumw1xnT1+nn^2*mu.n1^2*sumw1
 
  sumw2aTa<-sumw2xTx-2*mu.2*sumw2xT1+n*mu.2^2*sumw2
  sumw2aT1.sq<-sumw2xT1.sq-2*n*mu.2*sumw2xT1+n^2*mu.2^2*sumw2
 
  sumw3acTac<-sumw3xcTxc-2*mu.c3*sumw3xcT1+nc*mu.c3^2*sumw3
  sumw3acT1.sq<-sumw3xcT1.sq-2*nc*mu.c3*sumw3xcT1+nc^2*mu.c3^2*sumw3
 
  sumw3anTan<-sumw3xnTxn-2*mu.n3*sumw3xnT1+nn*mu.n3^2*sumw3
  sumw3anT1.sq<-sumw3xnT1.sq-2*nn*mu.n3*sumw3xnT1+nn^2*mu.n3^2*sumw3


  ##########
  dQdmu.c1<-(exp(-tau.c1)+exp(-r.c1-tau.c1))*(sumw1xcT1-nc*mu.c1*sumw1)/nc
  dQdmu.c1<-dQdmu.c1+(exp(-tau.n1)+exp(-r.n1-tau.n1))*(sumw1xnT1-nn*mu.n1*sumw1)/nn

  dQddelta.n1<- -exp(delta.n1)*(exp(-tau.n1)+exp(-r.n1-tau.n1))*(sumw1xnT1-nn*mu.n1*sumw1)/nn

  dQdmu.2<-(exp(-tau.2)+exp(-r.2-tau.2))*(sumw2xT1-n*mu.2*sumw2)/n

  dQdmu.c3<-(exp(-tau.c3)+exp(-r.c3-tau.c3))*(sumw3xcT1-nc*mu.c3*sumw3)/nc
  dQdmu.c3<-dQdmu.c3+(exp(-tau.n3)+exp(-r.n3-tau.n3))*(sumw3xnT1-nn*mu.n3*sumw3)/nn

  dQddelta.n3<- exp(delta.n3)*(exp(-tau.n3)+exp(-r.n3-tau.n3))*(sumw3xnT1-nn*mu.n3*sumw3)/nn

  #####
  dQdtau.c1<- -nc*sumw1/2+(nc-1)*(exp(-tau.c1)+exp(r.c1-tau.c1))*sumw1acTac/(2*nc)
  dQdtau.c1<- dQdtau.c1 - ((nc-2)*exp(-r.c1)+(nc-1)-exp(-2*r.c1))*sumw1acT1.sq/(2*nc^2)

  dQdtau.n1<- -nn*sumw1/2+(nn-1)*(exp(-tau.n1)+exp(r.n1-tau.n1))*sumw1acTac/(2*nn)
  dQdtau.n1<- dQdtau.n1 - ((nn-2)*exp(-r.n1)+(nn-1)-exp(-2*r.n1))*sumw1acT1.sq/(2*nn^2)

  dQdtau.2<- -n*sumw2/2+(n-1)*(exp(-tau.2)+exp(r.2-tau.2))*sumw2aTa/(2*n)
  tt<-exp(-r.2)
  if(tt!=Inf)
  {
    dQdtau.2<- dQdtau.2 - ((n-2)*exp(-r.2)+(n-1)-exp(-2*r.2))*sumw2aT1.sq/(2*n^2)
  } 

  dQdtau.c3<- -nc*sumw3/2+(nc-1)*(exp(-tau.c3)+exp(r.c3-tau.c3))*sumw3acTac/(2*nc)
  dQdtau.c3<- dQdtau.c3 - ((nc-2)*exp(-r.c3)+(nc-1)-exp(-2*r.c3))*sumw3acT1.sq/(2*nc^2)

  dQdtau.n3<- -nn*sumw3/2+(nn-1)*(exp(-tau.n3)+exp(r.n3-tau.n3))*sumw3acTac/(2*nn)
  dQdtau.n3<- dQdtau.n3 - ((nn-2)*exp(-r.n3)+(nn-1)-exp(-2*r.n3))*sumw3acT1.sq/(2*nn^2)

  #####
  tt<-exp(r.c1)
  if(tt!=Inf)
  {
    dQdr.c1<-(nc*exp(r.c1)/(2*(1+exp(r.c1)))-1/2)*sumw1
  } else {
    dQdr.c1<-(nc/(2*(1+exp(-r.c1)))-1/2)*sumw1
  }
  dQdr.c1<-dQdr.c1-(nc-1)*exp(r.c1-tau.c1)*sumw1acTac/(2*nc)
  dQdr.c1<-dQdr.c1+(exp(r.c1-tau.c1)*(nc-1)+exp(-r.c1-tau.c1))*sumw1acT1.sq/(2*nc^2)

  #####
  tt<-exp(r.n1)
  if(tt!=Inf)
  { 
    dQdr.n1<-(nn*exp(r.n1)/(2*(1+exp(r.n1)))-1/2)*sumw1
  } else {
    dQdr.n1<-(nn/(2*(exp(-r.n1)+1))-1/2)*sumw1
  }

  dQdr.n1<-dQdr.n1-(nn-1)*exp(r.n1-tau.n1)*sumw1anTan/(2*nn)
  dQdr.n1<-dQdr.n1+(exp(r.n1-tau.n1)*(nn-1)+exp(-r.n1-tau.n1))*sumw1anT1.sq/(2*nn^2*exp(tau.n1))

  #####
  tt<-exp(r.2)
  if(tt!=Inf)
  { 
    dQdr.2<-(n*exp(r.2)/(2*(1+exp(r.2)))-1/2)*sumw2
  } else {
    dQdr.2<-(n/(2*(exp(-r.2)+1))-1/2)*sumw2
  }
  dQdr.2<-dQdr.2-(n-1)*exp(r.2-tau.2)*sumw2aTa/(2*n)
  dQdr.2<-dQdr.2+(exp(r.2-tau.2)*(n-1)+exp(-r.2-tau.2))*sumw2aT1.sq/(2*n^2)

  #####
  tt<-exp(r.c3)
  if(tt!=Inf)
  {
    dQdr.c3<-(nc*exp(r.c3)/(2*(1+exp(r.c3)))-1/2)*sumw3
  } else {
    dQdr.c3<-(nc/(2*(1+exp(-r.c3)))-1/2)*sumw3
  }
  dQdr.c3<-dQdr.c3-(nc-1)*exp(r.c3-tau.c3)*sumw3acTac/(2*nc)
  dQdr.c3<-dQdr.c3+(exp(r.c3-tau.c3)*(nc-1)+exp(-r.c3-tau.c3))*sumw3acT1.sq/(2*nc^2)

  #####
  tt<-exp(r.n3)
  if(tt!=Inf)
  {
    dQdr.n3<-(nn*exp(r.n3)/(2*(1+exp(r.n3)))-1/2)*sumw3
  } else {
    dQdr.n3<-(nn/(2*(1+exp(-r.n3)))-1/2)*sumw3
  }
  dQdr.n3<-dQdr.n3-(nn-1)*exp(r.n3-tau.n3)*sumw3anTan/(2*nn)
  dQdr.n3<-dQdr.n3+(exp(r.n3-tau.n3)*(nn-1)+exp(-r.n3-tau.n3))*sumw3anT1.sq/(2*nn^2)

  #####
  dQfunc<-c(
    dQdmu.c1, dQdtau.c1, dQdr.c1, dQddelta.n1, dQdtau.n1, dQdr.n1, 
    dQdmu.2, dQdtau.2, dQdr.2,
    dQdmu.c3, dQdtau.c3, dQdr.c3, dQddelta.n3, dQdtau.n3, dQdr.n3)

  pos.Inf<-which(dQfunc==Inf)
  pos.nInf<-which(dQfunc== -Inf)
  pos.na<-which(is.na(dQfunc)== TRUE)
  pos<-unique(c(pos.Inf, pos.nInf, pos.na))
  if(length(pos))
  { dQfunc[pos]<- 0 }

  return(-dQfunc)
}


