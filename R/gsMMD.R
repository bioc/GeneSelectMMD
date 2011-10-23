# modified on Sept. 28, 2011
#  (1) added 'na.rm=TRUE' to function 'sum'
#
# define some constant
PI<-3.1415926
TNumPara<-18
paraNames<- c("pi.1", "pi.2", "pi.3",
              "mu.c1", "sigma2.c1", "rho.c1", 
              "mu.n1", "sigma2.n1", "rho.n1",
              "mu.2", "sigma2.2", "rho.2",
              "mu.c3", "sigma2.c3", "rho.c3", 
              "mu.n3", "sigma2.n3", "rho.n3")

# re-parametrization
TNumParaRP<-17
paraNamesRP<- c("pi.1", "pi.2", 
              "mu.c1", "tau.c1", "r.c1", 
              "delta.n1", "tau.n1", "r.n1",
              "mu.2", "tau.2", "r.2",
              "mu.c3", "tau.c3", "r.c3", 
              "delta.n3", "tau.n3", "r.n3")


"gsMMD"<-
function(obj.eSet, 
         memSubjects, 
         maxFlag=TRUE, 
         thrshPostProb=0.50, 
         geneNames=NULL, 
         alpha=0.05, 
         iniGeneMethod= "Ttest", 
         transformFlag=FALSE, 
         transformMethod="boxcox", 
         scaleFlag=TRUE, 
         criterion=c("cor", "skewness", "kurtosis"),
         minL=-10, 
         maxL=10, 
         stepL=0.1, 
         eps=1.0e-3, 
         ITMAX=100, 
         plotFlag=FALSE,
         quiet=TRUE)
{
  # get expression level matrix
  X<-exprs(obj.eSet)

  res<-gsMMD.default(X, 
         memSubjects, 
         maxFlag, 
         thrshPostProb, 
         geneNames, 
         alpha, 
         iniGeneMethod,
         transformFlag,
         transformMethod,
         scaleFlag,
         criterion,
         minL,
         maxL,
         stepL,
         eps,
         ITMAX,
         plotFlag,
         quiet)

  invisible(res)

}


"gsMMD.default"<-
function(X, 
         memSubjects, 
         maxFlag=TRUE, 
         thrshPostProb=0.50, 
         geneNames=NULL, 
         alpha=0.05, 
         iniGeneMethod= "Ttest", 
         transformFlag=FALSE, 
         transformMethod="boxcox", 
         scaleFlag=TRUE, 
         criterion=c("cor", "skewness", "kurtosis"),
         minL=-10, 
         maxL=10, 
         stepL=0.1, 
         eps=1.0e-3, 
         ITMAX=100, 
         plotFlag=FALSE,
         quiet=TRUE)
{
  transformMethod<-match.arg(transformMethod, choices=c("boxcox", "log2", "log10", "log", "none"))
  criterion<-match.arg(criterion, c("cor", "skewness", "kurtosis"))

  posMethod<-match(iniGeneMethod, c("Ttest","Wilcox"))
  tmppos<-which(is.na(posMethod)==TRUE)
  if(length(tmppos)>0)
  { msg<-paste("The initial gene partition method(s):", iniGeneMethod[tmppos], " not available!\n")
   stop(msg)
  }

  X<-as.matrix(X)
  nGenes<-nrow(X)
  nSubjects<-ncol(X)
  nMethods<-length(iniGeneMethod)
  nCases<-sum(memSubjects==1, na.rm=TRUE)
  nControls<-sum(memSubjects==0, na.rm=TRUE)

  if(sum(is.null(geneNames), na.rm=TRUE))
  { geneNames<-paste("gene", 1:nGenes, sep="") }

  cat("Programming is running. Please be patient...\n")
  lambda<-NA
  if(transformFlag)
  { 
    if(transformMethod!="none")
    {
      vec<-as.numeric(X)
      min.vec<-min(vec, na.rm=TRUE)
      if(min.vec<0)
      {
        cat("****** Begin Warning ******** \n")
        cat("Warning: Data contains non-positive values! To continue ",
          transformMethod, " transformation,\n")
        cat("We first perform the following transformation:\n")
        cat("x<-x+abs(min(x, na.rm=TRUE))+1\n")
        cat("****** End Warning ******** \n")
    
        X<-X+abs(min.vec)+1
      }
    }
    tmp<-transFunc(X, transformMethod, criterion, 
                   minL, maxL, stepL, eps, plotFlag, ITMAX=0) 
    if(transformMethod=="boxcox")
    { X<-tmp$dat 
      lambda<-tmp$lambda.avg
    }
    else {
     X<-tmp
    }
    if(!quiet)
    { cat(paste("Data transformation (", transformMethod, ") performed\n")) }
  }

  if(scaleFlag)
  {
    if(!quiet)
    { cat("Gene profiles are scaled so that they have mean zero and variance one!\n") }
    X<-t(apply(X, 1, scale, center=TRUE, scale=TRUE))

    # to avoid linear dependence of tissue samples after scaling
    # gene profiles, we delete a tissue sample.
    # We arbitrarily select the tissue sample, which has the biggest label number, 
    # from the tissue sample group that has larger size than the other 
    # tissue sample group. For example, if there are 6 cancer tissue samples 
    # and 10 normal tissue samples, we delete the 10-th normal tissue sample after scaling.

    if(nCases>nControls)
    { 
      pos<-which(memSubjects==1)
      pos2<-pos[nCases]
      X<-X[,-pos2]
      memSubjects<-memSubjects[-pos2]
    } else {
      pos<-which(memSubjects==0)
      pos2<-pos[nControls]
      X<-X[,-pos2]
      memSubjects<-memSubjects[-pos2]
    }
    nCases<-sum(memSubjects==1, na.rm=TRUE)
    nControls<-sum(memSubjects==0, na.rm=TRUE)
    nSubjects<-nCases+nControls
  }

  # records initial parameter estimates
  paraIniMatRP<-matrix(0, nrow=TNumParaRP, ncol=nMethods)
  # records initial gene-membership estimates
  memIniMat<-matrix(0, nrow=nGenes, ncol=nMethods)

  # records initial log-likelihood estimates
  llkhIniVec<-rep(0, nMethods)

  # records parameter estimates
  paraMatRP<-matrix(0, nrow=TNumParaRP, ncol=nMethods)
  # records gene-membership estimates
  memMat<-matrix(0, nrow=nGenes, ncol=nMethods)

  # records log-likelihood estimates
  llkhVec<-rep(0, nMethods)

  cat("Programming is running. Please be patient...\n")
  # records E(z_{ij} | x_i, Psi^{(m)})
  wiArray<-array(0, c(nGenes, 3, nMethods))
  for(i in 1:nMethods)
  {
    if(!quiet)
    { cat("******** initial parameter estimates method>>", 
        iniGeneMethod[i], " *******\n") }
    tmpIni<-getIniMemGenes(X, memSubjects, geneNames, 
                              iniGeneMethod[i], alpha, eps=eps)
    iniGeneMethod[i]<-tmpIni$iniGeneMethod
    memIniMat[,i]<-tmpIni$memGenes
    paraIniMatRP[,i]<-tmpIni$para
    llkhIniVec[i]<-tmpIni$llkh
    ttt<-paraIniMatRP[,i]
    names(ttt)<-paraNamesRP
    if(!quiet)
    { cat("paraIniMatRP[,i]>>\n"); print(round(ttt,3)); cat("\n"); }

    # Gene Selection based on EM algorithm
    res<- paraEst(X, tmpIni$para, memSubjects=memSubjects, 
                   maxFlag=maxFlag, thrshPostProb, geneNames=geneNames,
                   ITMAX=ITMAX, eps=eps, quiet=quiet)

    if(res$loop==0)
    {
      paraMatRP[,i]<-paraIniMatRP[,i]
      llkhVec[i]<-llkhIniVec[i]
      memMat[,i]<-memIniMat[,i]
      wiArray[,,i]<-res$memMat
    } else {
      paraMatRP[,i]<-res$para
      llkhVec[i]<-res$llkh
      memMat[,i]<-res$memGenes
      wiArray[,,i]<-res$memMat
    }
  }
  rownames(paraIniMatRP)<-paraNamesRP
  colnames(paraIniMatRP)<-iniGeneMethod

  rownames(memIniMat)<-geneNames
  colnames(memIniMat)<-iniGeneMethod
  names(llkhIniVec)<-iniGeneMethod
  rownames(paraMatRP)<-paraNamesRP

  colnames(paraMatRP)<-iniGeneMethod
  rownames(memMat)<-geneNames
  colnames(memMat)<-iniGeneMethod
  names(llkhVec)<-iniGeneMethod
  dimnames(wiArray)<-list(geneNames,
                          paste("cluster", 1:3, sep=""),
                          iniGeneMethod)                          


  # final results
  flagPi<-rep(0, nMethods)
  paraIniMat<-matrix(0, nrow=TNumPara, ncol=nMethods)
  rownames(paraIniMat)<-paraNames
  colnames(paraIniMat)<-iniGeneMethod

  paraMat<-matrix(0, nrow=TNumPara, ncol=nMethods)
  rownames(paraMat)<-paraNames
  colnames(paraMat)<-iniGeneMethod

  for(i in 1:nMethods)
  { 
    paraIniMat[,i]<-paraRPConverter(paraIniMatRP[,i], nCases, nControls)
    paraMat[,i]<-paraRPConverter(paraMatRP[,i], nCases, nControls)
    flagPi[i]<-sum(paraMat[2,i]>paraMat[1,i] & paraMat[2,i]>paraMat[3,i], na.rm=TRUE)
  }
  tmppos<-which(flagPi==0)
  if(length(tmppos))
  { llkhVec[tmppos]<- -Inf }
  if(!quiet)
  { cat("llkhVec>>\n"); print(llkhVec); cat("\n"); }
  pos<-which(llkhVec==max(llkhVec, na.rm=TRUE))
  len<-length(pos)
  tt<-sample(1:len, 1, rep=FALSE)
  pos<-pos[tt]

  memGenes<-as.vector(memMat[,pos])
  para<-paraMat[,pos]
  paraRP<-paraMatRP[,pos]
  llkh<-llkhVec[pos]
  wiMat<-wiArray[,,pos]

  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0 # non-differentially expressed genes

  if(!quiet)
  { cat("*******************************************************\n\n") 
    cat("Initial parameter estimates>>\n"); print(round(paraIniMat,3)); cat("\n");
    cat("Initial loglikelihood>>\n"); print(round(llkhIniVec,3)); cat("\n");
    #tmpMat<-matrix(0, nrow=TNumPara, ncol=nMethods)
    #rownames(tmpMat)<-rownames(paraMat)
    #colnames(tmpMat)<-colnames(paraMat)
    #for(i in 1:nMethods)
    #{
    #  #tmpMat[,i]<-paraMat[[i]][,1]
    #  tmpMat[,i]<-paraMat[,1]
    #}
    #cat("Final parameter estimates based on initial estimates>>\n"); print(round(tmpMat,3)); cat("\n");
    cat("Final loglikelihood based on initial estimates>>\n"); print(round(llkhVec,3)); cat("\n");
    cat("Final parameter estimates>>\n"); print(round(para,3)); cat("\n");
    cat("Final loglikelihood>>\n"); print(round(llkh,3)); cat("\n");
    cat("*******************************************************\n\n")
  }

  res<-list(dat=X, memSubjects=memSubjects, 
            memGenes=memGenes, memGenes2=memGenes2, 
            para=para, 
            llkh=llkh, wiMat=wiMat, wiArray=wiArray,
            memIniMat=memIniMat, paraIniMat=paraIniMat, llkhIniVec=llkhIniVec,
            memMat=memMat, paraMat=paraMat, llkhVec=llkhVec, lambda=lambda)
  invisible(res) 
}

getIniMemGenes<-function(X, memSubjects, geneNames, iniGeneMethod="Ttest",
                         alpha=0.05, eps=1e-6)
{
  iniGeneMethod<-match.arg(iniGeneMethod,  choices=c("Ttest","Wilcox"))

  if(iniGeneMethod=="Ttest")
  {
    # (1) two-sample t-test
    tmp<-iniMemGenesTestFunc(X, memSubjects=memSubjects, testFun=myTtest, 
                        geneNames=geneNames, alpha = alpha, eps=eps)
  } else { #iniGeneMethod=="Wilcox"
    # (2) two-sample wilcoxon test
    tmp<-iniMemGenesTestFunc(X, memSubjects=memSubjects, testFun=myWilcox, 
                        geneNames=geneNames, alpha = alpha, eps=eps)
  }
   
  res<-list(para=tmp$para, llkh=tmp$llkh, memGenes=tmp$memGenes, 
            memGenes2=tmp$memGenes2, iniGeneMethod=iniGeneMethod)
  return(res)
}

paraRPConverter<-function(paraRP, nCases, nControls)
{
  nc<-nCases
  nn<-nControls
  n<-nc+nn

  # mixture proportions
  pi.1<-paraRP[1]; pi.2<-paraRP[2]; 
  pi.3<-1-pi.1-pi.2

  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-paraRP[3]; 
  # variance of expression levels for cluster 1 for diseased subjects
  tau.c1<-paraRP[4]; 
  sigma2.c1<-exp(tau.c1)
  # modified logit of correlation among expression levels for cluster 1 for diseased subjects
  r.c1<-paraRP[5]; 
  rho.c1<-(exp(r.c1)-1/(nc-1))/(1+exp(r.c1))

  # mean expression level for cluster 1 for normal subjects
  delta.n1<-paraRP[6]; 
  mu.n1<-mu.c1-exp(delta.n1)
  # variance of expression levels for cluster 1 for normal subjects
  tau.n1<-paraRP[7]; 
  sigma2.n1<-exp(tau.n1)
  # modified logit of correlation among expression levels for cluster 1 for normal subjects
  r.n1<-paraRP[8]; 
  rho.n1<-(exp(r.n1)-1/(nn-1))/(1+exp(r.n1))

  # mean expression level for cluster 2
  mu.2<-paraRP[9]; 
  # variance of expression levels for cluster 2
  tau.2<-paraRP[10]; 
  sigma2.2<-exp(tau.2)
  # modified logit of correlation among expression levels for cluster 2
  r.2<-paraRP[11]; 
  rho.2<-(exp(r.2)-1/(n-1))/(1+exp(r.2))

  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-paraRP[12]; 
  # variance of expression levels for cluster 3 for diseased subjects
  tau.c3<-paraRP[13]; 
  sigma2.c3<-exp(tau.c3)
  # modified logit of correlation among expression levels for cluster 3 for diseased subjects
  r.c3<-paraRP[14]; 
  rho.c3<-(exp(r.c3)-1/(nc-1))/(1+exp(r.c3))

  # mean expression level for cluster 3 for normal subjects
  delta.n3<-paraRP[15]; 
  mu.n3<-mu.c3+exp(delta.n3)
  # variance of expression levels for cluster 3 for normal subjects
  tau.n3<-paraRP[16]; 
  sigma2.n3<-exp(tau.n3)
  # modified logit of correlation among expression levels for cluster 3 for normal subjects
  r.n3<-paraRP[17]; 
  rho.n3<-(exp(r.n3)-1/(nn-1))/(1+exp(r.n3))

  ##############################
  para<-c(pi.1, pi.2, pi.3,
               mu.c1, sigma2.c1, rho.c1,
               mu.n1, sigma2.n1, rho.n1,
               mu.2, sigma2.2, rho.2,
               mu.c3, sigma2.c3, rho.c3,
               mu.n3, sigma2.n3, rho.n3
  )
  names(para)<-paraNames

  return(para)
}


