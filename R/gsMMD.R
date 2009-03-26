# define some constant
PI<-3.1415926
TNumPara<-(18)
BigNum<-(1.0e+10)
paraNames<- c("pi.1", "pi.2", "pi.3", 
              "mu.c1", "sigma2.c1", "rho.c1", 
              "mu.n1", "sigma2.n1", "rho.n1",
              "mu.2", "sigma2.2", "rho.2",
              "mu.c3", "sigma2.c3", "rho.c3", 
              "mu.n3", "sigma2.n3", "rho.n3")

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
         scaleFlag=FALSE, 
         if.center=TRUE, 
         if.scale=TRUE,
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
         if.center,
         if.scale,
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
         scaleFlag=FALSE, 
         if.center=TRUE, 
         if.scale=TRUE,
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

  nGenes<-nrow(X)
  nSubjects<-ncol(X)
  nMethods<-length(iniGeneMethod)

  if(sum(is.null(geneNames)))
  { geneNames<-paste("gene", 1:nGenes, sep="") }

  cat("Programming is running. Please be patient...\n")
  lambda<-NA
  if(transformFlag)
  { 
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
    X<-t(apply(X, 1, scale, center=if.center, scale=if.scale))
  }

  # records initial parameter estimates
  paraIniMat<-matrix(0, nrow=TNumPara, ncol=nMethods)
  # records initial gene-membership estimates
  memIniMat<-matrix(0, nrow=nGenes, ncol=nMethods)

  # records initial log-likelihood estimates
  llkhIniVec<-rep(0, nMethods)

  # records parameter estimates
  paraMat<-matrix(0, nrow=TNumPara, ncol=nMethods)
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
    paraIniMat[,i]<-tmpIni$para
    llkhIniVec[i]<-tmpIni$llkh
    ttt<-paraIniMat[,i]
    names(ttt)<-paraNames
    if(!quiet)
    { cat("paraIniMat[,i]>>\n"); print(round(ttt,3)); cat("\n"); }

    # Gene Selection based on EM algorithm
    res<- paraEst(X, paraIniMat[,i], memSubjects=memSubjects, 
                   maxFlag=maxFlag, thrshPostProb, geneNames=geneNames,
                   ITMAX=ITMAX, eps=eps)

    if(res$loop==0)
    {
      paraMat[,i]<-paraIniMat[,i]
      llkhVec[i]<-llkhIniVec[i]
      memMat[,i]<-memIniMat[,i]
      wiArray[,,i]<-res$wiMat
    } else {
      paraMat[,i]<-res$para
      llkhVec[i]<-res$llkh
      memMat[,i]<-res$memGenes
      wiArray[,,i]<-res$wiMat
    }
  }
  rownames(paraIniMat)<-paraNames
  colnames(paraIniMat)<-iniGeneMethod

  rownames(memIniMat)<-geneNames
  colnames(memIniMat)<-iniGeneMethod
  names(llkhIniVec)<-iniGeneMethod
  rownames(paraMat)<-paraNames

  colnames(paraMat)<-iniGeneMethod
  rownames(memMat)<-geneNames
  colnames(memMat)<-iniGeneMethod
  names(llkhVec)<-iniGeneMethod
  dimnames(wiArray)<-list(geneNames,
                          paste("cluster", 1:3, sep=""),
                          iniGeneMethod)                          


  # final results
  flagPi<-rep(0, nMethods)
  for(i in 1:nMethods)
  { 
    flagPi[i]<-sum(paraMat[2,i]>paraMat[1,i] & paraMat[2,i]>paraMat[3,i])
  }
  tmppos<-which(flagPi==0)
  if(length(tmppos))
  { llkhVec[tmppos]<- -Inf }
  if(!quiet)
  { cat("llkhVec>>\n"); print(llkhVec); cat("\n"); }
  pos<-which(llkhVec==max(llkhVec))
  len<-length(pos)
  tt<-sample(1:len, 1, rep=FALSE)
  pos<-pos[tt]

  memGenes<-as.vector(memMat[,pos])
  para<-paraMat[,pos]
  llkh<-llkhVec[pos]
  wiMat<-wiArray[,,pos]

  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0 # non-differentially expressed genes

  if(!quiet)
  { cat("*******************************************************\n\n") 
    cat("Initial parameter estimates>>\n"); print(round(paraIniMat,3)); cat("\n");
    cat("Initial loglikelihood>>\n"); print(round(llkhIniVec,3)); cat("\n");
    cat("Final parameter estimates based on initial estimates>>\n"); print(round(paraMat,3)); cat("\n");
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
  } else { #if(iniGeneMethod=="Wilcox"){
    # (2) two-sample wilcoxon test
    tmp<-iniMemGenesTestFunc(X, memSubjects=memSubjects, testFun=myWilcox, 
                        geneNames=geneNames, alpha = alpha, eps=eps)
  }
   
  res<-list(para=tmp$para, llkh=tmp$llkh, memGenes=tmp$memGenes, 
            memGenes2=tmp$memGenes2, iniGeneMethod=iniGeneMethod)
  return(res)
}

