"gsMMD2"<-
function(obj.eSet, 
         memSubjects, 
         memIni,
         maxFlag=TRUE, 
         thrshPostProb=0.50, 
         geneNames=NULL, 
         alpha=0.05, 
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

  res<-gsMMD2.default(X, 
         memSubjects, 
         memIni,
         maxFlag, 
         thrshPostProb, 
         geneNames, 
         alpha, 
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


"gsMMD2.default"<-
function(X, 
         memSubjects, 
         memIni,
         maxFlag=TRUE, 
         thrshPostProb=0.50, 
         geneNames=NULL, 
         alpha=0.05, 
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

  nGenes<-nrow(X)
  nSubjects<-ncol(X)

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

  cat("Programming is running. Please be patient...\n")
  # records initial parameter estimates
  tmpIni<-getPara(X, memSubjects, memIni, eps)
  paraIni<-tmpIni$para

  # records initial log-likelihood estimates
  llkhIni<-tmpIni$llkh

  # records E(z_{ij} | x_i, Psi^{(m)})
  wiMat<-matrix(0, nrow=nGenes, ncol=3)
  {
    # Gene Selection based on EM algorithm
    res<- paraEst(X, paraIni, memSubjects=memSubjects, 
                   maxFlag=maxFlag, thrshPostProb, geneNames=geneNames,
                   ITMAX=ITMAX, eps=eps)

    if(res$loop==0)
    {
      para<-paraIni
      llkh<-llkhIni
      memGenes<-memIni
      wiMat<-res$wiMat
    } else {
      para<-res$para
      llkh<-res$llkh
      memGenes<-res$memGenes
      wiMat<-res$wiMat
    }
  }

  names(memGenes)<-geneNames
  rownames(wiMat)<-geneNames
  colnames(wiMat)<-paste("cluster", 1:3, sep="")

  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0 # non-differentially expressed genes

  if(!quiet)
  {
    cat("*******************************************************\n\n")
    cat("Initial parameter estimates>>\n"); print(round(paraIni,3)); cat("\n");
    cat("Initial loglikelihood>>\n"); print(round(llkhIni,3)); cat("\n");
    cat("Final parameter estimates>>\n"); print(round(para,3)); cat("\n");
    cat("Final loglikelihood>>\n"); print(round(llkh,3)); cat("\n");
    cat("*******************************************************\n\n")
  }
  
  res<-list(dat=X, memSubjects=memSubjects, 
            memGenes=memGenes, memGenes2=memGenes2, 
            para=para, 
            llkh=llkh, wiMat=wiMat, 
            memIni=memIni, paraIni=paraIni, llkhIni=llkhIni,
            lambda=lambda)
  invisible(res) 
}

