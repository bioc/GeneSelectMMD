
const.b1<-2
const.b2<-2
const.b3<-2
const.b.vec<-c(const.b1, const.b2, const.b3)


# para=(pi1, pi2, 
#       muc1, tau.c1, r.c1,
#       delta.n1, tau.n1, r.n1,
#       mu.2, tau.2, r.2,
#       mu.c3, tau.c3, r.c3,
#       delta.n3, tau.n3, r.n3)
# total  17 parameters
 
"paraEst" <-
function(X, paraIni, memSubjects, maxFlag=TRUE, thrshPostProb=0.50, 
         geneNames=NULL, ITMAX=100, eps=1.0e-03, quiet=TRUE)
{
  # check if parameters are in appropriate ranges
  checkPara(paraIni, eps)

  sumb<-sum(const.b.vec, na.rm=TRUE)
  const.b1.b2<-c(const.b1, const.b2)

  n<-ncol(X) # number of patients/subjects
  nc<-sum(memSubjects==1, na.rm=TRUE)  # number of cases
  nn<-sum(memSubjects==0, na.rm=TRUE)  # number of controls
  if(nc>=n){
    stop("Number of cases >= Total number of patients!\n")
  }
  if(nn>=n){
    stop("Number of controls >= Total number of patients!\n")
  }

  # start em algorithm
  Psi.m<-paraIni
  #cat("Psi.m>>\n"); print(Psi.m); cat("\n");
  nGenes<-nrow(X)

  ##############################

  xcMat<-X[,memSubjects==1,drop=FALSE]
  xnMat<-X[,memSubjects==0,drop=FALSE]

  xcTxc<-apply(xcMat, 1, function(x) {sum(x^2, na.rm=TRUE)})
  xcT1<-apply(xcMat, 1, function(x) {sum(x, na.rm=TRUE)})

  xnTxn<-apply(xnMat, 1, function(x) {sum(x^2, na.rm=TRUE)})
  xnT1<-apply(xnMat, 1, function(x) {sum(x, na.rm=TRUE)})

  xTx<-apply(X, 1, function(x) {sum(x^2, na.rm=TRUE)})
  xT1<-apply(X, 1, function(x) {sum(x, na.rm=TRUE)})


# initialize the wi arrays
  wi1<-rep(0,nrow(X))
  wi2<-rep(0,nrow(X))
  wi3<-rep(0,nrow(X))

  loop<-0

# call the paraestloop subroutine in paraEstLoop.f
# paraEstLoop.f includes the loop previously found in paraEst.R

# replace NAs with zeros and pass this matrix to paraestloop
  Xnna<-X
  Xnna[is.na(Xnna)] <- 0

  res <- .Fortran("paraestloop", w1=as.double(wi1), w2=as.double(wi2), w3=as.double(wi3),
          as.double(Xnna), Psi.m=as.double(Psi.m),
          as.integer(memSubjects),
          as.double(xcTxc), as.double(xcT1),
          as.double(xnTxn), as.double(xnT1),
          as.double(xTx), as.double(xT1),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          as.integer(sum(memSubjects==1, na.rm=TRUE)),
          as.integer(sum(memSubjects==0, na.rm=TRUE)),
          as.integer(ITMAX), as.double(eps),
          as.logical(quiet), loop=as.integer(loop)     )


# extract the parameters, w's and the number of loops required
  Psi.m <- res$Psi.m
  w1<-res$w1
  w2<-res$w2
  w3<-res$w3
  loop<-res$loop
  
  if(!quiet)
  { cat("Total iterations for EM algorithm=", loop, "\n") }


# gene cluster membership
#  memMat was changed to wiMat when Fortran code was included
#  memMat<-t(apply(X, 1, wiFun, Psi.m=Psi.m, memSubjects=memSubjects, eps=eps))

#   probably don't need to create these again because it is created at the start
  wi1<-rep(0,nrow(X))
  wi2<-rep(0,nrow(X))
  wi3<-rep(0,nrow(X))

# checkpara was called within wiFun.R, however this will be handled
# in Fortran, so checking the dimension of Psi.m will not be included
# in the future
#
# eps=1.0e-6 for wiFun

# replace NAs with zeros and pass this matrix to wifun
  Xnna<-X
  Xnna[is.na(Xnna)] <- 0

  res <- .Fortran("wifun", wi1=as.double(wi1), wi2=as.double(wi2), wi3=as.double(wi3), 
            as.double(Xnna), as.double(Psi.m),
            as.integer(memSubjects), as.double(1.0e-6),
            as.integer(nrow(X)), as.integer(ncol(X)),
            as.integer(sum(memSubjects==1, na.rm=TRUE)),
            as.integer(sum(memSubjects==0, na.rm=TRUE))     )


  wiMat<-cbind(res$wi1,res$wi2,res$wi3)

  colnames(wiMat)<-c("cluster1", "cluster2", "cluster3")
  rownames(wiMat)<-geneNames

#  update gene cluster membership
#  using Fortran routine maxposfun (converted maxPosFun in R)
#   memGenes<-apply(memMat, 1, maxPosFun, maxFlag=maxFlag, thrshPostProb=thrshPostProb)
#   memMat was changed to wiMat when Fortran code was included

# initialize memGenes array
  memGenes<-rep(0,as.integer(nrow(X)))

  res <- .Fortran("maxposfun", memGenes=as.integer(memGenes), as.double(wiMat), 
                  as.integer(nrow(X)), as.logical(maxFlag), as.double(thrshPostProb)  )

  memGenes<-res$memGenes


  sumw1<-mean(wiMat[,1], na.rm=TRUE)
  sumw2<-mean(wiMat[,2], na.rm=TRUE)
  sumw3<-mean(wiMat[,3], na.rm=TRUE)

  w1<-as.numeric(wiMat[,1])
  w2<-as.numeric(wiMat[,2])
  w3<-as.numeric(wiMat[,3])

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
 

# call the Fortran subroutine llkh
  llkh <- 0

  res <- .Fortran("llkhfun", llkh=as.double(llkh), 
            as.double(Xnna), as.double(Psi.m),
            as.integer(memSubjects), as.double(eps),
            as.integer(nrow(X)), as.integer(n),
            as.integer(nc),
            as.integer(nn)     )

  llkh <- res$llkh


  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0

  if(sum(is.null(geneNames), na.rm=TRUE))
  {
    geneNames<-paste("gene", 1:nGenes, sep="")
  } 

  names(Psi.m)<-paraNamesRP
  names(memGenes)<-geneNames
  names(memGenes2)<-geneNames

  invisible(list(para=Psi.m, llkh=llkh, memGenes=memGenes, 
                 memGenes2=memGenes2, wiMat=wiMat, loop=loop))
}



