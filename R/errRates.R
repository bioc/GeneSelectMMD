errRates<-function(obj.gsMMD)
{
  # number of genes
  p<-length(obj.gsMMD$memGenes2)
  # number of selected genes
  pr<-sum(obj.gsMMD$memGenes2==1)

  wiMat<-obj.gsMMD$wiMat
  # W2
  W2<-wiMat[,2]

  # calculate FDR
  flag<-wiMat[,2] < (apply(wiMat[,-2], 1, max))
  if(pr>0)
  { FDR<-sum(W2*flag)/pr }
  else
  { FDR<-0 }

  # calculate FNDR
  flag2<- wiMat[,2] == (apply(wiMat, 1, max))  
  pdiff<-p-pr 
  if(pdiff>0)
  { FNDR<-sum((1-W2)*flag2)/pdiff }
  else
  { FNDR<-0 }

  # calculate FPR
  denom<-sum(W2)
  if(denom>0)
  { FPR<-sum(W2*flag)/denom }
  else
  { FPR<-0 }

  # calculate FNR
  denom2<-sum(1-W2)
  if(denom2>0)
  { FNR<-sum((1-W2)*flag2)/denom2 }
  else
  { FNR<-0 }

  estRates<-c(FDR, FNDR, FPR, FNR)
  names(estRates)<-c("FDR", "FNDR", "FPR", "FNR")

  return(estRates)
}

