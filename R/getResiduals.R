# replace expression levels by the residuals of regression analysis 
# in which predictor of interest is not in the regression model.
# The purpose is to remove potential confounding factors
obtainResi<-function(es, fmla)
{
 
  dat<-exprs(es)
  pDat<-pData(es)
  pDat2 <- pDat
  rownames(pDat2) <- 1:nrow(pDat2)
  designMat <- model.matrix(fmla, pDat2)
  rn <- as.numeric(rownames(designMat))
  es2 <- es[, rn]
  dat2 <- dat[, rn, drop = FALSE]
  # use pDat instead of pDat2 to keep 
  # original row names 
  pDat3 <- pDat[rn, , drop = FALSE]
  rownames(designMat) <- colnames(dat2)

  fit = lmFit(dat2, designMat)
  ebFit = eBayes(fit)

  # obtain residual matrix
  resMat<-lapply(1:nrow(dat2), function(i) {
      yi<-dat2[i,]
      betai<-ebFit$coefficients[i,]
      ri<-yi - designMat%*%betai
      return(ri)
    }
  )
  
  resMat2<-t(sapply(resMat, function(x) {x}))
  rownames(resMat2)<-rownames(dat2)
  colnames(resMat2)<-colnames(dat2)
  
  exprs(es2)<-resMat2
  pData(es2)<-pDat3

  invisible(es2)
}


