transFunc<-function(x, transformMethod="boxcox", 
                    criterion=c("cor", "skewness", "kurtosis"), 
                    minL=-10, maxL=10, stepL=0.1, 
                    eps=1.0e-6, plotFlag=FALSE, ITMAX=0)
{
  transformMethod<-match.arg(transformMethod, choices=c("boxcox", "log2", "log10", "log", "none"))
  criterion<-match.arg(criterion, c("cor", "skewness", "kurtosis"))
  tmpx<-as.vector(x)
  minx<-min(tmpx)
  if(minx<0)
  { cat("****** Begin Warning ******** \n")
    cat("Warning: Data contains non-positive values! To continue box-cox transformation,\n")
    cat("We perform the following transformation:\n")
    cat("x<-x+abs(min(x))+1\n")
    cat("****** End Warning ******** \n")
    x<-x+abs(minx)+1 
  }

 
  if(transformMethod=="boxcox")
  { x2<-transBoxCoxMat(x, criterion, minL, maxL, stepL, eps, plotFlag, ITMAX)
    return(x2)
  }
  if(transformMethod=="log2")
  { 
    return(log2(x))
  }
  if(transformMethod=="log10")
  { 
    return(log10(x))
  }
  if(transformMethod=="log")
  { 
    return(log(x))
  }
  return(x)
}


