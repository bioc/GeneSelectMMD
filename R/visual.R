# v0.1.2
#  (1) obj.gsMMD is changed to contain 'memSubjects'
#  (2) fixed a bug relating to the 'ylim' of the plot.
#       new ylim is obtained from the ranges of both the estimated density
#       obtained from the model and the density obtained by 'hist' function
plotHistDensity<-function(obj.gsMMD,
                          plotFlag="case",
                          plotComponent=FALSE,
                          myxlab="expression level",
                          myylab="density",
                          mytitle="Histogram of gene expression levels\nimposed with estimated density (case)",
                          x.legend=NULL, y.legend=NULL,
                          numPoints=500,
                          mycol=1:4, mylty=1:4, mylwd=rep(3,4),
                          cex.main=2, cex.lab=1.5, cex.axis=1.5, cex=2,
                          bty="n")
{
  plotFlag<-match.arg(plotFlag, choices=c("case", "control"))

  dat<-obj.gsMMD$dat
  para<-obj.gsMMD$para
  memSubjects<-obj.gsMMD$memSubjects

  nCases<-sum(memSubjects==1)
  nControls<-sum(memSubjects==0)

  pi.1<-para[1]
  pi.2<-para[2]
  pi.3<-para[3]
  mu.c1<-para[4]
  sigma2.c1<-para[5]
  rho.c1<-para[6]
  mu.n1<-para[7]
  sigma2.n1<-para[8]
  rho.n1<-para[9]
  mu.2<-para[10]
  sigma2.2<-para[11]
  rho.2<-para[12]
  mu.c3<-para[13]
  sigma2.c3<-para[14]
  rho.c3<-para[15]
  mu.n3<-para[16]
  sigma2.n3<-para[17]
  rho.n3<-para[18]



  if(plotFlag=="case")
  { x<-as.vector(dat[,memSubjects==1, drop=FALSE])
  }
  else {
    x<-as.vector(dat[,memSubjects==0, drop=FALSE])
  }
  x<-sort(x)
  len.x<-length(x)
  delta<-floor(len.x/numPoints)
  x2<-x[seq(from=1,to=len.x, by=delta)]

  pi1<-pi.1
  pi2<-pi.2
  pi3<-pi.3
  muc1<-mu.c1
  mun1<-mu.n1
  sigmac12<-sigma2.c1
  sigman12<-sigma2.n1
  xi<-mu.2
  tau2<-sigma2.2
  muc3<-mu.c3
  mun3<-mu.n3
  sigmac32<-sigma2.c3
  sigman32<-sigma2.n3

  if(plotFlag=="case")
  { y1<-pi1*dnorm(x2, mean=muc1, sd=sqrt(sigmac12))
    y2<-pi2*dnorm(x2, mean=xi, sd=sqrt(tau2))
    y3<-pi3*dnorm(x2, mean=muc3, sd=sqrt(sigmac32))
    y<-y1+y2+y3

  } else {
    y1<-pi1*dnorm(x2, mean=mun1, sd=sqrt(sigman12))
    y2<-pi2*dnorm(x2, mean=xi, sd=sqrt(tau2))
    y3<-pi3*dnorm(x2, mean=mun3, sd=sqrt(sigman32))
    y<-y1+y2+y3
  }

  tmp<-hist(x, plot=FALSE)
  myylim<-range(c(y, tmp$density))
  hist(x, freq=FALSE, main=mytitle,xlab=myxlab, ylab=myylab, ylim=myylim,
    cex.main=cex.main, cex.lab=cex.lab)
  lines(x2, y, col=mycol[1], lty=mylty[1], lwd=mylwd[1])

  if(plotComponent)
  {
    lines(x2, y1, col=mycol[2], lty=mylty[2], lwd=mylwd[2])
    lines(x2, y2, col=mycol[3], lty=mylty[3], lwd=mylwd[2])
    lines(x2, y3, col=mycol[4], lty=mylty[4], lwd=mylwd[2])
  }
  if(is.null(x.legend))
  {
    x.max<-max(x)
    d<-max(x)-min(x)
    tmp1<-x.max-d/3
    x.legend<-c(tmp1, x.max)
  }
  if(is.null(y.legend))
  {
    y.max<-max(y)
    d<-max(y)-min(y)
    tmp1<-y.max-d/4
    y.legend<-c(tmp1, y.max)

  }

  if(plotComponent)
  { legend(x=x.legend, y=y.legend, legend=c("overall","component1", "component2"
, "component3"),
    lty=mylty, col=mycol, lwd=mylwd, cex=cex, bty=bty)
  }
  invisible(list(x=x,x2=x2, y=y, y1=y1, y2=y2, y3=y3))
}




