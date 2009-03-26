\name{plotHistDensity}
\alias{plotHistDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot of histogram and density estimate of the pooled gene expression levels. }
\description{
Plot of histogram of pooled gene expression levels, composited with density estimate based on the mixture of marginal distributions. The density estimate is based on the assumption that the marginal correlations between subjects are zero.
}
\usage{
plotHistDensity(obj.gsMMD,
                plotFlag="case",
                plotComponent=FALSE,
                myxlab="expression level",
                myylab="density",
                mytitle="Histogram of gene expression levels\nimposed with estimated density (case)",
                x.legend=NULL,
                y.legend=NULL,
                numPoints=500,
                mycol=1:4, 
                mylty=1:4, 
                mylwd=rep(3,4), 
                cex.main=2, 
                cex.lab=1.5, 
                cex.axis=1.5, 
                cex=2,
                bty="n")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.gsMMD}{an object returned by \code{gsMMD},
\code{gsMMD.default}, \code{gsMMD2},
or \code{gsMMD2.default}
}
  \item{plotFlag}{logical. Indicate the plot will based on which type of subjects.}
  \item{plotComponent}{logical. Indicate if components of the mixture of marginal distribution will be plotted.}
  \item{myxlab}{label for x-axis}
  \item{myylab}{label for y-axis}
  \item{mytitle}{title of the plot}
  \item{x.legend}{the x-corrdiates of the legend}
  \item{y.legend}{the y-corrdiates of the legend}
  \item{numPoints}{logical. Indicate how many genes will be plots.}
  \item{mycol}{color for the density estimates (overall and components)}
  \item{mylty}{line styles for the density estimates (overall and components)}
  \item{mylwd}{line width for the density estimates (overall and components)}
  \item{cex.main}{font for main title}
  \item{cex.lab}{font for x- and y-axis labels}
  \item{cex.axis}{font for x- and y-axis}
  \item{cex}{font for texts}
  \item{bty}{the type of box to be drawn around the legend.  The allowed
          values are '"o"' and '"n"' (the default).}
}
\details{
For a given type of subjects, we pool their expression levels together if
the marginal correlations among subjects are zero. We then draw a histogram
of the pooled expression levels. Next, we composite density estimates of gene
expression levels for the overal distribution and the 3 component distributions.
}
\value{
A list containing coordinates of the density estimates:
  \item{x }{sorted pooled gene expression levels for cases or controls.}
  \item{x2 }{a subset of \code{x} specified by the sequence:
\code{seq(from=1,to=len.x, by=delta)}, where \code{len.x} is the
length of the vector \code{x}, and \code{delta=floor(len.x/numpoints).}}
  \item{y }{density estimate corresponding to \code{x2}}
  \item{y1 }{weighted density estimate for gene cluster 1}
  \item{y2 }{weighted density estimate for gene cluster 2}
  \item{y3 }{weighted density estimate for gene cluster 3}
}
\references{ Qiu, W.-L., He, W., Wang, X.-G. and Lazarus, R. (2008). 
A Marginal Mixture Model for Selecting Differentially Expressed Genes across Two Types of Tissue Samples. \emph{The International Journal of Biostatistics. 4(1):Article 20.} \url{http://www.bepress.com/ijb/vol4/iss1/20}
}
\author{ 
Weiliang Qiu \email{stwxq@channing.harvard.edu},
Wenqing He \email{whe@stats.uwo.ca},
Xiaogang Wang \email{stevenw@mathstat.yorku.ca},
Ross Lazarus \email{ross.lazarus@channing.harvard.edu}
}
\note{ The density estimate is obtained based on the
assumption that the marginal correlation among
subjects is zero. If the estimated marginal correlation obtained by \code{gsMMD} is far from zero, then do not use this plot function. 
}
\examples{
  library(ALL)
  data(ALL)
  eSet1 <- ALL[1:100, ALL$BT == "B3" | ALL$BT == "T2"]
  
  mem.str <- as.character(eSet1$BT)
  nSubjects <- length(mem.str)
  memSubjects <- rep(0,nSubjects)
  # B3 coded as 0 (control), T2 coded as 1 (case)
  memSubjects[mem.str == "T2"] <- 1
  
  obj.gsMMD <- gsMMD(eSet1, memSubjects, transformFlag = TRUE, 
    transformMethod = "boxcox", scaleFlag = TRUE, quiet = FALSE)

  plotHistDensity(obj.gsMMD, plotFlag = "case", 
      mytitle = "Histogram of gene expression levels for T2\nimposed with estimated density (case)", 
      plotComponent = TRUE, 
      x.legend = c(0.8, 3), 
      y.legend = c(0.3, 0.4), 
      numPoints = 500)

}
\keyword{classif }% at least one, from doc/KEYWORDS
