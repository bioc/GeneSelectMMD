\name{errRates}
\alias{errRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculating FDR, FNDR, FPR, and FNR for a real microarray data set}
\description{
Calculating FDR, FNDR, FPR, and FNR for a real microarray data set based on the mixture of marginal distributions.
}
\usage{
errRates(obj.gsMMD)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.gsMMD}{an object returned by \code{gsMMD},
\code{gsMMD.default}, \code{gsMMD2},
or \code{gsMMD2.default}
}
}
\details{
We first fit the real microarray data set by the mixture of 
marginal distributions. Then we calculate the error rates
based on the posterior distributions of a gene belonging to
a gene cluster given its gene profiles. Please refer to Formula (7)
on the page 6 of the paper listed in the Reference section.
}
\value{
A vector of 4 elements:
  \item{FDR }{the percentage of nondifferentially expressed genes among selected genes.}
  \item{FNDR }{the percentage of differentially expressed genes among unselected genes.}
  \item{FPR }{the percentage of selected genes among nondifferentially expressed genes}
  \item{FNR }{the percentage of un-selected genes among differentially expressed genes}
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
\examples{
  \dontrun{
    library(ALL)
    data(ALL)
    eSet1 <- ALL[1:100, ALL$BT == "B3" | ALL$BT == "T2"]
    
    mem.str <- as.character(eSet1$BT)
    nSubjects <- length(mem.str)
    memSubjects <- rep(0,nSubjects)
    # B3 coded as 0, T2 coded as 1
    memSubjects[mem.str == "T2"] <- 1
    
    obj.gsMMD <- gsMMD(eSet1, memSubjects, transformFlag = TRUE, 
      transformMethod = "boxcox", scaleFlag = TRUE, quiet = FALSE)
    round(errRates(obj.gsMMD), 3)
  }
}
\keyword{classif }% at least one, from doc/KEYWORDS
