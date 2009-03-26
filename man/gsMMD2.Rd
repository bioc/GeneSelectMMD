\name{gsMMD2}
\alias{gsMMD2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gene selection based on a mixture of marginal distributions}
\description{
Gene selection based on the marginal distributions of gene profiles that characterized by a mixture of three-component multivariate distributions.
Input is an object derived from the class \code{ExpressionSet}. The user needs to provide initial gene cluster membership.
}
\usage{
gsMMD2(obj.eSet, 
       memSubjects, 
       memIni,
       maxFlag = TRUE, 
       thrshPostProb = 0.5, 
       geneNames = NULL, 
       alpha = 0.05, 
       transformFlag = FALSE, 
       transformMethod = "boxcox", 
       scaleFlag = FALSE, 
       if.center = TRUE, 
       if.scale = TRUE, 
       criterion = c("cor", "skewness", "kurtosis"), 
       minL = -10, 
       maxL = 10, 
       stepL = 0.1, 
       eps = 0.001, 
       ITMAX = 100, 
       plotFlag = FALSE,
       quiet=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{obj.eSet}{an object derived from the class \code{ExpressionSet} which contains the matrix of gene expression levels.
              The rows of the matrix are genes. The columns of the matrix
              are subjects.}
  \item{memSubjects}{a vector of membership of subjects. \code{memSubjects[i]=1} means that the \eqn{i}-th subject belongs to diseased group, \eqn{0} otherwise. 
  }
  \item{memIni}{a vector of user-provided gene cluster membership.}
  \item{maxFlag}{logical. Indicate how to assign gene class membership. \code{maxFlag}=TRUE means that a gene will be assigned
to a class in which the posterior probability of the gene belongs to this class is maximum. \code{maxFlag}=FALSE means that a gene will be assigned to class 1 if the posterior probability of the gene belongs to class 1 is
greater than \code{thrshPostProb}. Similarly, a gene will be assigned to class 1 if the posterior probability of the gene belongs to class 1 is
greater than \code{thrshPostProb}. If  the posterior probability is less than \code{thrshPostProb}, the gene
will be assigned to class 2 (non-differentially expressed gene group).}
  \item{thrshPostProb}{threshold for posterior probabilities. For example, if the posterior probability that a gene belongs to cluster 1 given its gene expression levels is larger than \code{thrshPostProb}, then this gene will be assigned to cluster 1.}
  \item{geneNames}{an optional character vector of gene names}
  \item{alpha}{significant level which is equal to \code{1-conf.level}, 
\code{conf.level} is the argument for the function \code{t.test}.
  }
  \item{transformFlag}{logical. Indicate if data transformation is needed}
  \item{transformMethod}{method for transforming data. Available methods include "boxcox", "log2", "log10", "log", "none".}
  \item{scaleFlag}{logical. Indicate if gene profiles are to be scaled. If \code{transformFlag=TRUE} and \code{scaleFlag=TRUE}, then scaling is performed after transformation.}
  \item{if.center}{logical. If \code{scaleFlag=TRUE} and \code{if.center=TRUE}, then each gene profile will be centered to have mean zero.}
  \item{if.scale}{logical. If \code{scaleFlag=TRUE} and \code{if.scale=TRUE}, then each gene profile will be scaled to have variance one.}
  \item{criterion}{if \code{transformFlag=TRUE}, \code{criterion} indicates what criterion to determine if data looks like normal.  \dQuote{cor} means using Pearson's correlation. The idea is that the observed quantiles after transformation should be close to theoretical normal quantiles. So we can use Pearson's correlation to check if the scatter plot of theoretical normal quantiles versus observed quantiles is a straightline.  \dQuote{skewness} means using skewness measure to check if the distribution of the transformed data are close to normal distribution; \dQuote{kurtosis} means using kurtosis measure to check normality.}
  \item{minL}{lower limit for the \code{lambda} parameter used in Box-Cox transformation}
  \item{maxL}{upper limit for the \code{lambda} parameter used in Box-Cox transformation}
  \item{stepL}{step increase when searching the optimal \code{lambda} parameter used in Box-Cox transformation}
  \item{eps}{a small positive value. If the absolute value of a value is smaller than \code{eps}, this value is regarded as zero.  }
  \item{ITMAX}{maximum iteration allowed for iterations in the EM algorithm}
  \item{plotFlag}{logical. Indicate if the Box-Cox normality plot should be output.}
  \item{quiet}{logical. Indicate if intermediate results should be printed out.}
}
\details{
We assume that the distribution of gene expression profiles is 
a mixture of 3-component multivariate normal distributions 
\eqn{\sum_{k=1}^{3} \pi_k f_k(x|\theta)}. Each component distribution \eqn{f_k} 
corresponds to a gene cluster. The 3 components correspond to 3 gene clusters:
(1) up-regulated gene cluster, (2) non-differentially expressed gene cluster, 
and (3) down-regulated gene cluster. 
The model parameter vector is \eqn{\theta=(\pi_1}, \eqn{\pi_2}, 
\eqn{\pi_3}, \eqn{\mu_{c1}}, 
\eqn{\sigma^2_{c1}}, \eqn{\rho_{c1}}, \eqn{\mu_{n1}}, \eqn{\sigma^2_{n1}}, 
\eqn{\rho_{n1}},
\eqn{\mu_2}, \eqn{\sigma^2_2}, \eqn{\rho_2}, 
\eqn{\mu_{c3}}, 
\eqn{\sigma^2_{c3}}, \eqn{\rho_{c3}}, \eqn{\mu_{n3}}, \eqn{\sigma^2_{n3}}, 
\eqn{\rho_{n3}}.
where \eqn{\pi_1}, \eqn{\pi_2}, and \eqn{\pi_3} are the mixing proportions; 
\eqn{\mu_{c1}}, \eqn{\sigma^2_{c1}}, and \eqn{\rho_{c1}} are 
the marginal mean, variance, and correlation of gene expression levels 
of cluster 1 (up-regulated genes) for diseased subjects; 
\eqn{\mu_{n1}}, \eqn{\sigma^2_{n1}}, and \eqn{\rho_{n1}} are 
the marginal mean, variance, and correlation of gene expression levels 
of cluster 1 (up-regulated genes) for non-diseased subjects; 
\eqn{\mu_2}, \eqn{\sigma^2_2}, and \eqn{\rho_2} are the marginal mean, 
variance, and correlation of gene 
expression levels of cluster 2 (non-differentially expressed genes); 
\eqn{\mu_{c3}}, \eqn{\sigma^2_{c3}}, and \eqn{\rho_{c3}} are 
the marginal mean, variance, and correlation of gene expression levels 
of cluster 3 (up-regulated genes) for diseased subjects; 
\eqn{\mu_{n3}}, \eqn{\sigma^2_{n3}}, and \eqn{\rho_{n3}} are 
the marginal mean, variance, and correlation of gene expression levels 
of cluster 3 (up-regulated genes) for non-diseased subjects. 

Note that genes in cluster 2 are non-differentially expressed across
abnormal and normal tissue samples. Hence there are only 3 parameters for
cluster 2.

We apply the EM algorithm to estimate the model parameters. 
We regard the cluster membership of genes as missing values.
}
\value{
A list contains 10 elements.
  \item{dat}{the (transformed) microarray data matrix. If tranformation
performed, then \code{dat} will be different from the input 
microarray data matrix.}
  \item{memSubjects}{the same as the input \code{memSubjects}.}
  \item{memGenes}{a vector of cluster membership of genes. \eqn{1} means up-regulated gene; \eqn{2} means non-differentially expressed gene; 
\eqn{3} means down-regulated gene.}
  \item{memGenes2}{an variant of the vector of cluster membership of genes. 
\eqn{1} means differentially expressed gene; \eqn{0} means non-differentially expressed gene.}
  \item{para}{parameter estimates (c.f. details).}
  \item{llkh}{value of the loglikelihood function.}
  \item{wiMat}{posterior probability that a gene belongs to a cluster given the expression levels of this gene. Column i is for cluster i.}
  \item{memIni}{the initial cluster membership of genes.}
  \item{paraIni}{the parameter estimates based on initial gene cluster membership.}
  \item{llkhIni}{the value of loglikelihood function.}
 \item{lambda}{the parameter used to do Box-Cox transformation}
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
\seealso{
\code{\link{gsMMD}},
\code{\link{gsMMD.default}},
\code{\link{gsMMD2.default}}
}
\note{
The speed of the program is slow for large data sets.
}
\examples{
  library(ALL)
  data(ALL)
  eSet1 <- ALL[1:100, ALL$BT == "B3" | ALL$BT == "T2"]
  
  mem.str <- as.character(eSet1$BT)
  nSubjects <- length(mem.str)
  memSubjects <- rep(0,nSubjects)
  # B3 coded as 0, T2 coded as 1
  memSubjects[mem.str == "T2"] <- 1
  
  myWilcox <-
  function(x, memSubjects, alpha = 0.05)
  {
    xc <- x[memSubjects == 1]
    xn <- x[memSubjects == 0]
  
    m <- sum(memSubjects == 1)
    res <- wilcox.test(x = xc, y = xn, conf.level = 1 - alpha)
    res2 <- c(res$p.value, res$statistic - m * (m + 1) / 2)
    names(res2) <- c("p.value", "statistic")
  
    return(res2)
  }
  
  mat <- exprs(eSet1)
  tmp <- t(apply(mat, 1, myWilcox, memSubjects = memSubjects))
  colnames(tmp) <- c("p.value", "statistic")
  memIni <- rep(2, nrow(mat))
  memIni[tmp[, 1] < 0.05 & tmp[, 2] > 0] <- 1
  memIni[tmp[, 1] < 0.05 & tmp[, 2] < 0] <- 3
  
  cat("initial gene cluster size>>\n"); print(table(memIni)); cat("\n");

  obj.gsMMD <- gsMMD2(eSet1, memSubjects, memIni, transformFlag = TRUE, 
       transformMethod = "boxcox", scaleFlag = TRUE, quiet = FALSE)
  round(obj.gsMMD$para, 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif  }
