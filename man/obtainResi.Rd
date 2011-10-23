\name{obtainResi}
\alias{obtainResi}
\title{
Replace expression levels by the residuals of
regression analysis to remove the confounding effects.
}
\description{
Replace expression levels by the residuals of regression analysis
in which predictor of interest is not in the regression model.
The purpose of this function is to remove potential confounding factors.
}
\usage{
obtainResi(es, fmla)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
An \code{ExpressionSet} object.
}
  \item{fmla}{
A formula object that specifies the covariates of the linear
regression model. The variable of interest should not be 
included. No response variable should be specified in \code{fmla} since
the response variable is always the expression level. See 
function \code{\link{lmFit}} of R Bioconductor package \code{limma}.
}
}
\details{
To remove confounding effects, we can replace the expression level by
the residuals
of  a linear regression model with
response variable the expression level and covariates the potential
confounders. The functions \code{lmFit}  and \code{eBayes} will
be used to obtain regression coefficients.
}
\value{
An ExpressionSet object with expression levels replaced by
residuals of linear regression analysis.
}
\author{
Weiliang Qiu \email{stwxq@channing.harvard.edu},
Wenqing He \email{whe@stats.uwo.ca},
Xiaogang Wang \email{stevenw@mathstat.yorku.ca},
Ross Lazarus \email{ross.lazarus@channing.harvard.edu}
}
\note{
The number of arrays of the returned ExpressionSet object
might be smaller than that of the original ExpressionSet object,
due to missing values in covariates.
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
