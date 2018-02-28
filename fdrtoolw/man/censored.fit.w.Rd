\name{censored.fit.w}
\alias{censored.fit.w}
\alias{fndr.cutoff.w}

\title{Fit Null Distribution To Censored Data by Maximum Likelihood}

\usage{
censored.fit.w(x, weight, cutoff, statistic=c("pvalue"))
fndr.cutoff.w(x, weight, statistic=c("pvalue"))
}

\description{
  \code{censored.fit.w} fits  a null distribution 
  to censored data.

  \code{fndr.cutoff.w} finds a suitable cutoff point based on the 
   (approximate) false non-discovery rate (FNDR).
}

\arguments{
  \item{x}{vector of test statistics.}
  \item{weight}{vector of weights.}
  \item{cutoff}{truncation point (this may a single value or a vector).}
  \item{statistic}{type of statistic - normal, correlation, or student t. Only pvalue is accepted.}
}
\details{

   As null model truncated normal, truncated student t or a truncated
   correlation density is assumed.  The truncation point is specified
   by the cutoff parameter.  All data points whose absolute value
   are large than the cutoff point are ignored when fitting the truncated
   null model via maximum likelihood.  The total number of data points is
   only used to estimate the fraction of null values eta0.
}
\value{
  \code{censored.fit.w} returns a matrix whose rows contain the estimated parameters and corresponding errors
  for each cutoff point. 

  \code{fndr.cutoff.w} returns a tentative cutoff point.
}


\seealso{
\code{\link{fdrtool.w}}

\code{censored.fit}, \code{fndr.cutoff} in \code{fdrtool}}

\examples{
# load library
library("fdrtoolw")

data(pvalues,package="fdrtool")
#Random weight
w = runif(length(pvalues), min=0.5, max=1)
co = fndr.cutoff.w(pvalues, w, statistic="pvalue")
co
censored.fit.w(pvalues, w, cutoff=co, statistic="pvalue")
}
\keyword{htest}
