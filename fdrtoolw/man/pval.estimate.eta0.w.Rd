\name{pval.estimate.eta0.w}
\alias{pval.estimate.eta0.w}
\title{Estimate the Proportion of Null p-Values}
\usage{
pval.estimate.eta0.w(p, weight, method=c("smoother","quantile"), lambda=seq(0,0.9,0.05), 
   diagnostic.plot=TRUE, q=0.1)
}
\description{
  \code{pval.estimate.eta0.w} estimates the proportion eta0 of null p-values in a given
  vector of p-values. Currently only (method=)smoother and quantile are accepted.
}

\arguments{
  \item{p}{vector of p-values}
  \item{weight}{vector of weights.}
  \item{method}{algorithm used to estimate the proportion of null p-values.
     Only smoother (default) and quantile are adapted for weights at the moment.}
   \item{lambda}{optional tuning parameter vector needed for "bootstrap"
     and "smoothing" methods (defaults to \code{seq(0,0.9,0.05)})
     - see Storey (2002) and Storey and Tibshirani (2003).}
   \item{diagnostic.plot}{if \code{TRUE} (the default) the histogram of the p-values together
          with the estimate of \code{eta0} null line is plotted.  This is useful to
	  visually check the fit of the estimated proportion of null p-values.}  
   \item{q}{quantile used for estimating eta0 - only if \code{method}="quantile"}
  }

\details{

This quantity \code{eta0}, i.e.  the proportion eta0 of null p-values in a given
  vector of p-values, is an important parameter 
  when controlling  the false discovery rate (FDR).  A conservative choice is
  eta0 = 1 but a choice closer to the true value will increase efficiency
  and power 
  - see   Benjamini and Hochberg (1995, 2000) and Storey (2002) for details.

The function \code{pval.estimate.eta0.w} provides five algorithms: the "conservative"
method always returns eta0 = 1  (Benjamini and Hochberg, 1995), "adaptive"
uses the approach suggested in Benjamini and Hochberg (2000), "bootstrap"
employs the method from Storey (2002), "smoother" uses the smoothing spline
approach in Storey and Tibshirani (2003), and "quantile" is a simplified version
of "smoother".
}
\value{
  The estimated proportion eta0 of null p-values. 
}

\references{
 
 \emph{"conservative" procedure:} Benjamini, Y., and Y. Hochberg (1995)  Controlling the false
 discovery rate: a practical and powerful approach to multiple testing.
 \emph{J. Roy. Statist. Soc. B}, \bold{57}, 289--300.

 \emph{"adaptive" procedure:} Benjamini, Y., and Y. Hochberg (2000) The adaptive control
 of the false discovery rate in multiple hypotheses testing with independent statistics.
 \emph{J. Behav. Educ. Statist.}, \bold{25}, 60--83.
 
 \emph{"bootstrap" procedure:} Storey, J. D. (2002) A direct approach to false
 discovery rates.
 \emph{J. Roy. Statist. Soc. B.}, \bold{64}, 479--498.
 
 \emph{"smoother" procedure:} Storey, J. D., and R. Tibshirani (2003)
  Statistical significance for genome-wide experiments. 
 \emph{Proc. Nat. Acad. Sci. USA}, \bold{100}, 9440-9445.

 \emph{"quantile" procedure:} similar to smoother, except that the lower q
 quantile of all eta0 computed in dependence of lambda is taken as overall estimate of eta0.

}
\author{
  Lingfei Wang <Lingfei.Wang.github@outlook.com>
  
  Adapted from code by K. Strimmer, Y. Benjamini, and J.D. Storey. 
}
\seealso{
\code{\link{fdrtool.w}}

\code{pval.estimate.eta0} in \code{fdrtool}}

\examples{
# load library and p-values
library("fdrtoolw")
data(pvalues,package="fdrtool")

w=runif(length(pvalues),min=0.5,max=1)
# Proportion of null p-values for different methods
pval.estimate.eta0.w(pvalues, w, method="smoother")
}
\keyword{htest}
