\name{estimateSizeFactorsForCounts}
\Rdversion{1.1}
\alias{estimateSizeFactorsForCounts}
\title{
   Low-level fucntion to estimate size factors with robust regression.
}
\description{
   Given a matrix or data frame of count data, this function estimates the size factors 
   by regressing each column on the row-wise geometric mean of the counts. The regression
   is done with the \code{rlm} function from the MASS package, with intercept forced to
   zero. The regression coefficients are the size factors.
   
   Typically, you will not call this function directly, but use \code{\link{estimateSizeFactors}}.
}
\usage{
estimateSizeFactorsForCounts(counts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
      a matrix or data frame of counts, i.e., non-negative integer values
}
}
\value{
   a vector with the estimates size factors, one element per column
}

\author{
   Simon Anders, EMBL Heidelberg, sanders@fs.tum.de
}
\seealso{
\code{\link{estimateSizeFactors}}
}
\examples{
# See the vignette for an example.
}
