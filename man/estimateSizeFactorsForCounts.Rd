\name{estimateSizeFactorsForCounts}
\Rdversion{1.1}
\alias{estimateSizeFactorsForCounts}
\title{
   Low-level function to estimate size factors with robust regression.
}
\description{
   Given a matrix or data frame of count data, this function estimates the size factors 
   as follows: Each column is divided by the geometric means of the rows. The median
   of these ratios (skipping the genes with a geometric mean of zero) is used as the
   size factor for this column.
   
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
