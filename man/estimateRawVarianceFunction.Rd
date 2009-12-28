\name{estimateRawVarianceFunction}
\Rdversion{1.1}
\alias{estimateRawVarianceFunction}
\title{
   Lower-level function to estimate a raw variance function.
}
\description{
   Usually, you should call \code{\link{estimateVarianceFunctions}} for a CountDataSet
   instead of calling this function directly. However, if you do not have your data
   in the form of a CountDataSet object, you may want to call this function directly.
}
\usage{
estimateRawVarianceFunction(counts, sizeFactors, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
     a matrix of data frame of count data. All the columns of this matrix will be
     considered as replicates of the same condition.
}
  \item{sizeFactors}{
     the size factors of the columns, as estimated e.g. with \code{\link{estimateSizeFactorsForCounts}}
}
  \item{\dots}{
      Options to be passed to the \code{lp} function of the locfit package. Use this
      to adjust the local fitting. For example, you may pass a value for \code{nn}
      different from the default (0.7) if the fit seems to smotth or to rough. See
      th documentation of the locfit package for details. Usually, you will not
      need to adjuste the fitting, as the defaults seem to work quite fine.  
}
}
\value{
   a raw variance function. This is a function that, given a base mean (i.e., the mean of counts 
   divided by size factors), returns a raw variance estimate (which needs to be multiplied with 
   the square of a size factor to get the biological variance on the count scale; to get the full
   variance, add the shot noise as well).
}
\author{
   Simon Anders, EMBL Heidelberg, sanders@fs.tum.de
}

