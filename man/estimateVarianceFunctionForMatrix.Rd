\name{estimateVarianceFunctionForMatrix}
\Rdversion{1.1}
\alias{estimateVarianceFunctionForMatrix}
\title{
   Lower-level function to estimate a raw variance function.
}
\description{
   Usually, you should call \code{\link{estimateVarianceFunctions}} for a CountDataSet
   instead of calling this function directly. However, if you do not have your data
   in the form of a CountDataSet object, you may want to call this function directly.
}
\usage{
estimateVarianceFunctionForMatrix(counts, sizeFactors,
   locfit_extra_args = list(), lp_extra_args = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
     a matrix of data frame of count data. All the columns of this matrix will be
     considered as replicates of the same condition.
}
  \item{sizeFactors}{
     the size factors of the columns, as estimated e.g. with \code{\link{estimateSizeFactorsForMatrix}}
}
  \item{locfit_extra_args, lp_extra_args}{
      Options to be passed to the \code{locfit} and to the \code{lp} function of the 
      locfit package. Use this  to adjust the local fitting. For example, you may pass 
      a value for \code{nn} different from the default (0.7) if the fit seems too smooth 
      or too rough by setting \code{lp_extra_agrs=list(nn=0.9)} or you can set
      \code{locfit_extra_args=list(maxk=200)} if you get the error that locfit ran
      out of nodes. See the documentation of the locfit package for details. Usually, 
      you will not need to adjuste the fitting parameters, as the defaults seem to 
      work quite fine.
}
}
\value{
   a raw variance function. This is a function that, given a base mean (i.e., the mean of counts 
   divided by size factors), returns a raw variance estimate (which needs to be multiplied with 
   the square of a size factor to get the biological variance on the count scale; to get the full
   variance, add the shot noise as well). For compatibility with \code{\link{rawVarFunc}}, the 
   function comes with an attribute, named \code{varAdjFactor}, which is always 1.
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
vf <- estimateVarianceFunctionForMatrix( counts(cds), sizeFactors(cds) )
vf( head( counts(cds)[,1] / sizeFactors(cds)[1] ) )
attr( vf, "varAdjFactor" )
}
