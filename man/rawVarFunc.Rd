\name{rawVarFunc}
\Rdversion{1.1}
\alias{rawVarFunc}
\title{
   Accessor for raw variance functions
}
\description{
   Given the name of a raw variance function or the name of a condition, get
   the function from the environment \code{rawVarFuncs} (which is a slot of
   CountDataSet).
}
\usage{
rawVarFunc(cds, condOrName, byName=FALSE)
}
\arguments{
  \item{cds}{
     a CountDataSet with variance fucntions (i.e., \code{\link{estimateVarianceFunctions}} has already
     been called)
}
  \item{condOrName}{
     For \code{byName==FALSE}, this argument is taken to be the name of a condition. \code{rawVarFunc}
     then looks up in \code{rawVarFuncTable(cds)} which raw variance function should be used for
     the given condition and returns this function.  For \code{byName==TRUE}, the look-up is skipped, i.e.,
     the argument is taken directly as the name under which the raw variance function is stored in
     the \code{cds}'s \code{rawVarFuncs} environment.
}
   \item{byName}{
     how to interpret the \code{condByName} argument, see above
}   
}
\value{
   a raw variance function. This is a function that, given a base mean (i.e., the mean of counts 
   divided by size factors), returns a raw variance estimate. To get the full variance, the raw variance
   needs to be scales up by multiplying it with the square of the size factor, and the shot noise
   has to be added. The function comes with an attribute, named \code{varAdjFactor}, which contains 
   (for \code{byName==FALSE}) the  variance adjustment factor (see \code{\link{varAdjFactors}}) set for 
   this condition. For \code{byName==TRUE}, the variance adjustment factor is always given as 1. Note
   that the raw variance fucntion does not multiply its result with this factor; you have to do this
   yourself if needed.
}

\author{
   Simon Anders, sanders@fs.tum.de
}
\seealso{
\code{\link{estimateSizeFactors}}
}
\examples{
# See example for estimateVarianceFunctions
}
