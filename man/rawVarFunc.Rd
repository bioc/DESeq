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
rawVarFunc(cds, condOrName)
}
\arguments{
  \item{cds}{
     a CountDataSet with variance fucntions (i.e., \code{\link{estimateVarianceFunctions}} has already
     been called)
}
  \item{condOrName}{
     The name of a condition or of a raw variance function. The name is first interpreted as a name 
     in the rawVarFuncs environment. If it is not present there, it is interpreted as the name of a  
     condition and translated into a name of a raw variance function by looking up in in the 
     \code{\link{rawVarFuncTable}}. See \code{\link{estimateVarianceFunctions}} for details.
}
}
\value{
   a raw variance function. This is a function that, given a base mean (i.e., the mean of counts 
   divided by size factors), returns a raw variance estimate. To get the full variance, the raw variance
   needs to be scales up by multiplying it with the square of the size factor, and the shot noise
   has to be added.
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
