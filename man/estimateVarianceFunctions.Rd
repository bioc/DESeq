\name{estimateVarianceFunctions}
\Rdversion{1.1}
\alias{estimateVarianceFunctions}
\title{
    Estimate the variance functions for a CountDataSet.
}
\description{
    This function calls, for each condition that has replicates, the 
    lower-level function \code{\link{estimateVarianceFunctionForMatrix}} 
    to estimate the raw variance function for this condition.     
}
\usage{
estimateVarianceFunctions(cds, pool = FALSE, ...)
}
\arguments{
  \item{cds}{
      a CountDataSet with size factors
}
  \item{pool}{
      Ignore the conditions and pool all data to estimate a single raw variance
      function. Use this option if you have no replicates but be sure to read the
      caveats in the vignette.
}
  \item{\dots}{
      Options to be passed to the \code{lp} function of the locfit package. Use this
      to adjust the local fitting. For example, you may pass a value for \code{nn}
      different from the default (0.7) if the fit seems to smotth or to rough. See
      th documentation of the locfit package for details. Usually, you will not
      need to adjuste the fitting, as the defaults seem to work quite fine.
}
}
\details{
   Behavior for pooled=FALSE: The estimated raw variance functions are placed in the
   environment rawVarFuncs, which is a slot in CountDataSet, using the condition labels
   as names. A further function, named "_max", is placed there as well, which always
   return the maximum of all the other functions.
   
   Then, the \code{\link{rawVarFuncTable}} (q.v.) is filled to assign to each replicated
   condition the raw variance fucntion estimated for it, and to each condition without
   replicates, the "_max" function.
   
   Behavior for pooled=TRUE: A single raw variance function is estimated from all the
   count data, ignoring the condition labels. It is stored in the rawVarFuncs slot
   under the name "_pooled". In the rawVarFuncTable, "_pooled" is assigned to all
   conditions.
}
\value{
   The CountDataSet cds, with the slots rawVarFuncs and rawVarFuncTable updated.
}
\seealso{
   \code{\link{scvPlot}} to visualize the result and \code{\link{varianceFitDiagnostics}} to check the fit
}
\author{
   Simon Anders, EMBL Heidelberg, sanders@fs.tum.de
}
\examples{
# See the vignette for an example.
}
