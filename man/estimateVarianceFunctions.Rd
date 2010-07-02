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
estimateVarianceFunctions(cds, 
   pool = NULL, method = c( "normal", "blind", "pooled" ),
   locfit_extra_args = list(), lp_extra_args = list())
}
\arguments{
  \item{cds}{
      a CountDataSet with size factors
}
  \item{method}{
     There are three ways how the variance functions can be estimated:
     
     \itemize{     
     \item{normal}{For each condition with replicates, estimate a variance function
        by considering the data from samples for this condition. Then, construct
        a variance function '_max' that takes the maximum over all other variance
        functions and assign this one to all samples of unreplicated conditions.}
     
     \item{blind}{Ignore the sample labels and pretend that all samples are 
        replicates of a single condition. This allows to get a variance estimate
        even if one does not have any biological replicates. However, this can 
        leed to drastic loss of power; see the vignette for details. The single
        estimated variance condition is the called "_blind" and assigned to all
        samples.}
     
     \item{pooled}{Use the samples from all conditions with replicates to estimate
        a single pooled variance function, to be called "_pool" and assign it to all
        samples.} }
}  
  \item{pool}{This argument is deprecated; do not use it. It is (for now) 
        retained for compatibility with code written for DESeq versions 1.1.6 or
        eralier. Setting 'pool=FALSE' is the same as 'method="normal"' and setting
        'pool=TRUE' is the same as 'method="blind"' (not 'method="pool"). Note that
        the deprecation also resolved the issue that calling the 'blind' estimation
        a 'pooled' was incorrect use of terminology.}
      
  \item{locfit_extra_args, lp_extra_args}{
      Options to be passed to the \code{locfit} and to the \code{lp} function of the 
      locfit package. Use this  to adjust the local fitting. For example, you may pass 
      a value for \code{nn} different from the default (0.7) if the fit seems too smooth 
      or too rough by setting \code{lp_extra_agrs=list(nn=0.9)} or you can set
      \code{locfit_extra_args=list(maxk=200)} if you get the error that locfit ran
      out of nodes. See the documentation of the locfit package for details. Usually, 
      you will not need to adjuste the fitting parameters, as the defaults seem to 
      work quite fine. }
      
}
\details{
   Behaviour for pooled=FALSE: The estimated raw variance functions are placed in the
   environment rawVarFuncs, which is a slot in CountDataSet, using the condition labels
   as names. A further function, named "_max", is placed there as well, which always
   return the maximum of all the other functions.
   
   Then, the \code{\link{rawVarFuncTable}} (q.v.) is filled to assign to each replicated
   condition the raw variance function estimated for it, and to each condition without
   replicates, the "_max" function. 
   
   Behaviour for pooled=TRUE: A single raw variance function is estimated from all the
   count data, ignoring the condition labels. It is stored in the rawVarFuncs slot
   under the name "_pooled". In the rawVarFuncTable, "_pooled" is assigned to all
   conditions.
   
   In either case, all the variance adjustment factors (see 
   \code{\link{varAdjFactors}}) are set to 1.
   
   It is advisable to always call \code{\link{residualsEcdfPlot}} afterwards to verify
   the fit.
}
\value{
   The CountDataSet cds, with the slots rawVarFuncs and rawVarFuncTable updated.
}
\seealso{
   \code{\link{scvPlot}} to visualize the result and \code{\link{varianceFitDiagnostics}} and 
   \code{\link{residualsEcdfPlot}} to check the fit.
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
vf <- rawVarFunc( cds, "A" )
vf( head( counts(cds)[,1] / sizeFactors(cds)[1] ) )
attr( vf, "size" )
}
