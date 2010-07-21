\name{residualsEcdfPlot}
\Rdversion{1.1}
\alias{residualsEcdfPlot}
\title{
Produce a diagnostic plot to check the fit of a raw variance function.
}
\description{
The residuals of the fit should follow a scaled chi-squared distribution. This function
calls \code{\link{varianceFitDiagnostics}} to get the cumulative chi-squared 
probabilities of the residuals and plots their ECDFs, stratified by base means.
}
\usage{
residualsEcdfPlot(cds, condition, ncuts = 7)
}
\arguments{
  \item{cds}{
     a CountDataSet with raw variance functions     
}
  \item{condition}{
     the name of a condition
}
  \item{ncuts}{
     the number of base mean strata (i.e, of curves)
}
}
\details{
   As the cumulative chi-square probabilities should be uniform, the ECDF curves
   should roughly follow the main diagonal (indicated by a green line). It is
   acceptable if the strata for very low counts deviate from this.
   
   If the ECDF curves deviate too strongly from the green line, your replicates are
   most likly heterogeneous (i.e., their some are more similar than others), and the
   variance estimation procedure gets mislead.
   If the ECDF curves are below the green line, variance is underestimated, which leads
   to too low p values (and hence, too many false positives). If the ECDF curves are above the green
   line, variance is overestimated, which leads to too high p values (and hence, DESeq might
   be more conservative in its test when necessary). 
   You can adjust this by setting variance adjustment factors with \code{\link{varAdjFactors}}.
}
\value{
None, but a plot is produced.
}
\seealso{
\code{\link{residualsEcdfPlotFromDiagnostics}}, \code{\link{varianceFitDiagnostics}}
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
residualsEcdfPlot( cds, "A" )
}
