\name{varianceFitDiagnosticsForMatrix}
\Rdversion{1.1}
\alias{varianceFitDiagnosticsForMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   A lower-level function to get fit diagnostics
}
\description{
Usually, you call \code{\link{varianceFitDiagnostics}}, which in turn counts 
\code{\link{varianceFitDiagnosticsForMatrix}}. If you do not have your data in a
CountDataSet object, you may want to call this function directly.
}
\usage{
varianceFitDiagnosticsForMatrix(counts, sizeFactors, rawVarFunc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
     a matrix of counts, with columns corresponding to replicates
}
  \item{sizeFactors}{
     the size factors for the columns of the counts matrix
}
  \item{rawVarFunc}{
     the raw variance function estimated from the counts data
}
}
\details{
See \code{\link{varianceFitDiagnostics}} for more explanation
}
\value{
The same return value as described for \code{\link{varianceFitDiagnostics}}.
}

\seealso{
\code{\link{varianceFitDiagnostics}},\code{\link{residualsEcdfPlotFromDiagnostics}}
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
colsA <- conditions(cds) == "A"
vfd <- varianceFitDiagnosticsForMatrix( 
          counts(cds)[,colsA], 
          sizeFactors(cds)[colsA],
          rawVarFunc( cds, "A" ) )
head( vfd )

}
