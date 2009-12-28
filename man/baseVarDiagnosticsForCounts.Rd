\name{baseVarDiagnosticsForCounts}
\Rdversion{1.1}
\alias{baseVarDiagnosticsForCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   A lower-level function to get fit diagnostics
}
\description{
Usually, you call \code{\link{baseVarDiagnostics}}, which in turn counts 
\code{\link{baseVarDiagnosticsForCounts}}. If you do not have your data in a
CountDataSet object, you may want to call this function directly.
}
\usage{
baseVarDiagnosticsForCounts(counts, sizeFactors, rawVarFunc)
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
See \code{\link{baseVarDiagnostics}} for more explanation
}
\value{
The same return value as described for \code{\link{baseVarDiagnostics}}.
}

\seealso{
\code{\link{baseVarDiagnostics}},\code{\link{residualsEcdfPlotFromDiagnostics}}
}
\author{
   Simon Anders, EMBL Heidelberg, sanders@fs.tum.de
}
\examples{
# See the vignette for an example.
}
