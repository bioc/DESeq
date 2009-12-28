\name{residualsEcdfPlotFromDiagnostics}
\Rdversion{1.1}
\alias{residualsEcdfPlotFromDiagnostics}
\title{
Produce a diagnostic plot to check the fit of a raw variance function.
}
\description{
This function produces the same plot as \code{\link{residualsEcdfPlot}}. While \code{\link{residualsEcdfPlot}}
takes a CountDataSet as argument and the calls \code{\link{baseVarDiagnostics}}, the present function
expects the output of a call to \code{\link{baseVarDiagnostics}} (or \code{\link{baseVarDiagnosticsForCounts}}) as
input.
}
\usage{
residualsEcdfPlotFromDiagnostics(fitdiag, ncuts = 7, plotTitle = "Residuals ECDF plot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitdiag}{
     a data frame as returned by \code{\link{baseVarDiagnostics}} or \code{\link{baseVarDiagnosticsForCounts}}
}
  \item{ncuts}{
     the number of baseMean strata
}
  \item{plotTitle}{
     the main title for the plot
}
}
\details{
   See \code{\link{residualsEcdfPlot}}.
}
\value{
None, but a plot is produced.
}
\seealso{
\code{\link{residualsEcdfPlot}}, \code{\link{baseVarDiagnostics}}, \code{\link{baseVarDiagnosticsForCounts}}
}
