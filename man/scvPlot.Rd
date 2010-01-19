\name{scvPlot}
\alias{scvPlot}
\title{
   Produces a diagnostic plot with the variance estimates, given as squared coefficients
   of variation (SCV)
}
\description{
   The produced plot shows the estimated variance functions of the given CountDataSet,
   in the form of the squared coefficient of variation (SCV), i.e., the variance divided
   by the squared mean.
   The solid lines are the raw SCV estimates, one per condition. The dashed lines are
   the full variance estimates for each sample, i.e., the vertical distance between a
   dashed line and its corresponding solid line (of the same colour) is the shot noise.
   As the x axis is scaled as base mean (size-adjusted mean), the amount of shot noise
   depends on the size factor.
   The solid black line is a density estimate of the base means. Only were a sufficient
   density of counts is present can a good estimate be expected.
}
\usage{
scvPlot(cds, xlim = NULL, ylim = NULL)
}
\arguments{
  \item{cds}{
      a CountDataSet with estimated variance functions
}
  \item{xlim, ylim}{
     the plot limits
}
}
\value{
   None; but a plot is produced.
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\note{
There is still a bug in this function, namely the colours between dashed and solid lines
fail to match correctly. Furthermore, the displayed SCV values are not bias corrected.}

\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
scvPlot( cds, ylim = c(0, 2) )
}
