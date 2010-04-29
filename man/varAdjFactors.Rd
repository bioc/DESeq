\name{varAdjFactors}
\Rdversion{1.1}
\alias{varAdjFactors}
\alias{varAdjFactors<-}
\title{
    Accessor functions for the varAdjFactors information in a CountDataSet
}
\description{
    In case that the estimated variance functions fail to fit well (see 
    \code{\link{residualsEcdfPlot}}, one may want to adjust the variance
    estimates. For each condition, a variance adjustment factor can be given,
    with which the raw variance estimates will be mutiplied.
    
    These accessor functions allow to get and set the vectors of variance
    adjustment factors.

}
\usage{
varAdjFactors(cds)
varAdjFactors(cds) <- value
}
\arguments{
  \item{cds}{
      a CountDataSet
}
  \item{value}{
      a vector of numbers, one variance adjustment factor for each condition
}
}

\seealso{
\code{\link{residualsEcdfPlot}}
}
\examples{
   cds <- makeExampleCountDataSet()
   cds <- estimateSizeFactors( cds )
   cds <- estimateVarianceFunctions( cds )
   varAdjFactors(cds)
   residualsEcdfPlot(cds, "A" )
   varAdjFactors(cds) <- c( A=.7, B=1 )
   residualsEcdfPlot(cds, "A" )
}
\author{
   Simon Anders, sanders@fs.tum.de
}

