\name{getVarianceStabilizedData}
\Rdversion{1.1}
\alias{getVarianceStabilizedData}
\title{
   Perform a variance stabilizing transformation (VST) on the count data
}
\description{
   This function calculates a variance stabilizing transformations (VST) from the
   raw variance functions and then transforms the count data, yielding a matrix
   of values which are now approximately homoskledastic. This is useful as input 
   to statistical analyses requiring homoskledasticity. 
}
\usage{
getVarianceStabilizedData(cds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cds}{
      a CountDataSet with estimated variance functions
}
}
\details{
   For each sample (i.e., column of \code{counts(cds)}), the full varaiance function
   is calculated from the raw variance (by scaling according to the size factor and adding 
   the shot noise). The reciprocal of the square root of the full variance is then numerically
   integrated up, and the integral (approximated by a spline function) evaluated for each
   count value in the column, yielding a tranformed value.
}
\value{
   A matrix of the same dimension as the count data, containing the transformed data.
}
\author{
   Simon Anders, EMBL Heidelberg, sanders@fs.tum.de
}

\examples{
# See vignette
}
