\name{getVarianceStabilizedData}
\Rdversion{1.1}
\alias{getVarianceStabilizedData}
\title{
   Perform a variance stabilising transformation (VST) on the count data
}
\description{
   This function calculates a variance stabilising transformations (VST) from the
   fitted dispersion-mean reltions and then transforms the count data (after normalization
   by division by the size factor), yielding a matrix
   of values which are now approximately homoskedastic. This is useful as input 
   to statistical analyses requiring homoskedasticity. 
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
   For each sample (i.e., column of \code{counts(cds)}), the full variance function
   is calculated from the raw variance (by scaling according to the size factor and adding 
   the shot noise). The function requires a blind estimate of the variance function, i.e.,
   one ignoring conditions, and hence, you need to call \code{\link{estimateDispersions}}
   with \code{method="blind"} before caling it. 

   The reciprocal of the square root of the variance of the normalized counts, as derived
   from the dispersion fit, is then numerically
   integrated up, and the integral (approximated by a spline function) evaluated for each
   count value in the column, yielding a transformed value.
   
   Limitations: In order to preserve normalization, the same transformation has to be used
   for all samples. This results in the variance stabilizition to be only approximate. The
   more the size factors differ, the more residual dependence of the variance on the mean
   you will find in the transformed data. (Compare the variance of the upper half of your
   transformed data with the lower half to see whether this is a problem in your case.)
}
\value{
   A matrix of the same dimension as the count data, containing the transformed data.
}
\author{
   Simon Anders, sanders@fs.tum.de
}

\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cds )
colsA <- conditions(cds) == "A"
plot( rank( rowMeans( vsd[,colsA] ) ), genefilter::rowVars( vsd[,colsA] ) )
}
