\name{getRawScvForSamplePair}
\alias{getRawScvForSamplePair}
\title{
for a pair of samples, calculate the average raw SCV between them, as a measure of noise
}
\description{
The average raw SCV between two samples is a good measure of the
sample-to-sample variability between two samples that does not include the shot
noise, i.e., only reflects the ``distance'' between the samples without
being affected by sequencing depth.
}
\usage{
getRawScvForSamplePair(cds, sample1, sample2)
}
\arguments{
  \item{cds}{A CountDataSet with estimated size factors.}
  \item{sample1}{ An index or column name identifying a sample. }
  \item{sample2}{ Ditto. }
}
\details{
This function calls \code{\link{estimateVarianceFunctionForMatrix}} with the 
data for just the two specified columns and so obtains a raw variance function,
which is the used to calculate the raw variance (i.e., the variance without
shot noise) for each gene. This is converted by dividing by the square of the
size-factors-normalized count for the gene, to get a raw squared coefficient of
variance (raw SCV), which is finally adjusted for the bias due to the division
with \code{\link{adjustScvForBias}}. The mean of the rawSCV values for all genes
is returned.

\emph{Interpretation}: The square root of the mean rawSCV is a measure of the
typical variation of the concentration of a gene's transcript between the
two samples. So, if you get, say, a mean raw SCV of 0.1 between two samples, it
means that, typicaly, the concentration of a genes transcript varies by 31%
(because sqrt(.1) is .31) from one sample to the other.

Check whether the value obtained between biological replicates seems accceptable
in the light of this interpretation. For technical replicates, the raw SCV value 
should be very low, maybe below 0.01, and large values indicate problems with
the library preparation.}

\value{
A single number, the mean raw SCV, as described above.
}
\author{
Simon Anders
}

\seealso{
\code{\link{getRawScvDistanceMatrix}}
}
\examples{
#cds <- makeExampleCountDataSet()
#cds <- estimateSizeFactors( cds )
#getRawScvForSamplePair( cds, "A1", "A2" )
}
