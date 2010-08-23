\name{getRawScvDistanceMatrix}
\alias{getRawScvDistanceMatrix}
\title{
Calculates a matrix of "distances" between samples
}
\description{
In order to judge the similarity of samples, a distance
matrix can be handy. This function calls 
\code{\link{getRawScvForSamplePair}} for each pair 
of samples and returns the a symmetric distance matrix.
}
\usage{
getRawScvDistanceMatrix(cds)
}
\arguments{
  \item{cds}{ a CountDataSet with estimates size factors}
}
\details{
See \code{\link{getRawScvForSamplePair}} for details on the metric.
}
\value{
A symmetric matrix, with one row and one column for each sample.
}
\author{
Simon Anders
}
\seealso{
\code{\link{getRawScvForSamplePair}}
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
dm <- getRawScvDistanceMatrix( cds )
heatmap( dm, symm=TRUE )
}
