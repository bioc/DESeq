\name{sizeFactors}
\Rdversion{1.1}
\alias{sizeFactors}
\alias{sizeFactors<-}
\title{
    Accessor functions for the sizeFactors information in a CountDataSet
}
\description{
    The sizeFactors vector is a factor that assigns to each column of the count data
    a value, the size factor, such that count values in teh columns can be brought to
    a common scale by dividing by the corresponding size factor.
}
\usage{
sizeFactors(cds)
sizeFactors(cds) <- value
}
\arguments{
  \item{cds}{
      a CountDataSet
}
  \item{value}{
      a vector of number, one size factor for each column in the count data
}
}

\seealso{
\code{\link{estimateSizeFactors}}
}
\examples{
   cds <- makeExampleCountDataSet()
   cds <- estimateSizeFactors( cds )
   sizeFactors(cds)
}
\author{
   Simon Anders, sanders@fs.tum.de
}

