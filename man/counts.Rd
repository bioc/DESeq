\name{counts}
\Rdversion{1.1}
\alias{counts}
\title{
   Accessor for the 'counts' slot of a CountDataSet object.
}
\description{
   The counts slot holds the count data as a matrix of non-negative integer
   count values, one row for each observational unit (gene or the like), and
   one column for each sample.
}
\usage{
counts(cds, normalized=FALSE)
}
\arguments{
  \item{cds}{a CountDataSet object.}
  \item{normalized}{logical indicating whether or not to divide the
    counts by the size factors before returning.}
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
   cds <- makeExampleCountDataSet()
   head( counts( cds ) )
}   
   
