\name{conditions}
\Rdversion{1.1}
\alias{conditions}
\alias{conditions,CountDataSet-method}
\title{
    Accessor function for the conditions information in a CountDataSet
}
\description{
    The conditions vector is a factor that assigns to each column of the count data
    a condition (or treatment, or phenotype, or the like). This information is stored
    in the CountDataSet's "phenoData" slot as a row named "condition".
}
\usage{
conditions(cds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cds}{
      a CountDataSet
}
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
   cds <- makeExampleCountDataSet()
   conditions( cds )
}   
   
