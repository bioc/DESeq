\name{estimateSizeFactors}
\Rdversion{1.1}
\alias{estimateSizeFactors}

\title{
   Estimate the size factors for a CountDataSet
}
\description{
   Estimate the size factors for a CountDataSet
}
\usage{
estimateSizeFactors(cds)
}
\arguments{
  \item{cds}{
     a CountDataSet
}
}
\details{
   You need to call this function right after \code{\link{newCountDataSet}} unless
   you have manually specified size factors.
   
   Typically, the function is called with the idiom
   
   \code{ cds <- estimateSizeFactors( cds ) }
   
   This estimates the size factors and stores the information in the object.

   Internally, the function calls \code{\link{estimateSizeFactorsForMatrix}}. See there
   for more details on the calculation.
}
\value{
   The CountDataSet passed as parameters, with the size factors filled in.
}

\author{
   Simon Anders, EMBL Heidelberg, sanders@fs.tum.de
}
\seealso{
\code{\link{estimateSizeFactorsForMatrix}}
}
\examples{
# See the vignette for an example.
}
