\name{nbinomTest}
\Rdversion{1.1}
\alias{nbinomTest}
\title{
   Perform row-wise tests for differences between the base means of two count matrices.
}
\description{
   This function is called by \code{\link{nbinomTestForContrast}}. Call it directly only
   if the S4 interface is unsuitable for your task.
}
\usage{
nbinomTest(countsA, countsB, sizeFactorsA, sizeFactorsB, rawVarFuncA, rawVarFuncB, eps = 0.003)
}
\arguments{
  \item{countsA}{ 
     A matrix of counts, where each column is a replicate
}
  \item{countsB}{
     Another matrix of counts, where each column is a replicate}
  \item{sizeFactorsA}{
     Size factors for the columns of the matrix 'countsA'
}
  \item{sizeFactorsB}{
     Size factors for the columns of the matrix 'countsB'
}
  \item{rawVarFuncA}{
     Raw variance function for 'countsA'
}
  \item{rawVarFuncB}{
     Raw variance function for 'countsB'
}
  \item{eps}{
     Precision goal for the p value. This is only a rough guidance with no guarantee
     of adherence.
}
}
\details{
   See the paper for an exact description of the null hypothesis tested.
}
\value{
   A vector of unadjusted p values, one for each row in the counts matrices.
}

\author{
   Simon Anders, EMBL Heidelberg, sanders@fs.tum.de
}
\examples{
# See the vignette for an example.
}
