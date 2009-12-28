\name{getBaseMeansAndVariances}
\Rdversion{1.1}
\alias{getBaseMeansAndVariances}
\title{
    Perform row-wise estimates of base-level means and variances for count data.
}
\description{
    This fucntion is called internally by a number of other functions. You will need
    to call it directly only in very special cases.
}
\usage{
getBaseMeansAndVariances(counts, sizeFactors)
}
\arguments{
  \item{counts}{
     a matrix of data frame of count data. All the columns of this matrix will be
     considered as replicates of the same condition.
}
  \item{sizeFactors}{
     the size factors of the columns, as estimated e.g. with \code{\link{estimateSizeFactorsForCounts}}
}
}
\value{
A data frame with one row for each row in 'counts' and two columns:
  \item{baseMean }{The base mean for each row. This is the mean of the counts after they have been divided
     by the size factors}
  \item{comp2 }{The base variance for each row. This is the variance of the counts after they have been divided
     by the size factors}
}

\author{
   Simon Anders, EMBL Heidelberg, sanders@fs.tum.de
}
