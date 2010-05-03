\name{nbinomTestForMatrices}
\Rdversion{1.1}
\alias{nbinomTestForMatrices}
\title{
   Perform row-wise tests for differences between the base means of two count matrices.
}
\description{
   This function is called by \code{\link{nbinomTest}}. Call it directly only
   if the S4 interface is unsuitable for your task.
}
\usage{
nbinomTestForMatrices(countsA, countsB, sizeFactorsA, sizeFactorsB, rawScvA, rawScvB, eps = 0.0001)
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
  \item{rawScvA}{
     Raw squared coefficient of variantion (SCV) for 'countsA', a vector with one value per gene
}
  \item{rawScvB}{
     The same for 'countsB'
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
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
colsA <- conditions(cds) == "A"
colsB <- conditions(cds) == "B"
bmvA <- getBaseMeansAndVariances( counts(cds)[,colsA], sizeFactors(cds)[colsA] )
bmvB <- getBaseMeansAndVariances( counts(cds)[,colsB], sizeFactors(cds)[colsB] )
pvals <- nbinomTestForMatrices( 
   counts(cds)[,colsA],
   counts(cds)[,colsB],
   sizeFactors(cds)[colsA],
   sizeFactors(cds)[colsB],
   adjustScvForBias( 
      rawVarFunc( cds, "A" )( bmvA$baseMean ) / bmvA$baseMean^2,
      length( colsA ) ),
   adjustScvForBias( 
      rawVarFunc( cds, "B" )( bmvB$baseMean ) / bmvB$baseMean^2,
      length( colsB ) ) )
names( pvals ) <- row.names( counts(cds) )
head( pvals )
}
