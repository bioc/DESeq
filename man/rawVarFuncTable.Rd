\name{rawVarFuncTable}
\Rdversion{1.1}
\alias{rawVarFuncTable}
\alias{rawVarFuncTable<-}
\title{
Accessor functions for the slot rawVarFuncTable of the class CountDataSet
}
\description{
Accesses the rawVarFuncTable.
}
\usage{
rawVarFuncTable(cds)
rawVarFuncTable(cds) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cds}{
     a CountDataSet object
}
\item{value}{
   The rawVarFuncTable is a named character vector with one element
   for each condition, and the condition labels used as names. The value of the element
   is the name, under which the raw variance function to be used for the corresponding
   condition can be found in the environment rawVarFuncs (which is a slot of CountDataSet).

   The function \code{\link{estimateVarianceFunctions}} fills in this table. The function
   \code{\link{rawVarFunc}}, and through it several other functions working on CountDataSets,
   use it to look up which raw variance function to use for which condition.
}
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\seealso{
\code{\link{estimateSizeFactors}}
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
rawVarFuncTable( cds )
}
