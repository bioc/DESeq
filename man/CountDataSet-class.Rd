\name{CountDataSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CountDataSet-class}

\title{Class "CountDataSet" -- a container for count data from HTS experiments }
\description{	This is the main class for the present package. }
\section{Objects from the Class}{
Objects should be created with calls to \code{\link{newCountDataSet}} (q.v.).
}
\section{Extends}{
Class \code{eSet} (package 'Biobase'), directly.
Class \code{VersionedBiobase} (package 'Biobase'), by class "eSet", distance 2.
Class \code{Versioned} (package 'Biobase'), by class "eSet", distance 3.
}
\note{ 
Note: This is a summary for reference. For an explanation of the actual usage, see the vignette.

A CountDataSet object stores counts from an HTS data set and offers
further slots which are populated during the analysis.

After creation with \code{\link{newCountDataSet}}, a CountDataSet typically contains 
a count table, i.e., a matrix of integer data, that is accessible with the accessor function \code{\link{counts}}.
Each row of the matrix corresponds to a gene (or binding region, or the like), and
each colum to an experimental sample. The experimental conditions of the samples
are stored in a factor (with one element for each row of the counts matrix), which
can be read with the accessor function \code{\link{conditions}}.

In the following analysis steps, further data slots are populated. First, the size factors
can be estimated with \code{\link{estimateSizeFactors}}, which are afterwards accessible 
via \code{\link{sizeFactors}}. Then, the variance functions are estimated with \code{\link{estimateVarianceFunctions}}. The resulting estimates are accessible via the function
\code{\link{rawVarFunc}}.

Internally, the mentioned data is stored in slots as follows: 

As \code{CountDataSet} is derived from \code{eSet}, it has a \code{phenoData} slots 
which allows to store sample annotation. This is used to store the factor with the conditions,
as a data frame column named \code{condition}, and to store the size factors, as an numeric data 
frame column named \code{sizeFactor}. The user may add further columns to the \code{phenoData} AnnotatedDataFrame.

The counts table is stored in the \code{eSet}'s \code{assayData} locked environment with 
the name \code{counts}.

The fits calculated by \code{\link{estimateVarianceFunctions}} are stored as closures in
an environment, which is part of the object as slot \code{rawVarFuncs}. The assignment of the
closures' names in the environment to the condition names is stored in a named character
vector, which forms the slot \code{rawVarFuncTable}. For further details on these two slots, see
\code{\link{estimateVarianceFunctions}}.

}

\examples{
# See the vignette
}

