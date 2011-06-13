\name{estimateDispersions}
\Rdversion{1.1}
\alias{estimateDispersions}
\title{
    Estimate and fit dispersions for a CountDataSet.
}
\description{
    This function estimates the dispersions for a count data set. For each condition
    (or collectively for all conditions, see 'method' argument below) is first estimates
    for each gene a dispersion value (a.k.a. a raw SCV value), then fits by regression
    a dispersion-mean variance and finally chooses for each gene a dispersion value to
    be used in tests from the initial and the fitted value according to the 'sharingMode'
    argument.
}
\usage{
estimateDispersions( cds, method = c( "normal", "blind", "pooled" ), 
   sharingMode = c( "fit-only", "maximum", "gene-est-only" ),
   locfit_extra_args=list(), lp_extra_args=list(), modelFrame = NULL )}
\arguments{
  \item{cds}{
      a CountDataSet with size factors
  }
  \item{method}{
     There are three ways how the variance functions can be estimated:
     
     \itemize{     
     \item normal -- For each condition with replicates, estimate the genes' dispersions
        by considering the data from samples for this condition. Then, construct
        a variance function '_max' that takes the maximum over all other variance
        functions and assign this one to all samples of unreplicated conditions.
        If the cds object has a multivariate design (i.e., a data frame was passed instead
        of a factor for the \code{condition} argument in \code{\link{newCountDataSet}}),
        this method is not available. 
     
     \item blind -- Ignore the sample labels and pretend that all samples are 
        replicates of a single condition. This allows to get a variance estimate
        even if one does not have any biological replicates. However, this typically
        leads to drastic loss of power; see the vignette for details. The single
        estimated variance condition is the called "blind" and assigned to all
        samples.
     
     \item pooled -- Use the samples from all conditions with replicates to estimate
        a single pooled variance function, to be called "pooled", and assign it to all
        samples. By "pooled variance", we mean a variance estimate where from each
        normalized count value, the mean for the value's condition is subtracted, for
        the calculation of the sum of squares (as it is typically done, e.g., in a standard
        t test). 
     }
  }  

  \item{sharingMode}{
     After dispersion values have been estimated for individual genes, a dispersion-mean
     relation is fitted for sharing information across genes in order to reduce variability
     of the dispersion estimates. Hence, for each gene, we have two values: the per-gene
     estimate (derived only from this gene's data), and the fitted value (i.e., the dispersion
     value typical for genes with an average expression similar to this gene's expression
     value). The \code{sharingMode} argument specifies which of these two values will be 
     written to the \code{featureData}'s \code{disp_} columns and hence be used 
     by the functions \code{\link{nbinomTest}} and \code{\link{fitNbinomGLMs}}.
     
     \itemize{     
     \item fit-only -- use only the fitted value, i.e., the per-gene value is used only as
         input to the fit. Use this only with very \emph{few} replicates
     \item maximum -- take the maximum of the two values. This is the conservative and hence
         prudent choice, recommended once you have at least three or four replicates and maybe even
         with only two replicates.
     \item gene-est-only -- no fitting or sharing, use only the per-gene value. This will lead
         to loss of type-I error control, unless you have many replicates, in which case it
         improves power over the 'maximum' mode, because information sharing is then unnecessary.
     }
  }

  \item{locfit_extra_args, lp_extra_args}{
      Options to be passed to the \code{locfit} and to the \code{lp} function of the 
      locfit package. Use this  to adjust the local fitting. For example, you may pass 
      a value for \code{nn} different from the default (0.7) if the fit seems too smooth 
      or too rough by setting \code{lp_extra_agrs=list(nn=0.9)} or you can set
      \code{locfit_extra_args=list(maxk=200)} if you get the error that locfit ran
      out of nodes. See the documentation of the locfit package for details. Usually, 
      you will not need to adjuste the fitting parameters, as the defaults seem to 
      work quite fine. }
      
   \item{modelFrame}{
      By default, the information in conditions(cds) or pData(cds) are used to 
      determine which samples are replicates (see \code{\link{newCountDataSet}}).
      Alternatively (and only for method "pooled"), a data frame can be passed 
      here, and all rows that are identical in this data frame are considered
      replicated.
   }
      
}
\details{
   Behaviour for \code{method} \code{normal}: For each replicated condition, a list, named
   with the condition's name, is placed in the slot \code{cds@fitInfo}. This list
   has three elements: The vector \code{perGeneDispEsts} contains the dispersion
   estimated for each gene that were used as input to the fit. The function \code{dispFunc}
   is the fit, i.e., it takes a normalized mean expression and returns the corresponding
   fitted dispersion. The values fitted according to this dispersion are in the third element,
   a vector named \code{fittedDispEst}.

   Behaviour for \code{method}s \code{blind} and \code{pooled}: Only one fitInfo object is produced,
   named \code{"blind"} or \code{"pooled"} and placed in \code{cds@fitInfo}.
   
   For each \code{fitInfo} object, the dispersion values to be used are derived according to
   the \code{sharingMode} and placed in the \code{featureData} data frame, in a column named
   with the name of the fitInfo object (i.e., the condition name or "blind" or "pooled"),
   prefixed with "\code{disp_}".
   
   Then, the \code{\link{dispTable}} (q.v.) is filled to assign to each
   condition the appropriate dispersion column in the phenoData frame.
   }
\value{
   The CountDataSet cds, with the slots rawVarFuncs and rawVarFuncTable updated.
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
dispTable( cds )
head( fData( cds ) )
}
