cat( "Note to DESeq users:\n\n",
   "You have just loaded the development version of DESeq, which\n",
   "is currently undergoing a major overhaul and hence will likely\n",
   "not work as expected.\n\n",
   "Until further notice, please only use the released version of DESeq\n",
   "(i.e., not this one, but version 1.4.x from Bioconductor release 2.8.)\n" )

setClass( "CountDataSet", 
   contains = "eSet",
   representation = representation( 
      fitInfo = "environment",
      dispTable = "character",
      multivariateConditions = "logical" ),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), CountDataSet = "1.1.0" ) ) )
)

newCountDataSet <- function( countData, conditions, sizeFactors=NULL,
      phenoData = NULL, featureData = NULL )
{
   countData <- as.matrix( countData )
   if( any( round( countData ) != countData ) )
      stop( "The countData is not integer." )
   mode( countData ) <- "integer"

   if( is.null( sizeFactors ) ) {
      sizeFactors <- rep( NA_real_, ncol(countData) )
   } else
      warning( "The 'sizeFactor' argument is deprecated. Use 'estimateSizeFactors'." )
   if( is.null( phenoData ) )
      phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
   if( is.null( featureData ) ) 
      featureData <- annotatedDataFrameFrom( countData, byrow=TRUE )
      
   phenoData$`sizeFactor` <- sizeFactors
   varMetadata( phenoData )[ "sizeFactor", "labelDescription" ] <-
      "size factor (relative estimate of sequencing depth)"
   
   if( is( conditions, "matrix" ) )
      conditions <- as.data.frame( conditions )
   
   if( is( conditions, "data.frame" ) || is( conditions, "AnnotatedDataFrame" ) ) {
      stopifnot( nrow( conditions ) == ncol( countData ) )
      conditions <- as( conditions, "AnnotatedDataFrame" )
      dimLabels( conditions ) <- dimLabels( phenoData )
      rownames( pData(conditions) ) <- rownames( pData(phenoData) )
         # TODO: What if the rownames were set?
      phenoData <- combine( phenoData, conditions )
      multivariateConditions <- TRUE
      rvft <- c( `_all` = NA_character_ )
   } else {
      conditions <- factor( conditions )
      stopifnot( length( conditions ) == ncol( countData ) )
      phenoData$`condition` <- factor( conditions )
      varMetadata( phenoData )[ "condition", "labelDescription" ] <-
         "experimental condition, treatment or phenotype"
      multivariateConditions <- FALSE
      rvft <- rep( NA_character_, length(levels(conditions)) )
   }
   
   cds <- new( "CountDataSet",
      assayData = assayDataNew( "environment", counts=countData ),
      phenoData = phenoData, 
      featureData = featureData,
      multivariateConditions = multivariateConditions,
      fitInfo = new.env( hash=TRUE ),
      dispTable = rvft )
            
   cds
}

setValidity( "CountDataSet", function( object ) {
   if( length( object@multivariateConditions ) != 1 )
      return( "multivariateConditions is not scalar." )
   if( ! "sizeFactor"  %in% names(pData(object)) )
      return( "phenoData does not contain a 'sizeFactor' columns.")
   if( ! is( pData(object)$`sizeFactor`, "numeric" ) )
      return( "The 'sizeFactor' column in phenoData is not numeric." )
   if( ! object@multivariateConditions ) {
      if( ! "condition"  %in% names(pData(object)) )
         return( "phenoData does not contain a 'condition' columns." )
      if( ! is( pData(object)$`condition`, "factor" ) )
         return( "The 'condition' column in phenoData is not a factor." )
   }
   if( !is.integer( counts(object) ) )
      return( "the count data is not in integer mode" )
   if( any( counts(object) < 0 ) )
      return( "the count data contains negative values" )
   TRUE
} )

setMethod("counts", signature(cds="CountDataSet"),
  function(cds, normalized=FALSE ) {
    if(!normalized) {
      assayData(cds)[["counts"]]
    } else {
      t( t( assayData(cds)[["counts"]] ) / sizeFactors(cds) )
    }
  })   
   
setMethod("sizeFactors", signature(cds="CountDataSet"),
  function(cds) {
   sf <- pData(cds)$sizeFactor
   names( sf ) <- colnames( counts(cds) )
   sf
 })   

setMethod("sizeFactors<-", signature(cds="CountDataSet"),
  function( cds, value ) {
   pData(cds)$sizeFactor <- value
   validObject( cds )
   cds
})   

setMethod("conditions", signature(cds="CountDataSet"),
  function( cds ) {
   if( cds@multivariateConditions )
      stop( "The 'conditions' accessor is only for simple single-factor conditions, but you have specified multivariate conditions. Access them via 'pData'." )
   conds <- pData(cds)$`condition`
   names( conds ) <- colnames( counts(cds) )
   conds
})   
   
setMethod("conditions<-", signature(cds="CountDataSet"),
  function( cds, value ) {
   if( cds@multivariateConditions )
      stop( "The 'conditions' accessor is only for simple single-factor conditions, but you have specified multivariate conditions. Access them via 'pData'." )
   pData(cds)$`condition` <- factor( value )
   validObject( cds )
   cds
})

setMethod("dispTable", signature(cds="CountDataSet"),
 function( cds ) {
    cds@dispTable
})   

setMethod("dispTable<-", signature(cds="CountDataSet"),
 function( cds, value ) {
    cds@dispTable <- value
    validObject( cds )
    cds
})   

