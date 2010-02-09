setClass( "CountDataSet", 
   contains = "eSet",
   representation = representation( 
      rawVarFuncs = "environment",
      rawVarFuncTable = "character"),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), CountDataSet = "1.0.0" ) ) )
)

newCountDataSet <- function( countData, conditions, sizeFactors=NULL,
      phenoData = NULL, featureData = NULL )
{
   countData <- as.matrix( countData )
   conditions <- factor( conditions )
   stopifnot( ncol( countData ) == length( conditions ) )

   if( is.null( sizeFactors ) )
      sizeFactors <- rep( NA_real_, length(conditions) )
   if( is.null( phenoData ) )
      phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
   if( is.null( featureData ) ) 
      featureData <- annotatedDataFrameFrom( countData, byrow=TRUE )
      
   phenoData$`condition` <- conditions
   phenoData$`sizeFactor` <- sizeFactors
   
   cds <- new( "CountDataSet",
      assayData = assayDataNew( "environment", counts=countData ),
      phenoData=phenoData, featureData=featureData,
      rawVarFuncs=new.env(hash=TRUE),
      rawVarFuncTable=rep( NA_character_, length(levels(conditions)) ) )
      
   names(cds@rawVarFuncTable) <- levels(conditions)
   validObject( cds )
   cds
}

setValidity( "CountDataSet", function( object ) {
   if( ! "condition"  %in% names(pData(object)) )
      return( "phenoData does not contain a 'condition' columns." )
   if( ! is( pData(object)$`condition`, "factor" ) )
      return( "The 'condition' column in phenoData is not a factor." )
   if( ! "sizeFactor"  %in% names(pData(object)) )
      return( "phenoData does not contain a 'sizeFactor' columns.")
   if( ! is( pData(object)$`sizeFactor`, "numeric" ) )
      return( "The 'sizeFactor' column in phenoData is not numeric." )
   if( length(object@rawVarFuncTable) != length(levels(conditions(object))) )
      return( "The rawVarFuncTable does not contain one element per condition." )
   if( any( names(object@rawVarFuncTable) != levels(conditions(object))) )
      return( "The names of the character vector 'rawVarFuncTable' are not identical to the levels of the factor 'conditions'." )
   for( cond in levels(conditions(object)) ) {
      if( ! is.na( object@rawVarFuncTable[cond] ) ) {
         bvfName <- object@rawVarFuncTable[cond]         
         if( is.null( object@rawVarFuncs[[bvfName]] ) )
            return( sprintf( "Condition '%s' has been assigned the rawVarFunction '%s' which is missing.",
               cond, bvfName ) )
      }
   }
   for( vmfName in ls(object@rawVarFuncs) ) {
      if( !is( object@rawVarFuncs[[vmfName]], "function" ) )
         return( sprintf( "rawVarFuncs contains a value, called '%s', which is not a function.", vmfName ) )
      if( length( formals( object@rawVarFuncs[[vmfName]] ) ) != 2 ) 
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does have the right argument list for a raw variance function.", vmfName ) )
      if( any( names( formals( object@rawVarFuncs[[vmfName]] ) ) != c( "q", "reportSize" ) ) )
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does have the right argument list for a raw variance function.", vmfName ) )
      testres <- object@rawVarFuncs[[vmfName]]( 1.5 ) 
      if( ! is( testres, "numeric" ) )
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does not return a numeric result.", vmfName ) )
      if( length(testres) != 1 )
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does not return a single number as result.", vmfName ) )
   }
   TRUE
} )

counts <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   assayData(cds)[["counts"]]
}   
   
sizeFactors <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   sf <- pData(cds)$`sizeFactor`
   names( sf ) <- colnames( counts(cds) )
   sf
}   
   
`sizeFactors<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   pData(cds)$`sizeFactor` <- value
   validObject( cds )
   cds
}   

conditions <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   conds <- pData(cds)$`condition`
   names( conds ) <- colnames( counts(cds) )
   conds
}   
   
`conditions<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   pData(cds)$`condition` <- factor( value )
   validObject( cds )
   cds
}   

rawVarFunc <- function( cds, condOrName ) {
   stopifnot( is( cds, "CountDataSet" ) )
   res <- cds@rawVarFuncs[[ as.character(condOrName) ]]
   if( is.null(res) ) {
      res <- cds@rawVarFuncs[[ cds@rawVarFuncTable[ as.character(condOrName) ] ]]
      if( is.null(res) )
         stop( sprintf( "No base variance function found for condition or with name '%s'.", condOrName ) )
   }
   res
}

rawVarFuncTable <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   cds@rawVarFuncTable
}   

`rawVarFuncTable<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( is.null( names(value) ) )
      names( value ) <- levels( conditions(cds) )
   cds@rawVarFuncTable <- value
   cds
}   


ensureHasVarFuncs <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( length(ls(cds@rawVarFuncs)) == 0 )
      stop( "CountDataSet object does not contain any base variance functions. Call 'estimateVarianceFunctions' first." )
   TRUE
}   
