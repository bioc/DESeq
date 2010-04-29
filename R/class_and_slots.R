setClass( "CountDataSet", 
   contains = "eSet",
   representation = representation( 
      rawVarFuncs = "environment",
      rawVarFuncTable = "data.frame"),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), CountDataSet = "1.0.0" ) ) )
)

newCountDataSet <- function( countData, conditions, sizeFactors=NULL,
      phenoData = NULL, featureData = NULL )
{
   countData <- as.matrix( countData )
   if( any( round( countData ) != countData ) )
      stop( "The countData is not integer." )
   mode( countData ) <- "integer"   
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
   
   rvft <- data.frame( 
      funcName = rep( NA_character_, length(levels(conditions)) ),
      varAdjFactor = rep( NA_real_, length(levels(conditions)) ),
      stringsAsFactors = FALSE )
   rownames(rvft) <- levels(conditions)
   
   cds <- new( "CountDataSet",
      assayData = assayDataNew( "environment", counts=countData ),
      phenoData=phenoData, featureData=featureData,
      rawVarFuncs=new.env(hash=TRUE),
      rawVarFuncTable=rvft )
            
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
   if( nrow(object@rawVarFuncTable) != length(levels(conditions(object))) )
      return( "The rawVarFuncTable does not contain one row per condition." )
   if( any( rownames(object@rawVarFuncTable) != levels(conditions(object))) )
      return( "The rownames of the character vector 'rawVarFuncTable' are not identical to the levels of the factor 'conditions'." )
   if( ncol(object@rawVarFuncTable) != 2 )
      return( "The rawVarFuncTable does not contain two columns." )
   if( all( colnames(object@rawVarFuncTable) != c( "funcName", "varAdjFactor" ) ) )
      return( "The rawVarFuncTable' columns are not named 'funcName' and 'varAdjFactor'." )
   if( ! is( object@rawVarFuncTable$funcName, "character" ) )
      return( "The rawVarFuncTable' column 'funcName' is not of type character." )
   if( ! is( object@rawVarFuncTable$varAdjFactor, "numeric" ) )
      return( "The rawVarFuncTable' column 'varAdjFactor' is not of type numeric." )
   for( cond in levels(conditions(object)) ) {
      if( ! is.na( object@rawVarFuncTable$funcName[cond] ) ) {
         bvfName <- object@rawVarFuncTable[cond]         
         if( is.null( object@rawVarFuncs[[bvfName]] ) )
            return( sprintf( "Condition '%s' has been assigned the rawVarFunction '%s' which is missing.",
               cond, bvfName ) )
      }
   }
   for( vmfName in ls(object@rawVarFuncs) ) {
      if( !is( object@rawVarFuncs[[vmfName]], "function" ) )
         return( sprintf( "rawVarFuncs contains a value, called '%s', which is not a function.", vmfName ) )
      if( length( formals( object@rawVarFuncs[[vmfName]] ) ) != 1 ) 
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does have the right argument list for a raw variance function.", vmfName ) )
      testres <- object@rawVarFuncs[[vmfName]]( 1.5 ) 
      if( ! is( testres, "numeric" ) )
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does not return a numeric result.", vmfName ) )
      if( length(testres) != 1 )
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does not return a single number as result.", vmfName ) )
      if( ! is( attr(testres, "size" ), "numeric" ) )
         return( sprintf( "rawVarFuncs contains a function, called '%s', which return a result with numeric 'size' attribute.", vmfName ) )
   }
   if( !is.integer( counts(object) ) )
      return( "the count data is not in integer mode" )
   if( any( counts(object) < 0 ) )
      return( "the count data contains negative values" )
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

rawVarFunc <- function( cds, condOrName, byName=FALSE ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( byName ) {
      res <- cds@rawVarFuncs[[ as.character(condOrName) ]]
      if( is.null(res) )
         stop( sprintf( "No raw variance function found with name '%s'.", condOrName ) )
      attr( res, "varAdjFactor" ) <- 1
   } else {
      res <- cds@rawVarFuncs[[ cds@rawVarFuncTable[ as.character(condOrName), "funcName" ] ]]
      if( is.null(res) )
         stop( sprintf( "No raw variance function found for condition '%s'.", condOrName ) )
      attr( res, "varAdjFactor" ) <- cds@rawVarFuncTable[ as.character(condOrName), "varAdjFactor" ]
   }
   res
}

rawVarFuncTable <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   cds@rawVarFuncTable
}   

`rawVarFuncTable<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( is.null( rownames(value) ) )
      rownames( value ) <- levels( conditions(cds) )
   cds@rawVarFuncTable <- value
   cds
}   

ensureHasVarFuncs <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( length(ls(cds@rawVarFuncs)) == 0 )
      stop( "CountDataSet object does not contain any variance functions. Call 'estimateVarianceFunctions' first." )
   TRUE
}   

varAdjFactors <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   ans <- cds@rawVarFuncTable$varAdjFactor
   names(ans) <- rownames( cds@rawVarFuncTable )
   ans
}

`varAdjFactors<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   validObject( cds )
   if( length( value ) != length( levels( conditions( cds ) ) ) )
      stop( "'varAdjFactors' has to be a vector with one value per condition." )
   if( !( is.null( names(value) ) || 
         all( names(value) == rownames( rawVarFuncTable(cds) ) ) ) )
      stop( "If 'varAdjFactors' is a named vector, the names being equal to the row names of rawVarFuncTable." )
   cds@rawVarFuncTable$varAdjFactor <- value
   cds
}
