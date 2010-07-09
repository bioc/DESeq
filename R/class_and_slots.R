setClass( "CountDataSet", 
   contains = "eSet",
   representation = representation( 
      rawVarFuncs = "environment",
      rawVarFuncTable = "data.frame",
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

   if( is.null( sizeFactors ) )
      sizeFactors <- rep( NA_real_, ncol(countData) )
   if( is.null( phenoData ) )
      phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
   if( is.null( featureData ) ) 
      featureData <- annotatedDataFrameFrom( countData, byrow=TRUE )
      
   phenoData$`sizeFactor` <- sizeFactors
   varMetadata( phenoData )[ "sizeFactor", "labelDescription" ] <-
      "size factor (relative estimate of sequencing depth)"
   
   if( is( conditions, "data.frame" ) || is( conditions, "AnnotatedDataFrame" ) ) {
      stopifnot( nrow( conditions ) == ncol( countData ) )
      conditions <- as( conditions, "AnnotatedDataFrame" )
      dimLabels( conditions ) <- dimLabels( phenoData )
      rownames( pData(conditions) ) <- rownames( pData(phenoData) )
         # TODO: What if the rownames were set?
      phenoData <- combine( phenoData, conditions )
      rvft <- data.frame( 
         row.names = "_pooled",
         funcName = NA_character_, 
         varAdjFactor = NA_real_,
         stringsAsFactors = FALSE )
      multivariateConditions <- TRUE
   } else {
      conditions <- factor( conditions )
      stopifnot( length( conditions ) == ncol( countData ) )
      phenoData$`condition` <- factor( conditions )
      varMetadata( phenoData )[ "condition", "labelDescription" ] <-
         "experimental condition, treatment or phenotype"
      rvft <- data.frame( 
         row.names = levels(conditions),
         funcName = rep( NA_character_, length(levels(conditions)) ),
         varAdjFactor = rep( NA_real_, length(levels(conditions)) ),
         stringsAsFactors = FALSE )
      multivariateConditions <- FALSE
   }
   
   cds <- new( "CountDataSet",
      assayData = assayDataNew( "environment", counts=countData ),
      phenoData=phenoData, featureData=featureData,
      multivariateConditions = multivariateConditions,
      rawVarFuncs=new.env(hash=TRUE),
      rawVarFuncTable=rvft )
            
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
      if( nrow(object@rawVarFuncTable) != length(levels(conditions(object))) )
         return( "The rawVarFuncTable does not contain one row per condition." )
      if( any( rownames(object@rawVarFuncTable) != levels(conditions(object))) )
         return( "The rownames of the character vector 'rawVarFuncTable' are not identical to the levels of the factor 'conditions'." )
   } else {
      if( any( dim(object@rawVarFuncTable) != c( 0, 0 ) ) )
         return( "The rawVarFuncTable must be empty if mutivariate conditions are used." )
   }
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
   if( cds@multivariateConditions )
      stop( "The 'conditions' accessor is only for simple single-factor conditions, but your have specified multivariate conditions. Access them via 'pData'." )
   conds <- pData(cds)$`condition`
   names( conds ) <- colnames( counts(cds) )
   conds
}   
   
`conditions<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( cds@multivariateConditions )
      stop( "The 'conditions' accessor is only for simple single-factor conditions, but your have specified multivariate conditions. Access them via 'pData'." )
   pData(cds)$`condition` <- factor( value )
   validObject( cds )
   cds
}   

rawVarFunc <- function( cds, condOrName=NULL, byName=FALSE ) {
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )   
   if( is.null( condOrName ) ) {
      if( length(cds@rawVarFuncs) == 1 )
         return( cds@rawVarFuncs[[ ls(cds@rawVarFuncs)[[1]] ]] )
      else
         stop( "There is more than one variance function. 'condOrName' may not be omitted." )
   }   
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
   validObject( cds )   
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
   if( cds@multivariateConditions )
      stop( "CountDataSets with multivariate conditions do not support variance adjustment factors." )
   ans <- cds@rawVarFuncTable$varAdjFactor
   names(ans) <- rownames( cds@rawVarFuncTable )
   ans
}

`varAdjFactors<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( cds@multivariateConditions )
      stop( "CountDataSets with multivariate conditions do not support variance adjustment factors." )
   if( length( value ) != length( levels( conditions( cds ) ) ) )
      stop( "'varAdjFactors' has to be a vector with one value per condition." )
   if( !( is.null( names(value) ) || 
         all( names(value) == rownames( rawVarFuncTable(cds) ) ) ) )
      stop( "If 'varAdjFactors' is a named vector, the names being equal to the row names of rawVarFuncTable." )
   cds@rawVarFuncTable$varAdjFactor <- value
   validObject( cds )
   cds
}
