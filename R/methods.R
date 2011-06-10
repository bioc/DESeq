estimateSizeFactors <- function( cds, locfunc=median )
{
   stopifnot( is( cds, "CountDataSet" ) )
   sizeFactors(cds) <- estimateSizeFactorsForMatrix( counts(cds), locfunc )
   cds
}

estimateDispersions <- function( cds, method = c( "normal", "blind", "pooled" ), 
   sharingMode = c( "fit-only", "maximum", "gene-est-only" ),
   locfit_extra_args=list(), lp_extra_args=list(), modelFrame = NULL )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   if( any( is.na( sizeFactors(cds) ) ) )
      stop( "NAs found in size factors. Have you called already 'estimateSizeFactors'?" )
   method <- match.arg( method )
   sharingMode <- match.arg( sharingMode )
   if( cds@multivariateConditions && ! method %in% c( "blind", "pooled" ) )
      stop( "You have specified multivariate conditions (i.e., passed a data frame with conditions). In this case, you need to specify method 'pooled' or 'blind'." )
   if( sharingMode == "gene-est-only" )
      warning( "sharingMode=='gene-est-only' will cause inflated numbers of false positives unless you have many replicates." )
   ## FIXME this warning should only be emitted when the number of replicates is indeed small.
   
   cds@fitInfo = new.env( hash=TRUE )
   
   switch( method,
     "blind" = {

      bmv <- getBaseMeansAndVariances( counts(cds), sizeFactors(cds) )
      dispFunc <- estimateDispersionFunctionFromBaseMeansAndVariances( bmv$baseMean,
         bmv$baseVar, sizeFactors(cds), locfit_extra_args, lp_extra_args )   
      cds@fitInfo[[ "blind" ]] <- list( 
         perGeneDispEsts = adjustScvForBias( 
            pmax( 0, ( bmv$baseVar - bmv$baseMean * mean(sizeFactors(cds)) ) / bmv$baseMean^2 ), ncol(counts(cds)) ),
         dispFunc = dispFunc,
         fittedDispEsts = dispFunc( bmv$baseMean ) )
      
      if( cds@multivariateConditions )
         dispTable(cds) <- c( "_all" = "blind" )
      else {
         a <- rep( "blind", length( levels( conditions(cds) ) ) )
         names(a) <- levels( conditions(cds) )
         cds@dispTable <- a }
      
   },
   "normal" = {
   
      replicated <- names( which( tapply( conditions(cds), conditions(cds), length ) > 1 ) )
      if( length( replicated ) < 1 )
         stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions." )
      nonreplicated <- names( which( tapply( conditions(cds), conditions(cds), length ) == 1 ) )
      overall_basemeans <- rowMeans( counts( cds, normalized=TRUE ) )

      for( cond in replicated ) {
         cols <- conditions(cds)==cond
         bmv <- getBaseMeansAndVariances( counts(cds)[ , cols ], sizeFactors(cds)[ cols ] )
         dispFunc <- estimateDispersionFunctionFromBaseMeansAndVariances( bmv$baseMean,
            bmv$baseVar, sizeFactors(cds)[cols], locfit_extra_args, lp_extra_args )   
         cds@fitInfo[[ cond ]] <- list( 
            perGeneDispEsts = adjustScvForBias( 
               pmax( 0, ( bmv$baseVar - bmv$baseMean * mean(sizeFactors(cds)[cols]) ) / bmv$baseMean^2 ), sum(cols) ),
            dispFunc = dispFunc,
            fittedDispEsts = dispFunc( overall_basemeans ) ) }   # Note that we do not use bmv$baseMean here
         
      cds@dispTable <- sapply( levels(conditions(cds)), function( cond )
            ifelse( cond %in% replicated, cond, "max" ) ) 
                        
   },
   "pooled"  = {
   
      if( cds@multivariateConditions ) {
         if( is.null( modelFrame ) )
            modelFrame <- pData(cds)[ , colnames(pData(cds)) != "sizeFactor" ]
         conds <- modelMatrixToConditionFactor( modelFrame ) }
      else
         conds <- conditions(cds)
   
      bmv <- getBaseMeansAndPooledVariances( counts(cds), sizeFactors(cds), conds )
      dispFunc <- estimateDispersionFunctionFromBaseMeansAndVariances( bmv$baseMean,
         bmv$baseVar, sizeFactors(cds), locfit_extra_args, lp_extra_args )   
      cds@fitInfo[[ "pooled" ]] <- list( 
         perGeneDispEsts = adjustScvForBias( 
            pmax( 0, ( bmv$baseVar - bmv$baseMean * mean(sizeFactors(cds)) ) / bmv$baseMean^2 ), length(sizeFactors(cds)) ),
         dispFunc = dispFunc,
         fittedDispEsts = dispFunc( bmv$baseMean ) )

      if( cds@multivariateConditions )
         dispTable(cds) <- c( "_all" = "pooled" )
      else {
         a <- rep( "pooled", length( levels( conditions(cds) ) ) )
         names(a) <- levels( conditions(cds) )
         cds@dispTable <- a }
   },
   stop(sprintf("Invalid method '%s'.", method))
   ) ## switch
   
   for( n in ls(cds@fitInfo) )
      fData(cds)[[ paste( "disp", n, sep="_" ) ]] <- 
         switch( sharingMode, 
            `fit-only`      = cds@fitInfo[[ n ]]$fittedDispEsts,
            `gene-est-only` = cds@fitInfo[[ n ]]$perGeneDispEsts,
            `maximum`       = pmax( cds@fitInfo[[ n ]]$fittedDispEsts, cds@fitInfo[[ n ]]$perGeneDispEsts, na.rm=TRUE ),
            stop(sprintf("Invalid sharingMode '%s'.", sharingMode))
         ) ## switch
        
   if( "max" %in% cds@dispTable )
      fData(cds)[["disp_max"]] <- do.call( pmax, 
         c( fData(cds)[ , colnames(fData(cds)) %in% paste( "disp", cds@dispTable, sep="_" ), drop=FALSE ], na.rm=TRUE ) )
        
        
   validObject( cds )
   cds
}


#varianceFitDiagnostics <- function( cds, cond=NULL, modelFrame=NULL )
#{
#   stopifnot( is( cds, "CountDataSet" ) )
#   ensureHasVarFuncs( cds )
#   
#   if( is.null( cond ) ) {
#      if( length(cds@rawVarFuncs) != 1 )
#         stop( "There is more than one variance function. 'cond' may not be omitted." )
#      cond <- ls( cds@rawVarFuncs )[[1]]
#      rvf <- rawVarFunc( cds, cond, byName=TRUE ) 
#   } 
#   else
#      rvf <- rawVarFunc( cds, cond ) 
#   
#   if( cond == "_pooled" ) {
#      if( cds@multivariateConditions ) {
#         if( is.null( modelFrame ) )
#            modelFrame <- pData(cds)[ , colnames(pData(cds)) != "sizeFactor" ]
#         poolconds <- modelMatrixToConditionFactor( modelFrame ) }
#      else
#         poolconds <- conditions(cds)

#      varianceFitDiagnosticsForMatrix(
#         counts(cds), sizeFactors(cds), rvf, poolconds )
#   }
#   
#   else if( cond == "_blind" )
#      varianceFitDiagnosticsForMatrix(
#         counts(cds), sizeFactors(cds), rvf )
#   
#   else 
#      varianceFitDiagnosticsForMatrix(
#         counts(cds)[,conditions(cds)==cond], 
#         sizeFactors(cds)[conditions(cds)==cond],
#         rvf ) 
#}

#residualsEcdfPlot <- function( cds, condition=NULL, ncuts=7 )
#{
#   stopifnot( is( cds, "CountDataSet" ) )   
#   ensureHasVarFuncs( cds )
#   fitdiag <- varianceFitDiagnostics( cds, condition )
#   residualsEcdfPlotFromDiagnostics( fitdiag, ncuts,
#      sprintf( "Residuals ECDF plot for condition '%s'", condition ) )
#}  

nbinomTest <- function( cds, condA, condB, pvals_only=FALSE, eps=1e-4 )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   if( cds@multivariateConditions )
      stop( "For CountDataSets with multivariate conditions, only the GLM-based test can be used." )
   stopifnot( condA %in% levels(conditions(cds)) )  
   stopifnot( condB %in% levels(conditions(cds)) )     
   
   colA <- conditions(cds)==condA
   colB <- conditions(cds)==condB

   bmv <- getBaseMeansAndVariances( counts(cds)[,colA|colB], 
      sizeFactors(cds)[colA|colB] )
   
   rawScvA <- fData(cds)[ , paste( "disp", dispTable(cds)[condA], sep="_" ) ]
   rawScvB <- fData(cds)[ , paste( "disp", dispTable(cds)[condB], sep="_" ) ]

   pval <- nbinomTestForMatrices( 
      counts(cds)[,colA], 
      counts(cds)[,colB], 
      sizeFactors(cds)[colA], 
      sizeFactors(cds)[colB], 
      rawScvA, 
      rawScvB,
      eps )
      
   if( pvals_only )
      pval
   else {
      bmvA <- getBaseMeansAndVariances( counts(cds)[,colA], sizeFactors(cds)[colA] )
      bmvB <- getBaseMeansAndVariances( counts(cds)[,colB], sizeFactors(cds)[colB] )
      data.frame( 
         id    = rownames( counts(cds) ),
         baseMean  = bmv$baseMean,
         baseMeanA = bmvA$baseMean,
         baseMeanB = bmvB$baseMean,
         foldChange = bmvB$baseMean / bmvA$baseMean,
         log2FoldChange = log2( bmvB$baseMean / bmvA$baseMean ), 
         pval = pval,
         padj = p.adjust( pval, method="BH" ), 
         # resVarA = cds@fitInfo[[ dispTable(cds)[condA] ]]$perGeneDispEsts / rawScvA,
         # resVarB = cds@fitInfo[[ dispTable(cds)[condB] ]]$perGeneDispEsts / rawScvB,
         stringsAsFactors = FALSE ) }
}

#scvPlot <- function( cds, xlim=NULL, ylim=c(0,.8), skipBiasCorrection = FALSE ) {
#   stopifnot( is( cds, "CountDataSet" ) )
#   ensureHasVarFuncs( cds )
#   
#   baseMeans <- getBaseMeansAndVariances( counts(cds), sizeFactors(cds) )$baseMean

#   xg <- exp( seq( log( max( min(baseMeans), 2/sum(sizeFactors(cds)) ) ), 
#      log( max(baseMeans) ), length.out = 1000 ) )
#   if( is.null( xlim ) )
#      xlim <- range( xg )
#   plot( NULL, xlim = xlim, ylim = ylim, yaxs="i",
#      xlab = "base mean", ylab = "squared coefficient of variation", log="x" )
#      
#   color <- 2
#   legendData <- data.frame( text=character(0), col=integer(0), lty=character(0), 
#      stringsAsFactors=FALSE ) 
#   for( funcName in unique( rawVarFuncTable(cds) ) ) {
#   
#      rawScv <- rawVarFunc( cds, funcName, byName=TRUE )( xg ) / xg^2      

#      lty <- if( funcName != "_max" ) "solid" else "dotted"
#      lines( xg, rawScv, col=color, lty=lty )
#      legendData <- rbind( legendData, data.frame( text=funcName, col=color, lty=lty,
#         stringsAsFactors=FALSE ) )
#	 
#      for( cond in names(rawVarFuncTable(cds))[ rawVarFuncTable(cds) == funcName ] ) {
#         if( cond != "_all" )
#	    cols <- which( conditions(cds) == cond )
#	 else
#	    cols <- 1:ncol(counts(cds))
#	 for( j in cols )
#	    lines( xg, rawScv + 1 / ( sizeFactors(cds)[[j]] * xg ), col=color, lty="dashed" )
#      }

#      color <- color + 1 
#   }
#   
#   dens <- density( log(baseMeans) )
#   lines( exp(dens$x), .7 * ylim[2] / max(dens$y) * dens$y, col=1, lty="solid" )
#   legendData <- rbind( legendData, data.frame( text="base mean density", col=1, lty="solid",
#      stringsAsFactors=FALSE ) )   

#   legend( "topright", 
#      legend = legendData$text,
#      col    = legendData$col,
#      lty    = legendData$lty )

#   invisible( NULL )
#}

getVarianceStabilizedData <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( ! "blind" %in% ls(cds@fitInfo) )
      stop( "Use 'estimateDispersions' with 'method=\"blind\"' before calling 'getVaraianceStabilizedData'" )
   ncounts <- t( t(counts(cds)) / sizeFactors(cds) )
   xg <- sinh( seq( asinh(0), asinh(max(ncounts)), length.out=1000 ) )[-1]
   xim <- mean( 1/sizeFactors(cds) )
   baseVarsAtGrid <- cds@fitInfo[["blind"]]$dispFunc( xg ) * xg^2 + xim * xg 
   integrand <- 1 / sqrt( baseVarsAtGrid )
   splf <- splinefun( 
      asinh( ( xg[-1] + xg[-length(xg)] )/2 ), 
      cumsum( 
         ( xg[-1] - xg[-length(xg)] ) * 
         ( integrand[-1] + integrand[-length(integrand)] )/2 ) )
   tc <- sapply( colnames(counts(cds)), function(clm)
      splf( asinh( ncounts[,clm] ) ) )
   rownames( tc ) <- rownames( counts(cds) )
   tc
}

makeExampleCountDataSet <- function( ) 
{
   ngenes <- 10000
   q0 <- rexp( ngenes, rate=1/250 )
   is_DE <- runif( ngenes ) < .3
   lfc <- rnorm( ngenes, sd=2 )
   q0A <- ifelse( is_DE, q0 * 2^(  lfc/2 ), q0 )
   q0B <- ifelse( is_DE, q0 * 2^( -lfc/2 ), q0 )
   true_sf <- c( 1., 1.3, .7, .9, 1.6 )   
   conds <- c( "A", "A", "B", "B", "B" )
   m <- t( sapply( 1:ngenes, function(i) 
      sapply( 1:5, function( j )
         rnbinom( 1, mu = true_sf[j] * ifelse( conds[j]=="A", q0A[i], q0B[i] ), 
            size = 1/.2 ) ) ) )
   colnames(m) <- c( "A1", "A2", "B1", "B2", "B3" )
   rownames(m) <- paste( "gene", 1:ngenes, 
      ifelse( is_DE, "T", "F" ), sep="_" )
   newCountDataSet( m, conds )
}

#nbinomFitGLM <- function( cds, modelFormula, glmControl=list() )
#{
#   stopifnot( is( cds, "CountDataSet" ) )
#   ensureHasVarFuncs( cds )
#   if( is.null( cds@rawVarFuncs[["_pooled"]] ) )
#      stop( "No pooled variance function found. Have you called 'estimateVarianceFunctions' with 'method=\"pooled\"'?" )
      
#   baseMeans <- colMeans(t(counts(cds))/sizeFactors(cds))
#   rawVars <- rawVarFunc( cds, "_pooled", TRUE )( baseMeans )
#   rawScv <- adjustScvForBias( rawVars/baseMeans^2, attr( rawVars, "size" ) )

#   nbinomGLMsForMatrix( counts(cds), sizeFactors(cds), rawScv, 
#      modelFormula, pData(cds), glmControl=glmControl )
#}

nbinomGLMTest <- function( resFull, resReduced )
   1 - pchisq( resReduced$deviance - resFull$deviance, 
   attr( resReduced, "df.residual" ) - attr( resFull, "df.residual" ) )

#getRawScvForSamplePair <- function( cds, sample1, sample2 ) {
#   if( any( is.na( sizeFactors(cds) ) ) )
#      stop( "Call 'estimateSizeFactors' first." )
#   data <- counts(cds)[ ,c( sample1, sample2 ) ]
#   sf <- sizeFactors(cds)[ c( sample1, sample2 ) ]
#   vf <- estimateVarianceFunctionForMatrix( data, sf )
#   bm <- t( t(data) / sf )
#   mean( adjustScvForBias( vf( bm ) / bm^2, 2 ), na.rm=TRUE ) }

#getRawScvDistanceMatrix <- function( cds ) {
#   res <- matrix( nrow=ncol(cds), ncol=ncol(cds) )
#   for( i in 1:ncol(cds) )
#      res[ i, i ] <- 0
#   for( i in 1:(ncol(cds)-1) )
#      for( j in (i+1):ncol(cds) ) {
#         res[ i, j ] <- getRawScvForSamplePair( cds, i, j )
#         res[ j, i ] <- res[ i, j ] }
#   colnames( res ) <- colnames( counts(cds) )
#   rownames( res ) <- colnames( counts(cds) )
#   res }
