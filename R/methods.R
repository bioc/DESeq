estimateSizeFactors <- function( cds )
{
   stopifnot( is( cds, "CountDataSet" ) )
   sizeFactors(cds) <- estimateSizeFactorsForMatrix( counts(cds) )
   cds
}

estimateVarianceFunctions <- function( cds, pool=NULL, 
   method = c( "normal", "blind", "pooled" ), 
   locfit_extra_args=list(), lp_extra_args=list(), modelFrame = NULL )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   if( any( is.na( sizeFactors(cds) ) ) )
      stop( "NAs found in size factors. Have you called already 'estimateSizeFactors'?" )
   
   if( length(method) != 3 & !is.null( pool ) )
      stop( "Do not specify both the 'pool' and the 'method' argument." )      
   if( !is.null( pool ) ) {
      if( pool == FALSE )
         method <- "normal"
      else if( pool == TRUE )
         method <- "blind"
      else
         stop( "Argument 'pool' used incorrectly." )
      warning( "The 'pool' argument to 'estimatevarianceFunction' is deprecated. Use the 'method' argument instead." )
   } else
      method <- match.arg( method )
   
   if( cds@multivariateConditions && method != "pooled" )
      stop( "You have specified multivariate conditions (i.e., passed a data frame with conditions). In this case, only the variance estimation is supported only with method='pooled'." )
   
   cds@rawVarFuncs <- new.env( hash=TRUE )
   
   if( method == "blind" ) {
      cds@rawVarFuncs[["_blind"]] <-
         estimateVarianceFunctionForMatrix( counts(cds), 
	    sizeFactors(cds), locfit_extra_args, lp_extra_args )
      rawVarFuncTable(cds) <- data.frame(
         row.names = levels(conditions(cds)),
         funcName = rep( "_blind", length( levels( conditions(cds) ) ) ),
         varAdjFactor = rep( 1, length( levels( conditions(cds) ) ) ),
         stringsAsFactors = FALSE ) }

   else if( method == "normal" ) {
      replicated <- names( which( tapply( conditions(cds), conditions(cds), length ) > 1 ) )
      if( length( replicated ) < 1 )
         stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions." )
      nonreplicated <- names( which( tapply( conditions(cds), conditions(cds), length ) == 1 ) )
      for( cond in replicated )
         cds@rawVarFuncs[[cond]] <- estimateVarianceFunctionForMatrix( 
            counts(cds)[ , conditions(cds)==cond ], sizeFactors(cds)[ conditions(cds)==cond ], 
	       locfit_extra_args, lp_extra_args )
      cds@rawVarFuncs[["_max"]] <- function( q ) {
         a <- lapply( replicated, function(cond) cds@rawVarFuncs[[cond]]( q ) )
         ans <- apply( array( unlist( a, recursive=FALSE ), dim = c( length(q), length(a) ) ), 1, max )
         rownames(ans) <- rownames(q)
         attr( ans, "size" ) <- min( sapply( a, attr, "size" ) )
         ans
      }
         
      rawVarFuncTable(cds) <- data.frame(
         funcName = sapply( levels(conditions(cds)), function( cond )
            ifelse( cond %in% replicated, cond, "_max" ) ),
         varAdjFactor = rep( 1, length( levels(conditions(cds)) ) ),
         stringsAsFactors = FALSE ) }

   else if( method == "pooled" ) {
   
      if( cds@multivariateConditions ) {
         if( is.null( modelFrame ) )
            modelFrame <- pData(cds)[ , colnames(pData(cds)) != "sizeFactor" ]
         conds <- modelMatrixToConditionFactor( modelFrame ) }
      else
         conds <- conditions(cds)

      cds@rawVarFuncs[["_pooled"]] <- estimatePooledVarianceFunctionForMatrix( 
         counts(cds), sizeFactors(cds), conds,
         locfit_extra_args, lp_extra_args )
         
      if( ! cds@multivariateConditions )
         rawVarFuncTable(cds) <- data.frame(
            funcName = sapply( levels(conditions(cds)), function( cond )
               "_pooled" ),
            varAdjFactor = rep( 1, length( levels(conditions(cds)) ) ),
            stringsAsFactors = FALSE ) }
      else {
         rawVarFuncTable(cds) <- data.frame( 
            row.names = "_pooled",
            funcName = "_pooled", 
            varAdjFactor = 1,
            stringsAsFactors = FALSE ) }
        
   validObject( cds )
   cds
}

varianceFitDiagnostics <- function( cds, cond, ignoreVarAdjFactor=FALSE )
{
   stopifnot( is( cds, "CountDataSet" ) )
   stopifnot(cond %in% levels(conditions(cds)) )  
   ensureHasVarFuncs( cds )
   
   rvf <-rawVarFunc( cds, cond ) 
   varianceFitDiagnosticsForMatrix(
      counts(cds)[,conditions(cds)==cond], 
      sizeFactors(cds)[conditions(cds)==cond],
      if( ignoreVarAdjFactor ) rvf else function(q) rvf(q) * attr(rvf,"varAdjFactor") )
}

residualsEcdfPlot <- function( cds, condition, ncuts=7, ignoreVarAdjFactor=FALSE )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   ensureHasVarFuncs( cds )
   fitdiag <- varianceFitDiagnostics( cds, condition,ignoreVarAdjFactor )
   residualsEcdfPlotFromDiagnostics( fitdiag, ncuts,
      sprintf( "Residuals ECDF plot for condition '%s'", condition ) )
}  

nbinomTest <- function( cds, condA, condB, pvals_only=FALSE, eps=1e-4 )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   ensureHasVarFuncs( cds )
   if( cds@multivariateConditions )
      stop( "For CountDataSets with multivariate conditions, only the GLM-based test can be used." )
   stopifnot( condA %in% levels(conditions(cds)) )  
   stopifnot( condB %in% levels(conditions(cds)) )     
   
   colA <- conditions(cds)==condA
   colB <- conditions(cds)==condB

   bmv <- getBaseMeansAndVariances( counts(cds)[,colA|colB], 
      sizeFactors(cds)[colA|colB] )

   rvfA <- rawVarFunc( cds, condA )
   rvfB <- rawVarFunc( cds, condB )
   
   rawScvA <- rvfA( bmv$baseMean ) * attr( rvfA, "varAdjFactor" ) / bmv$baseMean^2
   rawScvB <- rvfB( bmv$baseMean ) * attr( rvfB, "varAdjFactor" ) / bmv$baseMean^2
   
   rawScvA <- adjustScvForBias( rawScvA, attr( rawScvA, "size" ) )
   rawScvB <- adjustScvForBias( rawScvB, attr( rawScvB, "size" ) )

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
         resVarA = bmvA$baseVar / ( bmvA$baseMean * sum( 1/sizeFactors(cds)[colA] ) / length(condA) +
            rawVarFunc( cds, condA )( bmv$baseMean ) ),
         resVarB = bmvA$baseVar / ( bmvA$baseMean * sum( 1/sizeFactors(cds)[colB] ) / length(condB) +
            rawVarFunc( cds, condB )( bmv$baseMean ) ),
         stringsAsFactors = FALSE ) }
}

scvPlot <- function( cds, xlim=NULL, ylim=c(0,.8), ignoreVarAdjFactors = FALSE,
      skipBiasCorrection = FALSE ) {
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )
   
   baseMeans <- getBaseMeansAndVariances( counts(cds), sizeFactors(cds) )$baseMean
   xg <- exp( seq( log( max( min(baseMeans), 2/sum(sizeFactors(cds)) ) ), 
      log( max(baseMeans) ), length.out = 1000 ) )

   if( ! cds@multivariateConditions )
      conds <- conditions(cds)
   else
      conds <- factor( rep( "_pooled", ncol(counts(cds)) ) )
   
   rscv <- sapply( levels(conds), function(n) {
      if( !cds@multivariateConditions )
         rawScv <- rawVarFunc( cds, n )( xg ) / xg^2
      else
         rawScv <- rawVarFunc( cds )( xg ) / xg^2
      if( !ignoreVarAdjFactors )
         rawScv <- adjustScvForBias( rawScv, attr( rawScv, "size" ) )
      rawScv } )
      
   if( !ignoreVarAdjFactors & !cds@multivariateConditions )
      rscv <- sapply( levels(conds), function(n) 
         rscv[,n] * attr( rawVarFunc( cds, n ), "varAdjFactor" ) )
           
   bscv <- sapply( 1:ncol(counts(cds)), function(j)
      1 / ( sizeFactors(cds)[[j]] * xg ) + rscv[ , as.character(conds)[j] ] )
   colnames( bscv ) <- colnames( counts( cds ) )

   if( is.null( xlim ) )
      xlim <- range( xg )
   plot( NULL, xlim = xlim, ylim = ylim, yaxs="i",
      xlab = "base mean", ylab = "squared coefficient of variation", log="x" )      


   if( !cds@multivariateConditions )
      rvft <- rawVarFuncTable(cds)
   else
      rvft <- data.frame( row.names = "_pooled", funcName="_pooled", varAdjFactor=1 )

   cols <- 1 + ( 1 : length(unique(rvft$funcName)) )
   names(cols) <- unique(rvft$funcName)

   for( j in 1:ncol( bscv ) )
      lines( xg, bscv[,j], lty="dashed", 
	     col=cols[ rvft[ as.character(conds)[j], "funcName" ] ] )

   for( n in levels(conds) )
      lines( xg, rscv[,n], col=cols[ rvft[ n, "funcName" ] ], 
         lty = ifelse( rvft[ n, "funcName" ] != "_max", "solid", "dotted" ) )
   
   dens <- density( log(baseMeans) )
   lines( exp(dens$x), .7 * ylim[2] / max(dens$y) * dens$y, col=1, lty="solid" )

   legend( "topright", 
      legend = c( names(cols), "base mean density" ),
      col    = c( cols, 1 ),
      lty    = ifelse( c( names(cols), "" ) != "_max", "solid", "dotted" ) )
   invisible( NULL )
}

getVarianceStabilizedData <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   rvf <- estimateVarianceFunctionForMatrix( counts(cds), sizeFactors(cds) )
   tc <- sapply( colnames(counts(cds)), function(clm) {
      countcol <- counts(cds)[,clm]
      sf <- sizeFactors(cds)[clm]
      xg <- sinh( seq( asinh(0), asinh(max(countcol)), length.out=1000 ) )
      xim <- sum( 1/sizeFactors(cds) ) / ncol( counts(cds) )
      fullVarsAtGrid <- pmax( xg, xg + sf^2*( rvf( xg/sf ) ) )
      integrand <- 1 / ( sqrt( fullVarsAtGrid ) / sf )
      splf <- splinefun( asinh(xg[-1]/sf), cumsum( (xg[-1]-xg[-length(xg)])/sf * integrand[-1] ) )
      splf( asinh(countcol/sf) )
   } )
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
   colnames(m) <- c( "A1", "A2", "A3", "B1", "B2" )
   rownames(m) <- paste( "gene", 1:ngenes, 
      ifelse( is_DE, "T", "F" ), sep="_" )
   newCountDataSet( m, conds )
}

nbinomFitGLM <- function( cds, modelFormula )
{
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )
   if( is.null( cds@rawVarFuncs[["_pooled"]] ) )
      stop( "No pooled variance function found. Have you called 'estimateVarianceFunctions' with 'method=\"pooled\"'?" )
      
   baseMeans <- colMeans(t(counts(cds))/sizeFactors(cds))
   rawVars <- rawVarFunc( cds, "_pooled", TRUE )( baseMeans )
   rawScv <- adjustScvForBias( rawVars/baseMeans^2, attr( rawVars, "size" ) )

   nbinomGLMsForMatrix( counts(cds), sizeFactors(cds), rawScv, 
      modelFormula, pData(cds) )
}

nbinomGLMTest <- function( resFull, resReduced )
   1 - pchisq( resReduced$deviance - resFull$deviance, 
   attr( resReduced, "df.residual" ) - attr( resFull, "df.residual" ) )

