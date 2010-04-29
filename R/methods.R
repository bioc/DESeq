estimateSizeFactors <- function( cds )
{
   stopifnot( is( cds, "CountDataSet" ) )
   sizeFactors(cds) <- estimateSizeFactorsForMatrix( counts(cds) )
   cds
}

estimateVarianceFunctions <- function( cds, pool=FALSE, 
   locfit_extra_args=list(), lp_extra_args=list() )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   if( any( is.na( sizeFactors(cds) ) ) )
      stop( "NAs found in size factors. Have you called already 'estimateSizeFactors'?" )
   
   cds@rawVarFuncs <- new.env( hash=TRUE )
   
   if( pool ) {

      cds@rawVarFuncs[["_pooled"]] <-
         estimateVarianceFunctionForMatrix( counts(cds), 
	    sizeFactors(cds), locfit_extra_args, lp_extra_args )
      rawVarFuncTable(cds) <- data.frame(
         row.names = levels(conditions(cds)),
         funcName = rep( "_pooled", length( levels( conditions(cds) ) ) ),
         varAdjFactor = rep( 1, length( levels( conditions(cds) ) ) ),
         stringsAsFactors = FALSE )

   } else {
         
      replicated <- names( which( tapply( conditions(cds), conditions(cds), length ) > 1 ) )
      if( length( replicated ) < 1 )
         stop( "None of your conditions is replicated. Use pool=TRUE to pool across conditions." )
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
         stringsAsFactors = FALSE )
   }
        
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

nbinomTest <- function( cds, condA, condB, pvals_only=FALSE )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   ensureHasVarFuncs( cds )
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
      rawScvB )
      
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
         resVarA = bmvA$baseVar / ( bmvA$baseMean + rawVarFunc( cds, condA )( bmv$baseMean ) ),
         resVarB = bmvB$baseVar / ( bmvA$baseMean + rawVarFunc( cds, condB )( bmv$baseMean ) ),
         stringsAsFactors = FALSE ) }
}

scvPlot <- function( cds, xlim=NULL, ylim=c(0,.8), ignoreVarAdjFactors = FALSE ) {
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )
   
   baseMeans <- getBaseMeansAndVariances( counts(cds), sizeFactors(cds) )$baseMean
   xg <- exp( seq( log( max( min(baseMeans), 2/sum(sizeFactors(cds)) ) ), 
      log( max(baseMeans) ), length.out = 1000 ) )
   rscv <- sapply( levels(conditions(cds)), function(n) rawVarFunc( cds, n )( xg ) / xg^2 )
   if( !ignoreVarAdjFactors )
      rscv <- sapply( levels(conditions(cds)), function(n) 
         rscv[,n] * attr( rawVarFunc( cds, n ), "varAdjFactor" ) )
      
   bscv <- sapply( 1:ncol(counts(cds)), function(j)
      1 / ( sizeFactors(cds)[[j]] * xg ) + rscv[ , as.character(conditions(cds))[j] ] )
   colnames( bscv ) <- colnames( counts( cds ) )

   if( is.null( xlim ) )
      xlim <- range( xg )
   plot( NULL, xlim = xlim, ylim = ylim, yaxs="i",
      xlab = "base mean", ylab = "squared coefficient of variation", log="x" )      

   cols <- 1 + ( 1 : length(unique(rawVarFuncTable(cds)$funcName)) )
   names(cols) <- unique(rawVarFuncTable(cds)$funcName)

   for( j in 1:ncol( bscv ) )
      lines( xg, bscv[,j], lty="dashed", 
	     col=cols[ rawVarFuncTable(cds)[ as.character(conditions(cds))[j], "funcName" ] ] )

   for( n in levels(conditions(cds)) )
      lines( xg, rscv[,n], col=cols[ rawVarFuncTable(cds)[ n, "funcName" ] ], 
         lty = ifelse( rawVarFuncTable(cds)[ n, "funcName" ] != "_max", "solid", "dotted" ) )
   
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

