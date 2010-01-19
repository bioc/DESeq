estimateSizeFactors <- function( cds )
{
   stopifnot( is( cds, "CountDataSet" ) )
   sizeFactors(cds) <- estimateSizeFactorsForMatrix( counts(cds) )
   cds
}

estimateVarianceFunctions <- function( cds, pool=FALSE, ... )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   if( any( is.na( sizeFactors(cds) ) ) )
      stop( "NAs found in size factors. Have you called already 'estimateSizeFactors'?" )
   
   cds@rawVarFuncs <- new.env( hash=TRUE )
   
   if( pool ) {

      cds@rawVarFuncs[["_pooled"]] <-
         estimateVarianceFunctionForMatrix( counts(cds), sizeFactors(cds), ... )
      rawVarFuncTable(cds) <- rep( "_pooled", length( levels( conditions(cds) ) ) )

   } else {
         
      replicated <- names( which( tapply( conditions(cds), conditions(cds), length ) > 1 ) )
      if( length( replicated ) < 1 )
         stop( "None of your conditions is replicated. Use pool=TRUE to pool across conditions." )
      nonreplicated <- names( which( tapply( conditions(cds), conditions(cds), length ) == 1 ) )
      for( cond in replicated )
         cds@rawVarFuncs[[cond]] <- estimateVarianceFunctionForMatrix( 
            counts(cds)[ , conditions(cds)==cond ], sizeFactors(cds)[ conditions(cds)==cond ], ... )
      cds@rawVarFuncs[["_max"]] <- function( x, reportSize=FALSE )
         apply( rbind( sapply( replicated, function(cond) 
	    cds@rawVarFuncs[[cond]]( x, reportSize ) ) ), 1, max )
         
      rawVarFuncTable(cds) <- 
         sapply( levels(conditions(cds)), function( cond )
            ifelse( cond %in% replicated, cond, "_max" ) )
   }
        
   cds
}

varianceFitDiagnostics <- function( cds, cond )
{
   stopifnot( is( cds, "CountDataSet" ) )
   stopifnot(cond %in% levels(conditions(cds)) )  
   ensureHasVarFuncs( cds )
      
   varianceFitDiagnosticsForMatrix(
      counts(cds)[,conditions(cds)==cond], 
      sizeFactors(cds)[conditions(cds)==cond],
      rawVarFunc( cds, cond ) )
}

residualsEcdfPlot <- function( cds, condition, ncuts=7 )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   ensureHasVarFuncs( cds )
   fitdiag <- varianceFitDiagnostics( cds, condition )
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

   rawScvA <- rawVarFunc( cds, condA )( bmv$baseMean ) / bmv$baseMean^2
   rawScvB <- rawVarFunc( cds, condB )( bmv$baseMean ) / bmv$baseMean^2
   
   rawScvA <- adjustScvForBias( rawScvA, rawVarFunc( cds, condA )( reportSize=TRUE ) )
   rawScvB <- adjustScvForBias( rawScvB, rawVarFunc( cds, condB )( reportSize=TRUE ) )

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
         resVarB = bmvB$baseVar / rawVarFunc( cds, condB )( bmv$baseMean ),
         stringsAsFactors = FALSE ) }
}

scvPlot <- function( cds, xlim=NULL, ylim=NULL ) {
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )
   
   baseMeans <- getBaseMeansAndVariances( counts(cds), sizeFactors(cds) )$baseMean
   xg <- exp( seq( log( max( min(baseMeans), 2/sum(sizeFactors(cds)) ) ), 
      log( max(baseMeans) ), length.out = 1000 ) )
   rcv <- eapply( cds@rawVarFuncs, function( rvf ) rvf( xg ) / xg^2 )
   bcv <- sapply( 1:ncol(counts(cds)), function(j)
      1 / ( sizeFactors(cds)[[j]] * xg ) + rawVarFunc( cds, conditions(cds)[j] )( xg ) / xg^2 )
   colnames( bcv ) <- colnames( counts( cds ) )

   if( is.null( xlim ) )
      xlim <- range( xg )
   if( is.null( ylim ) )
      ylim <- c( 0, max( sapply( bcv, max ) ) )
   plot( NULL, xlim = xlim, ylim = ylim, yaxs="i",
      xlab = "base mean", ylab = "squared coefficient of variation", log="x" )
      

   nonmax <- names(rcv)[ names(rcv) != "_max" ]   
   cols <- 1 + (1:length(nonmax))
   names(cols) <- names(rcv[nonmax])
   if( "_max" %in% names(rcv) )
      cols <- c( cols, `_max`=1 )

   for( j in 1:ncol( bcv ) )
      lines( xg, bcv[,j], lty="dotted", 
         col=cols[ conditions(cds)[j] ] )

   for( j in nonmax ) 
      lines( xg, rcv[[j]], col=cols[j], lty="solid" )
      
   if( "_max" %in% names(rcv) )
      lines( xg, rcv[["_max"]], col=1, lty="dashed" )
   
   dens <- density( log(baseMeans) )
   lines( exp(dens$x), .7 * ylim[2] / max(dens$y) * dens$y, col=1, lty="solid" )

   legend( "topright", 
      legend = c( names(cols), "base mean density" ),
      col    = c( cols, 1 ),
      lty    = ifelse( c( names(cols), "" ) != "_max", "solid", "dashed" ) )
   invisible( NULL )
}

getVarianceStabilizedData <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )
   tc <- sapply( colnames(counts(cds)), function(clm) {   
      countcol <- counts(cds)[,clm]
      sf <- sizeFactors(cds)[clm]
      xg <- sinh( seq( asinh(0), asinh(max(countcol)), length.out=1000 ) )
      xim <- sum( 1/sizeFactors(cds) ) / ncol( counts(cds) )
      fullVarsAtGrid <- pmax( xg,
         xg + sf^2*( rawVarFunc( cds, as.character(conditions(cds)[clm]) )( xg/sf ) ) )
      integrand <- 1 / sqrt( fullVarsAtGrid )
      splf <- splinefun( asinh(xg[-1]), cumsum( (xg[-1]-xg[-length(xg)]) * integrand[-1] ) )
      splf( asinh(countcol) )
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
   m <- t( sapply( 1:ngenes, function(i) 
      sapply( true_sf, function( s )
         rnbinom( 1, mu = s * q0A[i], size = 1/.2 ) ) ) )
   colnames(m) <- c( "A1", "A2", "A3", "B1", "B2" )
   rownames(m) <- paste( "gene", 1:ngenes, 
      ifelse( is_DE, "T", "F" ), sep="_" )
   newCountDataSet( m, c( "A", "A", "B", "B", "B" ) )
}

