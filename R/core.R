estimateSizeFactorsForMatrix <- function( counts )
{
   geomeans <- exp( rowMeans( log(counts) ) )
   apply( counts, 2, function(cnts) 
      median( ( cnts / geomeans )[ geomeans>0 ] ) )
}


getBaseMeansAndVariances <- function( counts, sizeFactors ) {

   # Devides the counts by sizeFactors and calculates the estimates for
   # base means and variances for each gene.
   
   data.frame(
      baseMean = rowMeans( t( t(counts) / sizeFactors ) ),
      baseVar = rowVars( t( t(counts) / sizeFactors ) ) )
}   

estimateVarianceFunctionForMatrix <- function( counts, sizeFactors, ... ) {

   # This function should be called with a matrix of counts adjusted for
   # library size ratios, whose columns are replicates of an experimental
   # condition. It returns a function that, given a adjusted count value, 
   # returns an estimate for the variance at this count level.

   stopifnot( ncol( counts ) == length( sizeFactors ) )
   counts <- counts[ rowSums(counts) > 0, ]
   bmv <- getBaseMeansAndVariances( counts, sizeFactors ) 
   fit <- locfit( baseVar ~ lp( log(baseMean), ...), data = bmv, family="gamma" )
   xim <- sum( 1/sizeFactors ) / length( sizeFactors )
   rm( counts )
   rm( bmv )
   
   function( q, reportSize=FALSE ) 
      if( reportSize )
         length( sizeFactors )
      else
         pmax( safepredict( fit, log(q) ) - xim * q, 1e-8 * q )
   # Note: The 'pmax' construct above serves to limit the overdispersion to a minimum
   # of 10^-8, which should be indistinguishable from 0 but ensures numerical stability.
}   
   
safepredict <- function( fit, x )
{
   # A wrapper around predict to avoid the issue that predict.locfit cannot
   # propagate NAs and NaNs properly.

   res <- rep.int( NA_real_, length(x) )
   res[ is.finite(x) ] <- predict( fit, x[is.finite(x)] )
   res   
}

nbinomTestForMatricesRaw <- function( kA, kB, muA, vA, muB, vB, eps=0 )
{
   # Let kA and kB be two count observations from two random variables, for
   # which the null hypothesis assumes negative binomial distributions with
   # means muA, muB and vA and vB. Calculate the probability that kA and kB
   # or as or more extreme counts are observed, conditioned to the sum of the 
   # counts being kA+kB. "As or more extreme" means having conditional probability 
   # at most 
   #                     fNB( kA, muA, vA ) fNB( kA, muA, vA )
   #       -------------------------------------------------------------------
   #        sum of fNB( k, muA, vA ) fNB( kA+kB-k, muA, vA ) for k=0,..,kA+kB
   #
   # 'eps' is a roughly followed guidance on the required presision
   
   if( !all( is.finite( c( kA, kB, muA, vA, muB, vB, eps ) ) ) )
      return( NA )

   pobs <- dnbinom( kA, prob = muA/vA, size = muA^2/(vA-muA) ) * 
           dnbinom( kB, prob = muB/vB, size = muB^2/(vB-muB) )
           
   stopifnot( is.finite( pobs ) )
   
   pobs <- pobs * ( 1 + 1e-7 )
   # This is to avoid rounding errors in checking for p <= pobs

   totals <- 
      .Call( "add_from_middle_for_R", as.integer(kA+kB), pobs, muA, vA, muB, vB, FALSE, eps ) + 
      .Call( "add_from_middle_for_R", as.integer(kA+kB), pobs, muA, vA, muB, vB, TRUE, eps )
   unname( totals[2] / totals[1] )
}


nbinomTestForMatrices <- function( countsA, countsB, sizeFactorsA, sizeFactorsB, 
   rawScvA, rawScvB, eps=3e-3 )
{
   kAs <- rowSums( cbind(countsA) )
   kBs <- rowSums( cbind(countsB) )
   
   baseMeans <- rowMeans( cbind(      
      t( t( countsA ) / sizeFactorsA ),
      t( t( countsB ) / sizeFactorsB ) ) )      
   muAs <- baseMeans * sum( sizeFactorsA )
   muBs <- baseMeans * sum( sizeFactorsB )

   fullVarA <- muAs + rawScvA * baseMeans^2 * sum(sizeFactorsA^2)
   fullVarB <- muBs + rawScvB * baseMeans^2 * sum(sizeFactorsB^2)
   
   sapply( 1:nrow(cbind(countsA)), function(i) {
      nbinomTestForMatricesRaw( kAs[i], kBs[i], muAs[i], fullVarA[i], muBs[i], fullVarB[i], eps )
   } )
}


varianceFitDiagnosticsForMatrix <- function( counts, sizeFactors, rawVarFunc )
{
   res <- getBaseMeansAndVariances( counts, sizeFactors )
   res$fittedRawVar <- rawVarFunc( res$baseMean )
   res$fittedBaseVar <- res$fittedRawVar + 
      res$baseMean * sum( 1/sizeFactors ) / length( sizeFactors )
   df <- ncol( cbind(counts) ) - 1
   res$pchisq <- pchisq( df * res$baseVar / res$fittedBaseVar, df = df )
   res
}

multiecdfWithoutLegend <- function( x, ... )
{
   if( all( package_version( packageDescription( "geneplotter", fields="Version" ) ) 
          >= package_version( "1.21.2" ) ) )
      multiecdf( x, ..., legend=NULL )
   else
      multiecdf( x, ... ) 
}

residualsEcdfPlotFromDiagnostics <- function( fitdiag, ncuts=7, 
      plotTitle="Residuals ECDF plot" )
{
   ok <- !is.na(fitdiag$pchisq)
   cols <- colorRampPalette( c("red","blue") )( ncuts )
   cuts <- factor( cut( rank(fitdiag$baseMean[ok]), ncuts ) )

   multiecdfWithoutLegend( 
      fitdiag$pchisq[ok] ~ cuts,
      col = cols,
      xlab = "chi-squared probability of residual",
      ylab = "ECDF",
      main = plotTitle
   )
      
   segments( 0, 0, 1, 1, col="darkgreen" )
   legend( 0, 1, 
      c( sprintf( "%.1e .. %.1e", 
            tapply( fitdiag$baseMean[ok], cuts, min ), 
            tapply( fitdiag$baseMean[ok], cuts, max ) ), 
         "expected" ),
      col = c( cols, "darkgreen" ), lty="solid" )
}  

# Note: The following function is never called; it is here only for
# documentation purposes, as it has been used to produce the data object
# scvBiasCorrectionFits, which is stored in the file
# inst/scvBiasCorrectionFits.rda, gets loadewd by the line after this
# function and is used by the function adjustScvForBias

prepareScvBiasCorrectionFits <- function( maxnrepl=15, mu=100000, ngenes=10000,
      true_raw_scv = c( seq( 0, 2, length.out=100 )[-1], seq( 2, 10, length.out=20 )[-1] ) )
   lapply( 2:maxnrepl, function( m ) {
      est_raw_scv <- sapply( true_raw_scv, function( alpha ) {
	 k <- matrix( rnbinom( ngenes*m, mu=mu, size=1/alpha ), ncol=m )
	 k <- k[ rowSums(k)>0, ]
	 mean( rowVars(k) / rowMeans(k)^2 ) } )
      locfit( true_raw_scv ~ lp( est_raw_scv, nn=.2 ) ) } )

load( system.file ( "scvBiasCorrectionFits.rda", package="DESeq" ) )

adjustScvForBias <- function( scv, nsamples ) 
{
   if( nsamples - 1 > length( scvBiasCorrectionFits ) )
      scv
   else
      safepredict( scvBiasCorrectionFits[[ nsamples-1 ]], scv )
}      
