library( Biobase )
library( MASS )
library( genefilter )
library( locfit )

source( "R/class_and_slots.R" )
source( "R/methods.R" )
source( "R/core.R" )
dyn.load( "src/pval.so" )

load( "../RNA-Seq-stats/counts_naga.RData" )

cds <- newCountDataSet( 
   counts.naga, 
   c( "dT", "dTt", "dT", "RH", "RHt", "RH" ) )
   
cds <- estimateSizeFactors( cds )
print( sizeFactors(cds) )

cds <- estimateVarianceFunctions( cds )


df <- baseVarDiagnostics( cds, "dT" )
plot( baseVar ~ baseMean, df, log="xy", pch='.' ) 
lines( fittedBaseVar ~ baseMean, df[ order(df$baseMean), ], pch='.', col="red" )

hist( df$pchisq )

library(geneplotter)
residualsEcdfPlot( cds, "dT" )

p <- nbinomTestForContrast( cds, "dT", "RH" )

dnbinomMV <- function( k, mu, var ) 
   dnbinom( k, prob = mu/var, size = mu^2/(var-mu) )
   
rnbinomMV <- function( n, mu, var ) 
   rnbinom( n, prob = mu/var, size = mu^2/(var-mu) )   
   
