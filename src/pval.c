#include <math.h>
#include <limits.h>
#include <Rinternals.h>
#include <Rmath.h>

void add_from_middle( long int kS, double pobs, double muA, double vA, 
      double muB, double vB, int upwards, double eps,
      double * res_total, double * res_obstotal)
{
   long int k = kS * muA / (muA+muB);
   if( !upwards )
      k++;
   double total = 0;
   double esttotal = Rf_dnbinom( kS, (muA+muB)*(muA+muB) / (vA+vB-muA-muB), (muA+muB) / (vA+vB), 0 );
   double obstotal = 0;   
   long int step = 1;
   double val = 0;
   double sizeA = muA*muA / (vA-muA);
   double probA = muA/vA;
   double sizeB = muB*muB / (vB-muB);
   double probB = muB/vB;
   long int knew;
   double prevval;
   while( upwards ? (k < kS) : (k > 0) ) {
      prevval = val;
      while( 1 ) {
         if( upwards ) {
            if( k + step > kS )
               step = kS - k;
            knew = k + step; 
         } else {
            if( k - step < 0 )
               step = k;
            knew = k - step; 
         }
         val = Rf_dnbinom( knew, sizeA, probA, 0 ) * Rf_dnbinom( kS - knew, sizeB, probB, 0 );
         if( (step == 1) || (fabs( val - prevval ) < eps * esttotal/kS)  )
            break;
         step >>= 1;
         //Rprintf( "step decreased, k=%d, step=%d\n", k, step );         
      }
      k = knew;
      total += val * step;
      //Rprintf( "k=%d total=%g val=%g pobs=%g\n", k, total, val, pobs );
      if( val <= pobs ) {
         //Rprintf( "val<=pobs at k=%d\n", k );
         if( prevval <= pobs ) {
            obstotal += val * step;
            //Rprintf( "adding %g\n", val * step );
         } else {
            obstotal += val * ( 1 + (step-1) * ( pobs - val ) / ( prevval - val ) );
            //Rprintf( "adding %g (')\n", val * ( 1 + (step-1) * ( pobs - val ) / ( prevval - val ) ) );
            // ^^ TO DO: Double check whether this makes sense!!
         }
      }
      if( (step < INT_MAX) && (fabs( val - prevval ) < eps * esttotal/kS / 4) ) {
         step <<= 1;
         //Rprintf( "step increased, k=%d, step=%d\n", k, step );
      }
   }
   * res_total = total;
   * res_obstotal = obstotal;
}

SEXP add_from_middle_for_R( SEXP kS, SEXP pobs, SEXP muA, SEXP vA, 
   SEXP muB, SEXP vB, SEXP upwards, SEXP eps ) 
{
   SEXP res = Rf_allocVector( REALSXP, 2 );
   add_from_middle( INTEGER(kS)[0], REAL(pobs)[0], REAL(muA)[0], 
      REAL(vA)[0], REAL(muB)[0], REAL(vB)[0], LOGICAL(upwards)[0], 
      REAL(eps)[0], REAL(res), REAL(res) + 1 );
  return res;
}   
