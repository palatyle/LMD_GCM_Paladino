/* omega: correction to perfect gas compound */
/* 30 Oct 96 */

#include "titan.h"

double omega( ts, epsa, epsb )
double epsa,epsb,ts;
{
   double t;

   t = ts / sqrt( epsa * epsb );
   return 1.06036e0 * pow( t, -0.1561e0 ) + 0.193e0 * exp( -0.47635e0 * t )
          + 1.03587e0 * exp( -1.52996e0 * t ) + 1.76474e0 * exp( -3.89411e0 * t );
}
