/* Call to LAPACK matrix inversion dgesv */
/* ~60 times quicker than the Numerical Recipes standard solver */
/* Author : Jan Vatant d'Ollone 2018 */

#include "titan.h"

/* DGESV prototype */
extern void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv,
               double* b, int* ldb, int* info);

/*Main*/
void solve_lapack( double ***aa, int m, int n0, int n1 )
{
   int i,j,k;
   int n=n1-n0+1, ok;
   int ipiv[n];
   double a[n*n];
   double b[n*n];
   
   for( i = 0; i < n; i++ ) for( j = 0; j < n; j++ ) a[j+n*i]=aa[m][i+n0][j+n0];

   for( i = 0; i < n; i++ )
   {
      for( k = 0; k < n; k++ )
         b[k+n*i] = 0.0e0;
      b[i+n*i] = 1.0e0;
   }

   dgesv_(&n,&n,a,&n,ipiv,b,&n,&ok); /* Lapack subroutine */
   
   if( ok > 0 ) {
       printf( "The diagonal element of the triangular factor of A,\n" );
       printf( "U(%i,%i) is zero, so that A is singular;\n", ok, ok );
       printf( "the solution could not be computed.\n" );
       exit( 1 );
   }
   
   for( i = 0; i < n; i++ ) for( j = 0; j < n; j++ ) aa[m][i+n0][j+n0]=b[j+n*i];

}
