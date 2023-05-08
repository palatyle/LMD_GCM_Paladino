/* matrix inversion */
/* cf Numerical Recipes LU Method for equations numbers */
/* GCCM */
/* similaire a inv (GP) */
/* la matrice a est inversee seulement sur le bloc [n0;n1][n0;n1] */

#include "titan.h"

void solve_b( double ***aa, double **f, int m, int n0, int n1 )
{
   int    i,ii,imax,j,k,l,ll;
   double **a,aamax,dum,*indx,sum,*vv,tmp;
   FILE   *out;

   indx = dm1d( n0, n1 );
   vv   = dm1d( n0, n1 );
   a    = dm2d( n0, n1, n0, n1 );
   imax = n0;

   for( i = n0; i <= n1; i++ ) for( j = n0; j <= n1; j++ ) a[i][j]=aa[m][i][j];
   
   for( i = n0; i <= n1; i++ )
   {
      aamax = 0.0e0;
      for( k = n0; k <= n1; k++ )
         if( (tmp=fabs(a[i][k])) > aamax ) aamax = tmp;
      if( aamax < 1.0e-20 )
      {
         out = fopen( "err.log", "a" );
         fprintf( out, "Singular matrix. n0=%ld k=%ld aamax=%le\n",
                        n0,k,aamax);
         fclose( out );
         exit(0);   
/*         aamax = 1.e-30; */
      }
      vv[i] = 1.0e0 / aamax;   /* Save the scaling */
/*
      if( (aamax > 1.0e100)||(aamax < 1.0e-100) )
      {
         out = fopen( "err.log", "a" );
         fprintf( out, "ATTENTION aamax = %le\n", aamax );
         fclose( out );
         exit( 0 );
      }
*/
   }
   for( k = n0; k <= n1; k++ )
   {
      for( i = n0; i < k; i++ )      /* This is equation 2.3.12 except for i = j */
      {
         sum = a[i][k];
         for( l = n0; l < i; l++ )
            sum -= ( a[i][l] * a[l][k] );
         a[i][k] = sum;
      }
      aamax = 0.0e0;                 /* Initialize for the search for largest pivot element */
      for( i = k; i <= n1; i++ )     /* This is i = j of equation 2.3.12 and */
      {
         sum = a[i][k];              /* i = J + 1,...,N of equation 2.3.13 */
         for( l = n0; l < k; l++ )
            sum -= ( a[i][l] * a[l][k] );
         a[i][k] = sum;
         dum        = vv[i] * fabs(sum); /* Figure of merit for the pivot */
         if( dum >= aamax )          /* Is it better than the best so far ? */
         {
            imax  = i;
            aamax = dum;
         }
      }
      if( k != imax )                /* Do we need to interchange rows ? */
      {
         for( l = n0; l <= n1; l++ ) /* Yes, do so... */
         {
            dum        = a[imax][l];
            a[imax][l] = a[k][l];
            a[k][l]    = dum;
         }
         vv[imax] = vv[k];    /* Also interchange the scale factor */
      }
      indx[k] = imax;
      if( fabs(a[k][k]) < 1.0e-20 )
      {
         out = fopen( "err.log", "a" );
         fprintf( out, "Pivot too small. n0=%ld k=%ld fabs(a[k][k])=%le\n",
                        n0,k,fabs(a[k][k]) );
         fclose( out );
         exit(0);   
/*         a[k][k] = 1.e-20; */
      }
      if( k != n1 )                   /* If the pivot element is less than 1.0d-20 we */
      {                               /* assume that the matrix is singular */
          dum = a[k][k];              /* ( at least to the precision of the algorithm and the machine ) */
         for( i = k+1; i <= n1; i++ ) /* Now, finally devide by the pivot element */
            a[i][k] /= dum;
      }
   }                                  /* Go back to the next column in the reduction */

   for( i = n0; i <= n1; i++ ) vv[i]=f[i][m];

      ii = n0-1;
      for( i = n0; i <= n1; i++ )
      {
         ll    = indx[i];
         sum   = vv[ll];
         vv[ll] = vv[i];
         if( ii != (n0-1) )
            for( k = ii; k < i; k++ )
               sum -= ( a[i][k] * vv[k] );
         else if( sum != 0.0e0 ) ii = i;
         vv[i] = sum;
      }
      for( i = n1; i >= n0; i-- )
      {
         sum = vv[i];
         if( i < n1 )
            for( k = i+1; k <= n1; k++ )
               sum -= ( a[i][k] * vv[k] );
         vv[i] = sum / a[i][i];
      }
         
   for( i = n0; i <= n1; i++ ) f[i][m]=vv[i];

   fdm1d( indx, n0 );
   fdm1d(   vv, n0 );
   fdm2d(    a, n0, n1, n0 );
}
