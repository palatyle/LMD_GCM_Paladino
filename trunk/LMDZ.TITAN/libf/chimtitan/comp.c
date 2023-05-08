/* comp: Compounds characteristics. */
/* GCCM */

#include "titan.h"

void comp_(char CORPS[][10], double *CT, double *TEMP, 
           double *MASS, double MD[][NLEV])
{
   int   i,j;

   double m,ma,epsa,sig,siga,p;

   p         = 2.976e07;   /*9 10^7 R / 8 pi */

      /* WARNING BACKGROUND GAS IS N2 */

   ma   = 28.0134e0;               /* mass of background gas in g */
   siga = 3.798e0;         /* Lennard-Jones length of background gas 1/10 nm */
   epsa = 71.4e0;       /* Lennard-Jones energy of background gas */

   for( i = 0; i <= NC-1; i++)
   {
      for( j = 0; j <= NLEV-1; j++ ) MD[i][j] = 0.0e0;
   }

   for( i = 0; i <= NC-1; i++ )
   {
      if( strcmp(CORPS[i], "CH4") == 0 )
      {
         MASS[i] = 16.04e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.758e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,148.6e0) );
      }
      if( strcmp(CORPS[i], "H") == 0 )
      {
         MASS[i] = 1.01e0;
      }
      if( strcmp(CORPS[i], "H2") == 0 )
      {
         MASS[i] = 2.0158e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 2.827e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,59.7e0) );
      }
      if( strcmp(CORPS[i], "CH") == 0 )
      {
         MASS[i] = 13.02e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( ( strcmp( CORPS[i], "CH2" ) == 0 ) || ( strcmp( CORPS[i], "CH2s" ) == 0 ) )
      {
         MASS[i] = 14.03e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "CH3") == 0 )
      {
         MASS[i] = 15.03e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.7e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C") == 0 )
      {
         MASS[i] = 12.01e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 1.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C2") == 0 )
      {
         MASS[i] = 24.02e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.2e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C2H") == 0 )
      {
         MASS[i] = 25.03e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.5e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C2H3") == 0 )
      {
         MASS[i] = 27.05e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.8e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C2H4") == 0 )
      {
         MASS[i] = 28.05e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.163e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,224.7e0) );
      }
      if( strcmp(CORPS[i], "C2H2") == 0 )
      {
         MASS[i] = 26.04e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.033e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,231.8e0) );
      }
      if( strcmp(CORPS[i], "C2H5") == 0 )
      {
         MASS[i] = 29.06e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C2H6") == 0 )
      {
         MASS[i] = 30.07e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.443e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,215.7e0) );
      }
      if( strcmp(CORPS[i], "C3H2") == 0 )
      {
         MASS[i] = 38.05e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.6e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C3H3") == 0 )
      {
         MASS[i] = 39.06e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.7e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( ( strcmp(CORPS[i], "CH2CCH2") == 0 ) || ( strcmp(CORPS[i], "CH3CCH") == 0 ) )
      {
         MASS[i] = 40.07e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.761e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,251.8e0) );
      }
      if( strcmp(CORPS[i], "C3H5") == 0 )
      {
         MASS[i] = 41.07e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.78e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C3H6") == 0 )
      {
         MASS[i] = 42.08e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.807e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,248.9e0) );
      }
      if( strcmp(CORPS[i], "C3H7") == 0 )
      {
         MASS[i] = 43.09e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 5.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
       }
      if( strcmp(CORPS[i], "C3H8") == 0 )
      {
         MASS[i] = 44.11e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 5.118e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,237.1e0) );
      }
      if( strcmp(CORPS[i], "C4H") == 0 )
      {
         MASS[i] = 49.05e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.2e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( ( strcmp(CORPS[i], "C4H2") == 0 )||( strcmp(CORPS[i], "C4H2s") == 0 ) )
      {
         MASS[i] = 50.06e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.3e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C4H3") == 0 )
      {
         MASS[i] = 51.07e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C4H4") == 0 )
      {
         MASS[i] = 52.08e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.5e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C4H5") == 0 )
      {
         MASS[i] = 53.07e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.5e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C4H6") == 0 )
      {
         MASS[i] = 54.09e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.6e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C4H10") == 0 )
      {
         MASS[i] = 58.13e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.687e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,531.4e0) );
      }
      if( strcmp(CORPS[i], "C6H") == 0 )
      {
         MASS[i] = 73.07e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 5.2e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C6H2") == 0 )
      {
         MASS[i] = 74.08e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 5.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C8H2") == 0 )
      {
         MASS[i] = 98.10e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 6.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp( CORPS[i], "AC6H6" ) == 0 )
      {
         MASS[i] = 78.1136e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 5.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )   /* P. G. L. */
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( ( strcmp( CORPS[i], "C6H5" ) == 0 ) || ( strcmp( CORPS[i], "AC6H5" ) == 0 ) )
      {
         MASS[i] = 77.1136e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 5.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp( CORPS[i], "C6H6" ) == 0 )
      {
         MASS[i] = 78.1136e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 5.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "N2") == 0 )
      {
         MASS[i] = 28.0134e0;
      }
      if( strcmp(CORPS[i], "N4S") == 0 )
      {
         MASS[i] = 14.01e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 1.5e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "NH") == 0 )
      {
         MASS[i] = 15.01e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "CN") == 0 )
      {
         MASS[i] = 26.02e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.2e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "HCN") == 0 )
      {
         MASS[i] = 27.04e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.63e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,569.1e0) );
      }
      if( strcmp(CORPS[i], "H2CN") == 0 )
      {
         MASS[i] = 28.05e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.8e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "C2N") == 0 )         /* C2N */
      {
         MASS[i] = 39.05e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp( CORPS[i], "CHCN" ) == 0 )
      {
         MASS[i]   = 39.05e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp( CORPS[i], "CH2CN" ) == 0 )
      {
         MASS[i]   = 40.04e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp( CORPS[i], "CH3CN" ) == 0 )
      {
         MASS[i]   = 41.05e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp( CORPS[i], "C2H3CN" ) == 0 )
      {
         MASS[i]   = 53.06e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "NCCN") == 0 )        /* NCCN */
      {
         MASS[i] = 52.04e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.361e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,348.6e0) );
      }
      if( strcmp(CORPS[i], "C3N") == 0 )         /* C3N */
      {
         MASS[i] = 50.04e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "HC3N") == 0 )        /* HC3N */
      {
         MASS[i] = 51.05e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.5e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp( CORPS[i], "C4N2" ) == 0 )
      {
         MASS[i]   = 76.1e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "H2O") == 0 )       
      {
         MASS[i] = 18.02e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 2.641e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) 
           / ( CT[j] * sig * omega(TEMP[j],epsa,809.1e0) );
      }
      if( ( strcmp(CORPS[i], "O3P") == 0 ) || ( strcmp(CORPS[i], "O1D") == 0 ) )      
      {
         MASS[i] = 16.0e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 1.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "OH") == 0 )        
      {
         MASS[i] = 17.01e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.0e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "HO2") == 0 )
      {
         MASS[i] = 33.01e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.5e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "H2O2") == 0 )
      {
         MASS[i] = 33.01e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.5e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "O2") == 0 )
      {
         MASS[i] = 32.0e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.7e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "O3") == 0 )
      {
         MASS[i] = 32.0e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.7e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "CO") == 0 )          
      {
         MASS[i] = 28.01e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.69e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,91.7e0) );
      }
      if( strcmp(CORPS[i], "HCO") == 0 )        
      {
         MASS[i] = 29.02e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.7e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "CO2") == 0 )       
      {
         MASS[i] = 44.01e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.941e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) 
           / ( CT[j] * sig * omega(TEMP[j],epsa,195.2e0) );
      }
      if( strcmp(CORPS[i], "CH2CO") == 0 )    
      {
         MASS[i] = 42.04e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 4.5e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "CH2O") == 0 )    
      {
         MASS[i] = 30.03e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.75e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( ( strcmp(CORPS[i], "CH2OH") == 0 ) || ( strcmp(CORPS[i], "CH3O") == 0 ) )  
      {
         MASS[i] = 31.04e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.4e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m ) / ( CT[j] * sig );
      }
      if( strcmp(CORPS[i], "CH3OH") == 0 )       
      {
         MASS[i] = 32.042e0;
         m       = ( ma + MASS[i] ) / ( ma * MASS[i] );
         sig     = 1.0e-16 * pow( ( siga + 3.626e0 ), 2.0e0 );
         for( j = 0; j <= NLEV-1; j++ )
            MD[i][j] = sqrt( p * TEMP[j] * m )
           / ( CT[j] * sig * omega(TEMP[j],epsa,481.8e0) );
      }
   }
}
