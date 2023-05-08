/* disso: photodissociation rates */
/* GCCM */
/* correspond a chimie_simpnit (version 301105) */
/* !!! ATTENTION !!! */
/* Doit etre mis a jour en fonction de la chimie utilisee ! */

#include "titan.h"

void disso_( double KRPD[][NLRT][RDISS+1][15], int *NLAT )
{
   static double sH2[62] = {  /* incertain en dessous de 70 et en dessus de 85... */
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,1.000e-18,5.000e-18,1.000e-17,
          9.000e-18,6.500e-18,1.000e-18,1.000e-19 };
   static double sCH4[62] = {
          2.852e-19,7.816e-19,1.534e-18,2.069e-18,2.795e-18,4.088e-18,4.543e-18,
          4.223e-18,3.314e-18,1.565e-18,8.892e-19,8.760e-19,8.792e-19,9.163e-19,
          2.069e-18,9.378e-18,2.543e-17,3.785e-17,4.066e-17,3.302e-17,2.840e-17,
          1.800e-17,1.920e-17,1.820e-17,1.840e-17,1.140e-17,2.656e-18,1.256e-19,
          7.988e-22,1.366e-23,6.740e-24 };
   static double sCH3CN[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,8.000e-17,8.000e-17,7.500e-17,
          4.800e-17,3.700e-17,2.700e-17,4.100e-17,1.000e-17,8.600e-18,4.000e-18,
          1.800e-18,1.100e-18,6.400e-19,3.600e-19,2.000e-19,1.200e-19,5.200e-20};
   static double sC2H2[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,3.000e-17,
          3.000e-17,3.200e-17,3.640e-17,4.260e-17,1.040e-16,9.900e-18,4.800e-18,
          1.720e-17,1.782e-17,9.040e-19,9.600e-19,1.294e-18,1.352e-18,1.130e-18,
          6.680e-19,3.700e-19,3.900e-19,1.660e-19,5.400e-20 };
   static double sC2H4[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,.1992E-16,.2734E-16,.2594E-16,.1690E-16,.2258E-16,.8507E-17,
          .1583E-16,.2227E-16,.3056E-16,.3743E-16,.3788E-16,.2737E-16,.3171E-17,
          .6033E-18,.1223E-18,.7247E-19,1.000e-20 };
   static double sC2H3CN[62] = { 
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e-00,0.000e-00,0.000e-00,
          0.000e-00,0.000e-00,1.000e-17,6.600e-18,5.600e-18,3.600e-18,4.300e-18,
          4.000e-18,3.000e-18,2.900e-18,2.700e-18,2.700e-18,3.300e-18,4.500e-18,
          6.000e-18,7.100e-18,7.000e-18,5.000e-18,2.500e-18,6.600e-19,1.000e-19,
          1.500e-20,1.000e-21};
   static double sC2H6[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,2.000e-17,
          2.000e-17,2.000e-17,2.060e-17,2.160e-17,1.540e-17,8.060e-18,3.860e-18,
          1.484e-18,3.060e-19,9.600e-21 };
   static double sCH3C2H[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,1.000e-17,
          2.000e-17,4.580e-17,6.240e-17,5.880e-17,5.920e-17,1.940e-17,2.200e-17,
          2.820e-17,3.380e-17,1.142e-17,6.600e-18,8.100e-18,7.000e-18,2.800e-18,
          1.600e-18,2.300e-19,4.411e-19,2.119e-19,1.004e-19,2.934e-20,4.157e-21 };
   static double sCH2CCH2[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,2.200e-17,2.200e-17,1.700e-17,1.700e-17,1.800e-17,
          1.200e-17,1.500e-17,1.000e-17,2.100e-17,3.300e-17,4.000e-17,1.300e-17,
          5.000e-18,2.900e-18,2.601e-18,1.037e-18,9.046e-19,6.565e-19,4.672e-19,
          3.047e-19,1.579e-19,5.943e-20,2.261e-20 };
   static double sC3H6[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,3.400e-17,
          3.300e-17,2.400e-17,2.500e-17,4.000e-17,3.700e-17,2.300e-17,1.900e-17,
          2.000e-17,1.500e-17,2.200e-17,2.500e-17,4.400e-17,4.200e-17,2.700e-17,
          1.200e-17,5.800e-18,1.400e-19,9.100e-21 };
   static double sC3H8[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,4.000e-17,4.000e-17,4.000e-17,
          4.000e-17,4.000e-17,3.280e-17,3.100e-17,2.680e-17,2.200e-17,1.760e-17,
          6.440e-18,3.000e-18,9.140e-19,7.000e-20 };
   static double sC4H2[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,2.500e-17,6.100e-17,4.780e-17,5.780e-17,1.014e-16,
          5.200e-17,4.200e-17,1.030e-16,1.770e-16,9.100e-17,1.290e-17,2.380e-18,
          3.400e-19,2.800e-19,2.600e-19,1.513e-19,2.583e-19,3.353e-19,4.115e-19,
          4.755e-19,4.990e-19,4.399e-19,5.358e-19,2.485e-19,3.931e-19,1.067e-19,
          3.761e-20,2.370e-20,2.277e-20 };
   static double sC4H4[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,4.200e-17,5.000e-17,2.000e-17,6.700e-18,
          6.100e-18,8.000e-18,1.000e-17,1.400e-17,1.900e-17,2.200e-17,2.600e-17,
          1.600e-17,1.600e-17,4.000e-18,5.700e-19,1.400e-19 };
   static double sC4H6[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,6.000e-18,1.200e-17,1.500e-17,1.300e-17,
          1.800e-17,3.100e-17,3.800e-17,6.600e-17,9.600e-17,1.060e-16,8.500e-17,
          6.700e-17,1.100e-17,1.700e-18,3.300e-19,6.500e-20 };
   static double sC4H10[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,6.000e-17,6.000e-17,6.000e-17,
          6.000e-17,6.000e-17,5.500e-17,4.400e-17,4.400e-17,3.800e-17,3.100e-17,
          1.900e-17,4.000e-18,1.300e-18,3.200e-19,2.000e-20 };
   static double sC6H6[62] = { 
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,7.000e-17,
          8.000e-17,5.000e-17,3.500e-17,3.500e-17,3.500e-17,2.000e-17,1.500e-17,
          1.500e-17,1.500e-17,2.000e-17,2.500e-17,4.000e-17,9.500e-17,2.200e-16,
          1.000e-16,2.000e-17,2.000e-17,2.000e-17,2.000e-17,5.000e-18,1.000e-20,
          1.000e-20,1.000e-20,1.000e-20,1.000e-20,1.000e-19,2.000e-19,2.500e-19,
          4.000e-19,2.000e-19,2.000e-19,1.000e-19 };
   static double sN2[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,4.898e-18,1.097e-17,
          2.192e-17,2.214e-17,2.336e-17,1.679e-17,1.893e-17 };
   static double sHCN[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,2.800e-17,3.300e-17,2.800e-17,3.500e-17,
          4.800e-17,2.800e-17,1.700e-17,3.600e-18,2.500e-18,4.700e-18,2.700e-18,
          8.600e-19,4.400e-19,2.000e-19,1.429e-19,1.145e-19,7.482e-20,3.852e-20,
          1.009e-20 };
   static double sHC3N[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,3.402e-17,
          3.647e-17,4.317e-17,3.759e-17,7.927e-17,3.796e-17,9.565e-17,1.716e-16,
          1.247e-16,2.360e-17,8.411e-18,4.400e-18,8.600e-19,7.400e-19,6.200e-19,
          4.899e-19,3.307e-19,2.128e-19,2.561e-19,2.621e-19,2.737e-19,3.601e-19,
          1.564e-19,1.816e-19,8.427e-20 };
   static double sC2N2[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,2.635e-17,
          3.126e-17,4.392e-17,9.788e-17,1.682e-16,1.064e-16,6.513e-18,2.039e-18,
          2.828e-18,5.136e-18,8.188e-18,8.857e-18,1.489e-18,9.000e-20,6.200e-20,
          3.800e-20,4.483e-20,8.618e-20,1.008e-19,7.579e-20,6.666e-20,2.907e-20,
          2.476e-20,1.142e-20 };
   static double sH2O[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,1.927e-17,9.633e-18,7.440e-18,
          1.087e-18,8.600e-18,8.850e-18,7.300e-18,4.100e-18,1.240e-18,5.600e-19,
          7.977e-19,1.888e-18,3.333e-18,4.729e-18,4.963e-18,3.513e-18,1.440e-18,
          1.563e-19 };
   static double sCO[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,7.300e-17,2.050e-16,9.000e-17,7.180e-18,3.700e-17};
   static double sCO2[62] = {
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
          0.000e+00,1.022e-19,1.070e-19,4.220e-19,8.220e-19,6.640e-19,5.650e-19,
          5.850e-19,4.400e-19,2.100e-19,8.700e-20,4.000e-20,2.130e-20,3.700e-21,
          6.930e-22,4.040e-22,2.470e-23,4.800e-24,9.500e-25 };
   static double sol[62] = {
          9.45e+08,2.84e+08,4.24e+09,3.11e+09,1.22e+10,6.83e+09,2.54e+09,8.20e+08,
          1.82e+09,1.30e+09,1.05e+09,1.37e+09,3.87e+08,4.29e+08,2.45e+09,4.09e+09,
          1.06e+10,1.03e+10,4.08e+09,2.98e+10,7.00e+09,1.94e+09,7.40e+09,4.12e+11,
          1.36e+10,4.92e+10,3.04e+10,3.81e+10,5.49e+10,1.00e+11,1.26e+11,1.67e+11,
          2.97e+11,4.92e+11,7.17e+11,9.43e+11,1.34e+12,1.96e+12,2.96e+12,4.34e+12,
          7.25e+12,1.76e+13,2.21e+13,3.10e+13,2.80e+13,2.97e+13,3.03e+13,3.77e+13,
          3.44e+13,3.42e+13,7.08e+13,9.58e+13,1.63e+14,1.38e+14,1.26e+14,1.60e+14,
          2.23e+14,3.70e+14,3.31e+14,3.26e+14,3.59e+14,4.05e+14};

   int   i,j,l,s,lat,x;
   double f,**flux;
   double flact;
   char   name[60],dir[24];
   FILE   *fp,*out;

   flux = dm2d(0,NLRT,0,14);

/* lecture des flux actiniques:
   - suppose que l'executable est dans $LMDGCM/RUN/xxx/
   - et les moyennes dans $LMDGCM/INPUT/PHOT(NLAT)/Moy_(lat:1 a NLAT)
*/
   strcpy( dir, "../../INPUT/PHOT" );
   if( (*NLAT) < 10 )
   {
     strcat( dir, "0"  );
     strcat( dir, (const char *)ecvt((float)(*NLAT),1,&x,&x) );
   }
   else 
     strcat( dir, (const char *)ecvt((float)(*NLAT),2,&x,&x) );
   strcat( dir, "x/Moy_"  );
   printf( "Directories for actinic fluxes: %s \n", dir );
 
   for( lat = 0; lat <= (*NLAT)-1; lat++ )   /* Old array is set equal to 0. */
     for( j = 0; j <= NLRT-1; j++ )
       for( i = 0; i <= RDISS; i++ )                  
         for( s = 0; s <= 14; s++ )
            KRPD[lat][j][i][s] = 0.0e0;
            
   for( i = 0; i <= 13; i++ ) sCH4[i]  = 0.0e0;
   for( i = 0; i <= 13; i++ ) sC2H2[i] = 0.0e0;
   for( i = 0; i <= 16; i++ ) sC2H4[i] = 0.0e0;
   for( i = 0; i <= 16; i++ ) sC2H6[i] = 0.0e0;

   for( lat = 0; lat <= (*NLAT)-1; lat++ ) /*     Main loop on latitude */
   for( i = 10; i <= 310; i += 5 )             /* Main loop on wavelength. */
   {
      strcpy( name, dir ); 
      if( (lat+1) < 10 )
      {
       strcat( name, "0"  );
       strcat( name, (const char *)ecvt((float)(lat+1),1,&x,&x) );
      }
      else 
       strcat( name, (const char *)ecvt((float)(lat+1),2,&x,&x) );
      if( i < 160 ) strcat( name, "/photmoy3a" );
      else          strcat( name, "/photmoy3a" );
      if( i < 100 )
      {
         strcat( name, ".0" );
         strcat( name, (const char *)ecvt((float)i,2,&x,&x) );
      }
      else
      {
         strcat( name, "." );
         strcat( name, (const char *)ecvt((float)i,3,&x,&x) );
      }
      if( !( fp = fopen( name, "r" ) ) )
      {
         out = fopen( "err.log", "a" );
         fprintf( out, "I cannot open %s\n", name );
         fclose( out );
         exit(0);
      }
      for( j = 0; j <= NLRT-1; j++ )
      {
         fscanf( fp,"%d ",&l );
         for( s = 0; s < 15; s++ )
         {
           fscanf( fp,"%lg ", &flact );
           flux[j][s] = flact;
         }
      }
      fclose(fp);

      l = i / 5 - 2;                  /* Pointer on wavelength. */
      
/* taux de photodissociations */

      for( s = 0; s <= 14; s++ )
       for( j = 0; j <= NLRT-1; j++ )
       {
         f = flux[j][s] * sol[l] / ( 9.5e0 * 9.5e0 );   /* !! # de reac de 0 a RDISS-1 !! */
         if( i == 220 ) KRPD[lat][j][1][s] += 4.4e-17 * f;  /* CH3 -> 1CH2 + H */
         KRPD[lat][j][ 0][s] += sH2[l]     * f * 1.00;      /* H2  -> H + H */
         KRPD[lat][j][ 7][s] += sC2H4[l]   * f * 0.51;      /* C2H4 -> C2H2 + H2 */
         KRPD[lat][j][ 8][s] += sC2H4[l]   * f * 0.49;      /* C2H4 -> C2H2 + 2H */
         KRPD[lat][j][15][s] += sCH2CCH2[l]* f * 0.89;      /* CH2CCH2 -> C3H3 + H    Jackson 91 */
         KRPD[lat][j][16][s] += sCH2CCH2[l]* f * 0.11;      /* CH2CCH2 -> C3H2 + H2   Jackson 91 */
         KRPD[lat][j][17][s] += sCH3C2H[l] * f * 0.89;      /* CH3C2H -> C3H3 + H     Jackson 91 */
         KRPD[lat][j][18][s] += sCH3C2H[l] * f * 0.11;      /* CH3C2H -> C3H2 + H2    Jackson 91 */
         KRPD[lat][j][19][s] += sC3H6[l]   * f * 0.33;      /* C3H6 -> CH2CCH2 + H2 */
         KRPD[lat][j][20][s] += sC3H6[l]   * f * 0.17;      /* C3H6 -> CH3CCH + H2 */
         KRPD[lat][j][21][s] += sC3H6[l]   * f * 0.03;      /* C3H6 -> C2H4 + 3CH2 */
         KRPD[lat][j][22][s] += sC3H6[l]   * f * 0.35;      /* C3H6 -> C2H3 + CH3 */
         KRPD[lat][j][23][s] += sC3H6[l]   * f * 0.05;      /* C3H6 -> C2H2 + CH4 */
         KRPD[lat][j][32][s] += sC4H4[l]   * f * 0.80;      /* C4H4 -> C4H2 + H2     Gladstone 96 */
         KRPD[lat][j][33][s] += sC4H4[l]   * f * 0.20;      /* C4H4 -> C2H2 + C2H2   Gladstone 96  */
         KRPD[lat][j][34][s] += sC4H6[l]   * f * 0.04;      /* C4H6 -> C4H4 + H2 */
         KRPD[lat][j][35][s] += sC4H6[l]   * f * 0.27;      /* C4H6 -> C2H4 + C2H2 */
         KRPD[lat][j][36][s] += sC4H6[l]   * f * 0.69;      /* C4H6 -> CH3 + C3H3 */
         KRPD[lat][j][45][s] += sC6H6[l]   * f * 0.04;      /* AC6H6 -> C5H3 (prod...) + CH3 */
         KRPD[lat][j][46][s] += sC6H6[l]   * f * 0.96;      /* AC6H6 -> AC6H5 + H */
         KRPD[lat][j][47][s] += sN2[l]     * f;             /* N2   -> 2N2d */
         KRPD[lat][j][48][s] += sHCN[l]    * f;             /* HCN  -> H + CN */
         KRPD[lat][j][51][s] += sC2N2[l]   * f * 0.3;       /* C2N2 -> 2CN */
         KRPD[lat][j][52][s] += sCH3CN[l]  * f * 1.0;       /* CH3CN -> CH3 + CN */
         KRPD[lat][j][53][s] += sC2N2[l]   * f * 0.3;       /* C4N2 -> C3N + CN */

         if( i != 125 )          /* Not Lyman alpha */
         {
            KRPD[lat][j][ 2][s] += sCH4[l]  * f;            /* CH4 -> 1CH2 + H2 */
            KRPD[lat][j][ 9][s] += sC2H6[l] * f * 0.56;     /* C2H6 -> C2H4 + H2 */
            KRPD[lat][j][10][s] += sC2H6[l] * f * 0.14;     /* C2H6 -> C2H4 + 2H */
            KRPD[lat][j][11][s] += sC2H6[l] * f * 0.27;     /* C2H6 -> C2H2 + 2H2 */
            KRPD[lat][j][12][s] += sC2H6[l] * f * 0.02;     /* C2H6 -> CH4 + 3CH2 */
            KRPD[lat][j][13][s] += sC2H6[l] * f * 0.01;     /* C2H6 -> 2CH3 */
            KRPD[lat][j][24][s] += sC3H8[l] * f * 0.94;     /* C3H8 -> C3H6 + H2 */
            KRPD[lat][j][27][s] += sC3H8[l] * f * 0.06;     /* C3H8 -> C2H4 + CH4 */
         }
         else                       /* Lyman alpha */
         {
            KRPD[lat][j][ 2][s] += sCH4[l]  * f * 0.64;     /* CH4 -> 1CH2 + H2 */
            KRPD[lat][j][ 3][s] += sCH4[l]  * f * 0.07;     /* CH4 -> CH + H2 + H */
            KRPD[lat][j][ 4][s] += sCH4[l]  * f * 0.29;     /* CH4 -> CH3 + H */
            KRPD[lat][j][ 9][s] += sC2H6[l] * f * 0.13;     /* C2H6 -> C2H4 + H2 */
            KRPD[lat][j][10][s] += sC2H6[l] * f * 0.3;      /* C2H6 -> C2H4 + 2H */
            KRPD[lat][j][11][s] += sC2H6[l] * f * 0.25;     /* C2H6 -> C2H2 + 2H2 */
            KRPD[lat][j][12][s] += sC2H6[l] * f * 0.25;     /* C2H6 -> CH4 + 3CH2 */
            KRPD[lat][j][13][s] += sC2H6[l] * f * 0.08;     /* C2H6 -> 2CH3 */
            KRPD[lat][j][24][s] += sC3H8[l] * f * 0.33;     /* C3H8 -> C3H6 + H2 */
            KRPD[lat][j][25][s] += sC3H8[l] * f * 0.09;     /* C3H8 -> C2H6 + 3CH2 */
            KRPD[lat][j][26][s] += sC3H8[l] * f * 0.39;     /* C3H8 -> C2H5 + CH3 */
            KRPD[lat][j][27][s] += sC3H8[l] * f * 0.2;      /* C3H8 -> C2H4 + CH4 */
         }
         if( i < 145 )   /* C4H10: a revoir avec Jackson & Lias, 1974... */
         {
            KRPD[lat][j][37][s] += sC4H10[l]* f * 0.18;      /* C4H10 -> C4H8(ieC3H5+CH3)+H2 */
            KRPD[lat][j][38][s] += sC4H10[l]* f * 0.20;      /* C4H10 -> 2 C2H4 + H2 */
            KRPD[lat][j][39][s] += sC4H10[l]* f * 0.03;      /* C4H10 -> C3H6 + CH4 */
            KRPD[lat][j][40][s] += sC4H10[l]* f * 0.07;      /* C4H10 -> C3H6 + CH3 + H */
            KRPD[lat][j][41][s] += sC4H10[l]* f * 0.00;      /* C4H10 -> C2H6 + C2H4 */
            KRPD[lat][j][42][s] += sC4H10[l]* f * 0.15;      /* C4H10 -> C2H6 + C2H2 + H2 */
            KRPD[lat][j][43][s] += sC4H10[l]* f * 0.27;      /* C4H10 -> CH3 + C3H7 */
            KRPD[lat][j][44][s] += sC4H10[l]* f * 0.10;      /* C4H10 -> 2 C2H5 */
         }
         else
         {
            KRPD[lat][j][37][s] += sC4H10[l]* f * 0.41;      /* C4H10 -> C4H8(ieC3H5+CH3)+H2 */
            KRPD[lat][j][38][s] += sC4H10[l]* f * 0.12;      /* C4H10 -> 2 C2H4 + H2 */
            KRPD[lat][j][39][s] += sC4H10[l]* f * 0.01;      /* C4H10 -> C3H6 + CH4 */
            KRPD[lat][j][40][s] += sC4H10[l]* f * 0.07;      /* C4H10 -> C3H6 + CH3 + H */
            KRPD[lat][j][41][s] += sC4H10[l]* f * 0.02;      /* C4H10 -> C2H6 + C2H4 */
            KRPD[lat][j][42][s] += sC4H10[l]* f * 0.06;      /* C4H10 -> C2H6 + C2H2 + H2 */
            KRPD[lat][j][43][s] += sC4H10[l]* f * 0.24;      /* C4H10 -> CH3 + C3H7 */
            KRPD[lat][j][44][s] += sC4H10[l]* f * 0.07;      /* C4H10 -> 2 C2H5 */
         }
         if( i < 150 )
         {
            KRPD[lat][j][ 5][s] += sC2H2[l] * f * 0.3;      /* C2H2 -> C2H + H */
            KRPD[lat][j][ 6][s] += sC2H2[l] * f * 0.1;      /* C2H2 -> C2 + H2 */
            KRPD[lat][j][49][s] += sHC3N[l] * f * 0.3;      /* HC3N -> C2H + CN */
            KRPD[lat][j][50][s] += sHC3N[l] * f * 0.09;     /* HC3N -> H + C3N */
         }
         else if( i < 205 )
         {
            KRPD[lat][j][ 5][s] += sC2H2[l] * f * 0.08;     /* C2H2 -> C2H + H */
            KRPD[lat][j][ 6][s] += sC2H2[l] * f * 0.1;      /* C2H2 -> C2 + H2 */
            KRPD[lat][j][49][s] += sHC3N[l] * f * 0.05;     /* HC3N -> C2H + CN */
            KRPD[lat][j][50][s] += sHC3N[l] * f * 0.09;     /* HC3N -> H + C3N */
         }
         else if( i < 245 )
         {
            KRPD[lat][j][50][s] += sHC3N[l] * f * 0.09;     /* HC3N -> H + C3N */
         }
         if( i < 165 )
         {
            KRPD[lat][j][28][s] += sC4H2[l] * f * 0.2;      /* C4H2 -> C4H + H */
            KRPD[lat][j][29][s] += sC4H2[l] * f * 0.03;     /* C4H2 -> 2C2H */
            KRPD[lat][j][30][s] += sC4H2[l] * f * 0.1;      /* C4H2 -> C2H2 + C2 */
            KRPD[lat][j][31][s] += sC4H2[l] * f * 0.67;     /* C4H2 -> C4H2* */
         }
         else if( i < 205 )
         {
            KRPD[lat][j][29][s] += sC4H2[l] * f * 0.01;     /* C4H2 -> 2C2H */
            KRPD[lat][j][30][s] += sC4H2[l] * f * 0.06;     /* C4H2 -> C2H2 + C2 */
            KRPD[lat][j][31][s] += sC4H2[l] * f * 0.93;     /* C4H2 -> C4H2* */
         }
         else
         {
            KRPD[lat][j][31][s] += sC4H2[l] * f * 1.00;     /* C4H2 -> C4H2* */
         }
         if( i < 190 ) KRPD[lat][j][14][s] += 4.e-17 * f;   /* C3H3 -> C3H2 + H */
      }
   }

/* taux de dissociation de N2 ( e- et GCR ) */

   for( lat = 0; lat <= (*NLAT)-1; lat++ ) 
    for( s = 0; s <= 14; s++ )
    {
     for( j = 99; j <= NLRT-1; j++ )  /* level 100 = 200 km */
      KRPD[lat][j][RDISS][s] = 1.0e-16;
     for( j = 49; j <= 98; j++ )      /* level 50 = 100 km */
      KRPD[lat][j][RDISS][s] = 1.0e-17+1.8e-18*(j-49);
     for( j = 34; j <= 48; j++ )      /* level 35 = 70 km */
      KRPD[lat][j][RDISS][s] = pow(10.,(-23+0.4*(j-34)));
     for( j = 0; j <= 33; j++ )
      KRPD[lat][j][RDISS][s] = 0.0e0;
    }

   fdm2d( flux, 0, NLEV, 0 );
}