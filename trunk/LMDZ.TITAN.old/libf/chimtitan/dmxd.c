#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

                /* dm1d: memory allocation for a vector of double */
                /* ---------------------------------------------- */
                
double *dm1d( int nl, int nh )
{
   double *m;
   FILE  *out;

   m = (double *)calloc( (unsigned)(nh-nl+1), sizeof(double) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in dm1d" );
      fclose( out );
      exit(0);
   }
   return m-nl;
}

                /* dm2d: memory allocation for a matrix of double */
                /* ---------------------------------------------- */

double **dm2d( int nrl, int nrh, int ncl, int nch )
{
   int   i;
   double **m;
   FILE  *out;

   m = (double **)calloc( (unsigned)(nrh-nrl+1), sizeof(double *) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in dm2d" );
      fclose( out );
      exit(0);
   }
   m -= nrl;
   for( i = nrl; i <= nrh; i++ )
   {
      m[i] = (double *)calloc( (unsigned)(nch-ncl+1), sizeof(double) );
      if(!m[i])
      {
         out = fopen( "err.log", "a" );
         fprintf( out, "Memory allocation error in dm2d" );
         fclose( out );
         exit(0);
      }
      m[i] -= ncl;
   }
   return m;
}

                /* dm3d: memory allocation for a 3D vector of double */
                /* ------------------------------------------------- */

double ***dm3d( int nrl, int nrh, int ncl, int nch, int nal, int nah )
{
   int   i;
   double ***m;
   FILE  *out;

   m = (double ***)calloc( (unsigned)(nrh-nrl+1), sizeof(double **) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in dm3d" );
      fclose( out );
      exit(0);
   }
   m -= nrl;
   for( i = nrl; i <= nrh; i++ )
      m[i] = dm2d( ncl, nch, nal, nah );
   return m;
}

                /* dm4d: memory allocation for a 4D vector of double */
                /* ------------------------------------------------- */

double ****dm4d( int nrl, int nrh, int ncl, int nch, int nal, int nah, int nll, int nlh )
{
   int   i;
   double ****m;
   FILE  *out;

   m = (double ****)calloc( (unsigned)(nrh-nrl+1), sizeof(double ***) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in dm4d" );
      fclose( out );
      exit(0);
   }
   m -= nrl;
   for( i = nrl; i <= nrh; i++ )
      m[i] = dm3d( ncl, nch, nal, nah, nll, nlh );
   return m;
}

                     /* fdm1d: release a vector of double */
                     /* --------------------------------- */

void fdm1d(v,nl)
double *v;
int nl;
{
   free((char *)(v+nl));
}

                     /* fdm2d: release a matrix of double */
                     /* --------------------------------- */

void fdm2d(m,nrl,nrh,ncl)
double **m;
int ncl,nrh,nrl;
{
   int i;

   for( i = nrh; i >= nrl; i-- ) free((char *)(m[i]+ncl));
   free((char *)(m+nrl));
}

                     /* fdm3d: release a 3D vector of double */
                     /* ------------------------------------ */

void fdm3d(m,nrl,nrh,ncl,nch,nal)
double ***m;
int nch,ncl,nrh,nrl,nal;
{
   int i;

   for( i = nrh; i >= nrl; i-- ) fdm2d(m[i],ncl,nch,nal);
   free((char *)(m+nrl));
}

                     /* frm1d: release a vector of float */
                     /* -------------------------------- */

void frm1d(v,nl)
float *v;
int nl;
{
   free((char *)(v+nl));
}

                     /* frm2d: release a matrix of float */
                     /* -------------------------------- */

void frm2d(m,nrl,nrh,ncl)
float **m;
int ncl,nrh,nrl;
{
   int i;

   for( i = nrh; i >= nrl; i-- ) free((char *)(m[i]+ncl));
   free((char *)(m+nrl));
}

                     /* frm3d: release a 3D vector of float */
                     /* ----------------------------------- */

void frm3d(m,nrl,nrh,ncl,nch,nal)
float ***m;
int nch,ncl,nrh,nrl,nal;
{
   int i;

   for( i = nrh; i >= nrl; i-- ) frm2d(m[i],ncl,nch,nal);
   free((char *)(m+nrl));
}

                /* rm1d: memory allocation for a vector of float */
                /* --------------------------------------------- */
                
float *rm1d( int nl, int nh )
{
   float *m;
   FILE  *out;

   m = (float *)calloc( (unsigned)(nh-nl+1), sizeof(float) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in rm1d" );
      fclose( out );
      exit(0);
   }
   return m-nl;
}

                /* rm2d: memory allocation for a matrix of float */
                /* --------------------------------------------- */

float **rm2d( int nrl, int nrh, int ncl, int nch )
{
   int   i;
   float **m;
   FILE  *out;

   m = (float **)calloc( (unsigned)(nrh-nrl+1), sizeof(float *) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in rm2d" );
      fclose( out );
      exit(0);
   }
   m -= nrl;
   for( i = nrl; i <= nrh; i++ )
   {
      m[i] = (float *)calloc( (unsigned)(nch-ncl+1), sizeof(float) );
      if(!m[i])
      {
         out = fopen( "err.log", "a" );
         fprintf( out, "Memory allocation error in rm2d" );
         fclose( out );
         exit(0);
      }
      m[i] -= ncl;
   }
   return m;
}

                /* rm3d: memory allocation for a 3D vector of float */
                /* ------------------------------------------------ */

float ***rm3d( int nrl, int nrh, int ncl, int nch, int nal, int nah )
{
   int   i;
   float ***m;
   FILE  *out;

   m = (float ***)calloc( (unsigned)(nrh-nrl+1), sizeof(float **) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in rm3d" );
      fclose( out );
      exit(0);
   }
   m -= nrl;
   for( i = nrl; i <= nrh; i++ )
      m[i] = rm2d( ncl, nch, nal, nah );
   return m;
}

                /* rm4d: memory allocation for a 4D vector of float */
                /* ------------------------------------------------ */

float ****rm4d( int nrl, int nrh, int ncl, int nch, int nal, int nah, int nll, int nlh )
{
   int   i;
   float ****m;
   FILE  *out;

   m = (float ****)calloc( (unsigned)(nrh-nrl+1), sizeof(float ***) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in rm4d" );
      fclose( out );
      exit(0);
   }
   m -= nrl;
   for( i = nrl; i <= nrh; i++ )
      m[i] = rm3d( ncl, nch, nal, nah, nll, nlh );
   return m;
}

                /* im1d: memory allocation for a vector of integer */
                /* ----------------------------------------------- */

int *im1d( int nl, int nh )
{
   int *m;
   FILE *out;

   m = (int *)calloc( (unsigned)(nh-nl+1), sizeof(int ) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in im1d" );
      fclose( out );
      exit(0);
   }
   return m-nl;
}

                /* im2d: memory allocation for a matrix of integer */
                /* ----------------------------------------------- */

int **im2d( int nrl, int nrh, int ncl, int nch )
{
   int i,**m;
   FILE *out;

   m = (int **)calloc( (unsigned)(nrh-nrl+1), sizeof(int *) );
   if( !m )
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in im2d" );
      fclose( out );
      exit(0);
   }
   m -= nrl;
   for( i = nrl; i <= nrh; i++ )
   {
      m[i] = (int *)calloc( (unsigned)(nch-ncl+1), sizeof(int) );
      if(!m[i])
      {
         out = fopen( "err.log", "a" );
         fprintf( out, "Memory allocation error in im2d" );
         fclose( out );
         exit(0);
      }
      m[i] -= ncl;
   }
   return m;
}

                /* im3d: memory allocation for a 3D vector of integer */
                /* -------------------------------------------------- */

int ***im3d( int nrl, int nrh, int ncl, int nch, int nal, int nah )
{
   int i,***m;
   FILE *out;

   m = (int ***)calloc( (unsigned)(nrh-nrl+1), sizeof(int **) );
   if(!m)
   {
      out = fopen( "err.log", "a" );
      fprintf( out, "Memory allocation error in im3d" );
      fclose( out );
      exit(0);
   }
   m -= nrl;
   for( i = nrl; i <= nrh; i++ )
      m[i] = im2d( ncl, nch, nal, nah );
   return m;
}
