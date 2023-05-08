/* titan.h: parameters for gptitan.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define R0    (double)(2575.0) /* Titan's radius */
#define NLEV  (int)(133)  /* Nbre de niv verticaux -> Need to be coherent with the vertical grid used !! */
#define NLD   (int)(40)   /* Nbre de niv verticaux faits sans diff -> Need to be coherent with the vertical grid used !! */
#define NLRT  (int)(650)  /* Nbre de niv verticaux dans table fmoy - aussi dans common_mod */

/* fluxes at 1300 km : upward is +, downward is - */
#define top_H   (double)(+1.1e4)
#define top_H2  (double)(+3.7e3)
#define top_N4S (double)(-1.1e8) /* = -2.5e8/2.27 ... */

/* DEPEND DE LA VERSION CHIMIE: */
#define VERCHIM "chimie_simpnit_051006_bis"
#define NREAC (int)(377)    /* nombre de reactions - aussi dans common_mod */
#define RDISS (int)(54)     /* nombre de photodiss - aussi dans common_mod */
#define NC    (int)(44)     /* nb de composes      - aussi dans common_mod */
#define ST    (int)(NC)     /* nb de composes inverses */
#define NHC   (int)(32)     /* nb hydrocarbons */

#define THETA (double)(0.501)
#ifndef M_PI
#define M_PI  (double)(3.14159265358979323846e0)
#endif
#define RAD   (double)(M_PI / 180.0e0)
#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<=(b)?(a):(b))
#endif

void  chimie_(char (*)[10], double *, double *, double (*)[NLEV], 
              int (*)[5], int *, int *, int (*)[200][2], int (*)[200]);
void  comp_(char (*)[10], double *, double *, double *, double (*)[NLEV]);
void  disso_(double (*)[15][NLRT][RDISS+1], int *);
double omega( double, double, double );
void  solve( double ***, int, int, int );
void  solve_lapack( double ***, int, int, int );
void  solve_b( double ***, double **, int, int, int );
float *rm1d( int, int );
float **rm2d( int, int, int, int );
float ***rm3d( int, int, int, int, int, int );
float ****rm4d( int, int, int, int, int, int, int, int );
double *dm1d( int, int );
double **dm2d( int, int, int, int );
double ***dm3d( int, int, int, int, int, int );
double ****dm4d( int, int, int, int, int, int, int, int );
void  frm1d( float *, int );
void  frm2d( float **, int, int, int );
void  frm3d( float ***, int, int, int, int, int );
void  fdm1d( double *, int );
void  fdm2d( double **, int, int, int );
void  fdm3d( double ***, int, int, int, int, int );
int   *im1d( int, int );
int   **im2d( int, int, int, int );
int   ***im3d( int, int, int, int, int, int );
