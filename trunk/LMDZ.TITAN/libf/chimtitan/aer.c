/* aer: production of aerosol precursors by photochemistry */

#include "titan.h"

void aer( char corps[][10], double *tp, double *nb, double y[][NLEV],
          int *zj, double **k_dep, double **f,
          double *r1, double *r2, double *r3, double *r4, double *r5, double *r6,
          int *utilaer, double *m, double *prod, double *csn, double *csh )
{
  int   re,x,z;
  int   i,j,c2h2,c2,c2h,hc3n,c3n,c4h2,hcn,h2cn,c4h3,ac6h5,ac6h6;
  int   nccn,ch3cn,c2h3cn,NMAX;
  double v1,v2,v3,alpha,beta1,beta2,gamma,denom;
  double temp,ct,dy_c2h2,dy_hc3n,dy_hcn;
  double dy_nccn,dy_ch3cn,dy_c2h3cn;
  char  sor1[20],sor2[20],sor3[20],outlog[10];
  FILE  *fp, *out;

  NMAX = 20;
  z    = (*zj);
  temp = tp[z];
  ct   = nb[z];

/* debug
      out = fopen( "aer.log", "a" );
      fprintf( out, "" );
      fclose( out );
*/

/* composes interessants */
/* --------------------- */
/* !! decalage de 1 par rapport a calchim !! */

     c2h2 = utilaer[2];
       c2 = utilaer[3];
      c2h = utilaer[4];
     hc3n = utilaer[5];
      hcn = utilaer[6];
      c3n = utilaer[7];
     h2cn = utilaer[8];
     c4h2 = utilaer[9];
     c4h3 = utilaer[10];
    ac6h5 = utilaer[11];
    ac6h6 = utilaer[12];
     nccn = utilaer[13];
    ch3cn = utilaer[14];
   c2h3cn = utilaer[15];
 

/* si C2H3CN n'est pas pris en compte: attribue a CH3CN, mais reaction nulle */

/* vitesses de reactions */
/* --------------------- */

  /* v (et K) en cm3.s-1 */
  v1 = 2.e-16;
  v2 = 3.72e-13*exp(-1561./temp);
  v3 = 1.0e-12*exp(-900./temp);

  alpha = 6.;
  beta1 = 1.;
  beta2 = 1.;
  gamma = 1.;

  /* k en s-1 */
  k_dep[1][1] = v1  *ct*y[c2h2][z]*y[c4h3][z];         /* C2H2 + C4H3 */
  k_dep[2][1] = v1  *ct*y[hc3n][z]*y[c4h3][z] *alpha;  /* HC3N + C4H3 */
  for( i = 3; i <= 5; i++ ) k_dep[i][1] = 0.;

  k_dep[1][2] = v2  *ct*y[c2h2][z]*y[ac6h5][z];        /* C2H2 + AC6H5 */
  k_dep[2][2] = v2  *ct* y[hcn][z]*y[ac6h5][z]*beta1;  /*  HCN + AC6H5 */
  k_dep[3][2] = v2  *ct*y[hc3n][z]*y[ac6h5][z]*beta2;  /* HC3N + AC6H5 */
/*  for( i = 4; i <= 5; i++ ) k_dep[i][2] = 0.; */
  k_dep[4][2] = 5.2e-10 *ct* y[c2][z]*y[ac6h6][z];  /* C2  + AC6H6 */
  k_dep[5][2] = 8.2e-11 *ct*y[c2h][z]*y[ac6h6][z];  /* C2H + AC6H6 */

  k_dep[1][3] = v3  *ct* y[hcn][z]*y[h2cn][z];         /* HCNH +    HCN */
/* !! debug !!
  k_dep[2][3] = v3  *ct*y[hc3n][z]*y[h2cn][z]*gamma;      HCNH +   HC3N */
  k_dep[2][3] = v3  *ct*y[hc3n][z]*y[h2cn][z]*1.e-10;  /* HCNH +   HC3N */
  k_dep[3][3] = v3  *ct*y[nccn][z]*y[h2cn][z]*gamma;   /* HCNH +   NCCN */
  k_dep[4][3] = v3 *ct*y[ch3cn][z]*y[h2cn][z]*gamma;   /* HCNH +  CH3CN */
 if(c2h3cn!=ch3cn)
  k_dep[5][3] = v3*ct*y[c2h3cn][z]*y[h2cn][z]*gamma;   /* HCNH + C2H3CN */
 else 
  k_dep[5][3] = 0.;

/* Fractions de chaque compose dans les polymeres */
/* ---------------------------------------------- */

/* polyC2H2 */

   denom   = y[c2h2][z] + alpha*y[hc3n][z];
   f[1][1] = y[c2h2][z] / denom;       /* C2H2 */
   f[2][1] = alpha*y[hc3n][z] / denom; /* HC3N */
   for( i = 3; i <= 5; i++ ) f[i][1] = 0.;

/* PAHs */

   denom   = y[c2h2][z] + beta1*y[hcn][z] + beta2*y[hc3n][z];
   f[1][2] = y[c2h2][z] / denom;       /* C2H2 */
   f[2][2] = beta1*y[hcn][z] / denom;  /*  HCN */
   f[3][2] = beta2*y[hc3n][z] / denom; /* HC3N */
   for( i = 4; i <= 5; i++ ) f[i][2] = 0.;

/* polyHCN */

  if(c2h3cn!=ch3cn)
   denom   = y[hcn][z] + gamma*(y[hc3n][z]+y[nccn][z]+y[ch3cn][z]+y[c2h3cn][z]);
  else
   denom   = y[hcn][z] + gamma*(y[hc3n][z]+y[nccn][z]+y[ch3cn][z]);
   f[1][3] =          y[hcn][z] / denom;         /*    HCN */
   f[2][3] = gamma*  y[hc3n][z] / denom;         /*   HC3N */
   f[3][3] = gamma*  y[nccn][z] / denom;         /*   NCCN */
   f[4][3] = gamma* y[ch3cn][z] / denom;         /*  CH3CN */
  if(c2h3cn!=ch3cn)
   f[5][3] = gamma*y[c2h3cn][z] / denom;         /* C2H3CN */
  else
   f[5][3] = 0.;         /* C2H3CN */

/* Masse molaire et Rapports C/N et C/H */
/* Taux de production en masse          */
/* Taux de destruction molecules        */
/* ------------------------------------ */

/* polyC2H2 */

    m[1]    = NMAX*(f[1][1]*26 + f[2][1]*51);  /* g.mol-1 */

    prod[1] = 0.;
    for( i = 1; i <= 5; i++ ) prod[1] += k_dep[i][1]; /* s-1 */

    dy_c2h2 = prod[1] * NMAX * f[1][1];
    dy_hc3n = prod[1] * NMAX * f[2][1];

    if( f[2][1] != 0.0e0 ) csn[1] = (2*f[1][1] + 3*f[2][1]) / f[2][1];
    else                   csn[1] = 1.0e30;
                           csh[1] = (2*f[1][1] + 3*f[2][1]) / (2*f[1][1] + f[2][1]);

/* PAHs */

    m[2]    = NMAX*(f[1][2]*26 + f[2][2]*27 + f[3][2]*51);  /* g.mol-1 */

    prod[2] = 0.;
    for( i = 1; i <= 5; i++ ) prod[2] += k_dep[i][2]; /* s-1 */

    dy_c2h2 += prod[2] * NMAX * f[1][2];
    dy_hcn   = prod[2] * NMAX * f[2][2];
    dy_hc3n += prod[2] * NMAX * f[3][2];

    if( (f[2][2]+f[3][2]) != 0.0e0 )
           csn[2] = (2*f[1][2] + f[2][2]+ 3*f[3][2]) / (f[2][2] + f[3][2]);
    else   csn[2] = 1.0e30;
           csh[2] = 1.;         /* probleme du nombre exact de H */

/* polyHCN */

    m[3]    = NMAX*(f[1][3]*27+f[2][3]*51+f[3][3]*52+f[4][3]*41+f[5][3]*53);  /* g.mol-1 */

    prod[3] = 0.;
    for( i = 1; i <= 5; i++ ) prod[3] += k_dep[i][3]; /* s-1 */

    dy_hcn  += prod[3] * NMAX * f[1][3];
    dy_hc3n += prod[3] * NMAX * f[2][3];
    dy_nccn  = prod[3] * NMAX * f[3][3];
    dy_ch3cn = prod[3] * NMAX * f[4][3];
    dy_c2h3cn= prod[3] * NMAX * f[5][3];

    csn[3] = (f[1][3]+3*f[2][3]+2*f[3][3]+2*f[4][3]+3*f[5][3])
           / (f[1][3]+  f[2][3]+2*f[3][3]+  f[4][3]+  f[5][3]);
    csh[3] = (f[1][3]+3*f[2][3]+2*f[3][3]+2*f[4][3]+3*f[5][3])
           / (f[1][3]+  f[2][3]          +3*f[4][3]+3*f[5][3]);

/* melange */

    csn[0] =  (  prod[1]*(2*f[1][1]+3*f[2][1])
               + prod[2]*(2*f[1][2]+  f[2][2]+3*f[3][2])
               + prod[3]*(  f[1][3]+3*f[2][3]+2*f[3][3]+2*f[4][3]+3*f[5][3]) )
            / (  prod[1]*(            f[2][1])
               + prod[2]*(            f[2][2]+  f[3][2])
               + prod[3]*(  f[1][3]+  f[2][3]+2*f[3][3]+  f[4][3]+  f[5][3]) );
    csh[0] =  (  prod[1]*(2*f[1][1]+3*f[2][1])
               + prod[2]*(2*f[1][2]+  f[2][2]+3*f[3][2])
               + prod[3]*(  f[1][3]+3*f[2][3]+2*f[3][3]+2*f[4][3]+3*f[5][3]) )
            / (  prod[1]*(2*f[1][1]+  f[2][1])
               + prod[2]*(2*f[1][2]+  f[2][2]+  f[3][2])
               + prod[3]*(  f[1][3]+  f[2][3]          +3*f[4][3]+3*f[5][3]) );

/* mass production rates (in kg m-3 s-1) */

    prod[0] = 0.;
    for( i = 1; i <= 3; i++ )
    {
       prod[i]  =  prod[i] * ct * m[i] / 6.022e23 *1.e3;
       prod[0] += prod[i];
    }

    *r1   = dy_c2h2;
    *r2   = dy_hc3n;
    *r3   = dy_hcn;
    *r4   = dy_nccn;
    *r5   = dy_ch3cn;
    *r6   = dy_c2h3cn;

}

