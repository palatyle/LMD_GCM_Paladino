/* htoh2: production of H2 from heterogenous recombination of H
   on the haze particles */

#include "titan.h"

void heterohtoh2( char corps[][10], double *tp, double *nb, double y[][NLEV],
            double *sh, int *zj,
            double *out1, double *out2, int *utilaer )
{
  int   z;
  int   i,j,h,h2;
  double dy_h2,dy_h,nbCsites;
  double surfhaze,temp,ct;

  z        = (*zj);
  temp     = tp[z];
  ct       = nb[z];
  surfhaze = sh[z];

/* composes interessants */
/* --------------------- */
/* !! decalage de 1 par rapport a calchim !! */

  h  = utilaer[0];
  h2 = utilaer[1];

/* nbCsites: total nb of C sites */
/* ----------------------------- */

/* HYPOTHESE POUR LA TAILLE DU SITE D'UN C */

/* 2e-9*4pi = 2.5e-8 = surface (um2) d'1 C
   en supposant un disque PAH (correspond a un rayon de 0.9AA) */

  nbCsites = surfhaze / 2.5e-8; 

/* taux de recombinaison */
/* --------------------- */

/* H + bounded H -> H2 */
  
  dy_h2 = y[h][z]
         * 1.58e4 * sqrt(temp)      /* kinetic speed of H atoms (cm s-1) */
         * nbCsites                 /* haze: total nb of C sites (cm-3) */
         * 1.8e-18*exp(-300/temp);  /* X-section Y.Sekine (cm2) */
//         * 1.e-15*exp(-1700/temp);  /* X-section for bounded H atoms (cm2) */

  dy_h  = -dy_h2;

/* H + surface -> bounded H */

  if(1==1)  // si faux, surface saturee
  dy_h  = dy_h - y[h][z]
         * 1.58e4 * sqrt(temp)      /* kinetic speed of H atoms (cm s-1) */
	 * nbCsites                 /* haze: total nb of C sites (cm-3) */
	 * 8.8e-16*exp(-1100/temp); /* Xsection Y.Sekine (cm2) */
//         * 1.e-15*exp(-1700/temp);  /* X-section for bounded H atoms (cm2) */

  *out1  = dy_h;
  *out2  = dy_h2;
}

