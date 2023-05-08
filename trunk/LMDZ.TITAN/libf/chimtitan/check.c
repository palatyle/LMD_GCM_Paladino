/* Sanity check routine about coherence between C and Fortran parameters for chemistry */
/* Author : Jan Vatant d'Ollone - 2018 */

#include "titan.h"

void check_(
     int *NLEV_fort, int *NLD_fort,
     int *NLRT_fort, int *NC_fort   )
{

   if( (*NLEV_fort) != NLEV )
   {
    printf("\n");
    printf("!! Coherence problem between C and Fortran !!\n");
    printf("In C, NLEV=%d whereas in Fortran nlaykim_tot=%d\n",NLEV,*(NLEV_fort));
    printf("You surely didn't modified titan.h according to your startfile ! \n");
    printf("I abort ...\n");
    exit(0);
   }

   if( (*NLD_fort) != NLD )
   {
    printf("\n");
    printf("!! Coherence problem between C and Fortran !!\n");
    printf("The number of levels with diffusion must be Numbers-of-levels-of-GCM minus 15 !\n");
    printf("In C, NLD=%d whereas in Fortran klev-15=%d\n",NLD,*(NLD_fort));
    printf("You surely didn't modified titan.h according to your startfile ! \n");
    printf("I abort ...\n");
    exit(0);
   }

   if( (*NLRT_fort) != NLRT )
   {
    printf("\n");
    printf("!! Coherence problem between C and Fortran for the number of UV levels !!\n");
    printf("In C in titan.h, NLRT=%d whereas in Fortran comchem_h, nlrt_kim=%d\n",NLRT,*(NLRT_fort));
    printf("I abort ...\n");
    exit(0);
   }

   if( (*NC_fort) != NC )
   {
    printf("\n");
    printf("!! Coherence problem between C and Fortran for number of compounds !!\n");
    printf("In C, NC=%d whereas in Fortran nkim=%d\n",NC,*(NC_fort));
    printf("You have maybe changed the chemistry but running with old startfiles ! \n");
    printf("I abort ...\n");
    exit(0);
   }

}
