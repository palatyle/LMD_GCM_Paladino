/* chimie: Chemistry fabric */
/* GCCM */

/* ko: low pressure limit:  3 body reaction (1st) */
/* ki: high pressure limit: 2 body reaction (2nd) */

#include "titan.h"

void chimie_( char CORPS[][10], double *NB, double *TEMP, double KRATE[][NLEV], int REACTIF[][5],
             int *NOM_PERTE, int *NOM_PROD, int PERTE[][200][2], int PROD[][200] )
{
   int    dep,i,j,k,l,lat;
   double  ai,ao,ei,eo,ki,ko,m,ti,to;
   FILE   *out;
   static char   reaction[NREAC+1][12][10]={
#include VERCHIM
   "",       "",       "",       "",     "",  "","","","","","",""};

   for( i = 0; i <= NC-1; i++ )
   {
      NOM_PERTE[i] = NOM_PROD[i] = 0;
      for( j = 0; j < 200; j++ )
         PERTE[i][j][0] = PERTE[i][j][1] = PROD[i][j] = 0;
   }
   dep = 0;
   for( i = dep; i <= NREAC-1; i++ )            /* Number of reactions */
   {
      for( j = 0; j <= 4; j++ )                 /* Number of compouds in each reactions */
      {
         k = 0;
         if( (strcmp(reaction[i][j],"prod"))   /* prod and soot are not */
          && (strcmp(reaction[i][j],"soot"))   /* considered in the compounds */
          && (strcmp(reaction[i][j],""))     ) /* Which compound ? */
         {
            while( strcmp(reaction[i][j],CORPS[k]) )     /* Compound in reaction j, column i */
            {
               if( k == NC+1 )
               {
                  out = fopen( "err.log", "a" );
                  fprintf( out, "I cannot find %s\n", reaction[i][j] );
                  fclose( out );
                  exit(0);
               }
               k++;
            }
            REACTIF[i][j] = k;                         /* is number k */
         }
         else REACTIF[i][j] = NC;
      }
   }
   for( i = dep; i <= NREAC-1; i++ )     /* Total loss and production */
   {
      j = REACTIF[i][0];                 /* First compound reaction i */
      k = NOM_PERTE[j] + 1;              /* Loss one more times */
      if( k >= 200 ) exit(0);
      NOM_PERTE[j] = k;
      PERTE[j][k-1][0] = i;              /* Compound j loss in reaction i for the kth time */
      PERTE[j][k-1][1] = REACTIF[i][1];  /* Compound j reacts with number 2 in reaction i */
      j = REACTIF[i][1];                 /* Second compound reaction i */
      if( j != NC )                      /* Neither photodissociation nor desexcitation */
      {
       k = NOM_PERTE[j] + 1;             /* Loss one more times */
       if( k >= 200 ) exit(0);
       NOM_PERTE[j] = k;
       PERTE[j][k-1][0] = i;             /* Compound j is loss in reaction i for the kth times */
       PERTE[j][k-1][1] = REACTIF[i][0]; /* Compound j reacts with number 1 reaction i */
      }
      for( j = 2; j < 5; j++ )           /* Compounds from 3 to 5 */
      {
         k = REACTIF[i][j];              /* Number of compound */
         if( k != NC )
         {
            l = NOM_PROD[k] + 1;         /* One more time produced */
            if( l >= 200 ) exit(0);
            NOM_PROD[k]  = l;
            PROD[k][l-1] = i;             /* at reaction i */
         }
      }
   }
   for( i = RDISS+1; i <= NREAC-1; i++ )   /* 0 a RDISS-1: photodiss, RDISS: disso N2 */
   {
      ao = strtod(reaction[i][5], NULL);
      to = strtod(reaction[i][6], NULL);
      eo = strtod(reaction[i][7], NULL);
      ai = strtod(reaction[i][8], NULL);
      ti = strtod(reaction[i][9], NULL);
      ei = strtod(reaction[i][10],NULL);
      m  = strtod(reaction[i][11],NULL);
      for( j = 0; j <= NLEV-1; j++ )
      {
            ko = ao * pow( TEMP[j], to) * exp( eo / TEMP[j] );
            if( m == 1.0e0 )
               KRATE[i][j] = ko;
            else if( m == 2.0e0 )
            {
               KRATE[i][j] = ko * pow( NB[j], 2.0e0 );
            }
            else if( m == 3.0e0 )
            {
               if( ai == 0.0e0 )
                  KRATE[i][j] = ko * pow( NB[j], 3.0e0 );
               else
               {
                  ki = ai * pow( TEMP[j], ti ) * exp( ei / TEMP[j] );
                  KRATE[i][j] = pow( NB[j], 3.0e0 ) * ko * ki / ( ko * NB[j] + ki );
               }
            }
      }
   }
}
