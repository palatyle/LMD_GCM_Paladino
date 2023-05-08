/* tractitan: suivi de traceurs avec constantes de temps de rappel */
/* GCCM */

#include "titan.h"

void tractitan_( double *RB, char CORPS[][10], double Y[][NLEV], 
     double Y0[][NLEV], double *FIN )
{
   char   outlog[100];
   int    i,j,k,l;
   double  annee,**tau,**ym1;
   double  cm,conv,cp,delta,deltamax,deltao;
   double  test,time,ts;
   char   str2[15];
   FILE   *out;

   time     = ts = 0.0e0;
   annee    = 9.46728e8;
   strcpy( outlog, "chimietitan" );
   strcat( outlog, ".log" );
   deltamax = 2.e5;
   
   deltao   = delta = 1.e5;

   ym1       = dm2d( 0,   NC-1, 0, NLEV-1 );
   tau       = dm2d( 0,   NC-1, 0, NLEV-1 );

/* debug */
/*
            out = fopen( "err.log", "a" );
            fprintf( out,"%s\n", );
            fclose( out );
*/

/* Composition pour le rappel (identique a inichim): Y0 */
   
/* initialisation ym1 */
     
   for( j = 0; j <= NLEV-1; j++ )
      for( i = 0; i <= NC-1; i++ ) ym1[i][j] = Y[i][j];
       
/* initialisation tau sans dependance en lat */
    
   for( i = 0; i <= NC-1; i++ )
   { 
    for( j = NLEV-1; j >= 0; j-- ) 
    { 
        tau[i][j] = 1.e6;   /* autres corps = 1.e6 s, donc rappel tres fort */
        
        if( strcmp(CORPS[i],"C2H2") == 0 ) 
           tau[i][j] = annee*pow( 10., 2.+1.*(100.-(RB[j]-R0))/200. );
        if( strcmp(CORPS[i],"C2H6") == 0 ) 
           tau[i][j] = annee*pow( 10., 1.+1.*(200.-(RB[j]-R0))/300. );  
        if( strcmp(CORPS[i],"HCN") == 0 )
        { 
          if( (RB[j]-R0) >= 350. ) 
             tau[i][j] = annee*pow( 10., 1.+1.3*((RB[j]-R0)-350.)/150. );
          else
             tau[i][j] = annee*10.;  
        }
        if( strcmp(CORPS[i],"C4H2") == 0 )
        {
          if( (RB[j]-R0) >= 300. ) 
             tau[i][j] = annee*pow( 10.,-1.+0.3*(300.-(RB[j]-R0))/200. );
          else
             tau[i][j] = annee*pow( 10., 0.+1.0*(100.-(RB[j]-R0))/200. );
        }
    }
/* COUCHES HAUTES: RAPPEL FORCE PLUS GRAND */
    tau[i][NLEV-1] = min(tau[i][NLEV-1],annee/100.);
    tau[i][NLEV-2] = min(tau[i][NLEV-2],annee/50.);
    tau[i][NLEV-3] = min(tau[i][NLEV-3],annee/10.);
/*    tau[i][NLEV-4] = min(tau[i][NLEV-4],annee/100.); */
   }

/*   out = fopen( outlog, "a" ); */
/* vu la rapidite, on laisse le fichier ouvert pendant toute la boucle */ 

/* ***************** */
/*  Main time loop.  */
/* ***************** */

   while( time < (*FIN) )             
   {
     for( i = 0; i <= NC-1; i++ )
     {
/* rappel */   
/* ------ */
       for( j = NLEV-2; j >= 0; j-- )  
          Y[i][j] += delta * ( Y0[i][j] - Y[i][j] ) / tau[i][j];
/* on laisse fixe la couche la plus haute */
       Y[i][NLEV-1] = Y0[i][NLEV-1]; 
     }
     
/* test evolution delta */
/* -------------------- */
     for( j = 0; j <= NLEV-2; j++ ) if( (RB[j]-R0) >= 90. ) 
         for( i = 0; i <= NC-1; i++ )
         {
            test = 1.0e-15; 
            if( ( Y[i][j] > test ) && ( ym1[i][j] > test ) )
            {
               conv = fabs( Y[i][j] - ym1[i][j] ) / ym1[i][j];
               if( conv > ts )
               {
/*
                  if( conv > 0.1 )
                  {
                     fprintf(out, "%d %s %e %e\n",j,CORPS[i],ym1[i][j],Y[i][j]); 
                  }
*/
                  ts = conv;
               }
            }
         }
/*
     fprintf(out, "%e %e %e\n",time,delta,ts); 
*/
     if( ts < 0.1e0 )
     {
         for( i = 0; i <= NC-1; i++ )
               for( j = 0; j <= NLEV-1; j++ )
                 if( Y[i][j] >= 1.0e0 )
                 {
/*
                  fprintf( out, "WARNING %s mixing ratio is %e %e at %d",
                           CORPS[i], ym1[i][j], Y[i][j], j );
                  fclose( out );  
*/
                  exit(0);
                 }
         for( j = 0; j <= NLEV-1; j++ ) 
               for( i = 0; i <= NC-1; i++ ) ym1[i][j] = Y[i][j];
         time += ( deltao = delta );
         if(   ts < 1.00e-5 )                      delta *= 10.e0;
         if( ( ts > 1.00e-5 ) && ( ts < 1.0e-4 ) ) delta *= 5.0e0;
         if( ( ts > 1.00e-4 ) && ( ts < 1.0e-3 ) ) delta *= 2.0e0;
         if( ( ts > 0.001e0 ) && ( ts < 0.01e0 ) ) delta *= 1.5e0;
         if( ( ts > 0.010e0 ) && ( ts < 0.05e0 ) ) delta *= 1.1e0;
         
         delta = min( deltamax, delta );
         
         if( ( time + delta ) > (*FIN) )
         {
            delta = (*FIN) - time;
            time = (*FIN);
         }
     }
     else
     {
         for( j = 0; j <= NLEV-1; j++ ) 
               for( i = 0; i <= NC-1; i++ ) Y[i][j] = ym1[i][j];
         delta *= 0.3e0;
     }
     ts = 0.0e0;
   }               

/* **************** */        
/* end of main loop */
/* **************** */        
     
/*
   fprintf( out, "%e\n", time ); 
   fclose( out ); 
*/   
   fdm2d(  ym1, 0,   NC-1, 0 );
   fdm2d(  tau, 0,   NC-1, 0 );
}
