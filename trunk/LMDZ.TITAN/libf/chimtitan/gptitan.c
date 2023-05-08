/* gptitan: photochimie */
/* GCCM */

/* nitriles et hydrocarbures separes pour l'inversion */

/* flux variable au sommet */

#include "titan.h"

void gptitan_(
	      double *RA, double *TEMP, double *NB, 
	      char CORPS[][10], double Y[][NLEV],
	      double *FIN, double *LAT, double *MASS, double MD[][NLEV],
	      double *KEDD, double KRATE[][NLEV],
	      int reactif[][5], int *nom_prod, int *nom_perte, 
	      int prod[][200], int perte[][200][2], int *aerprod, int *utilaer, 
	      double MAER[][NLEV], double PRODAER[][NLEV], 
	      double CSN[][NLEV], double CSH[][NLEV],
	      int *htoh2, double *surfhaze)
{
  char   outlog[100];
  int    i,j,k,l;
  int    ireac,ncom1,ncom2;
  int    in2, ih, ih2, in4s;
  double  ***a,***b,**c;
  double  *fl,*fp,*mu,**jac,**ym1,**f;
  double  conv,delta,deltamax;
  double  cm,cp,dim,dip,dm,dp,dym,dyp,km,kp,r,dra,dram,drap;
  double  np,nm,s,test,time,ts,v,dv;
  char   str2[15];
  FILE   *out;

  /* va avec htoh2 */
  double  dyh,dyh2;

  /* va avec aer */
  double  dyc2h2,dyhc3n,dyhcn,dynccn,dych3cn,dyc2h3cn;
  double  **k_dep,**faer;
  double  *productaer,*csurn,*csurh,*mmolaer;

  if( (*aerprod) == 1 )
    {
      k_dep = dm2d( 1, 5, 1, 3 );     /* k en s-1, reactions d'initiation */
      faer  = dm2d( 1, 5, 1, 3 );     /* fraction de chaque compose */
      productaer  = dm1d( 0, 3 );     /* local production rate by pathways */
      mmolaer     = dm1d( 0, 3 );     /* local molar mass by pathways */
      csurn = dm1d( 0, 3 );           /* local C/N by pathways */
      csurh = dm1d( 0, 3 );           /* local C/H by pathways */
    }

  /* DEBUG */
  printf("CHIMIE: lat=%g\n",(*LAT));
  /**/

  /* Warning about computational time */
#ifndef LAPACK
  printf("Hey there !! LAPACK key is turned off, the chemistry could be way quicker if you were activating it !\n");
#endif
  /**/
   
  strcpy( outlog, "chimietitan" );
  strcat( outlog, ".log" );
  out = fopen( outlog, "w" );
  fprintf(out,"CHIMIE: lat=%g\n",(*LAT));
  fclose( out );

  deltamax = 1.e5;
  test = 1.0e-15; 

  /* valeur de r:
     r = g0 R0^2 / R * 2 * 1E-3
     avec g0 en cm/s2, R0 en km, mu et mass en g
  */
  r        = 21.595656e0;

  /*  DEBUG 
      out = fopen( outlog, "a" );
      fprintf(out,"CHIMIE: lat=%g\n",(*LAT));
      fclose( out );
  */
  fl      = dm1d( 0,   NC-1 );
  fp      = dm1d( 0,   NC-1 );
  mu      = dm1d( 0, NLEV-1 );
  ym1     = dm2d( 0,   NC-1, 0, NLEV-1 );
  f       = dm2d( 0,   NC-1, 0, NLEV-1 );
  jac     = dm2d( 0,   NC-1, 0,   NC-1 );
  c       = dm2d( 0, NLEV-1, 0,   NC-1 );
  a       = dm3d( 0, NLEV-1, 0,   NC-1, 0,   NC-1 );
  b       = dm3d( 0, NLEV-1, 0,   NC-1, 1,    2 );

  /* DEBUG */
  /* 
     out = fopen( "err.log", "a" );
     fprintf( out,"%s\n", );
     fclose( out );
  */

  /* init indexes of special species */
  for( i = 0; i <= ST-1; i++ ) 
    {
      if( strcmp(CORPS[i], "H"  ) == 0 ) ih   = i;
      if( strcmp(CORPS[i], "H2" ) == 0 ) ih2  = i;
      if( strcmp(CORPS[i], "N2" ) == 0 ) in2  = i;
      if( strcmp(CORPS[i], "N4S") == 0 ) in4s = i;
    }

  /* initialisation mu */
     
  for( j = 0; j <= NLEV-1; j++ )
    {
      mu[j] = 0.0e0;
      for( i = 0; i <= ST-1; i++ ) 
	{
	  mu[j] += ( MASS[i] * Y[i][j] );
	}
    }

  /* initialisation compo avant calcul */
  for( j = NLEV-1; j >= 0; j-- )
    for( i = 0; i <= ST-1; i++ ) ym1[i][j] = max(Y[i][j],1.e-30);

  /* 
     ==========================================================================
     STRATEGIE:
     INVERSION COMPLETE AVEC DIFFUSION ENTRE NLEV-1 et NLD 
     PUIS INVERSION LOCALE PAR BLOC ENTRE NLD ET LA SURFACE
     ==========================================================================

     PREMIERE ETAPE: 
     ===============
     INVERSION COMPLETE AVEC DIFFUSION ENTRE NLEV-1 et NLD 
     ===============
  */

  /* ****************** */
  /*  Time loop:        */
  /* ****************** */

  time     = ts = 0.0e0;
  delta    = 1.e-3;

  while( time < (*FIN) )     
    {


      /* DEBUG 
	 for( j = NLEV-1; j >= NLD; j-- )
	 {
	 out = fopen( outlog, "a" );
	 fprintf(out,"j=%d z=%e nb=%e T=%e\n",j,(RA[j]-R0),NB[j],TEMP[j]);
	 fclose( out );

	 out = fopen( "profils.log", "a" );
	 fprintf(out,"%d %e %e %e\n",j,(RA[j]-R0),NB[j],TEMP[j]);
	 for (i=0;i<=NREAC-1;i++) fprintf(out,"%d %e\n",i,KRATE[i][j]);
	 for (i=0;i<=ST-1;i++) fprintf(out,"%10s %e\n",CORPS[i],Y[i][j]);
	 fclose( out );
	 }
	 exit(0);
      */


      /* ------------------------------ */
      /* Calculs variations et jacobien */
      /* ------------------------------ */

      for( j = NLEV-1; j >= NLD; j-- )
	{

	  /* init of step */
	  /* ------------ */
	  for( i = 0; i <= ST-1; i++ ) 
	    {
	      fp[i] = fl[i] = 0.0e0; 
	      for( l = 0; l <= ST-1; l++ ) jac[i][l] = 0.0e0;
	    }

	  /* Chimie */
	  /* ------ */

	  /* productions et pertes chimiques */
	  for( i = 0; i <= ST-1; i++ )    
	    {
	      Y[i][j] = max(Y[i][j],1.e-30);                /* minimum */

	      for( l = 0; l <= nom_prod[i]-1; l++ )    /* Production term */
		{
		  ireac = prod[i][l];                  /* Number of the reaction involves. */
		  ncom1 = reactif[ireac][0];           /* First compound which reacts. */
		  if( reactif[ireac][1] == NC )        /* Photodissociation or relaxation */
		    {
		      jac[i][ncom1] += ( KRATE[ireac][j] * NB[j] );
		      fp[i]         += ( KRATE[ireac][j] * NB[j] * Y[ncom1][j] );
		    }
		  else                                 /* General case. */
		    {
		      ncom2          = reactif[ireac][1];                       /* Second compound which reacts. */
		      jac[i][ncom1] += ( KRATE[ireac][j] * Y[ncom2][j] );       /* Jacobian compound #1. */
		      jac[i][ncom2] += ( KRATE[ireac][j] * Y[ncom1][j] );       /* Jacobian compound #2. */
		      fp[i] += ( KRATE[ireac][j] * Y[ncom1][j] * Y[ncom2][j] ); /* Production term. */
		    }
		}
            
	      for( l = 0; l <= nom_perte[i]-1; l++ )   /* Loss term. */
		{
		  ireac = perte[i][l][0];              /* Reaction number. */
		  ncom2 = perte[i][l][1];              /* Compound #2 reacts. */
		  if( reactif[ireac][1] == NC )        /* Photodissociation or relaxation */
		    {
		      jac[i][i] -= ( KRATE[ireac][j] * NB[j] );
		      fl[i]     += ( KRATE[ireac][j] * NB[j] );
		    }
		  else                                 /* General case. */
		    {
		      jac[i][ncom2] -= ( KRATE[ireac][j] * Y[i][j] );       /* Jacobian compound #1. */
		      jac[i][i]     -= ( KRATE[ireac][j] * Y[ncom2][j] );   /* Jacobien compound #2. */
		      fl[i]         += ( KRATE[ireac][j] * Y[ncom2][j] );   /* Loss term. */
		    }
		}
	    }


	  /* Aerosols */
	  /* -------- */
	  if( (*aerprod) == 1 )
	    {
	      aer(CORPS,TEMP,NB,Y,&j,k_dep,faer,
		  &dyc2h2,&dyhc3n,&dyhcn,&dynccn,&dych3cn,&dyc2h3cn,utilaer,
		  mmolaer,productaer,csurn,csurh);

	      for( i = 0; i <= 3; i++ )
		{
		  PRODAER[i][j] = productaer[i];
                  MAER[i][j] = mmolaer[i];
		  CSN[i][j] = csurn[i];
		  CSH[i][j] = csurh[i];
		}
	      /* DEBUG
		 printf("AERPROD : LAT = %g - J = %d\n",(*LAT),j);
		 if(fabs(dyc2h2*NB[j])>fabs(fp[utilaer[2]]/10.))
		 printf("fp(%s) =%e; dyc2h2 =%e\n",CORPS[utilaer[2]],
		 fp[utilaer[2]],dyc2h2*NB[j]);
		 if(fabs(dyhcn*NB[j])>fabs(fp[utilaer[5]]/10.))
		 printf("fp(%s) =%e; dyhcn  =%e\n",CORPS[utilaer[5]],
		 fp[utilaer[5]],dyhcn*NB[j]);
		 if(fabs(dyhc3n*NB[j])>fabs(fp[utilaer[6]]/10.))
		 printf("fp(%s) =%e; dyhc3n =%e\n",CORPS[utilaer[6]],
		 fp[utilaer[6]],dyhc3n*NB[j]);
		 if(fabs(dynccn*NB[j])>fabs(fp[utilaer[13]]/10.))
		 printf("fp(%s) =%e; dynccn =%e\n",CORPS[utilaer[13]],
		 fp[utilaer[13]],dynccn*NB[j]);
		 if(fabs(dych3cn*NB[j])>fabs(fp[utilaer[14]]/10.))
		 printf("fp(%s) =%e; dych3cn=%e\n",CORPS[utilaer[14]],
		 fp[utilaer[14]],dych3cn*NB[j]);
		 if(fabs(dyc2h3cn*NB[j])>fabs(fp[utilaer[15]]/10.))
		 printf("fp(%s) =%e; dyc2h3cn=%e\n",CORPS[utilaer[15]],
		 fp[utilaer[15]],dyc2h3cn*NB[j]);
	      */

	      fp[utilaer[2]] -= (   dyc2h2 * NB[j] );
	      fp[utilaer[5]] -= (    dyhcn * NB[j] );
	      fp[utilaer[6]] -= (   dyhc3n * NB[j] );
	      fp[utilaer[13]]-= (   dynccn * NB[j] );
	      fp[utilaer[14]]-= (  dych3cn * NB[j] );
	      fp[utilaer[15]]-= ( dyc2h3cn * NB[j] );
	      if( Y[utilaer[2]][j]  != 0.0 )
		jac[utilaer[2]][utilaer[2]] -= (  dyc2h2 * NB[j] / Y[utilaer[2]][j] );
	      if( Y[utilaer[5]][j]  != 0.0 )
		jac[utilaer[5]][utilaer[5]] -= (   dyhcn * NB[j] / Y[utilaer[5]][j] );
	      if( Y[utilaer[6]][j]  != 0.0 )
		jac[utilaer[6]][utilaer[6]] -= (  dyhc3n * NB[j] / Y[utilaer[6]][j] );
	      if( Y[utilaer[13]][j] != 0.0 )
		jac[utilaer[13]][utilaer[13]] -= (  dynccn * NB[j] / Y[utilaer[13]][j] );
	      if( Y[utilaer[14]][j] != 0.0 )
		jac[utilaer[14]][utilaer[14]] -= ( dych3cn * NB[j] / Y[utilaer[14]][j] );
	      if( Y[utilaer[15]][j] != 0.0 )
		jac[utilaer[15]][utilaer[15]] -= (dyc2h3cn * NB[j] / Y[utilaer[15]][j] );
	    }
     
       
	  /* H -> H2 on haze particles */
	  /* ------------------------- */
	  if( (*htoh2) == 1 )
	    {
              heterohtoh2(CORPS,TEMP,NB,Y,surfhaze,&j,&dyh,&dyh2,utilaer);
	      /* dyh <= 0 / 1.0 en adsor., 1 en reac. */

	      /* DEBUG 
		 printf("HTOH2 : LAT = %g - J = %d\n",(*LAT),j);
		 if(fabs(dyh*NB[j])>fabs(fp[utilaer[0]]/10.))
		 printf("fp(%s) = %e; dyh  = %e\n",CORPS[utilaer[0]],fp[utilaer[0]],dyh*NB[j]);
		 if(fabs(dyh2*NB[j])>fabs(fp[utilaer[1]]/10.))
		 printf("fp(%s) = %e; dyh2 = %e\n",CORPS[utilaer[1]],fp[utilaer[1]],dyh2*NB[j]);
	      */

              fp[utilaer[0]] += ( dyh  * NB[j] ); 
	      /* pourquoi pas *2 ?? cf gptit dans 2da... */

              fp[utilaer[1]] += ( dyh2 * NB[j] );
              if( Y[utilaer[0]][j] != 0.0 )
		jac[utilaer[0]][utilaer[0]] += ( dyh  * NB[j] / Y[utilaer[0]][j] );
	      /* pourquoi pas *2 ?? cf gptit dans 2da... */
	    }


	  /* Backup jacobian level j. */
	  /* ------------------------ */
	  for( i = 0; i <= ST-1; i++ )
            for( k = 0; k <= ST-1; k++ )
	      a[j][i][k] = jac[i][k];     


	  /* Diffusion verticale et flux exterieurs */
	  /* -------------------------------------- */

	  /*
	    pour dy/dr, dr doit etre en cm...
	    pareil pour dphi/dr
	  */
	  for( i = 0; i <= ST-1; i++ )
	    {

	      /* First level. */
	      if( j == NLD )
		{
		  v = dv = 0.0e0;
		  dra = RA[j+1]-RA[j];

		  cp  = (NB[j+1]+NB[j])/2.;  /* Mean total concentration. */
		  dip = r * (MASS[i]-(mu[j+1]+mu[j])/2.) / (TEMP[j+1]+TEMP[j]) /
		    pow( RA[j+1], 2.0e0 );    /* Delta i,j level +1. */
		  dp  = (MD[i][j]+MD[i][j+1])/2.;     /* Mean molecular diffusion. */
		  dyp = (Y[i][j+1]-Y[i][j])/(RA[j+2]-RA[j])*2.e-5; /* Delta y level +1. */
		  kp  = (KEDD[j+1]+KEDD[j])/2.;       /* Mean eddy diffusion. */
		  /* div phi. */
		  f[i][j] = cp * ( dp * ( (Y[i][j+1]+Y[i][j])/2. * dip + dyp ) 
				   + kp * dyp ) 
		    * (4.e-5/dra/pow((1.+RA[j]/RA[j+1]),2.)) 
		    + fp[i] - Y[i][j]*fl[i] + v;
		  /* dphi / dy this level. */
		  a[j][i][i] += ( cp * ( dp * 0.5e0 * dip 
					 - 2.e-5/(RA[j+2]-RA[j]) * (dp + kp) ) 
				  * (4.e-5/dra/pow((1.+RA[j]/RA[j+1]),2.)) + dv );
		  /* dphi / dy level +1. */
		  c[j][i]     = -THETA * delta 
		    * cp * ( dp * 0.5e0 * dip 
			     + 2.e-5/(RA[j+2]-RA[j]) * (dp + kp) ) 
		    * (4.e-5/dra/pow((1.+RA[j]/RA[j+1]),2.));
		}
	      /* Last level. */
	      else if( j == NLEV-1 )
		{
		  v = dv = 0.0e0;
		  dra = RA[NLEV-1]-RA[NLEV-2];

		  /* Jeans escape */
		  if( i == ih )  
		    {
		      dv = top_H  * NB[NLEV-1] 
                        * (4.e-5/dra/pow((2.-dra/(RA[NLEV-1]+dra)),2.)); 
		      v  = dv * Y[i][NLEV-1];
		    }
		  else if( i == ih2 )
		    {
		      dv = top_H2 * NB[NLEV-1]
                        * (4.e-5/dra/pow((2.-dra/(RA[NLEV-1]+dra)),2.)); 
		      v  = dv * Y[i][NLEV-1];
		    }
		  /* Input flux for N(4S) */
		  else if( i == in4s )
		    v  = top_N4S
		      * (4.e-5/dra/pow((2.-dra/(RA[NLEV-1]+dra)),2.)); 

		  cm  = (NB[NLEV-1]+NB[NLEV-2])/2.;  /* Mean total concentration. */
		  dim = r * (MASS[i]-(mu[NLEV-1]+mu[NLEV-2])/2.)
		    / (TEMP[NLEV-1]+TEMP[NLEV-2]) 
		    / pow( RA[NLEV-1],   2.0e0 );  /* Delta i,j level -1. */
		  dm  = (MD[i][NLEV-1]+MD[i][NLEV-2])/2.;    /* Mean molecular diffusion. */
		  dym = (Y[i][NLEV-1]-Y[i][NLEV-2])/dra*1.e-5; /* Delta y level -1. */
		  km  = (KEDD[NLEV-1]+KEDD[NLEV-2])/2.;      /* Mean eddy diffusion. */
		  /* div phi. */
		  f[i][NLEV-1] = fp[i] - Y[i][NLEV-1]*fl[i] - v 
		    - cm * ( dm * ( (Y[i][NLEV-1]+Y[i][NLEV-2])/2. * dim + dym ) 
			     + km * dym ) 
		    * (4.e-5/dra/pow((2.+dra/RA[NLEV-1]),2.)); 
		  /* dphi / dy this level */
		  a[NLEV-1][i][i] -= ( cm * ( dm * 0.5e0 * dim 
					      + 1.e-5/dra * (dm + km ) )
				       * (4.e-5/dra/pow((2.+dra/RA[NLEV-1]),2.)) + dv );
		  /* dphi / dy level -1. */
		  b[NLEV-1][i][2]  =  THETA * delta 
		    * cm * ( dm * 0.5e0 * dim 
			     - 1.e-5/dra * (dm + km ) )
		    * (4.e-5/dra/pow((2.+dra/RA[NLEV-1]),2.));
		}
	      else
		{
		  v = dv = 0.0e0;
		  dram=(RA[j+1]-RA[j-1])/2.;
		  if (j<NLEV-2)
		    drap=(RA[j+1]-RA[j-1])/2.;
		  else
		    drap=dram;

		  cm  = (NB[j]+NB[j-1])/2.;       /* Mean concentration level -1. */
		  cp  = (NB[j]+NB[j+1])/2.;       /* Mean concentration level +1. */
		  dip = r * (MASS[i]-(mu[j+1]+mu[j])/2.) / (TEMP[j+1]+TEMP[j]) /
		    pow( RA[j+1], 2.0e0 );    /* Delta i,j level +1. */
		  dim = r * (MASS[i]-(mu[j]+mu[j-1])/2.) / (TEMP[j]+TEMP[j-1]) /
		    pow( RA[j],   2.0e0 );    /* Delta i,j level -1. */
		  dm  = (MD[i][j-1]+MD[i][j])/2.;    /* Mean molecular diffusion level -1. */
		  dp  = (MD[i][j+1]+MD[i][j])/2.;    /* Mean molecular diffusion level +1. */
		  dym = (Y[i][j]-Y[i][j-1])/dram*1.e-5; /* Delta y level -1. */
		  dyp = (Y[i][j+1]-Y[i][j])/drap*1.e-5; /* Delta y level +1. */
		  km  = (KEDD[j]+KEDD[j-1])/2.;      /* Mean eddy diffusion level -1. */
		  kp  = (KEDD[j]+KEDD[j+1])/2.;      /* Mean eddy diffusion level +1. */
		  /* div phi. */
		  f[i][j] = cp * ( dp * ( (Y[i][j+1]+Y[i][j])/2. * dip + dyp ) 
				   + kp * dyp ) 
		    * (4.e-5/(RA[j+1]-RA[j])/pow((1.+RA[j]/RA[j+1]),2.)) 
		    - cm * ( dm * ( (Y[i][j]+Y[i][j-1])/2. * dim + dym ) 
			     + km * dym ) 
		    * (4.e-5/(RA[j+1]-RA[j])/pow((1.+RA[j+1]/RA[j]),2.)) 
		    + fp[i] - fl[i] * Y[i][j] + v;
		  /* dphi / dy this level */
		  a[j][i][i] += ( cp * ( dp * 0.5e0 * dip 
					 - 1.e-5/drap * (dp + kp) ) 
				  * (4.e-5/(RA[j+1]-RA[j])/pow((1.+RA[j]/RA[j+1]),2.))
				  - cm * ( dm * 0.5e0 * dim 
					   + 1.e-5/dram * (dm + km ) )
				  * (4.e-5/(RA[j+1]-RA[j])/pow((1.+RA[j+1]/RA[j]),2.)) );
		  /* dphi / dy level -1. */
		  b[j][i][2]  =  THETA * delta 
		    * cm * ( dm * 0.5e0 * dim
			     - 1.e-5/dram * (dm + km ) )
		    * (4.e-5/(RA[j+1]-RA[j])/pow((1.+RA[j+1]/RA[j]),2.));
		  /* dphi / dy level +1. */
		  c[j][i]     = -THETA * delta
		    * cp * ( dp * 0.5e0 * dip
			     + 1.e-5/drap * (dp + kp) ) 
		    * (4.e-5/(RA[j+1]-RA[j])/pow((1.+RA[j]/RA[j+1]),2.));
		}
	    }



	  /* finition pour inversion */
	  /* ----------------------- */

	  for( i = 0; i <= ST-1; i++ )
	    {
	      for( k = 0; k <= ST-1; k++ )
		{
		  a[j][i][k] *= ( -THETA * delta );  /* Correction time step. */
		  if( k == i ) a[j][k][k] += NB[j];  /* Correction diagonal. */
		}
	      f[i][j] *= delta;
	    }

	}


      /* -------------------------------- */
      /* Inversion of matrix cf method LU */
      /* -------------------------------- */

      for( j = NLD+1; j <= NLEV-1; j++ )
	{
#ifdef LAPACK
	  solve_lapack( a, j-1, 0, ST-1 );
#else
	  solve( a, j-1, 0, ST-1 );
#endif
	  for( i = 0; i <= ST-1; i++ )
	    {
	      s = 0.0e0;
	      for( k = 0; k <= ST-1; k++ )
		{
		  a[j][i][k] -= ( b[j][i][2] * c[j-1][k] * a[j-1][i][k] );
		  s          += ( b[j][i][2] * f[k][j-1] * a[j-1][i][k] );
		}
	      f[i][j] -= s;
	    }
	}
#ifdef LAPACK
      solve_lapack( a, NLEV-1, 0, ST-1 );
#else
      solve( a, NLEV-1, 0, ST-1 );
#endif
      for( j = NLEV-1; j >= NLD; j-- )     
	{
	  if( j != NLEV-1 )
            for( i = 0; i <= ST-1; i++ ) f[i][j] -= ( c[j][i] * b[j+1][i][1] );
	  for( i = 0; i <= ST-1; i++ )
	    {
	      s = 0.0e0;
	      for( k = 0; k <= ST-1; k++ ) s += ( a[j][i][k] * f[k][j] );
	      b[j][i][1]  = s;
	      Y[i][j]    += s;
	      if( Y[i][j] <= 1.0e-30 ) Y[i][j] = 0.0e0;
	    }
	}

      /* ------------------ */
      /* Tests et evolution */
      /* ------------------ */

      /* Calcul deviation */
      /* ---------------- */

      for( j = NLD; j <= NLEV-1; j++ )
	for( i = 0; i <= ST-1; i++ )
	  if( ( Y[i][j] > test ) && ( ym1[i][j] > test ) )
	    {
	      conv = fabs( Y[i][j] - ym1[i][j] ) / ym1[i][j];
	      if( conv > ts )
		{
		  /*
		    if( conv >= 0.1 )
		    {
		    out = fopen( outlog, "a" );
		    fprintf( out, "Latitude %g;", (*LAT));
		    fprintf(out, " alt:%e; %s %e %e ; %e %e\n",(RA[j]-R0),CORPS[i],ym1[i],Y[i][j],time,delta);
		    fclose( out );
		    }
		  */
                  ts = conv;
		}
	    }

      /* test deviation */
      /* -------------- */

      if( ts < 0.1e0 )
	{
	  for( i = 0; i <= ST-1; i++ )
	    for( j = NLD; j <= NLEV-1; j++ )
	      if( (Y[i][j] >= 0.5e0) && (i != in2) )
		{
                  out = fopen( outlog, "a" );
                  fprintf( out, "WARNING %s mixing ratio is %e %e at %d\n",
                           CORPS[i], ym1[i], Y[i][j], j );
                  for( k = 0; k <= NLEV-1; k++ ) fprintf( out, "%d %e %e\n",k,ym1[i],Y[i][k] );
                  fclose( out );
                  exit(0); 
		  //                  Y[i][j] = 1.e-20;
		}
	  for( j = NLD; j <= NLEV-1; j++ )
	    for( i = 0; i <= NC-1; i++ ) ym1[i][j] = max(Y[i][j],1.e-30);
	  time += delta;
	  if(   ts < 1.00e-5 )                      delta *= 1.0e2;
	  if( ( ts > 1.00e-5 ) && ( ts < 1.0e-4 ) ) delta *= 1.0e1;
	  if( ( ts > 1.00e-4 ) && ( ts < 1.0e-3 ) ) delta *= 5.0e0;
	  if( ( ts > 1.00e-3 ) && ( ts < 5.0e-3 ) ) delta *= 3.0e0;
	  if( ( ts > 5.00e-3 ) && ( ts < 0.01e0 ) ) delta *= 1.5e0;
	  if( ( ts > 0.010e0 ) && ( ts < 0.03e0 ) ) delta *= 1.2e0;
	  if( ( ts > 0.030e0 ) && ( ts < 0.05e0 ) ) delta *= 1.1e0;

	  //            if( ( ts > 0.001e0 ) && ( ts < 0.01e0 ) ) delta *= 3.0e0;
	  //            if( ( ts > 0.010e0 ) && ( ts < 0.05e0 ) ) delta *= 1.5e0;
         
	  delta = min( deltamax, delta );
	}
      else
	{
	  for( j = NLD; j <= NLEV-1; j++ )
	    for( i = 0; i <= NC-1; i++ ) Y[i][j] = ym1[i][j];

	  if(   ts > 0.8 )                    delta *= 1.e-6;
	  if( ( ts > 0.6 ) && ( ts <= 0.8 ) ) delta *= 1.e-4;
	  if( ( ts > 0.4 ) && ( ts <= 0.6 ) ) delta *= 1.e-2;
	  if( ( ts > 0.3 ) && ( ts <= 0.4 ) ) delta *= 0.1;
	  if( ( ts > 0.2 ) && ( ts <= 0.3 ) ) delta *= 0.2;
	  if( ( ts > 0.1 ) && ( ts <= 0.2 ) ) delta *= 0.3;
	}
      ts = 0.0e0;

      out = fopen( outlog, "a" );
      fprintf(out, "delta:%e; time:%e; fin:%e\n",delta,time,(*FIN));
      fclose( out );

    }
  /* **************** */        
  /* end of time loop */
  /* **************** */        

  /*
    ==========================================================================

    SECONDE ETAPE: 
    ===============
    INVERSION LOCALE PAR BLOC ENTRE NLD ET LA SURFACE
    ===============
  */
  if( NLD != 0 ) 
    for( j = NLD-1; j >= 0; j-- )
      {
	time     = ts = 0.0e0;
	delta    = 1.e-3;

	/* ++++++++++++ */
	/*  time loop.  */
	/* ++++++++++++ */

	while( time < (*FIN) )     
	  {

	    /* init of step */
	    /* ------------ */
	    for( i = 0; i <= ST-1; i++ ) 
	      {
		fp[i] = fl[i] = 0.0e0; 
		for( l = 0; l <= ST-1; l++ ) jac[i][l] = 0.0e0;
	      }

	    /* Chimie */
	    /* ------ */

	    /* productions et pertes chimiques */
	    for( i = 0; i <= ST-1; i++ )    
	      {
		Y[i][j] = max(Y[i][j],1.e-30);                /* minimum */

		for( l = 0; l <= nom_prod[i]-1; l++ )    /* Production term */
		  {
		    ireac = prod[i][l];                  /* Number of the reaction involves. */
		    ncom1 = reactif[ireac][0];           /* First compound which reacts. */
		    if( reactif[ireac][1] == NC )        /* Photodissociation or relaxation */
		      {
			jac[i][ncom1] += ( KRATE[ireac][j] * NB[j] );
			fp[i]         += ( KRATE[ireac][j] * NB[j] * Y[ncom1][j] );
		      }
		    else                                 /* General case. */
		      {
			ncom2          = reactif[ireac][1];                       /* Second compound which reacts. */
			jac[i][ncom1] += ( KRATE[ireac][j] * Y[ncom2][j] );       /* Jacobian compound #1. */
			jac[i][ncom2] += ( KRATE[ireac][j] * Y[ncom1][j] );       /* Jacobian compound #2. */
			fp[i] += ( KRATE[ireac][j] * Y[ncom1][j] * Y[ncom2][j] ); /* Production term. */
		      }
		  }
            
		for( l = 0; l <= nom_perte[i]-1; l++ )   /* Loss term. */
		  {
		    ireac = perte[i][l][0];              /* Reaction number. */
		    ncom2 = perte[i][l][1];              /* Compound #2 reacts. */
		    if( reactif[ireac][1] == NC )        /* Photodissociation or relaxation */
		      {
			jac[i][i] -= ( KRATE[ireac][j] * NB[j] );
			fl[i]     += ( KRATE[ireac][j] * NB[j] );
		      }
		    else                                 /* General case. */
		      {
			jac[i][ncom2] -= ( KRATE[ireac][j] * Y[i][j] );       /* Jacobian compound #1. */
			jac[i][i]     -= ( KRATE[ireac][j] * Y[ncom2][j] );   /* Jacobien compound #2. */
			fl[i]         += ( KRATE[ireac][j] * Y[ncom2][j] );   /* Loss term. */
		      }
		  }
	      }


	    /* Aerosols */
	    /* -------- */
	    if( (*aerprod) == 1 )
	      {
		aer(CORPS,TEMP,NB,Y,&j,k_dep,faer,
		    &dyc2h2,&dyhc3n,&dyhcn,&dynccn,&dych3cn,&dyc2h3cn,utilaer,
		    mmolaer,productaer,csurn,csurh);

		for( i = 0; i <= 3; i++ )
		  {
		    PRODAER[i][j] = productaer[i];
		    MAER[i][j] = mmolaer[i];
		    CSN[i][j] = csurn[i];
		    CSH[i][j] = csurh[i];
		  }
		/* DEBUG
		   printf("AERPROD : LAT = %g - J = %d\n",(*LAT),j);
		   if(fabs(dyc2h2*NB[j])>fabs(fp[utilaer[2]]/10.))
		   printf("fp(%s) =%e; dyc2h2 =%e\n",CORPS[utilaer[2]],
		   fp[utilaer[2]],dyc2h2*NB[j]);
		   if(fabs(dyhcn*NB[j])>fabs(fp[utilaer[5]]/10.))
		   printf("fp(%s) =%e; dyhcn  =%e\n",CORPS[utilaer[5]],
		   fp[utilaer[5]],dyhcn*NB[j]);
		   if(fabs(dyhc3n*NB[j])>fabs(fp[utilaer[6]]/10.))
		   printf("fp(%s) =%e; dyhc3n =%e\n",CORPS[utilaer[6]],
		   fp[utilaer[6]],dyhc3n*NB[j]);
		   if(fabs(dynccn*NB[j])>fabs(fp[utilaer[13]]/10.))
		   printf("fp(%s) =%e; dynccn =%e\n",CORPS[utilaer[13]],
		   fp[utilaer[13]],dynccn*NB[j]);
		   if(fabs(dych3cn*NB[j])>fabs(fp[utilaer[14]]/10.))
		   printf("fp(%s) =%e; dych3cn=%e\n",CORPS[utilaer[14]],
		   fp[utilaer[14]],dych3cn*NB[j]);
		   if(fabs(dyc2h3cn*NB[j])>fabs(fp[utilaer[15]]/10.))
		   printf("fp(%s) =%e; dyc2h3cn=%e\n",CORPS[utilaer[15]],
		   fp[utilaer[15]],dyc2h3cn*NB[j]);
		*/

		fp[utilaer[2]] -= (   dyc2h2 * NB[j] );
		fp[utilaer[5]] -= (    dyhcn * NB[j] );
		fp[utilaer[6]] -= (   dyhc3n * NB[j] );
		fp[utilaer[13]]-= (   dynccn * NB[j] );
		fp[utilaer[14]]-= (  dych3cn * NB[j] );
		fp[utilaer[15]]-= ( dyc2h3cn * NB[j] );
		if( Y[utilaer[2]][j]  != 0.0 )
		  jac[utilaer[2]][utilaer[2]] -= (  dyc2h2 * NB[j] / Y[utilaer[2]][j] );
		if( Y[utilaer[5]][j]  != 0.0 )
		  jac[utilaer[5]][utilaer[5]] -= (   dyhcn * NB[j] / Y[utilaer[5]][j] );
		if( Y[utilaer[6]][j]  != 0.0 )
		  jac[utilaer[6]][utilaer[6]] -= (  dyhc3n * NB[j] / Y[utilaer[6]][j] );
		if( Y[utilaer[13]][j] != 0.0 )
		  jac[utilaer[13]][utilaer[13]] -= (  dynccn * NB[j] / Y[utilaer[13]][j] );
		if( Y[utilaer[14]][j] != 0.0 )
		  jac[utilaer[14]][utilaer[14]] -= ( dych3cn * NB[j] / Y[utilaer[14]][j] );
		if( Y[utilaer[15]][j] != 0.0 )
		  jac[utilaer[15]][utilaer[15]] -= (dyc2h3cn * NB[j] / Y[utilaer[15]][j] );
	      }
     
       
	    /* H -> H2 on haze particles */
	    /* ------------------------- */
	    if( (*htoh2) == 1 )
	      {
		heterohtoh2(CORPS,TEMP,NB,Y,surfhaze,&j,&dyh,&dyh2,utilaer);
		/* dyh <= 0 / 1.0 en adsor., 1 en reac. */

		/* DEBUG 
		   printf("HTOH2 : LAT = %g - J = %d\n",(*LAT),j);
		   if(fabs(dyh*NB[j])>fabs(fp[utilaer[0]]/10.))
		   printf("fp(%s) = %e; dyh  = %e\n",CORPS[utilaer[0]],fp[utilaer[0]],dyh*NB[j]);
		   if(fabs(dyh2*NB[j])>fabs(fp[utilaer[1]]/10.))
		   printf("fp(%s) = %e; dyh2 = %e\n",CORPS[utilaer[1]],fp[utilaer[1]],dyh2*NB[j]);
		*/

		fp[utilaer[0]] += ( dyh  * NB[j] ); 
		/* pourquoi pas *2 ?? cf gptit dans 2da... */

		fp[utilaer[1]] += ( dyh2 * NB[j] );
		if( Y[utilaer[0]][j] != 0.0 )
		  jac[utilaer[0]][utilaer[0]] += ( dyh  * NB[j] / Y[utilaer[0]][j] );
		/* pourquoi pas *2 ?? cf gptit dans 2da... */
	      }


	    /* Backup jacobian level j. */
	    /* ------------------------ */
	    for( i = 0; i <= ST-1; i++ )
	      {
		for( k = 0; k <= ST-1; k++ )
		  a[j][i][k] = jac[i][k];     
		f[i][j] = fp[i] - fl[i] * Y[i][j];
	      }


	    /* finition pour inversion */
	    /* ----------------------- */

	    for( i = 0; i <= ST-1; i++ )
	      {
		for( k = 0; k <= ST-1; k++ )
		  {
		    a[j][i][k] *= ( -THETA * delta );  /* Correction time step. */
		    if( k == i ) a[j][k][k] += NB[j];  /* Correction diagonal. */
		  }
		f[i][j] *= delta;
	      }


	    /* Inversion of matrix cf method LU */
	    /* -------------------------------- */

	    /* inversion by blocs: */
	    /* Hydrocarbons */

	    solve_b( a, f, j, 0, NHC-1 );             
	    for( i = 0; i <= NHC-1; i++ )
	      {
		Y[i][j] += f[i][j];
		if( Y[i][j] <= 1.0e-30 ) Y[i][j] = 0.0e0;
	      }

	    /* Nitriles */

	    solve_b( a, f, j, NHC, ST-1 );             
	    for( i = NHC+1; i <= ST-1; i++ )
	      {
		Y[i][j] += f[i][j];
		if( Y[i][j] <= 1.0e-30 ) Y[i][j] = 0.0e0;
	      }

	    /* end inversion by blocs: */

	    /* ------------------ */
	    /* Tests et evolution */
	    /* ------------------ */

	    /* Calcul deviation */
	    /* ---------------- */

	    for( i = 0; i <= ST-1; i++ )
	      {
		test = 1.0e-15; 
		if( ( Y[i][j] > test ) && ( ym1[i][j] > test ) )
		  {
		    conv = fabs( Y[i][j] - ym1[i][j] ) / ym1[i][j];

		    if( conv > ts )
		      {
			/*
			  if( conv >= 0.1 )
			  {
			  out = fopen( outlog, "a" );
			  fprintf( out, "Latitude: %g; declin:%e;", (*LAT), (*DECLIN) );
			  fprintf(out, " alt:%e; %s %e %e ; %e %e\n",(RA[j]-R0),CORPS[i],ym1[i],Y[i][j],time,delta);
			  fclose( out );
			  }
			*/
			ts = conv;
		      }
		  }
	      }

	    /* test deviation */
	    /* -------------- */

	    if( ts < 0.1e0 )
	      {
		for( i = 0; i <= ST-1; i++ )
		  if( (Y[i][j] >= 0.5e0) && (i != in2) )
		    {
		      out = fopen( outlog, "a" );
		      fprintf( out, "WARNING %s mixing ratio is %e %e at %d\n",
			       CORPS[i], ym1[i][j], Y[i][j], j );
		      for( k = 0; k <= NLEV-1; k++ ) fprintf( out, "%d %e %e\n",k,ym1[i][j],Y[i][k] );
		      fclose( out );
		      //     exit(0); 
		      Y[i][j] = 1.e-20;
		    }
		for( i = 0; i <= NC-1; i++ ) ym1[i][j] = max(Y[i][j],1.e-30);
		time += delta;
		if     (   ts < 1.00e-5 )                      delta *= 1.0e2;
		else if( ( ts > 1.00e-5 ) && ( ts < 1.0e-4 ) ) delta *= 1.0e1;
		else if( ( ts > 1.00e-4 ) && ( ts < 1.0e-3 ) ) delta *= 5.0e0;
		else if( ( ts > 0.001e0 ) && ( ts < 0.01e0 ) ) delta *= 3.0e0;
		else if( ( ts > 0.010e0 ) && ( ts < 0.05e0 ) ) delta *= 1.5e0;
         
		delta = min( deltamax, delta );
	      }
	    else
	      {
		for( i = 0; i <= NC-1; i++ ) Y[i][j] = ym1[i][j];

		if     (   ts > 0.8 )                    delta *= 1.e-6;
		else if( ( ts > 0.6 ) && ( ts <= 0.8 ) ) delta *= 1.e-4;
		else if( ( ts > 0.4 ) && ( ts <= 0.6 ) ) delta *= 1.e-2;
		else if( ( ts > 0.3 ) && ( ts <= 0.4 ) ) delta *= 0.1;
		else if( ( ts > 0.2 ) && ( ts <= 0.3 ) ) delta *= 0.2;
		else if( ( ts > 0.1 ) && ( ts <= 0.2 ) ) delta *= 0.3;
	      }
	    ts = 0.0e0;
	    /*
	      out = fopen( outlog, "a" );
	      fprintf(out, " alt:%e; delta:%e; time:%e; fin:%e\n",(RA[j]-R0),delta,time,(*FIN));
	      fclose( out );
	    */
	  }               

	/* +++++++++++++++++++ */
	/*  end of time loop.  */
	/* +++++++++++++++++++ */

      }  /*  boucle j */


  /*
    ==========================================================================

    FINALISATION: 
    ===============

  /* Niveau de N2 */
  /* ------------ */

  for( j = 0; j <= NLEV-1; j++ ) 
    {
      conv = 0.0e0;
      for( i = 0; i <= ST-1; i++ ) 
	if( i != in2 ) conv += Y[i][j];
      Y[in2][j] = 1. - conv;
    }

  if( (*aerprod) == 1 )
    {
      fdm2d( k_dep, 1, 5, 1 );
      fdm2d(  faer, 1, 5, 1 );
      fdm1d( productaer,  0 );
      fdm1d( mmolaer,  0 );
      fdm1d( csurn, 0 );
      fdm1d( csurh, 0 );
    }

  fdm1d(      fl, 0 );
  fdm1d(      fp, 0 );
  fdm1d(      mu, 0 );
  fdm2d(     ym1, 0,   NC-1, 0 );
  fdm2d(       f, 0,   NC-1, 0 );
  fdm2d(     jac, 0,   NC-1, 0 );
  fdm2d(       c, 0, NLEV-1, 0 );
  fdm3d(       a, 0, NLEV-1, 0,   NC-1, 0 );
  fdm3d(       b, 0, NLEV-1, 0,   NC-1, 1 );
}
