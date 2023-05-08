!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/radlwsw.F,v 1.2 2004/10/27 10:14:46 lmdzadmin Exp $
!
      SUBROUTINE radlwsw_newtoncool(dist, rmu0, fract, zzlev,
     .                  paprs, pplay,tsol, pt)
      
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19960719
c Objet: interface entre le modele et les rayonnements
c Arguments:
c dist-----input-R- distance astronomique terre-soleil
c rmu0-----input-R- cosinus de l'angle zenithal
c fract----input-R- duree d'ensoleillement normalisee
c paprs----input-R- pression a inter-couche (Pa)
c pplay----input-R- pression au milieu de couche (Pa)
c tsol-----input-R- temperature du sol (en K)
c pt-------input-R- temperature (K)
c
      
c   S. Lebonnois    12/04/2007
c  VERSION NEWTONIAN COOLING pour Venus (no diurnal cycle)
c  update 01/2014

c======================================================================
      use dimphy
      USE geometry_mod, ONLY: latitude ! in radians
      USE phys_state_var_mod, only: heat,cool,radsol,
     .  topsw,toplw,solsw,sollw,sollwdown,lwnet,swnet,zt_eq
      USE write_field_phy
      IMPLICIT none
#include "YOMCST.h"
#include "clesphys.h" 

c ARGUMENTS
      real rmu0(klon), fract(klon), dist
 
      real zzlev(klon,klev+1),paprs(klon,klev+1), pplay(klon,klev)
      real tsol(klon)
      real pt(klon,klev)
 
c LOCAL VARIABLES
      INTEGER i,j,k
      integer    nlevCLee,level
      parameter (nlevCLee=30)
      REAL   pressCLee(nlevCLee+1),tempCLee(nlevCLee+1)
      real   dt_epCLee(nlevCLee+1),etaCLee(nlevCLee+1)
      real   tauCLee
      parameter (tauCLee=25*86400) ! en s
      real   ztemp,zdt,fact
      real   dTsdt(klev)
      
      data     etaCLee/9.602e-1,8.679e-1,7.577e-1,6.420e-1,5.299e-1, 
     .                 4.273e-1,3.373e-1,2.610e-1,1.979e-1,1.472e-1,
     .                 1.074e-1,7.672e-2,5.361e-2,3.657e-2,2.430e-2,
     .                 1.569e-2,9.814e-3,5.929e-3,3.454e-3,1.934e-3,
     .                 1.043e-3,5.400e-4,2.710e-4,1.324e-4,6.355e-5,
     .                 3.070e-5,1.525e-5,7.950e-6,4.500e-6,2.925e-6,
     .                 2.265e-6/
      data   tempCLee/728.187,715.129,697.876,677.284,654.078,628.885,
     .                602.225,574.542,546.104,517.339,488.560,459.932,
     .                431.741,404.202,377.555,352.042,327.887,305.313,
     .                284.556,265.697,248.844,233.771,220.368,208.247,
     .                197.127,187.104,178.489,171.800,167.598,165.899,
     .                165.676/
      data   dt_epCLee/6.101 , 6.136 , 6.176 , 6.410 , 6.634 , 6.678 ,
     .                 6.719 , 6.762 , 7.167 , 7.524 , 9.840 ,14.948 ,
     .                21.370 ,28.746 ,36.373 ,43.315 ,48.534 ,51.175 ,
     .                50.757 ,47.342 ,41.536 ,34.295 ,26.758 ,19.807 ,
     .                14.001 , 9.599 , 6.504 , 4.439 , 3.126 , 2.370 , 
     .                2.000/
c

      logical firstcall
      data    firstcall/.true./
      save    firstcall
      
c  Initialisations
c-----------------

      if (firstcall) then

        PRINT*,"******* ATTENTION, NEWTONIAN COOLING ********"

        pressCLee = etaCLee * 9.2e6

        DO i = 1, klon
	
	  do k = 1,klev
	  
            level = 1
	    do j=1,nlevCLee
	      if (pressCLee(j).gt.pplay(i,k)) level = j
	    enddo
	    
	    fact  = (log10(pplay(i,k))-log10(pressCLee(level)))
     .        /(log10(pressCLee(level+1))-log10(pressCLee(level)))
	    ztemp = tempCLee(level)*(1-fact)+tempCLee(level+1)*fact
	    zdt   = dt_epCLee(level)*(1-fact)+dt_epCLee(level+1)*fact
c           zt_eq(i,k) = ztemp + zdt*(cos(latitude(i))-2./RPI)
	    zt_eq(i,k) = ztemp + zdt*(cos(latitude(i))-RPI/4.)
	    
	  enddo
	  
        ENDDO !i

      endif ! firstcall
      
c+++++++ BOUCLE SUR LA GRILLE +++++++++++++++++++++++++
      DO j = 1,klon
 
	  do k = 1,klev
	     dTsdt(k) = -(pt(j,k)-zt_eq(j,k))/tauCLee   ! en K/s
	  enddo
       
         radsol(j) = 0.           ! + vers bas
         topsw(j) = 0.            ! + vers bas
         toplw(j) = 0.            ! + vers haut
         solsw(j) = 0.            ! + vers bas
         sollw(j) = 0.            ! + vers bas
         sollwdown(j) = 0.        ! + vers bas

        DO k = 1, klev+1
         lwnet  (j,k)   = 0.
         swnet  (j,k)   = 0.
        ENDDO

        DO k = 1, klev
         heat (j,k) = dTsdt(k)    ! K/s
         cool (j,k) = 0.
        ENDDO
c
      ENDDO !j
c+++++++ FIN BOUCLE SUR LA GRILLE +++++++++++++++++++++++++

      firstcall = .false.
      RETURN
      END

