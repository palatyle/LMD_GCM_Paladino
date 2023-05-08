      SUBROUTINE nirco2abs(ngrid,nlayer,pplay,dist_sol,
     $     mu0,fract,declin,pdtnirco2)
                                                   
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing heating rate due to
c   absorption by CO2 in the near-infrared
c   This version includes NLTE effects
c
c   (Scheme to be described in Forget et al., JGR, 2003)
c   (old Scheme described in Forget et al., JGR, 1999)
c
c   This version updated with a new functional fit,
c   see NLTE correction-factor of Lopez-Valverde et al (1998)
c   Stephen Lewis 2000
c  
c   08/2002 : correction for bug when running with diurnal=F
c
c   author:  Frederic Hourdin 1996 
c   ------
c            Francois Forget 1999 
c
c   input:
c   ----- 
c   ngrid                 number of gridpoint of horizontal grid
c   nlayer                Number of layer
c   dist_sol              sun-Mars distance (AU)
c   mu0(ngridmx)          
c   fract(ngridmx)        day fraction of the time interval
c   declin                latitude of subslar point
c
c   output:
c   -------
c
c   pdtnirco2(ngrid,nlayer)      Heating rate (K/s)
c
c
c=======================================================================
c
c    0.  Declarations :
c    ------------------
c
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "comdiurn.h"


c-----------------------------------------------------------------------
c    Input/Output
c    ------------
      INTEGER ngrid,nlayer

      REAL pplay(ngrid,nlayer)
      REAL dist_sol,mu0(ngridmx),fract(ngridmx),declin

      REAL pdtnirco2(ngrid,nlayer)
c
c    Local variables :
c    -----------------
      INTEGER l,ig, n, nstep
      REAL co2heat0, zmu(ngridmx)

c     special diurnal=F 
      real mu0_int(ngridmx),fract_int(ngridmx),zday_int
      real ztim1,ztim2,ztim3,step

c
c   local saved variables
c   ---------------------

c     p0noonlte is a pressure below which non LTE effects are significant.
c     REAL p0nonlte
c     DATA p0nonlte/7.5e-3/
c     SAVE p0nonlte

c     parameters for CO2 heating fit
      real n_a, n_p0, n_b
      parameter (n_a=1.1956475)
      parameter (n_b=1.9628251)
      parameter (n_p0=0.0015888279)

c----------------------------------------------------------------------

c     Initialisation
c     --------------
c     co2heat is the heating by CO2 at 700Pa for a zero zenithal angle.
      co2heat0=n_a*(1.52/dist_sol)**2/daysec

c     Simple calcul for a given sun incident angle (if diurnal=T)
c     --------------------------------------------

      IF (diurnal) THEN  
         do ig=1,ngrid
            zmu(ig)=sqrt(1224.*mu0(ig)*mu0(ig)+1.)/35.
         enddo
         do l=1,nlayer
           do ig=1,ngrid
             if(fract(ig).gt.0.) pdtnirco2(ig,l)=
     &             co2heat0*sqrt((700.*zmu(ig))/pplay(ig,l))
     &             /(1.+n_p0/pplay(ig,l))**n_b

c           OLD SCHEME (forget et al. 1999)
c    s           co2heat0*sqrt((700.*zmu(ig))/pplay(ig,l))
c    s          / (1+p0nonlte/pplay(ig,l))
           enddo
         enddo

c     Averaging over diurnal cycle (if diurnal=F)
c     -------------------------------------------
c     NIR CO2 abs is slightly non linear. To remove the diurnal
c     cycle, it is better to average the heating rate over 1 day rather
c     than using the mean mu0 computed by mucorr in physiq.F (FF, 1998)

      ELSE      ! if (.not.diurnal) then

         nstep = 20   ! number of integration step /sol
         do n=1,nstep
            zday_int = (n-1)/float(nstep)
            ztim2=COS(declin)*COS(2.*pi*(zday_int-.5))
            ztim3=-COS(declin)*SIN(2.*pi*(zday_int-.5))
            CALL solang(ngrid,sinlon,coslon,sinlat,coslat,
     s             ztim1,ztim2,ztim3,
     s             mu0_int,fract_int)
            do ig=1,ngrid
               zmu(ig)=sqrt(1224.*mu0_int(ig)*mu0_int(ig)+1.)/35.
            enddo
            do l=1,nlayer
               do ig=1,ngrid
                  if(fract_int(ig).gt.0.) pdtnirco2(ig,l)=
     &                 pdtnirco2(ig,l) + (1/float(nstep))*
     &                 co2heat0*sqrt((700.*zmu(ig))/pplay(ig,l))
     &                 /(1.+n_p0/pplay(ig,l))**n_b
               enddo
            enddo
         end do
      END IF

      return
      end 

