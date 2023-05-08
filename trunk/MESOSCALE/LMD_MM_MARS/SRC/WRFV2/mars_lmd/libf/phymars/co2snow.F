      SUBROUTINE co2snow (ngrid,nlayer,ptimestep,
     &          emisref,condsub,pplev,
     &        pcondicea,pcondices,pfallice,pemisurf)

       IMPLICIT NONE


c*****WRF: ligne trop longue

c=======================================================================
c     Program for simulate the impact of the CO2 snow fall on
c     the surface infrared emission (emissivity)  and on 
c     the airborne dust
c     F.Forget 1996
c=======================================================================

c Declarations
c ------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "surfdat.h"
#include "callkeys.h"

c     input
c     -----

      INTEGER ngrid,nlayer
      REAL ptimestep , emisref(ngrid) ! grd or ice  emissivity without snow
      logical condsub(ngrid)

      REAL pplev(ngrid,nlayer+1)
      REAL pcondicea(ngrid,nlayer), pcondices(ngrid)
      REAL pfallice(ngrid,nlayer+1)

c     output
c     ------

      REAL pemisurf(ngrid)

c     local
c     -----
      integer ig , l , icap

      REAL zdemisurf ,dtemis
      REAL sumdaer

c     saved
c     -----
      REAL Kscat(2), scaveng
      LOGICAL firstcall
      SAVE firstcall
      save Kscat , scaveng
      DATA firstcall/.true./

c --------------
c Initialisation
c --------------

      if (firstcall) then

c   Kscat : coefficient for decreasing the surface emissivity
c           =(0.001/3.)*alpha/iceradius ,
c           with 0.3< alpha < 0.6, set to 0.45 (coeff from emis = f (tau))
c           and iceradius the mean radius of the
c           scaterring particles  (200.e-6<iceradius<10.e-6)

            Kscat(1)=(0.001/3.)*0.45/iceradius(1)
            Kscat(2)=(0.001/3.)*0.45/iceradius(2)

c          Scavenging Ratio (dust concentration in the air / in the snow)
           scaveng = 100.0
           
c          Collision Scavenging coefficient  (m2.kg-1)
c          Csca  = 2.3  ! not used yet !!!!!!!!!!!
           firstcall = .false.

      end if


c     LOOP on grid points
c     -------------------
      do ig=1,ngrid 
         if (condsub(ig)) then

c         IF (scavenging) then
c          Airborne Dust 
c          -------------
c          sumdaer=0.
c          do l=nlayer, 1, -1
c             pdaerosol(ig,l)= -paerosol(ig,l,1)*
c    &              (1-exp(-scaveng*pcondicea(ig,l)*ptimestep*g/
c    &               (pplev(ig,l)-pplev(ig,l+1))))/ptimestep  

c    &           - Csca*paerosol(ig,l,1) ! Scavenging by collision
c    &           * 0.5*(pfallice(ig,l)) ! not included

c             test to avoid releasing to much dust when subliming:
c             if(pdaerosol(ig,l).gt.-sumdaer)pdaerosol(ig,l)=-sumdaer 
c             sumdaer=sumdaer + pdaerosol(ig,l)

c            if (-pdaerosol(ig,l)*ptimestep.gt.paerosol(ig,l,1)) then
c                write(*,*) 'ds co2snow: aerosol < 0.0 !!!'
c                write(*,*) 'ig =' , ig
c            end if
c          end do
c         END IF

c          Surface emissivity
c          ------------------
c   dtemis: Time scale for increasing the ice emissivity

           IF(ig.GT.ngrid/2+1) THEN
              icap=2
           ELSE
              icap=1
           ENDIF

           zdemisurf = 
     &    (emisref(ig)-pemisurf(ig))/(dtemisice(icap)*daysec) 
c Using directly the diferential equation:
c    &       -Kscat(icap)*emisref(ig)*
c    &      (pemisurf(ig)/emisref(ig))**4 *pfallice(ig,1) 
c Using an integrated form for numerical safety instead
     & +(emisref(ig)* ((pemisurf(ig)/emisref(ig))**(-3)+3.*Kscat(icap)*
     & pfallice(ig,1)*ptimestep)**(-1/3.) -pemisurf(ig))/ptimestep   


           pemisurf(ig) = pemisurf(ig) + zdemisurf*ptimestep 

           if (pemisurf(ig).lt.0.1) then
                 write(*,*) 'ds co2snow: emis < 0.1 !!!'
                 write(*,*) 'ig =' , ig
                 write(*,*)'pemisurf(ig)',pemisurf(ig)
                 write(*,*) 'zdemisurf*ptimestep',zdemisurf*ptimestep
           end if
         else
           pemisurf(ig) = emisref(ig)
         end if
      end do

      return
      end 
