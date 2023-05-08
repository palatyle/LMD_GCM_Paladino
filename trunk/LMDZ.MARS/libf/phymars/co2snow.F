      SUBROUTINE co2snow (ngrid,nlayer,ptimestep,emisref,condsub
     &         ,pplev,pcondicea,pcondices,pfallice,pemisurf)

      use surfdat_h, only: iceradius, dtemisice
      use geometry_mod, only: latitude ! grid point latitudes (rad)
      use time_phylmdz_mod, only: daysec
      IMPLICIT NONE

c=======================================================================
c     Program for simulate the impact of the CO2 snow fall on
c     the surface infrared emission (emissivity)  and on 
c     the airborne dust
c     F.Forget 1996
c=======================================================================

c Declarations
c ------------

#include "callkeys.h"

c     Arguments
c     ---------

      INTEGER,INTENT(IN) :: ngrid ! number of atmospheric columns
      INTEGER,INTENT(IN) :: nlayer ! number of atmospheric layers
      REAL,INTENT(IN) :: ptimestep ! timestep of the physics (s)
      REAL,INTENT(IN) :: emisref(ngrid) ! grd or ice  emissivity without snow
      logical,intent(in) :: condsub(ngrid) ! true if there is CO2 condensation
                                           ! or sublimation in the column
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1) ! interlayer pressure (Pa)
      REAL,INTENT(IN) :: pcondicea(ngrid,nlayer) ! CO2 condensation rate (kg/m2/s)
      REAL,INTENT(IN) :: pcondices(ngrid) ! CO2 condensation rate (kg/m2/s)
                                          ! on the surface
      REAL,INTENT(IN) :: pfallice(ngrid,nlayer+1) ! falling CO2 ice (kg/m2/s)

      REAL,INTENT(OUT) :: pemisurf(ngrid) ! surface emissivity

c     local
c     -----
      integer ig , l , icap

      REAL zdemisurf ,dtemis
      REAL sumdaer

c     saved
c     -----
      REAL,save :: Kscat(2), scaveng
      LOGICAL,SAVE :: firstcall=.true.

c --------------
c Initialisation
c --------------

      ! AS: firstcall OK absolute
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

           IF(latitude(ig).LT. 0.) THEN
              icap=2 ! Southern hemisphere
           ELSE
              icap=1 ! Northern Hemisphere
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
