      SUBROUTINE SW_venus_dc_1Dglobave(PRMU0, PFRAC, 
     S              PPB, pt,
     S              PHEAT, 
     S              PTOPSW,PSOLSW,ZFSNET)
      
      use dimphy
      use cpdet_phy_mod, only: cpdet
      IMPLICIT none

#include "YOMCST.h"
C
C     ------------------------------------------------------------------
C
C     PURPOSE.
C     --------
C
c      this routine loads and interpolates the shortwave radiation
c     fluxes taken from Dave Crisp calculations for Venus.
c     Ref: Crisp 1986.
C
C     AUTHOR.
C     -------
C        Sebastien Lebonnois
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 27/07/2005
c        L.Salmi  : june 2013 astuce to reduce the excess of  NIR 
c                   in the transition region LTE/LTE
c
c        G.Gilli  : feb  2014          
C     ------------------------------------------------------------------
C
C* ARGUMENTS:
C
c inputs

      REAL   PRMU0  ! COSINE OF ZENITHAL ANGLE
      REAL   PFRAC  ! fraction de la journee
      REAL   PPB(klev+1)  ! inter-couches PRESSURE (bar)
      REAL   pt(klev)     ! mid-layer temperature
C
c output

      REAL   PHEAT(klev)  ! SHORTWAVE HEATING (K/s) within each layer
      REAL   PTOPSW       ! SHORTWAVE FLUX AT T.O.A. (net)
      REAL   PSOLSW       ! SHORTWAVE FLUX AT SURFACE (net)
      REAL   ZFSNET(klev+1) ! net solar flux at ppb levels

C
C* LOCAL VARIABLES:
C
      integer nldc,nszadc
      parameter (nldc=49)  ! fichiers Crisp
      parameter (nszadc=8) ! fichiers Crisp
      
      integer i,j,nsza,nsza0,nl0
      real   solarrate               ! solar heating rate (K/earthday)
      real   zsnet(nldc+1,nszadc)    ! net solar flux (W/m**2) (+ vers bas)
      real   zsdn,zsup               ! downward/upward solar flux (W/m**2)
      real   solza(nszadc)           ! solar zenith angles in table
      real   presdc(nldc+1)          ! pressure levels in table (bar)
      real   tempdc(nldc+1)          ! temperature in table (K)
      real   altdc(nldc+1)           ! altitude in table (km)
      real   coolrate                ! IR heating rate (K/earthday) ?
      real   totalrate               ! total rate (K/earthday)
      real   zldn                    ! downward IR flux (W/m**2) ?
      real   zlup                    !   upward IR flux (W/m**2) ?
      real   zsolnet(nldc+1)         ! for testing mean net solar flux in DC
      character*22 nullchar
      real   deltasza
      real   sza0,factflux
      logical firstcall
      data    firstcall/.true./
      save   solza,zsnet,presdc,tempdc,altdc,zsolnet
      save   firstcall
      
c ------------------------
c Loading the file
c ------------------------

      if (firstcall) then

       open(11,file='dataDCrisp.dat')
       read(11,*) nullchar
      
       do nsza=1,nszadc
        read(11,*) nullchar
        read(11,*) nullchar
        read(11,*) nullchar
        read(11,'(22x,F11.5)') solza(nsza)
        read(11,*) nullchar
        read(11,*) nullchar
        read(11,*) nullchar
        read(11,'(3(2x,F10.4),36x,4(2x,F11.5))')
     .          presdc(nldc+1),tempdc(nldc+1), altdc(nldc+1),
     .          zsdn,zsup,zldn,zlup
        zsnet(nldc+1,nsza)=zsdn-zsup
        do i=1,nldc
           j = nldc+1-i        ! changing: vectors from surface to top
           read(11,'(6(2x,F10.4),4(2x,F11.5))') 
     .          presdc(j),tempdc(j),altdc(j),
     .          solarrate,coolrate,totalrate,
     .          zsdn,zsup,zldn,zlup
           zsnet(j,nsza)=zsdn-zsup
        enddo
       enddo

       close(11)

c ----------- TEST ------------
c      Moyenne planetaire
c -----------------------------
      
      deltasza=(solza(2)-solza(1))*RPI/180.

      do j=1,nldc+1
        zsolnet(j) = zsnet(j,1)*deltasza*deltasza/16.
        do nsza=2,nszadc
        zsolnet(j) = zsolnet(j)+zsnet(j,nsza)*0.5*deltasza*
     .             sin(solza(nsza)*RPI/180.)
        enddo
c overestimation:
        zsolnet(j) = zsolnet(j)*0.84  
c        print*,j,altdc(j),zsolnet(j)
      enddo
c      stop
c -----------------------------
c --------  FIN TEST ----------

       firstcall=.false.
      endif

c --------------------------------------
c Interpolation in the GCM vertical grid
c --------------------------------------

c Pressure levels
c ---------------

      do j=1,klev+1
        nl0 = 2
        do i=1,nldc
           if (presdc(i).ge.PPB(j)) then
                nl0 = i+1
           endif
        enddo
        
        factflux = (log10(max(PPB(j),presdc(nldc+1)))
     .                          -log10(presdc(nl0-1)))
     .            /(log10(presdc(nl0))-log10(presdc(nl0-1)))
        ZFSNET(j) =  factflux     *zsolnet(nl0)
     .             + (1.-factflux)*zsolnet(nl0-1)
        
      enddo

      PTOPSW = ZFSNET(klev+1)
      PSOLSW = ZFSNET(1) 
      
c Heating rates
c -------------
c On utilise le gradient du flux pour calculer le taux de chauffage:
c   heat(K/s) = d(fluxnet)  (W/m2)
c              *g           (m/s2)
c              /(-dp)  (epaisseur couche, en Pa=kg/m/s2)
c              /cp  (J/kg/K) 

      do j=1,klev
! ADAPTATION GCM POUR CP(T)
        PHEAT(j) = (ZFSNET(j+1)-ZFSNET(j))
     .            *RG/cpdet(pt(j)) / ((PPB(j)-PPB(j+1))*1.e5)
c--------------
c BIDOUILLE POUR AJUSTEMENT ET TEST
c       if (PPB(j).lt.1.e-2) then
c         PHEAT(j) = PHEAT(j)*0.3
c       endif
c       if ((PPB(j).gt.1.e-2).and.(PPB(j).lt.2.e-1)) then
c         PHEAT(j) = PHEAT(j)*0.7
c       endif
c       if (PPB(j).gt.1.) then
c         PHEAT(j) = PHEAT(j)*2.
c       endif
c--------------
      enddo

      return
      end

