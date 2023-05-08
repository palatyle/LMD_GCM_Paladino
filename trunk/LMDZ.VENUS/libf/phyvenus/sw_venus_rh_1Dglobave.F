      SUBROUTINE SW_venus_rh_1Dglobave(PRMU0, PFRAC,
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
c     fluxes taken from Rainer Haus calculations for Venus.
c     Ref: Haus et al. 2016
C
C     AUTHOR.
C     -------
C        Sebastien Lebonnois
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 5/2016
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
      integer nlrh,nszarh,nlatrh
      parameter (nlrh=118)  ! fichiers Rainer Haus
      parameter (nszarh=7) ! fichiers Rainer Haus
      parameter (nlatrh=19) ! fichiers Rainer Haus
      
      integer i,j,lat,nsza,nsza0,nl0
      real   zsnetmoy(nlrh+1,nlatrh)    ! net solar flux (W/m**2) (+ vers bas)
      real   presrh(nlrh+1)             ! pressure in table (bar)
      real   altrh(nlrh+1)              ! altitude in table (km)
      real   latrh(nlatrh)              ! latitude in table (degrees)
      real   zsolnet(nlrh+1)            ! for mean net solar flux in RH
      character*22 nullchar
      real   factflux
      real   zsnet(nszarh)    ! net solar flux (W/m**2) (+ vers bas)
      real   deltalat
      logical firstcall
      data    firstcall/.true./
      save   zsolnet,altrh,presrh
      save   firstcall
      
c ------------------------
c Loading the file
c ------------------------

      if (firstcall) then

       open(11,file='SolarNetFlux_RH.dat')

       do i=1,nlrh+1
          read(11,'(E5.1,4x,F8.2)') altrh(i),presrh(i)
       enddo

       do lat=1,nlatrh
         latrh(lat)=5.*(lat-1)
         read(11,*) nullchar
         read(11,*) nullchar
         read(11,*) nullchar
         read(11,*) nullchar

         do i=1,nlrh+1
          read(11,'(E6.1,7(2x,F11.5),7x,F11.5)')
     .          altrh(i),zsnet,zsnetmoy(i,lat)
         enddo
         read(11,*) nullchar
       enddo
       latrh(nlatrh)=89.

       close(11)

c ----------- TEST ------------
c      Moyenne planetaire
c -----------------------------

      zsolnet=0.
      do lat=1,nlatrh-1
        deltalat=(latrh(lat+1)-latrh(lat))*RPI/180.
        do j=1,nlrh+1
        zsolnet(j) = zsolnet(j)+
     .         (zsnetmoy(j,lat+1)+zsnetmoy(j,lat))/2.*
     .       deltalat*cos((latrh(lat+1)+latrh(lat))*RPI/360.)
        enddo
      enddo
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
        nl0 = nlrh
        do i=nlrh+1,2,-1
           if (presrh(i).ge.PPB(j)) then
                nl0 = i-1
           endif
        enddo
        
        factflux = (log10(max(PPB(j),presrh(1)))-log10(presrh(nl0+1)))
     .            /(log10(presrh(nl0))-log10(presrh(nl0+1)))

        ZFSNET(j) =  factflux     *zsolnet(nl0)
     .             + (1.-factflux)*zsolnet(nl0+1)
        
c-----TEST-------
c tayloring the solar flux...
c        if ((PPB(j).gt.0.236).and.(PPB(j).le.22.52)) then
c         ZFSNET(j) = ZFSNET(j)+2.5*(1.+cos((log10(PPB(j)/3.5)/
c     .                                  log10(0.236/3.5))*RPI))
c        endif
c----------------
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
c-----TEST-------
c tayloring the solar flux...
        if ((PPB(j).gt.1.4).and.(PPB(j).le.10.)) then
          PHEAT(j) = PHEAT(j)*3
c       elseif ((PPB(j).gt.10.).and.(PPB(j).le.30.)) then
c         PHEAT(j) = PHEAT(j)*1.5
        endif
c----------------
c       print*,PPB(j),ZFSNET(j),PHEAT(j)
      enddo
c     stop
      return
      end

