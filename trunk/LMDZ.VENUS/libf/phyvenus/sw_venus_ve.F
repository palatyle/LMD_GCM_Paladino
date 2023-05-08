      SUBROUTINE SW_venus_ve( PRMU0, PFRAC, 
     S              PPB, pt, pz,
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
c     fluxes and heating rates computed from Vincent Eymet 3D MC code
C
C     AUTHOR.
C     -------
C        Sebastien Lebonnois
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 06/2014
C     ------------------------------------------------------------------
C
C* ARGUMENTS:
C
c inputs

      REAL   PRMU0  ! COSINE OF ZENITHAL ANGLE
      REAL   PFRAC  ! fraction de la journee
      REAL   PPB(klev+1)  ! inter-couches PRESSURE (bar)
      REAL   pt(klev)     ! mid-layer temperature
      REAL   pz(klev+1)   ! inter-couches altitude (m)
C
c output

      REAL   PHEAT(klev) ! SHORTWAVE HEATING (K/VENUSDAY) within each layer
      REAL   PTOPSW       ! SHORTWAVE FLUX AT T.O.A. (net)
      REAL   PSOLSW       ! SHORTWAVE FLUX AT SURFACE (net)
      REAL   ZFSNET(klev+1) ! net solar flux at ppb levels

C
C* LOCAL VARIABLES:
C
      integer nlve,nszave
      parameter (nlve=78)   ! fichiers planet_EMC
      parameter (nszave=20) ! fichiers planet_EMC
      
      integer i,j,nsza,nsza0,nl0
      real   solarrate               ! solar heating rate (K/earthday)
      real   zsnet(nlve,nszave)      ! net solar flux (W/m**2) (+ vers bas)
      real   zheat(nlve-1,nszave)    ! rad budget (W/m**2)
      real   zsdn,zsup               ! downward/upward solar flux (W/m**2)
      real   solza(nszave)           ! solar zenith angles in table (rad)
      real   altve(nlve)             ! altitude in table (m)
      real   zsolnet(nlve)           ! for testing mean net solar flux
      character*22 nullchar
      real   sza0,factsza,factflux,alt
      logical firstcall
      data    firstcall/.true./
      save   solza,zsnet,altve,zheat
      save   firstcall
      
c ------------------------
c Loading the files
c ------------------------

      if (firstcall) then

! FLUXES (W/m2)

       open(11,file='solar_fluxes_GCM.dat')
       read(11,*) nullchar
       read(11,*) nullchar
       read(11,*) nullchar
       read(11,*) nullchar
      
       do nsza=1,nszave
        read(11,*) nullchar
        read(11,*) solza(nsza)
        read(11,*) nullchar
        read(11,*) nullchar
        do j=1,nlve
           read(11,'(4(2x,F12.5))') 
     .          altve(j),zsdn,zsup,zsnet(j,nsza)
        enddo
       enddo

       close(11)

! HEATING RATES (W/m2)

       open(12,file='solar_budgets_GCM.dat')
       read(12,*) nullchar
       read(12,*) nullchar
       read(12,*) nullchar
       read(12,*) nullchar
      
       do nsza=1,nszave
        read(12,*) nullchar
        read(12,*) solza(nsza)
        read(12,*) nullchar
        read(12,*) nullchar
        do j=1,nlve-1
           read(12,'(2(2x,F12.5))') 
     .          alt,zheat(j,nsza)
        enddo
       enddo

       close(12)

       firstcall=.false.
      endif

c --------------------------------------
c Interpolation in the GCM vertical grid
c --------------------------------------

c Zenith angle
c ------------
      
      sza0 = acos(PRMU0)  ! in radians
c        print*,'Angle Zenithal =',sza0,' PFRAC=',PFRAC

      nsza0=1
      do nsza=1,nszave
         if (solza(nsza).le.sza0) then
              nsza0 = nsza+1
         endif
      enddo
      
      if ((nsza0.ne.1).and.(nsza0.ne.nszave+1)) then
          factsza = (sza0-solza(nsza0-1))/(solza(nsza0)-solza(nsza0-1))
      endif

c Pressure levels
c ---------------

      do j=1,klev+1
        nl0 = 2
        do i=1,nlve-1
           if (altve(i).le.pz(j)) then
                nl0 = i+1
           endif
        enddo
        
        factflux = (min(pz(j),altve(nlve))
     .                          -altve(nl0-1))
     .            /(altve(nl0)-altve(nl0-1))

! FLUXES

        ZFSNET(j) = 0.
        if ((nsza0.ne.1).and.(nsza0.ne.nszave+1)) then
          ZFSNET(j) =  factflux   *  factsza   *zsnet(nl0,nsza0)
     .             +   factflux   *(1.-factsza)*zsnet(nl0,nsza0-1)
     .             + (1.-factflux)*  factsza   *zsnet(nl0-1,nsza0)
     .             + (1.-factflux)*(1.-factsza)*zsnet(nl0-1,nsza0-1)
        else if (nsza0.eq.1) then
          ZFSNET(j) =  factflux   *zsnet(nl0,1)
     .             + (1.-factflux)*zsnet(nl0-1,1)
        endif
        ZFSNET(j) = ZFSNET(j)*PFRAC

! HEATING RATES

        if (j.ne.klev+1) then
          PHEAT(j) = 0.

          if ((nsza0.ne.1).and.(nsza0.ne.nszave+1)) then
            PHEAT(j) =  factflux   *  factsza   *zheat(nl0,nsza0)
     .             +   factflux   *(1.-factsza)*zheat(nl0,nsza0-1)
     .             + (1.-factflux)*  factsza   *zheat(nl0-1,nsza0)
     .             + (1.-factflux)*(1.-factsza)*zheat(nl0-1,nsza0-1)
          else if (nsza0.eq.1) then
            PHEAT(j) =  factflux   *zheat(nl0,1)
     .              + (1.-factflux)*zheat(nl0-1,1)
          endif
          PHEAT(j) = PHEAT(j)*PFRAC
        endif

      enddo

      PTOPSW = ZFSNET(klev+1)
      PSOLSW = ZFSNET(1) 
      
c Heating rates
c -------------
c Conversion from W/m2 to K/s:
c   heat(K/s) = d(fluxnet)  (W/m2)
c              *g           (m/s2)
c              /(-dp)  (epaisseur couche, en Pa=kg/m/s2)
c              /cp  (J/kg/K) 

      do j=1,klev
! ADAPTATION GCM POUR CP(T)
        PHEAT(j) = PHEAT(j)
     .            *RG/cpdet(pt(j)) / ((PPB(j)-PPB(j+1))*1.e5)
      enddo

      return
      end

