      SUBROUTINE SW_venus_dc(PRMU0, PFRAC, 
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
      character*22 nullchar
      real   sza0,factsza,factflux
      logical firstcall
      data    firstcall/.true./
      save   solza,zsnet,presdc,tempdc,altdc
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

       firstcall=.false.
      endif

c --------------------------------------
c Interpolation in the GCM vertical grid
c --------------------------------------

c Zenith angle
c ------------
      
      sza0 = acos(PRMU0)/3.1416*180.
c        print*,'Angle Zenithal =',sza0,' PFRAC=',PFRAC

      do nsza=1,nszadc
         if (solza(nsza).le.sza0) then
              nsza0 = nsza+1
         endif
      enddo
      
      if (nsza0.ne.nszadc+1) then
          factsza = (sza0-solza(nsza0-1))/(solza(nsza0)-solza(nsza0-1))
      else
          factsza = min((sza0-solza(nszadc))/(90.-solza(nszadc)), 1.)
      endif

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
c       factflux = (max(PPB(j),presdc(nldc+1))-presdc(nl0-1))
c    .            /(presdc(nl0)-presdc(nl0-1))
        if (nsza0.ne.nszadc+1) then
          ZFSNET(j) =  factflux   *  factsza   *zsnet(nl0,nsza0)
     .             +   factflux   *(1.-factsza)*zsnet(nl0,nsza0-1)
     .             + (1.-factflux)*  factsza   *zsnet(nl0-1,nsza0)
     .             + (1.-factflux)*(1.-factsza)*zsnet(nl0-1,nsza0-1)
        else
          ZFSNET(j) =  factflux   *(1.-factsza)*zsnet(nl0,nsza0-1)
     .             + (1.-factflux)*(1.-factsza)*zsnet(nl0-1,nsza0-1)
        endif
        
        ZFSNET(j) = ZFSNET(j)*PFRAC

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
      enddo

      return
      end

