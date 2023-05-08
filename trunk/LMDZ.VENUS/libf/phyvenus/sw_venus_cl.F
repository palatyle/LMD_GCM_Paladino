      SUBROUTINE SW_venus_cl(PRMU0, PFRAC, 
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
c     fluxes taken from Chris Lee calculations for Venus.
c     Ref: Lee and Richardson 2011
C
C     AUTHOR.
C     -------
C        Sebastien Lebonnois
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 11/2014
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
      integer nlcl,nszacl
      parameter (nlcl=80)  ! fichiers Crisp
      parameter (nszacl=18) ! fichiers Crisp
      
      integer i,j,nsza,nsza0,nl0
      real   solarrate               ! solar heating rate (K/earthday)
      real   zsnet(nlcl+1,nszacl)    ! net solar flux (W/m**2) (+ vers bas)
      real   zsdn,zsup               ! downward/upward solar flux (W/m**2)
      real   solza(nszacl)           ! solar zenith angles in table
      real   prescl(nlcl+1)          ! pressure levels in table (bar)
      real   tempcl(nlcl+1)          ! temperature in table (K)
      real   altcl(nlcl+1)           ! altitude in table (km)
      real   coolrate                ! IR heating rate (K/earthday) ?
      real   totalrate               ! total rate (K/earthday)
      character*22 nullchar
      real   sza0,factsza,factflux
      real   zlnet,tmpzsnet(nszacl)
      logical firstcall
      data    firstcall/.true./
      save   solza,zsnet,prescl,tempcl,altcl
      save   firstcall
      
c ------------------------
c Loading the file
c ------------------------

      if (firstcall) then

       do nsza=1,nszacl
          solza(nsza)=(nsza-1)*5.
       enddo
       
       open(11,file='CLee-SW.dat')
       read(11,*) nullchar
       
       do i=1,nlcl+1
        read(11,'(4(F10.4,1x),18(F11.4,1x))')
     .          altcl(i),prescl(i),tempcl(i),zlnet,tmpzsnet
c change of sign convention:
        zsnet(i,:)=tmpzsnet*(-1.)
	prescl(i)=prescl(i)*1.e-5 ! conversion to bars...
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

      do nsza=1,nszacl
         if (solza(nsza).le.sza0) then
              nsza0 = nsza+1
         endif
      enddo
      
      if (nsza0.ne.nszacl+1) then
          factsza = (sza0-solza(nsza0-1))/(solza(nsza0)-solza(nsza0-1))
      else
          factsza = min((sza0-solza(nszacl))/(90.-solza(nszacl)), 1.)
      endif

c Pressure levels
c ---------------

      do j=1,klev+1
        nl0 = 2
        do i=1,nlcl
           if (prescl(i).ge.PPB(j)) then
                nl0 = i+1
           endif
        enddo
        
        factflux = (log10(max(PPB(j),prescl(nlcl+1)))
     .                          -log10(prescl(nl0-1)))
     .            /(log10(prescl(nl0))-log10(prescl(nl0-1)))
        if (nsza0.ne.nszacl+1) then
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

