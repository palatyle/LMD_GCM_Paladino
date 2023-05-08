      SUBROUTINE SW_venus_rh(PRMU0, PFRAC, latdeg,
     S              PPA, PPB, pt,
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
      REAL   latdeg ! |latitude| (in degrees)
      REAL   PPB(klev+1)  ! inter-couches PRESSURE (bar)
      REAL   PPA(klev)
      REAL   pt(klev)     ! mid-layer temperature
C
c output

      REAL   PHEAT(klev)  ! SHORTWAVE HEATING (K/s) within each layer
      REAL   PHEATPPA(klev)
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
      
      integer i,j,k,lat,nsza,nsza0(2),nl0,nlat0
      real   zsnet(nlrh+1,nszarh+1,nlatrh+1)! net solar flux (W/m**2) (+ vers bas)
      real   solza(nszarh,nlatrh)       ! solar zenith angles in table
      real   presrh(nlrh+1)             ! pressure in table (bar)
      real   logplayrh(nlrh)
      real   altrh(nlrh+1)              ! altitude in table (km)
      real   latrh(nlatrh)              ! latitude in table (degrees)
      character*22 nullchar
      real   sza0,factsza(2),factflux,factlat
      real   zsnetmoy
      logical firstcall
      data    firstcall/.true./
      save   solza,zsnet,altrh,latrh,presrh
      save   firstcall
      real   Tplay(nlrh)
      real   Qrh1(nlrh)
      real   Qrh2(nlrh)
      real   Qrh3(nlrh)
      real   Qrh4(nlrh)
      
c ------------------------
c Loading the file
c ------------------------
      if (firstcall) then

       zsnet=0.

       open(11,file='SolarNetFlux_RH.dat')

       do i=1,nlrh+1
          read(11,'(E5.1,4x,F8.2)') altrh(i),presrh(i)
       enddo

       do lat=1,nlatrh
         latrh(lat)=5.*(lat-1)
         read(11,*) nullchar
         read(11,*) nullchar
         read(11,'(3x,7(5x,E8.5))') solza(:,lat)
         read(11,*) nullchar

         do i=1,nlrh+1
          read(11,'(E6.1,7(2x,F11.5),7x,F11.5)')
     .          altrh(i),zsnet(i,1:nszarh,lat),zsnetmoy
         enddo
         read(11,*) nullchar
       enddo
       latrh(nlatrh)=89.

c Correction of factor 2 in the table...
       zsnet=zsnet*2.

       close(11)

       firstcall=.false.
      endif

c --------------------------------------
c Interpolation in the GCM vertical grid
c --------------------------------------

c Latitude
c ---------
      
      do lat=1,nlatrh
         if (latrh(lat).le.latdeg) then
              nlat0 = lat+1
         endif
      enddo

      if (nlat0.ne.nlatrh+1) then
        factlat = (latdeg-latrh(nlat0-1))/(latrh(nlat0)-latrh(nlat0-1))
      else
        factlat = min((latdeg-latrh(nlatrh))/(90.-latrh(nlatrh)), 1.)
      endif

c Zenith angle
c ------------
      
      sza0 = acos(PRMU0)/3.1416*180.
      nsza0(:)=2

      do nsza=1,nszarh
         if (solza(nsza,nlat0-1).le.sza0) then
              nsza0(1) = nsza+1
         endif
      enddo
      if (nsza0(1).ne.nszarh+1) then
          factsza(1) = (sza0-solza(nsza0(1)-1,nlat0-1))/
     .        (solza(nsza0(1),nlat0-1)-solza(nsza0(1)-1,nlat0-1))
      else
          factsza(1) = min((sza0-solza(nszarh,nlat0-1))/
     .         (90.-solza(nszarh,nlat0-1)), 1.)
      endif
      if (nlat0.ne.nlatrh+1) then
       do nsza=1,nszarh
         if (solza(nsza,nlat0).le.sza0) then
              nsza0(2) = nsza+1
         endif
       enddo
       if (nsza0(2).eq.nszarh+1) then
          factsza(2) = min((sza0-solza(nszarh,nlat0))/
     .         (90.-solza(nszarh,nlat0)), 1.)
       elseif ((nsza0(2).eq.2).and.(solza(1,nlat0).gt.sza0)) then
          factsza(2) = 0.
       else
          factsza(2) = (sza0-solza(nsza0(2)-1,nlat0))/
     .        (solza(nsza0(2),nlat0)-solza(nsza0(2)-1,nlat0))
       endif
      else
        nsza0(2)   = nszarh+1
        factsza(2) = 1.
      endif
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

        ZFSNET(j) =  factlat*(
     .      factflux   *  factsza(2)   *zsnet(nl0,nsza0(2),nlat0)
     . +   factflux   *(1.-factsza(2))*zsnet(nl0,nsza0(2)-1,nlat0)
     . + (1.-factflux)*  factsza(2)   *zsnet(nl0+1,nsza0(2),nlat0)
     . + (1.-factflux)*(1.-factsza(2))*zsnet(nl0+1,nsza0(2)-1,nlat0) )
     .            + (1.-factlat)*(
     .      factflux   *  factsza(1)   *zsnet(nl0,nsza0(1),nlat0-1)
     . +   factflux   *(1.-factsza(1))*zsnet(nl0,nsza0(1)-1,nlat0-1)
     . + (1.-factflux)*  factsza(1)   *zsnet(nl0+1,nsza0(1),nlat0-1)
     . + (1.-factflux)*(1.-factsza(1))*zsnet(nl0+1,nsza0(1)-1,nlat0-1) )

        ZFSNET(j) = ZFSNET(j)*PFRAC
        
      enddo
      PTOPSW = ZFSNET(klev+1)
      PSOLSW = ZFSNET(1) 

#ifdef MESOSCALE
! extrapolation play RH pressure
      do j=1,nlrh
         logplayrh(j)=(log(presrh(j+1))+log(presrh(j)))/2.
      enddo
! Extrapolation of temperature over RH play pressure
      do i=nlrh,2,-1
        nl0 = 2
        do j=1,klev-1
           if (exp(logplayrh(i)).le.PPA(j)) then
                nl0 = j+1
           endif
        enddo
        factflux = (log10(max(exp(logplayrh(i)),PPA(klev)))
     .                         -log10(PPA(nl0-1)))
     .       /(log10(PPA(nl0))-log10(PPA(nl0-1)))
        Tplay(i)=factflux*pt(nl0)
     .             + (1.-factflux)*pt(nl0-1)

      ENDDO
! RH PHEAT over RH play pressure
      DO k=1,nlrh
c
       Qrh1(k)=((RG/cpdet(Tplay(k)))
     .     *((zsnet(k+1,nsza0(1),nlat0-1)-zsnet(k,nsza0(1),nlat0-1))
     .         *PFRAC))
     .      /((presrh(k)-presrh(k+1))*1.e5)
       Qrh2(k)=((RG/cpdet(Tplay(k)))
     . *((zsnet(k+1,nsza0(1)-1,nlat0-1)-zsnet(k,nsza0(1)-1,nlat0-1))
     .         *PFRAC))
     .      /((presrh(k)-presrh(k+1))*1.e5)
       Qrh3(k)=((RG/cpdet(Tplay(k)))
     .       *((zsnet(k+1,nsza0(2),nlat0)-zsnet(k,nsza0(2),nlat0))
     .         *PFRAC))
     .      /((presrh(k)-presrh(k+1))*1.e5)
       Qrh4(k)=((RG/cpdet(Tplay(k)))
     .       *((zsnet(k+1,nsza0(2)-1,nlat0)-zsnet(k,nsza0(2)-1,nlat0))
     .        *PFRAC))
     .      /((presrh(k)-presrh(k+1))*1.e5)
      ENDDO
! Interapolation of PHEAT over GCM/MESOSCALE play levels
      do j=1,klev
        nl0 = nlrh-1
        do i=nlrh,2,-1
           if (exp(logplayrh(i)).ge.PPA(j)) then
                nl0 = i-1
           endif
        enddo
c        factflux = (log10(max(PPB(j),presrh(1)))-log10(presrh(nl0+1)))
c     .            /(log10(presrh(nl0))-log10(presrh(nl0+1)))
        factflux = (log10(max(PPA(j),exp(logplayrh(1))))
     .                         -log10(exp(logplayrh(nl0+1))))
     .     /(log10(exp(logplayrh(nl0)))-log10(exp(logplayrh(nl0+1))))
        PHEATPPA(j)=factlat*(
     .      factflux   *  factsza(2)  *Qrh3(nl0)
     . +   factflux   *(1.-factsza(2))*Qrh4(nl0)
     . + (1.-factflux)*  factsza(2)   *Qrh3(nl0+1)
     . + (1.-factflux)*(1.-factsza(2))*Qrh4(nl0+1))
     .            + (1.-factlat)*(
     .      factflux   *  factsza(1)  *Qrh1(nl0)
     . +   factflux   *(1.-factsza(1))*Qrh2(nl0)
     . + (1.-factflux)*  factsza(1)   *Qrh1(nl0+1)
     . + (1.-factflux)*(1.-factsza(1))*Qrh2(nl0+1) )
        PHEAT(j)=PHEATPPA(j)
      ENDDO


#else      
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
c        if ((PPB(j).gt.0.04).and.(PPB(j).le.0.1)) then
c         PHEAT(j) = PHEAT(j)*1.5
c        endif
c        if ((PPB(j).gt.0.1).and.(PPB(j).le.0.5)) then
c         PHEAT(j) = PHEAT(j)*2.
c        endif
        if ((PPB(j).gt.1.4).and.(PPB(j).le.10.)) then
c        if ((PPB(j).gt.1.4).and.(PPB(j).le.100.)) then
         PHEAT(j) = PHEAT(j)*3
        endif
c        if ((PPB(j).gt.10.).and.(PPB(j).le.120.)) then
c         PHEAT(j) = PHEAT(j)*2
c        endif
c----------------
      enddo
#endif
     

      return
      end

