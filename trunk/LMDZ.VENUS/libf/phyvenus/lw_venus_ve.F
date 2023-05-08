      SUBROUTINE LW_venus_ve( 
     S              PPB, pt, psi, deltapsi,
     S              PCOOL, 
     S              PTOPLW,PSOLLW,PSOLLWDN,
     S              ZFLNET)
      
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
c     This routine uses the NER matrix
c     (computed for a given cell and temp profile in radlwsw)
c     to compute cooling rates and radiative fluxes.
c
C     AUTHOR.
C     -------
C        Sebastien Lebonnois
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 27/07/2005
C        version multimatrice (topographie, sommet nuages): 20/12/2006
C     ------------------------------------------------------------------
C
C* ARGUMENTS:
C
c inputs

      REAL    PPB(klev+1)  ! inter-couches PRESSURE (bar)
      REAL    pt(klev)     ! mid-layer temperature
      real    psi(0:klev+1,0:klev+1)      ! NER in W/m**2
      real    deltapsi(0:klev+1,0:klev+1) ! D NER / DT in W/m**2/K
C
c output

      REAL   PCOOL(klev)  ! LONGWAVE COOLING (K/s) within each layer
      REAL   PTOPLW       ! LONGWAVE FLUX AT T.O.A. (net, + vers le haut)
      REAL   PSOLLW       ! LONGWAVE FLUX AT SURFACE (net, + vers le haut)
      REAL   PSOLLWDN     ! LONGWAVE FLUX AT SURFACE (down, + vers le bas)
      REAL   ZFLNET(klev+1) ! net thermal flux at ppb levels (+ vers le haut)

C
C* LOCAL VARIABLES:
C
      integer i,j,p
      real   zlnet(klev+1)    ! net thermal flux (W/m**2)
      real   dzlnet(0:klev)   ! Radiative budget (W/m**2)
      real   pdp(klev)        ! epaisseur de la couche en pression (Pa)

c --------------------------
c Calculation of the fluxes
c --------------------------

c flux aux intercouches:
c zlnet(i+1) est le flux net traversant le plafond de la couche i (+ vers le haut)
      do p=0,klev ! numero de la couche
        zlnet(p+1) = 0.0
        do j=p+1,klev+1
         do i=0,p
           zlnet(p+1) = zlnet(p+1)+ psi(i,j)
         enddo
        enddo
      enddo

c flux net au sol, + vers le haut:
      PSOLLW = zlnet(1)
c flux vers le bas au sol, + vers le bas:
      PSOLLWDN = 0.0
      do i=1,klev+1
        PSOLLWDN = PSOLLWDN+max(psi(i,0),0.0)
      enddo

c dfluxnet = radiative budget (W m-2)
      do p=0,klev ! numero de la couche
        dzlnet(p) = 0.0
        do j=0,klev+1
           dzlnet(p) = dzlnet(p)+psi(p,j)
        enddo
      enddo
      
c --------------------------------------
c Interpolation in the GCM vertical grid
c --------------------------------------

c Flux net
c --------
      
      do j=1,klev+1
        ZFLNET(j) =  zlnet(j)
      enddo
      PTOPLW   = ZFLNET(klev+1)
      
c Heating rates
c -------------

c  cool (K/s) = dfluxnet (W/m2)    ! positif quand nrj sort de la couche
c              *g        (m/s2)
c              /(-dp)  (epaisseur couche, en Pa=kg/m/s2)
c              /cp  (J/kg/K) 
      
      do j=1,klev
        pdp(j)=(PPB(j)-PPB(j+1))*1.e5
      enddo

c calcul direct OU calcul par schema implicit
      if (1.eq.1) then 
        do j=1,klev
! ADAPTATION GCM POUR CP(T)
         PCOOL(j) = dzlnet(j) *RG/cpdet(pt(j)) / pdp(j)
        enddo
      else
        call lwi(klev,dzlnet,deltapsi,pdp,pt,PCOOL)
      endif
c     print*,dzlnet
c     print*,pdp
c     print*,PCOOL

      return
      end

