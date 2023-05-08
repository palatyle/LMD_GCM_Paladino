
C  SUBROUTINE DE PARAMETRISATION DES MONTAGNES D ECHELLE SOUS MAILLE

      SUBROUTINE drag_noro (nlon,nlev,dtime,paprs,pplay,pgeop,pn2,
     e                   pmea,pstd, psig, pgam, pthe,ppic,pval,
     e                   kgwd,kdx,ktest,
     e                   t, u, v,
     s                   pulow, pvlow, pustr, pvstr,
     s                   d_t, d_u, d_v,
     s                   blustr,blvstr,pnlow,zeff,zbl,
     s                   ptau,tau0,knu2,kbreak)
c
      use dimphy
      IMPLICIT none

c======================================================================
c Auteur(s): F.Lott (LMD/CNRS) date: 19950201
c Object: Mountain drag interface. Made necessary because:
C 1. in the LMD-GCM Layers are from bottom to top,
C    contrary to most European GCM.
c 2. the altitude above ground of each model layers
c    needs to be known (variable zgeom)
c======================================================================
c Explicit Arguments:
c ==================
c nlon----input-I-Total number of horizontal points that get into physics
c nlev----input-I-Number of vertical levels
c dtime---input-R-Time-step (s)
c paprs---input-R-Pressure in semi layers    (Pa)
c pplay---input-R-Pressure model-layers      (Pa)
c pgeop---input-R-Geopotential model layers (reference to ground)
c pn2-----input-R-Brunt-Vaisala freq.^2 at 1/2 layers
c t-------input-R-temperature (K)
c u-------input-R-Horizontal wind (m/s)
c v-------input-R-Meridional wind (m/s)
c pmea----input-R-Mean Orography (m)
C pstd----input-R-SSO standard deviation (m)
c psig----input-R-SSO slope
c pgam----input-R-SSO Anisotropy
c pthe----input-R-SSO Angle
c ppic----input-R-SSO Peacks elevation (m)
c pval----input-R-SSO Valleys elevation (m)
c
c kgwd- -input-I: Total nb of points where the orography schemes are active
c ktest--input-I: Flags to indicate active points
c kdx----input-I: Locate the physical location of an active point.

c pulow, pvlow -output-R: Low-level wind
c pustr, pvstr -output-R: Surface stress due to SSO drag      (Pa)
c
c d_t-----output-R: T increment            
c d_u-----output-R: U increment              
c d_v-----output-R: V increment              
c
c Implicit Arguments:
c ===================
c
c iim--common-I: Number of longitude intervals
c jjm--common-I: Number of latitude intervals
c klon-common-I: Number of points seen by the physics
c                (iim+1)*(jjm+1) for instance
c klev-common-I: Number of vertical layers
c======================================================================
c Local Variables:
c ================
c
c zgeom-----R: Altitude (m) of layer above ground (from top to bottom)
c pt, pu, pv --R: t u v from top to bottom
c pdtdt, pdudt, pdvdt --R: t u v tendencies (from top to bottom) 
c papmf: pressure at model layer (from top to bottom)
c papmh: pressure at model 1/2 layer (from top to bottom)
c 
c======================================================================

#include "YOMCST.h"
#include "YOEGWD.h"

c  ARGUMENTS
c
      INTEGER nlon,nlev
      REAL dtime
      REAL paprs(nlon,nlev+1)
      REAL pplay(nlon,nlev)
      REAL pgeop(nlon,nlev),pn2(nlon,nlev)
      REAL pmea(nlon),pstd(nlon),psig(nlon),pgam(nlon),pthe(nlon)
      REAL ppic(nlon),pval(nlon)
      REAL pulow(nlon),pvlow(nlon),pustr(nlon),pvstr(nlon)
      REAL t(nlon,nlev), u(nlon,nlev), v(nlon,nlev)
      REAL d_t(nlon,nlev), d_u(nlon,nlev), d_v(nlon,nlev)
      REAL blustr(nlon),blvstr(nlon),pnlow(nlon),zeff(nlon),zbl(nlon)
      REAL knu2(nlon),kbreak(nlon)
      REAL ztau(klon,klev+1), ptau(klon,klev), tau0(klon)
c
      INTEGER i, k, kgwd,  kdx(nlon), ktest(nlon)
      INTEGER ikenvh(nlon)
      INTEGER iknu2(nlon)
      INTEGER ikbreak(nlon)
c
c LOCAL VARIABLES:
c
      REAL zgeom(klon,klev),zn2(klon,klev)
      REAL pdtdt(klon,klev), pdudt(klon,klev), pdvdt(klon,klev)
      REAL pt(klon,klev), pu(klon,klev), pv(klon,klev)
      REAL papmf(klon,klev),papmh(klon,klev+1)
c
c INITIALIZE OUTPUT VARIABLES 
c
      DO i = 1,klon
         pulow(i) = 0.0
         pvlow(i) = 0.0
         pustr(i) = 0.0
         pvstr(i) = 0.0
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_u(i,k) = 0.0
         d_v(i,k) = 0.0
         pdudt(i,k)=0.0
         pdvdt(i,k)=0.0
         pdtdt(i,k)=0.0
      ENDDO
      ENDDO
c
c PREPARE INPUT VARIABLES FOR ORODRAG (i.e., ORDERED FROM TOP TO BOTTOM)
C CALCULATE LAYERS HEIGHT ABOVE GROUND)
c
      DO k = 1, klev
      DO i = 1, klon
         pt(i,k) = t(i,klev-k+1) 
         pu(i,k) = u(i,klev-k+1)
         pv(i,k) = v(i,klev-k+1)
         papmf(i,k) = pplay(i,klev-k+1)
      ENDDO
      ENDDO
      DO k = 1, klev+1
      DO i = 1, klon
         papmh(i,k) = paprs(i,klev-k+2)
      ENDDO
      ENDDO

      DO k = klev, 1, -1
      DO i = 1, klon
         zgeom(i,k) = pgeop(i,klev-k+1)/RG
         zn2(i,k)   = pn2(i,klev-k+1)
      ENDDO
      ENDDO

c CALL SSO DRAG ROUTINES        
c
      CALL orodrag(klon,klev,kgwd,kdx,ktest,
     .            dtime,
     .            papmh, papmf, zgeom, zn2,
     .            pt, pu, pv,
     .            pmea, pstd, psig, pgam, pthe, ppic,pval,
     .            pulow,pvlow,
     .            pdudt,pdvdt,pdtdt,
     .            blustr,blvstr,pnlow,zeff,ikenvh,
     .            ztau,iknu2,ikbreak)

      zbl(:) = real(klev-ikenvh(:))
      knu2(:) = real(klev-iknu2(:))
      kbreak(:) = real(klev-ikbreak(:))
      tau0 = ztau(:,klev+1)

C
C COMPUTE INCREMENTS AND STRESS FROM TENDENCIES

      DO k = 1, klev
      DO i = 1, klon
         ptau(i,klev+1-k) = ztau(i,k)
         d_u(i,klev+1-k) = dtime*pdudt(i,k)
         d_v(i,klev+1-k) = dtime*pdvdt(i,k)
         d_t(i,klev+1-k) = dtime*pdtdt(i,k)
         pustr(i)        = pustr(i)
     .                    +pdudt(i,k)*(papmh(i,k+1)-papmh(i,k))/rg
         pvstr(i)        = pvstr(i)
     .                    +pdvdt(i,k)*(papmh(i,k+1)-papmh(i,k))/rg
      ENDDO
      ENDDO
c
      RETURN
      END

