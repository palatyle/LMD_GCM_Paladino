      SUBROUTINE drag_noro (klon,klev,dtime,pplay,pplev,
     e                   pvar, psig, pgam, pthe,
     e                   kgwd,kgwdim,kdx,ktest,
     e                   t, u, v,
     s                   pulow, pvlow, pustr, pvstr,
     s                   d_t, d_u, d_v)
C**** *DRAG_NORO* INTERFACE FOR SUB-GRID SCALE OROGRAPHIC SCHEME
C
C     PURPOSE.
C     --------
C           ZEROS TENDENCIES, COMPUTES GEOPOTENTIAL HEIGHT AND UPDATES THE
C           TENDENCIES AFTER THE SCHEME HAS BEEN CALLED.
C
C     EXPLICIT ARGUMENTS :
C     --------------------
C
C     INPUT :
C
C     NLON               : NUMBER OF HORIZONTAL GRID POINTS
C     NLEV               : NUMBER OF LEVELS
C     DTIME              : LENGTH OF TIME STEP
C     PPLAY(NLON,NLEV+1) : PRESSURE AT MIDDLE LEVELS
C     PPLEV(NLON,NLEV)   : PRESSURE ON MODEL LEVELS
C     PVAR(NLON)         : SUB-GRID SCALE STANDARD DEVIATION
C     PSIG(NLON)         : SUB-GRID SCALE SLOPE
C     PGAM(NLON)         : SUB-GRID SCALE ANISOTROPY
C     PTHE(NLON)         : SUB-GRID SCALE PRINCIPAL AXES ANGLE
C     KGWD               : NUMBER OF POINTS AT WHICH THE SCHEME IS CALLED
C     KGWDIM             : NUMBER OF POINTS AT WHICH THE SCHEME IS CALLED
C     KDX(NLON)          : POINTS AT WHICH TO CALL THE SCHEME
C     KTEST(NLON)        : MAP OF CALLING POINTS
C     T(NLON,NLEV)       : TEMPERATURE
C     U(NLON,NLEV)       : ZONAL WIND
C     V(NLON,NLEV)       : MERIDIONAL WIND
C
C     OUTPUT :
C
C     PULOW(NLON)        : LOW LEVEL ZONAL WIND
C     PVLOW(NLON)        : LOW LEVEL MERIDIONAL WIND
C     PUSTR(NLON)        : LOW LEVEL ZONAL STRESS
C     PVSTR(NLON)        : LOW LEVEL MERIDIONAL STRESS
C     D_T(NLON,NLEV)     : TEMPERATURE TENDENCY
C     D_U(NLON,NLEV)     : ZONAL WIND TENDENCY
C     D_V(NLON,NLEV)     : MERIDIONAL WIND TENDENCY
C
C     IMPLICIT ARGUMENTS :
C     --------------------
C
C     comcstfi.h
C     dimphys.h
C
c
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li F.Lott (LMD/CNRS) date: 19950201
c Objet: Frottement de la montagne Interface
c======================================================================
c Arguments:
c dtime---input-R- pas d'integration (s)
c s-------input-R-la valeur "s" pour chaque couche
c pplay--input-R- pression au milieu des couches en Pa
c pplev--input-R-pression au bords des couches en Pa
c t-------input-R-temperature (K)
c u-------input-R-vitesse horizontale (m/s)
c v-------input-R-vitesse horizontale (m/s)
c
c d_t-----output-R-increment de la temperature t
c d_u-----output-R-increment de la vitesse u
c d_v-----output-R-increment de la vitesse v
c======================================================================
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "comcstfi.h"
c
c ARGUMENTS
c
      REAL dtime
      INTEGER klon,klev
      real pplay(NDLO2,klev),pplev(NDLO2,klev+1)
      REAL pvar(NDLO2),psig(NDLO2),pgam(NDLO2),pthe(NDLO2)
      REAL pulow(NDLO2),pvlow(NDLO2),pustr(NDLO2),pvstr(NDLO2)
      REAL u(NDLO2,klev), v(NDLO2,klev),t(NDLO2,klev)
      REAL d_t(NDLO2,klev), d_u(NDLO2,klev), d_v(NDLO2,klev)
c
      INTEGER i, k, kgwd, kgwdim, kdx(NDLO2), ktest(NDLO2)
c
c Variables locales:
c
      REAL paprs(NDLO2,nlayermx+1)
      REAL paprsf(NDLO2,nlayermx)
      REAL zgeom(NDLO2,nlayermx)
      REAL pdtdt(NDLO2,nlayermx)
      REAL pdudt(NDLO2,nlayermx), pdvdt(NDLO2,nlayermx)
      REAL pt(NDLO2,nlayermx), pu(NDLO2,nlayermx)
      REAL pv(NDLO2,nlayermx)
c
c initialiser les variables de sortie (pour securite)
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
c preparer les variables d'entree (attention: l'ordre des niveaux 
c verticaux augmente du haut vers le bas)
c
      DO k = 1, klev
      DO i = 1, klon
         pt(i,k) = t(i,klev-k+1)
         pu(i,k) = u(i,klev-k+1)
         pv(i,k) = v(i,klev-k+1)
         paprsf(i,k) = pplay(i,klev-k+1)
         paprs(i,k) = pplev(i,klev+1-k+1)
      ENDDO
      ENDDO
      DO i = 1, klon
         paprs(i,klev+1) = pplev(i,1)
      ENDDO
      DO i = 1, klon
         zgeom(i,klev) = r * pt(i,klev)
     .                  * LOG(paprs(i,klev+1)/paprsf(i,klev))
      ENDDO
      DO k = klev-1, 1, -1
      DO i = 1, klon
         zgeom(i,k) = zgeom(i,k+1) + r * (pt(i,k)+pt(i,k+1))/2.0
     .               * LOG(paprsf(i,k+1)/paprsf(i,k))
      ENDDO
      ENDDO
c
c appeler la routine principale
c

      CALL ORODRAG(klon,klev,kgwd,kgwdim,kdx,ktest,
     .            dtime,
     .            paprs, paprsf, zgeom,
     .            pt, pu, pv, pvar, psig, pgam, pthe,
     .            pulow,pvlow,
     .            pdudt,pdvdt,pdtdt)
C
      DO k = 1, klev
      DO i = 1, klon
         d_u(i,klev+1-k) = dtime*pdudt(i,k)
         d_v(i,klev+1-k) = dtime*pdvdt(i,k)
         d_t(i,klev+1-k) = dtime*pdtdt(i,k)
         pustr(i)        = pustr(i)
     .                    +g*pdudt(i,k)*(paprs(i,k+1)-paprs(i,k))
         pvstr(i)        = pvstr(i)
     .                    +g*pdvdt(i,k)*(paprs(i,k+1)-paprs(i,k))
      ENDDO
      ENDDO
c
      RETURN
      END
