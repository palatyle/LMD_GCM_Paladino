      SUBROUTINE ORODRAG( klon,klev 
     I                 , KGWD, KGWDIM, KDX, KTEST
     R                 , PTSPHY
     R                 , PAPHM1,PAPM1,PGEOM1,PTM1,PUM1
     R                 , PVM1, PVAROR, PSIG, PGAMMA, PTHETA
C OUTPUTS
     R                 , PULOW,PVLOW
     R                 , PVOM,PVOL,PTE )
C
C
C**** *ORODRAG* - DOES THE GRAVITY WAVE PARAMETRIZATION.
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
C          PROGNOSTIC VARIABLES U,V  AND T DUE TO  VERTICAL TRANSPORTS BY
C          SUBGRIDSCALE OROGRAPHICALLY EXCITED GRAVITY WAVES
C
C     EXPLICIT ARGUMENTS :
C     --------------------
C
C     INPUT :
C
C     NLON               : NUMBER OF HORIZONTAL GRID POINTS
C     NLEV               : NUMBER OF LEVELS
C     KGWD               : NUMBER OF POINTS AT WHICH THE SCHEME IS CALLED
C     KGWDIM             : NUMBER OF POINTS AT WHICH THE SCHEME IS CALLED
C     KDX(NLON)          : POINTS AT WHICH TO CALL THE SCHEME
C     KTEST(NLON)        : MAP OF CALLING POINTS
C     PTSPHY             : LENGTH OF TIME STEP
C     PAPHM1(NLON,NLEV+1): PRESSURE AT MIDDLE LEVELS
C     PAPM1(NLON,NLEV)   : PRESSURE ON MODEL LEVELS
C     PGEOM1(NLON,NLEV)  : GEOPOTENTIAL HIEGHT OF MODEL LEVELS
C     PTM1(NLON,NLEV)    : TEMPERATURE
C     PUM1(NLON,NLEV)    : ZONAL WIND
C     PVM1(NLON,NLEV)    : MERIDIONAL WIND
C     PVAROR(NLON)       : SUB-GRID SCALE STANDARD DEVIATION
C     PSIG(NLON)         : SUB-GRID SCALE SLOPE
C     PGAMMA(NLON)       : SUB-GRID SCALE ANISOTROPY
C     PTHETA(NLON)       : SUB-GRID SCALE PRINCIPAL AXES ANGLE
C
C     OUTPUT :
C
C     PULOW(NLON)        : LOW LEVEL ZONAL WIND
C     PVLOW(NLON)        : LOW LEVEL MERIDIONAL WIND
C     PVOM(NLON,NLEV)    : ZONAL WIND TENDENCY
C     PVOL(NLON,NLEV)    : MERIDIONAL WIND TENDENCY
C     PTE(NLON,NLEV)     : TEMPERATURE TENDENCY
C
C     IMPLICIT ARGUMENTS :
C     --------------------
C
C     comcstfi.h
C     dimphys.h
C     yoegwd.h
C
C     METHOD.
C     -------
C
C     EXTERNALS.
C     ----------
C
C     REFERENCE.
C     ----------
C
C     AUTHOR.
C     -------
C     M.MILLER + B.RITTER   E.C.M.W.F.     15/06/86.
C
C     F.LOTT + M. MILLER    E.C.M.W.F.     22/11/94
C-----------------------------------------------------------------------
      implicit none
C
C
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
      integer klon,klev,kidia,kfdia
      parameter(kidia=1,kfdia=NDLO2)

#include "comcstfi.h"
#include "yoegwd.h"
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
C
      REAL  PTE(NDLO2,klev),
     *      PVOL(NDLO2,klev),
     *      PVOM(NDLO2,klev),
     *      PULOW(NDLO2),
     *      PVLOW(NDLO2)
      REAL  PUM1(NDLO2,klev),
     *      PVM1(NDLO2,klev),
     *      PTM1(NDLO2,klev),
     *      PVAROR(NDLO2),PSIG(NDLO2),PGAMMA(NDLO2),PTHETA(NDLO2),
     *      PGEOM1(NDLO2,klev),
     *      PAPM1(NDLO2,klev),
     *      PAPHM1(NDLO2,klev+1)
C
      integer kgwd,kgwdim
      real ptsphy
      INTEGER  KDX(NDLO2),KTEST(NDLO2)
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
      INTEGER  ISECT(NDLO2),
     *         ICRIT(NDLO2),
     *         IKCRITH(NDLO2),
     *         IKenvh(NDLO2),
     *         IKNU(NDLO2),
     *         IKNU2(NDLO2),
     *         IKCRIT(NDLO2),
     *         IKHLIM(NDLO2)
      integer ji,jk,jl,klevm1,ilevp1
C      real gkwake
      real ztmst,pvar,ztauf,zrtmst,zdelp,zb,zc,zbet
      real zconb,zabsv,zzd1,ratio,zust,zvst,zdis,ztemp
C
      REAL   ZTAU(NDLO2,nlayermx+1),
     *       ZSTAB(NDLO2,nlayermx+1),
     *       ZVPH(NDLO2,nlayermx+1),
     *       ZRHO(NDLO2,nlayermx+1),
     *       ZRI(NDLO2,nlayermx+1),
     *       ZpsI(NDLO2,nlayermx+1),
     *       Zzdep(NDLO2,nlayermx)
      REAL   ZDUDT(NDLO2),
     *       ZDVDT(NDLO2),
     *       ZDTDT(NDLO2),
     *       ZDEDT(NDLO2),
     *       ZVIDIS(NDLO2),
     *       ZTFR(NDLO2),
     *       Znu(NDLO2),
     *       Zd1(NDLO2),
     *       Zd2(NDLO2),
     *       Zdmod(NDLO2)
C
C------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         1.1   COMPUTATIONAL CONSTANTS
C                -----------------------
C
 110  CONTINUE
C
c     ZTMST=TWODT
c     IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      KLEVM1=KLEV-1
      ZTMST=PTSPHY
      ZRTMST=1./ZTMST
C     ------------------------------------------------------------------
C
 120  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         1.3   CHECK WHETHER ROW CONTAINS POINT FOR PRINTING
C                ---------------------------------------------
C
 130  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.     PRECOMPUTE BASIC STATE VARIABLES.
C*                ---------- ----- ----- ----------
C*                DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
C*                LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
C*                THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.
C
  200 CONTINUE
C
C
C
      CALL OROSETUP
     *     ( klon, klev , KTEST 
     *     , IKCRIT, IKCRITH, ICRIT, ISECT, IKHLIM, ikenvh,IKNU,iknu2
     *     , PAPHM1, PAPM1 , PUM1   , PVM1 , PTM1 , PGEOM1, pvaror
     *     , ZRHO  , ZRI   , ZSTAB  , ZTAU , ZVPH , zpsi, zzdep
     *     , PULOW, PVLOW 
     *     , ptheta,pgamma,znu  ,zd1,  zd2,  zdmod                )
C
C
C
C***********************************************************
C
C
C*         3.      COMPUTE LOW LEVEL STRESSES USING SUBCRITICAL AND
C*                 SUPERCRITICAL FORMS.COMPUTES ANISOTROPY COEFFICIENT
C*                 AS MEASURE OF OROGRAPHIC TWODIMENSIONALITY.
C
  300 CONTINUE
C
      CALL GWSTRESS
     *    ( klon  , klev
     *    , IKCRIT, ISECT, IKHLIM, KTEST, IKCRITH, ICRIT, ikenvh, IKNU
     *    , ZRHO  , ZSTAB, ZVPH  , PVAR , pvaror,  psig
     *    , ZTFR   , ZTAU 
     *    , pgeom1,pgamma,zd1,zd2,zdmod,znu)
C
C*         4.      COMPUTE STRESS PROFILE.
C*                 ------- ------ --------
C
  400 CONTINUE
C
C
      CALL GWPROFIL
     *       (  klon , klev
     *       , kgwd   , kdx  , KTEST
     *       , IKCRIT, IKCRITH, ICRIT  , ikenvh, IKNU
     *       ,iknu2 , pAPHM1, ZRHO   , ZSTAB , ZTFR   , ZVPH
     *       , ZRI   , ZTAU   , ztauf
     *       , zdmod , znu    , psig  , pgamma , pvaror   )
C
C
C*         5.      COMPUTE TENDENCIES.
C*                 -------------------
C
  500 CONTINUE
C
C  EXPLICIT SOLUTION AT ALL LEVELS FOR THE GRAVITY WAVE
C  IMPLICIT SOLUTION FOR THE BLOCKED LEVELS

      DO 510 JL=KIDIA,KFDIA
      ZVIDIS(JL)=0.0
      ZDUDT(JL)=0.0
      ZDVDT(JL)=0.0
      ZDTDT(JL)=0.0
  510 CONTINUE
C
      ILEVP1=KLEV+1
C
C
      DO 524 JK=1,klev
C
CDIR$ IVDEP
C
C      GKWAKE=0.5
C
C     NOW SET IN SUGWD.F
C
      DO 523 JL=1,KGWD
      JI=KDX(JL)
      ZDELP=pAPHM1(Ji,JK+1)-pAPHM1(Ji,JK)
      ZTEMP=-g*(ZTAU(Ji,JK+1)-ZTAU(Ji,JK))/(ZVPH(Ji,ILEVP1)*ZDELP)
      ZDUDT(JI)=(PULOW(JI)*Zd1(ji)-pvlow(ji)*zd2(ji))*ztemp/zdmod(ji)
      ZDVDT(JI)=(pvLOW(JI)*Zd1(ji)+pulow(ji)*zd2(ji))*ztemp/zdmod(ji)
      if(jk.ge.ikenvh(ji)) then
         zb=1.0-0.18*pgamma(ji)-0.04*pgamma(ji)**2
         zc=0.48*pgamma(ji)+0.3*pgamma(ji)**2
         zconb=2.*ztmst*GKWAKE*psig(ji)/(4.*pvaror(ji))
         zabsv=sqrt(PUM1(JI,JK)**2+PVM1(JI,JK)**2)/2.
         zzd1=zb*cos(zpsi(ji,jk))**2+zc*sin(zpsi(ji,jk))**2
         ratio=(cos(zpsi(ji,jk))**2+pgamma(ji)*sin(zpsi(ji,jk))**2)/
     *   (pgamma(ji)*cos(zpsi(ji,jk))**2+sin(zpsi(ji,jk))**2)
         zbet=max(0.,2.-1./ratio)*zconb*zzdep(ji,jk)*zzd1*zabsv
         zdudt(ji)=-pum1(ji,jk)/ztmst
         zdvdt(ji)=-pvm1(ji,jk)/ztmst
         zdudt(ji)=zdudt(ji)*(zbet/(1.+zbet))
         zdvdt(ji)=zdvdt(ji)*(zbet/(1.+zbet))
      end if
      PVOM(JI,JK)=ZDUDT(JI)
      PVOL(JI,JK)=ZDVDT(JI)
      ZUST=PUM1(JI,JK)+ZTMST*ZDUDT(JI)
      ZVST=PVM1(JI,JK)+ZTMST*ZDVDT(JI)
      ZDIS=0.5*(PUM1(JI,JK)**2+PVM1(JI,JK)**2-ZUST**2-ZVST**2)
      ZDEDT(JI)=ZDIS/ZTMST
      ZVIDIS(JI)=ZVIDIS(JI)+ZDIS*ZDELP
      ZDTDT(JI)=ZDEDT(JI)/cpp
      PTE(JI,JK)=ZDTDT(JI)

  523 CONTINUE

  524 CONTINUE
C
C
      RETURN
      END
