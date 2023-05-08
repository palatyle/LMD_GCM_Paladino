      SUBROUTINE SUGWD(NLON,NLEV,paprs,pplay)
C     
C
C**** *SUGWD* INITIALIZE COMMON YOEGWD CONTROLLING GRAVITY WAVE DRAG
C
C     PURPOSE.
C     --------
C           INITIALIZE YOEGWD, THE COMMON THAT CONTROLS THE
C           GRAVITY WAVE DRAG PARAMETRIZATION.
C    VERY IMPORTANT:
C    ______________
C           THIS ROUTINE SET_UP THE "TUNABLE PARAMETERS" OF THE
C           VARIOUS SSO SCHEMES
C
C**   INTERFACE.
C     ----------
C        CALL *SUGWD* FROM *SUPHEC*
C              -----        ------
C (called not from suphec but from first call of physiq.F)
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C        PAPRS,PPLAY : Pressure at semi and full model levels
C        NLEV        : number of model levels
c        NLON        : number of points treated in the physics
C
C        IMPLICIT ARGUMENTS :
C        --------------------
C        COMMON YOEGWD
C-GFRCRIT-R:  Critical Non-dimensional mountain Height
C             (HNC in (1),    LOTT 1999)
C-GKWAKE--R:  Bluff-body drag coefficient for low level wake
C             (Cd in (2),     LOTT 1999)
C-GRCRIT--R:  Critical Richardson Number 
C             (Ric, End of first column p791 of LOTT 1999) 
C-GKDRAG--R:  Gravity wave drag coefficient
C             (G in (3),      LOTT 1999)
C-GKLIFT--R:  Mountain Lift coefficient
C             (Cl in (4),     LOTT 1999)
C-GHMAX---R:  Not used
C-GRAHILO-R:  Set-up the trapped waves fraction
C             (Beta , End of first column,  LOTT 1999)
C
C-GSIGCR--R:  Security value for blocked flow depth
C-NKTOPG--I:  Security value for blocked flow level
C-NTOP----I:  An estimate to qualify the upper levels of
C             the model where one wants to impose strees
C             profiles
C-GSSECC--R:  Security min value for low-level B-V frequency
C-GTSEC---R:  Security min value for anisotropy and GW stress.
C-GVSEC---R:  Security min value for ulow
C         
C
C     METHOD.
C     -------
C        SEE DOCUMENTATION
C
C     EXTERNALS.
C     ----------
C        NONE
C
C     REFERENCE.
C     ----------
C     Lott, 1999: Alleviation of stationary biases in a GCM through...
C                 Monthly Weather Review, 127, pp 788-801.
C
C     AUTHOR.
C     -------
C        FRANCOIS LOTT        *LMD*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 90-01-01 (MARTIN MILLER, ECMWF)
C        LAST:  99-07-09     (FRANCOIS LOTT,LMD)
C     ------------------------------------------------------------------
      use dimphy
      IMPLICIT NONE

#include "YOEGWD.h"
C
C  ARGUMENTS
      integer nlon,nlev
      REAL paprs(nlon,nlev+1)
      REAL pplay(nlon,nlev)
C
      INTEGER JK
      REAL ZPR,ZTOP,ZSIGT,ZPM1R

C
C*       1.    SET THE VALUES OF THE PARAMETERS
C              --------------------------------
C
 100  CONTINUE
C
      PRINT *,' DANS SUGWD NLEV=',NLEV
      GHMAX=10000.
C
      ZPR=100000.
      ZTOP=0.001 
c valeurs dans la dernière routine de FLott
c      ZSIGT=0.94
c valeurs dans les routines Mars
      ZSIGT=0.85
C
Coff  CALL gather(pplay,pplay_glo)
Coff  CALL bcast(pplay_glo)
Coff  CALL gather(paprs,paprs_glo)
Coff  CALL bcast(paprs_glo)

      DO 110 JK=1,NLEV
Coff  ZPM1R=pplay_glo(klon_glo/2,jk)/paprs_glo(klon_glo/2,1) 
      ZPM1R=pplay(klon/2,jk)/paprs(klon/2,1) 
      IF(ZPM1R.GE.ZSIGT)THEN
         nktopg=JK
      ENDIF
Coff  ZPM1R=pplay_glo(klon_glo/2,jk)/paprs_glo(klon_glo/2,1) 
      ZPM1R=pplay(klon/2,jk)/paprs(klon/2,1) 
      IF(ZPM1R.GE.ZTOP)THEN
         ntop=JK
      ENDIF
  110 CONTINUE
c
c  inversion car dans orodrag on compte les niveaux a l'envers
      nktopg=nlev-nktopg+1
      ntop=nlev-ntop
      print *,' DANS SUGWD nktopg=', nktopg
      print *,' DANS SUGWD ntop=', ntop
C
      GSIGCR=0.80
C
c valeurs dans la dernière routine de FLott
c      GKDRAG=0.1875
c      GRAHILO=0.1   
c      GRCRIT=1.00 
c      GFRCRIT=1.00
c      GKWAKE=0.50
C
c      GKLIFT=0.25
c      GVCRIT =0.1

c valeurs dans les routines Mars
c      GKDRAG=0.1
c      GRAHILO=1.0   
c      GRCRIT=0.25 
c      GFRCRIT=1.00
c      GKWAKE=1.0
C
C VENUS
      GKDRAG=0.5      ! G
      GRAHILO=1.0     ! beta - useless
      GRCRIT=0.25     ! Ric  - useless
      GFRCRIT=1.0     ! Hnc
      GKWAKE=1.0      ! Cd
      TAUBS=2.0       ! VENUS: stress threshold is 2 Pa
      !TAUBS=1.0       ! VENUS: stress threshold is 1 Pa
      !TAUBS=0.5       ! VENUS: stress threshold is 0.5 Pa
      LEVBS=nlev-9   ! VENUS: level release is 9
      !LEVBS=nlev-19   ! VENUS: level release is 19
      !LEVBS=nlev-13   ! VENUS: level release is 13

      GKLIFT=0.25
      GVCRIT =0.0

      WRITE(UNIT=6,FMT='('' *** SSO essential constants ***'')')
      WRITE(UNIT=6,FMT='('' *** SPECIFIED IN SUGWD ***'')')
      WRITE(UNIT=6,FMT='('' Gravity wave ct '',E14.7,'' '')')GKDRAG
      WRITE(UNIT=6,FMT='('' Trapped/total wave dag '',E14.7,'' '')')
     S                  GRAHILO
      WRITE(UNIT=6,FMT='('' Critical Richardson   = '',E14.7,'' '')')
     S                  GRCRIT
      WRITE(UNIT=6,FMT='('' Critical Froude'',e14.7)') GFRCRIT
      WRITE(UNIT=6,FMT='('' Low level Wake bluff cte'',e14.7)') GKWAKE
      WRITE(UNIT=6,FMT='('' Low level lift  cte'',e14.7)') GKLIFT

      WRITE(UNIT=6,FMT='('' VENUS: Mountain stress threshold'',E14.7)')
     S                  TAUBS
      WRITE(UNIT=6,FMT='('' VENUS: Level release'',I5)') nlev - LEVBS
C
C
C      ----------------------------------------------------------------
C
C*       2.    SET VALUES OF SECURITY PARAMETERS
C              ---------------------------------
C
 200  CONTINUE
C
c valeurs dans la dernière routine de FLott
c      GVSEC=0.10
c      GSSEC=0.0001
C
c      GTSEC=0.00001
C
c valeurs dans les routines Mars
      GVSEC=0.10
      GSSEC=1.e-12
C
      GTSEC=1.e-7
C
      RETURN
      END

