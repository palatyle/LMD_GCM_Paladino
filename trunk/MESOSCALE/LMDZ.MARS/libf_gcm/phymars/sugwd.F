      SUBROUTINE SUGWD(NLEV,PSIG)
C     
C**** *SUGWD* INITIALIZE COMMON YOEGWD CONTROLLING GRAVITY WAVE DRAG
C     
C     PURPOSE.
C     --------
C     INITIALIZE YOEGWD, THE COMMON THAT CONTROLS THE
C     GRAVITY WAVE DRAG PARAMETRIZATION.
C     
C     EXPLICIT ARGUMENTS :
C     --------------------
C     
C     INPUT :
C     
C     PSIG        : VERTICAL COORDINATE TABLE
C     NLEV        : NUMBER OF MODEL LEVELS
C     
C     OUTPUT :
C     
C     NONE
C     
C     IMPLICIT ARGUMENTS :
C     --------------------
C     COMMON YOEGWD
C     
C     METHOD.
C     -------
C     SEE DOCUMENTATION
C     
C     EXTERNALS.
C     ----------
C     NONE
C     
C     REFERENCE.
C     ----------
C     ECMWF Research Department documentation of the IFS
C     
C     AUTHOR.
C     -------
C     MARTIN MILLER             *ECMWF*
C     
C     MODIFICATIONS.
C     --------------
C     ORIGINAL : 90-01-01
C     ------------------------------------------------------------------
      implicit none
C     
C     -----------------------------------------------------------------
c     !-*-      include 'yoegwd.h'
#include "yoegwd.h"
c     !-*-
C     ----------------------------------------------------------------
C     
      integer nlev
      REAL PSIG(NLEV+1)

      real zsigt,zpr,zpm1r
      integer jk
C     
C     *       1.    SET THE VALUES OF THE PARAMETERS
C     --------------------------------
C     
 100  CONTINUE
C     
c     PRINT *,' Dans sugwd nlev=',NLEV,' SIG=',PSIG
      GHMAX=10000.
C     
c     old  ZSIGT=0.94
c     old  ZPR=80000.
      ZSIGT=0.85
      ZPR=100000.
C     
      DO 110 JK=1,NLEV-1
         ZPM1R=0.5*ZPR*(PSIG(JK)+PSIG(JK+1))
         IF((ZPM1R/ZPR).GE.ZSIGT)THEN
            NKTOPG=JK
         ENDIF
 110  CONTINUE
      WRITE(*,*) 'In sugwd NKTOPG=',NKTOPG
C     
      GSIGCR=0.80
C     
      GKDRAG= 0.1  ! used to be 0.1 for mcd Version 1 and 2 (before 10/2000)
      GRCRIT=0.25
      GFRCRIT=1.0
      GKWAKE=1.0

C     
      GKDRAGL=4.*GKDRAG
      GRAHILO=1.
      GVCRIT =0.0
C     
C     
C     ----------------------------------------------------------------
C     
C     *       2.    SET VALUES OF SECURITY PARAMETERS
C     ---------------------------------
C     
 200  CONTINUE
C     
      GVSEC=0.10
      GSSEC=1.E-12
C     
      GTSEC=1.E-07
C     
C     ----------------------------------------------------------------
C     
      RETURN
      END
