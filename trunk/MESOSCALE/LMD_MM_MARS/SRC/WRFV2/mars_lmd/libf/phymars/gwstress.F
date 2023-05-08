      SUBROUTINE GWSTRESS
     *         (  klon  , klev
     *         , KKCRIT, KSECT, KKHLIM, KTEST, KKCRITH, KCRIT, kkenvh
     *         , kknu
     *         , PRHO  , PSTAB, PVPH  , PVAR ,PVARor, psig  
     *         , PTFR  , PTAU  
     *         ,pgeom1 , pgamma, pd1  , pd2   ,pdmod ,pnu )
C
C**** *GWSTRESS*
C
C     PURPOSE.
C     --------
C
C**   INTERFACE.
C     ----------
C     CALL *GWSTRESS*  FROM *GWDRAG*
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C     ==== OUTPUTS ===
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C
C     EXTERNALS.
C     ----------
C
C
C     REFERENCE.
C     ----------
C
C        SEE ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S."
C
C     AUTHOR.
C     -------
C
C     MODIFICATIONS.
C     --------------
C     F. LOTT PUT THE NEW GWD ON IFS      22/11/93
C
C-----------------------------------------------------------------------
      implicit none
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
      integer klon,klev,kidia,kfdia

#include "yoegwd.h"

C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER KKCRIT(NDLO2),KKCRITH(NDLO2),KCRIT(NDLO2),KSECT(NDLO2),
     *        KKHLIM(NDLO2),KTEST(NDLO2),KKENVH(NDLO2),KKNU(NDLO2)
C
      REAL PRHO(NDLO2,klev+1),PSTAB(NDLO2,klev+1),PTAU(NDLO2,klev+1),
     *     PVPH(NDLO2,klev+1),PVAR(NDLO2,4),PTFR(NDLO2),
     *     pgeom1(NDLO2,klev),PVARor(NDLO2)
C
      real pd1(NDLO2),pd2(NDLO2),pnu(NDLO2),psig(NDLO2),pgamma(NDLO2)
      real pdmod(NDLO2)
C
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
      integer jl
      real zblock,zvar,zeff
      logical lo

C
C-----------------------------------------------------------------------
C
C*       0.3   FUNCTIONS
C              ---------
C     ------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------


      kidia=1
      kfdia=klon


C
 100  CONTINUE
C
C*         3.1     GRAVITY WAVE STRESS.
C
  300 CONTINUE
C
C
      DO 301 JL=kidia,kfdia
      IF(KTEST(JL).EQ.1) THEN
      
C  EFFECTIVE MOUNTAIN HEIGHT ABOVE THE BLOCKED FLOW
  
c        IF(KKENVH(JL).EQ.KLEV)THEN
         ZBLOCK=0.0 
c        ELSE
c         ZBLOCK=(PGEOM1(JL,KKENVH(JL))+PGEOM1(JL,KKENVH(JL)+1))/2./RG          
c        ENDIF
      
        ZVAR=PVAROR(JL)
        ZEFF=AMAX1(0.,2.*ZVAR-ZBLOCK)

        PTAU(JL,KLEV+1)=PRHO(JL,KLEV+1)*GKDRAG*psig(jl)*ZEFF**2
     *    /4./ZVAR*PVPH(JL,KLEV+1)*pdmod(jl)*sqrt(pstab(jl,klev+1))

C  TOO SMALL VALUE OF STRESS OR  LOW LEVEL FLOW INCLUDE CRITICAL LEVEL
C  OR LOW LEVEL FLOW:  GRAVITY WAVE STRESS NUL.
                
        LO=(PTAU(JL,KLEV+1).LT.GTSEC).OR.(KCRIT(JL).GE.KKNU(JL))
     *      .OR.(PVPH(JL,KLEV+1).LT.GVCRIT)
c       IF(LO) PTAU(JL,KLEV+1)=0.0
      
      ELSE
      
          PTAU(JL,KLEV+1)=0.0
          
      ENDIF
      
  301 CONTINUE
C
      RETURN
      END
