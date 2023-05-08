      SUBROUTINE GWPROFIL
     *         ( klon, klev
     *         , kgwd ,kdx  , ktest
     *         , KKCRIT, KKCRITH, KCRIT ,  kkenvh, kknu,kknu2
     *         , PAPHM1, PRHO   , PSTAB , PTFR , PVPH , PRI , PTAU
     *         , ptauf ,pdmod   , pnu   , psig ,pgamma, pvar      )

C**** *GWPROFIL*
C
C     PURPOSE.
C     --------
C
C**   INTERFACE.
C     ----------
C          FROM *GWDRAG*
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C     ==== OUTPUTS ===
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD:
C     -------
C     THE STRESS PROFILE FOR GRAVITY WAVES IS COMPUTED AS FOLLOWS:
C     IT IS CONSTANT (NO GWD) AT THE LEVELS BETWEEN THE GROUND
C     AND THE TOP OF THE BLOCKED LAYER (KKENVH).
C     IT DECREASES LINEARLY WITH HEIGHTS FROM THE TOP OF THE 
C     BLOCKED LAYER TO 3*VAROR (kKNU), TO SIMULATES LEE WAVES OR 
C     NONLINEAR GRAVITY WAVE BREAKING.
C     ABOVE IT IS CONSTANT, EXCEPT WHEN THE WAVE ENCOUNTERS A CRITICAL
C     LEVEL (KCRIT) OR WHEN IT BREAKS.
C     
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
C     PASSAGE OF THE NEW GWDRAG TO I.F.S. (F. LOTT, 22/11/93)
C-----------------------------------------------------------------------
      implicit none
C

C

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
      integer kgwd
      INTEGER KKCRIT(NDLO2),KKCRITH(NDLO2),KCRIT(NDLO2)
     *       ,kdx(NDLO2),ktest(NDLO2)
     *       ,kkenvh(NDLO2),kknu(NDLO2),kknu2(NDLO2)
C
      REAL PAPHM1(NDLO2,klev+1), PSTAB(NDLO2,klev+1),
     *     PRHO  (NDLO2,klev+1), PVPH (NDLO2,klev+1),
     *     PRI   (NDLO2,klev+1), PTFR (NDLO2), PTAU(NDLO2,klev+1),
     *     ptauf (NDLO2,klev+1)
     
      REAL pdmod (NDLO2) , pnu (NDLO2) , psig(NDLO2),
     *     pgamma(NDLO2) , pvar(NDLO2)
     
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
c   declarations pour 'implicit none"
      real zsqr,zalfa,zriw,zalpha,zb,zdel,zdz2n,zdelp,zdelpt

      integer ji,jk,jl,ilevh
      REAL ZDZ2 (NDLO2,nlayermx) , ZNORM(NDLO2) , zoro(NDLO2)
      REAL ZTAU (NDLO2,nlayermx+1)
C
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------


      kidia=1
      kfdia=klon

 100  CONTINUE
C
C
C*    COMPUTATIONAL CONSTANTS.
C     ------------- ----------
C
      ilevh=KLEV/3
C
      DO 400 ji=1,kgwd
      jl=kdx(ji)
      Zoro(JL)=Psig(JL)*Pdmod(JL)/4./max(pvar(jl),1.0)
      ZTAU(JL,KKNU(JL)+1)=PTAU(JL,KKNU(JL)+1)
      ZTAU(JL,KLEV+1)=PTAU(JL,KLEV+1)
  400 CONTINUE
C
      DO 430 JK=KLEV,2,-1
C
C
C*         4.1    CONSTANT WAVE STRESS UNTIL TOP OF THE
C                 BLOCKING LAYER.
  410 CONTINUE
C
      DO 411 ji=1,kgwd
      jl=kdx(ji)
           IF(JK.GE.KKNU2(JL)) THEN
           PTAU(JL,JK)=ZTAU(JL,KLEV+1)
           ENDIF
 411  CONTINUE             
C
C*         4.15   CONSTANT SHEAR STRESS UNTIL THE TOP OF THE
C                 LOW LEVEL FLOW LAYER.
 415  CONTINUE
C        
C
C*         4.2    WAVE DISPLACEMENT AT NEXT LEVEL.
C
  420 CONTINUE
C
      DO 421 ji=1,kgwd
      jl=kdx(ji)
      IF(JK.LT.KKNU2(JL)) THEN
      ZNORM(JL)=gkdrag*PRHO(JL,JK)*SQRT(PSTAB(JL,JK))*PVPH(JL,JK)
     *                                                    *zoro(jl)
      ZDZ2(JL,JK)=PTAU(JL,JK+1)/max(ZNORM(JL),gssec)
      ENDIF
  421 CONTINUE
C
C*         4.3    WAVE RICHARDSON NUMBER, NEW WAVE DISPLACEMENT
C*                AND STRESS:  BREAKING EVALUATION AND CRITICAL 
C                 LEVEL
C
      DO 431 ji=1,kgwd
      jl=kdx(ji)
          IF(JK.LT.KKNU2(JL)) THEN
          IF((PTAU(JL,JK+1).LT.GTSEC).OR.(JK.LE.KCRIT(JL))) THEN
            PTAU(JL,JK)=0.0
          ELSE
               ZSQR=SQRT(PRI(JL,JK))
               ZALFA=SQRT(PSTAB(JL,JK)*ZDZ2(JL,JK))/PVPH(JL,JK)
               ZRIW=PRI(JL,JK)*(1.-ZALFA)/(1+ZALFA*ZSQR)**2
               IF(ZRIW.LT.GRCRIT) THEN
                 ZDEL=4./ZSQR/GRCRIT+1./GRCRIT**2+4./GRCRIT
                 ZB=1./GRCRIT+2./ZSQR
                 ZALPHA=0.5*(-ZB+SQRT(ZDEL))
                 ZDZ2N=(PVPH(JL,JK)*ZALPHA)**2/PSTAB(JL,JK)
                 PTAU(JL,JK)=ZNORM(JL)*ZDZ2N
               ELSE
                 PTAU(JL,JK)=ZNORM(JL)*ZDZ2(JL,JK)
               ENDIF
            PTAU(JL,JK)=MIN(PTAU(JL,JK),PTAU(JL,JK+1))
          ENDIF
          ENDIF
  431 CONTINUE
  
  430 CONTINUE
  440 CONTINUE
  
c     write(*,*) 'ptau'
c     write(*,99) ((ji,ilevh,ptau(ji,ilevh),ji=1,NDLO2),
c    .                  ilevh=1,nlayermx+1)
 99   FORMAT(i3,i3,f15.5)


C  REORGANISATION OF THE STRESS PROFILE
C  IF BREAKING OCCURS AT LOW LEVEL:

      DO 530 ji=1,kgwd
      jl=kdx(ji)
      ZTAU(JL,KKENVH(JL))=PTAU(JL,KKENVH(JL))
      ZTAU(JL,KKCRITH(JL))=PTAU(JL,KKCRITH(JL))
 530  CONTINUE      

      DO 531 JK=1,KLEV
      
      DO 532 ji=1,kgwd
      jl=kdx(ji)
                
         IF(JK.GT.KKCRITH(JL).AND.JK.LT.KKENVH(JL))THEN

          ZDELP=PAPHM1(JL,JK)-PAPHM1(JL,KKENVH(JL))
          ZDELPT=PAPHM1(JL,KKCRITH(JL))-PAPHM1(JL,KKENVH(JL))
          PTAU(JL,JK)=ZTAU(JL,KKENVH(JL)) +
     .                (ZTAU(JL,KKCRITH(JL))-ZTAU(JL,KKENVH(JL)) )*
     .                ZDELP/ZDELPT
     
        ENDIF
            
 532  CONTINUE    
 
 531  CONTINUE        

      RETURN
      END
