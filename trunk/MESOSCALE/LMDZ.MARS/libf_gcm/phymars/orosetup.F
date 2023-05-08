      SUBROUTINE OROSETUP
     *         ( klon   , klev  , KTEST
     *         , KKCRIT, KKCRITH, KCRIT, KSECT , KKHLIM
     *         , kkenvh, kknu  , kknu2
     *         , PAPHM1, PAPM1 , PUM1   , PVM1 , PTM1  , PGEOM1, pvaror
     *         , PRHO  , PRI   , PSTAB  , PTAU , PVPH  ,ppsi, pzdep
     *         , PULOW , PVLOW  
     *         , Ptheta, pgamma,  pnu  ,  pd1  ,  pd2  ,pdmod  )
C
C**** *GWSETUP*
C
C     PURPOSE.
C     --------
C
C**   INTERFACE.
C     ----------
C          FROM *ORODRAG*
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
C     F.LOTT  FOR THE NEW-GWDRAG SCHEME NOVEMBER 1993
C
C-----------------------------------------------------------------------
      implicit none
C

#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
      integer klon,klev,kidia,kfdia

#include "comcstfi.h"
#include "yoegwd.h"

C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER KKCRIT(NDLO2),KKCRITH(NDLO2),KCRIT(NDLO2),KSECT(NDLO2),
     *        KKHLIM(NDLO2),KTEST(NDLO2),kkenvh(NDLO2)

C
      REAL PAPHM1(NDLO2,KLEV+1),PAPM1(NDLO2,KLEV),PUM1(NDLO2,KLEV),
     *     PVM1(NDLO2,KLEV),PTM1(NDLO2,KLEV),PGEOM1(NDLO2,KLEV),
     *     PRHO(NDLO2,KLEV+1),PRI(NDLO2,KLEV+1),PSTAB(NDLO2,KLEV+1),
     *     PTAU(NDLO2,KLEV+1),PVPH(NDLO2,KLEV+1),ppsi(NDLO2,klev+1),
     *     pzdep(NDLO2,klev)
       REAL PULOW(NDLO2),PVLOW(NDLO2),ptheta(NDLO2),pgamma(NDLO2),
     *     pnu(NDLO2),
     *     pd1(NDLO2),pd2(NDLO2),pdmod(NDLO2)
      real pvaror(NDLO2)
C
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
C
      LOGICAL LL1(NDLO2,nlayermx+1)
      integer kknu(NDLO2),kknu2(NDLO2),kknub(NDLO2),kknul(NDLO2),
     *        kentp(NDLO2),ncount(NDLO2)  
C
      REAL ZHCRIT(NDLO2,nlayermx),ZNCRIT(NDLO2,nlayermx),
     *     ZVPF(NDLO2,nlayermx), ZDP(NDLO2,nlayermx)
      REAL ZNORM(NDLO2),zpsi(NDLO2),zb(NDLO2),zc(NDLO2),
     *      zulow(NDLO2),zvlow(NDLO2),znup(NDLO2),znum(NDLO2)
C
c   declarations pour "implicit none"
      integer jk,jl,ilevm1,ilevm2,ilevh
      real zu,zphi,zcons1,zcons2,zcons3,zwind,zdwind,zhgeo
      real zvt1,zvt2,zst,zvar,zdelp,zstabm,zstabp,zrhom,zrhop
      real alpha,zggeenv,zggeom1,zgvar
      logical lo

C     ------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
c     print *,' entree gwsetup'
 100  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         1.1   COMPUTATIONAL CONSTANTS
C                -----------------------
C

      kidia=1
	  kfdia=klon

 110  CONTINUE
C
      ILEVM1=KLEV-1
      ILEVM2=KLEV-2
      ILEVH =KLEV/3
C
      ZCONS1=1./r
cold  ZCONS2=G**2/CPD
      ZCONS2=g**2/cpp
cold  ZCONS3=1.5*API
      ZCONS3=1.5*PI
C
C
C     ------------------------------------------------------------------
C
C*         2.
C                --------------
C
 200  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
C*                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
C*                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.
C
C
C
      DO 2001 JL=kidia,kfdia
      kknu(JL)    =klev
      kknu2(JL)   =klev
      kknub(JL)   =klev
      kknul(JL)   =klev
      pgamma(JL) =max(pgamma(jl),gtsec)
      ll1(jl,klev+1)=.false.
 2001 CONTINUE
C
C*      DEFINE TOP OF LOW LEVEL FLOW
C       ----------------------------
      DO 2002 JK=KLEV,ilevh,-1
      DO 2003 JL=kidia,kfdia
      LO=(PAPHM1(JL,JK)/PAPHM1(JL,KLEV+1)).GE.GSIGCR
      IF(LO) THEN
        KKCRIT(JL)=JK
      ENDIF
      ZHCRIT(JL,JK)=4.*pvaror(JL)
      ZHGEO=PGEOM1(JL,JK)/g
      ll1(JL,JK)=(ZHGEO.GT.ZHCRIT(JL,JK))
      IF(ll1(JL,JK).NEQV.ll1(JL,JK+1)) THEN
        kknu(JL)=JK
      ENDIF
 2003 CONTINUE
 2002 CONTINUE
      DO 2004 JK=KLEV,ilevh,-1
      DO 2005 JL=kidia,kfdia
      ZHCRIT(JL,JK)=3.*pvaror(JL)
      ZHGEO=PGEOM1(JL,JK)/g
      ll1(JL,JK)=(ZHGEO.GT.ZHCRIT(JL,JK))
      IF(ll1(JL,JK).NEQV.ll1(JL,JK+1)) THEN
        kknu2(JL)=JK
      ENDIF
 2005 CONTINUE
 2004 CONTINUE
      DO 2006 JK=KLEV,ilevh,-1
      DO 2007 JL=kidia,kfdia
      ZHCRIT(JL,JK)=2.*pvaror(JL)
      ZHGEO=PGEOM1(JL,JK)/g
      ll1(JL,JK)=(ZHGEO.GT.ZHCRIT(JL,JK))
      IF(ll1(JL,JK).NEQV.ll1(JL,JK+1)) THEN
        kknub(JL)=JK
      ENDIF
 2007 CONTINUE
 2006 CONTINUE
      DO 2008 JK=KLEV,ilevh,-1
      DO 2009 JL=kidia,kfdia
      ZHCRIT(JL,JK)=pvaror(JL)
      ZHGEO=PGEOM1(JL,JK)/g
      ll1(JL,JK)=(ZHGEO.GT.ZHCRIT(JL,JK))
      IF(ll1(JL,JK).NEQV.ll1(JL,JK+1)) THEN
        kknul(JL)=JK
      ENDIF
 2009 CONTINUE
 2008 CONTINUE
C
      do 2010 jl=kidia,kfdia  
      kknu(jl)=min(kknu(jl),nktopg)
      kknub(jl)=min(kknub(jl),nktopg)
      if(kknub(jl).eq.nktopg) kknul(jl)=klev
C
C     CHANGE IN HERE TO STOP KKNUL=KKNUB
C
      if(kknul(jl).le.kknub(jl)) kknul(jl)=nktopg
 2010 continue      
C

 210  CONTINUE
C
C
CC*     INITIALIZE VARIOUS ARRAYS
C
      DO 2107 JL=kidia,kfdia
      PRHO(JL,KLEV+1)  =0.0
      PSTAB(JL,KLEV+1) =0.0
      PSTAB(JL,1)      =0.0
      PRI(JL,KLEV+1)   =9999.0
      ppsi(JL,KLEV+1)  =0.0
      PRI(JL,1)        =0.0
      PVPH(JL,1)       =0.0
      PULOW(JL)        =0.0
      PVLOW(JL)        =0.0
      zulow(JL)        =0.0
      zvlow(JL)        =0.0
      KKCRITH(JL)      =KLEV
      KKenvH(JL)       =KLEV
      Kentp(JL)        =KLEV
      KCRIT(JL)        =1
      ncount(JL)       =0
      ll1(JL,klev+1)   =.false.
 2107 CONTINUE
C
C*     DEFINE LOW-LEVEL FLOW
C      ---------------------
C
      DO 223 JK=KLEV,2,-1
      DO 222 JL=kidia,kfdia
      IF(KTEST(JL).EQ.1) THEN
        ZDP(JL,JK)=PAPM1(JL,JK)-PAPM1(JL,JK-1)
        PRHO(JL,JK)=2.*PAPHM1(JL,JK)*ZCONS1/(PTM1(JL,JK)+PTM1(JL,JK-1))
        PSTAB(JL,JK)=2.*ZCONS2/(PTM1(JL,JK)+PTM1(JL,JK-1))*
     *  (1.-cpp*PRHO(JL,JK)*(PTM1(JL,JK)-PTM1(JL,JK-1))/ZDP(JL,JK))
        PSTAB(JL,JK)=MAX(PSTAB(JL,JK),GSSEC)
      ENDIF
  222 CONTINUE
  223 CONTINUE
C
C********************************************************************
C
C*     DEFINE blocked FLOW
C      -------------------
      DO 2115 JK=klev,ilevh,-1
      DO 2116 JL=kidia,kfdia
      if(jk.ge.kknub(jl).and.jk.le.kknul(jl)) then
        pulow(JL)=pulow(JL)+PUM1(JL,JK)*(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
        pvlow(JL)=pvlow(JL)+PVM1(JL,JK)*(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
      end if
 2116 CONTINUE
 2115 CONTINUE
      DO 2110 JL=kidia,kfdia
      pulow(JL)=pulow(JL)/(PAPHM1(JL,Kknul(jl)+1)-PAPHM1(JL,kknub(jl)))
      pvlow(JL)=pvlow(JL)/(PAPHM1(JL,Kknul(jl)+1)-PAPHM1(JL,kknub(jl)))
      ZNORM(JL)=MAX(SQRT(PULOW(JL)**2+PVLOW(JL)**2),GVSEC)
      PVPH(JL,KLEV+1)=ZNORM(JL)
 2110 CONTINUE
C
C*******  SETUP OROGRAPHY AXES AND DEFINE PLANE OF PROFILES  *******
C
      DO 2112 JL=kidia,kfdia
      LO=(PULOW(JL).LT.GVSEC).AND.(PULOW(JL).GE.-GVSEC)
      IF(LO) THEN
        ZU=PULOW(JL)+2.*GVSEC
      ELSE
        ZU=PULOW(JL)
      ENDIF
      Zphi=ATAN(PVLOW(JL)/ZU)
      ppsi(jl,klev+1)=ptheta(jl)*pi/180.-zphi
      zb(jl)=1.-0.18*pgamma(jl)-0.04*pgamma(jl)**2
      zc(jl)=0.48*pgamma(jl)+0.3*pgamma(jl)**2
      pd1(jl)=zb(jl)-(zb(jl)-zc(jl))*(sin(ppsi(jl,klev+1))**2)
      pd2(jl)=(zb(jl)-zc(jl))*sin(ppsi(jl,klev+1))*cos(ppsi(jl,klev+1))
      pdmod(jl)=sqrt(pd1(jl)**2+pd2(jl)**2)
 2112 CONTINUE
C
C  ************ DEFINE FLOW IN PLANE OF LOWLEVEL STRESS *************
C
      DO 213 JK=1,KLEV
      DO 212 JL=kidia,kfdia
      IF(KTEST(JL).EQ.1)  THEN
        ZVt1       =PULOW(JL)*PUM1(JL,JK)+PVLOW(JL)*PVM1(JL,JK)
        ZVt2       =-PvLOW(JL)*PUM1(JL,JK)+PuLOW(JL)*PVM1(JL,JK)
        ZVPF(JL,JK)=(zvt1*pd1(jl)+zvt2*pd2(JL))/(znorm(jl)*pdmod(jl))
      ENDIF
      PTAU(JL,JK)  =0.0
      Pzdep(JL,JK) =0.0
      Ppsi(JL,JK)  =0.0
      ll1(JL,JK)   =.FALSE.
  212 CONTINUE
  213 CONTINUE
      DO 215 JK=2,KLEV
      DO 214 JL=kidia,kfdia
      IF(KTEST(JL).EQ.1) THEN
        ZDP(JL,JK)=PAPM1(JL,JK)-PAPM1(JL,JK-1)
        PVPH(JL,JK)=((PAPHM1(JL,JK)-PAPM1(JL,JK-1))*ZVPF(JL,JK)+
     *            (PAPM1(JL,JK)-PAPHM1(JL,JK))*ZVPF(JL,JK-1))
     *            /ZDP(JL,JK)
        IF(PVPH(JL,JK).LT.GVSEC) THEN
          PVPH(JL,JK)=GVSEC
          KCRIT(JL)=JK
        ENDIF
      ENDIF
  214 CONTINUE
  215 CONTINUE
C
C
C*         2.2     BRUNT-VAISALA FREQUENCY AND DENSITY AT HALF LEVELS.
C
  220 CONTINUE
C
      DO 2211 JK=ilevh,KLEV
      DO 221 JL=kidia,kfdia
      IF(KTEST(JL).EQ.1) THEN
      IF(jk.ge.(kknub(jl)+1).and.jk.le.kknul(jl)) THEN
           ZST=ZCONS2/PTM1(JL,JK)*(1.-cpp*PRHO(JL,JK)*
     *                   (PTM1(JL,JK)-PTM1(JL,JK-1))/ZDP(JL,JK))
           PSTAB(JL,KLEV+1)=PSTAB(JL,KLEV+1)+ZST*ZDP(JL,JK)
           PSTAB(JL,KLEV+1)=MAX(PSTAB(JL,KLEV+1),GSSEC)
           PRHO(JL,KLEV+1)=PRHO(JL,KLEV+1)+PAPHM1(JL,JK)*2.*ZDP(JL,JK)
     *                   *ZCONS1/(PTM1(JL,JK)+PTM1(JL,JK-1))
      ENDIF
      ENDIF
  221 CONTINUE
 2211 CONTINUE
C
      DO 2212 JL=kidia,kfdia
C*****************************************************************************
C
C     O.K. THERE IS A POSSIBLE PROBLEM HERE. IF KKNUL=KKNUB THEN
C     DIVISION BY ZERO OCCURS. I HAVE PUT A FIX IN HERE BUT WILL ASK FRANCOIS 
C     LOTT ABOUT IT IN PARIS.
C
C     MAT COLLINS 30.1.96
C
C     ALSO IF THIS IS THE CASE PSTAB AND PRHO ARE NOT DEFINED AT KLEV+1
C     SO I HAVE ADDED THE ELSE
C
C*****************************************************************************
      IF (KKNUL(JL).NE.KKNUB(JL)) THEN
        PSTAB(JL,KLEV+1)=PSTAB(JL,KLEV+1)/(PAPM1(JL,Kknul(jl))
     *                                          -PAPM1(JL,kknub(jl)))
        PRHO(JL,KLEV+1)=PRHO(JL,KLEV+1)/(PAPM1(JL,Kknul(jl))
     *                                          -PAPM1(JL,kknub(jl)))
      ELSE
      WRITE(*,*) 'OROSETUP: KKNUB=KKNUL= ',KKNUB(JL),' AT JL= ',JL
        PSTAB(JL,KLEV+1)=PSTAB(JL,KLEV)
        PRHO(JL,KLEV+1)=PRHO(JL,KLEV)
      ENDIF
        ZVAR=PVARor(JL)
 2212 CONTINUE
C
C*         2.3     MEAN FLOW RICHARDSON NUMBER.
C*                 AND CRITICAL HEIGHT FOR FROUDE LAYER
C
  230 CONTINUE
C
      DO 232 JK=2,KLEV
      DO 231 JL=kidia,kfdia
      IF(KTEST(JL).EQ.1) THEN
        ZDWIND=MAX(ABS(ZVPF(JL,JK)-ZVPF(JL,JK-1)),GVSEC)
        PRI(JL,JK)=PSTAB(JL,JK)*(ZDP(JL,JK)
     *          /(g*PRHO(JL,JK)*ZDWIND))**2
        PRI(JL,JK)=MAX(PRI(JL,JK),GRCRIT)
      ENDIF
  231 CONTINUE
  232 CONTINUE
C
C
C*      DEFINE TOP OF 'envelope' layer
C       ----------------------------

      DO 233 JL=kidia,kfdia
      pnu (jl)=0.0
      znum(jl)=0.0
 233  CONTINUE
      
      DO 234 JK=2,KLEV-1
      DO 234 JL=kidia,kfdia
      
      IF(KTEST(JL).EQ.1) THEN
       
      IF (JK.GE.KKNU2(JL)) THEN
          
            ZNUM(JL)=PNU(JL)
            ZWIND=(pulow(JL)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/
     *            max(sqrt(pulow(jl)**2+pvlow(jl)**2),gvsec)
            ZWIND=max(sqrt(zwind**2),gvsec)
            ZDELP=PAPHM1(JL,JK+1)-PAPHM1(JL,JK)
            ZSTABM=SQRT(MAX(PSTAB(JL,JK  ),GSSEC))
            ZSTABP=SQRT(MAX(PSTAB(JL,JK+1),GSSEC))
            ZRHOM=PRHO(JL,JK  )
            ZRHOP=PRHO(JL,JK+1)
            PNU(JL) = PNU(JL) + (ZDELP/g)*
     *            ((zstabp/zrhop+zstabm/zrhom)/2.)/ZWIND     
            IF((ZNUM(JL).LE.GFRCRIT).AND.(PNU(JL).GT.GFRCRIT)
     *                          .AND.(KKENVH(JL).EQ.KLEV))
     *      KKENVH(JL)=JK
     
      ENDIF    

      ENDIF
      
 234  CONTINUE

C  CALCULATION OF A DYNAMICAL MIXING HEIGHT FOR THE BREAKING
C  OF GRAVITY WAVES:

              
      DO 235 JL=kidia,kfdia
      znup(jl)=0.0
      znum(jl)=0.0
 235  CONTINUE

      DO 236 JK=KLEV-1,2,-1
      DO 236 JL=kidia,kfdia
      
      IF(KTEST(JL).EQ.1) THEN
        
      IF (JK.LT.KKENVH(JL)) THEN

            ZNUM(JL)=ZNUP(JL)
            ZWIND=(pulow(JL)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/
     *            max(sqrt(pulow(jl)**2+pvlow(jl)**2),gvsec)
            ZWIND=max(sqrt(zwind**2),gvsec)
            ZDELP=PAPHM1(JL,JK+1)-PAPHM1(JL,JK)
            ZSTABM=SQRT(MAX(PSTAB(JL,JK  ),GSSEC))
            ZSTABP=SQRT(MAX(PSTAB(JL,JK+1),GSSEC))
            ZRHOM=PRHO(JL,JK  )
            ZRHOP=PRHO(JL,JK+1)
            ZNUP(JL) = ZNUP(JL) + (ZDELP/g)*
     *            ((zstabp/zrhop+zstabm/zrhom)/2.)/ZWIND     
            IF((ZNUM(JL).LE.1.5).AND.(ZNUP(JL).GT.1.5)
     *                          .AND.(KKCRITH(JL).EQ.KLEV))
     *      KKCRITH(JL)=JK
     
      ENDIF
      
      ENDIF
      
 236  CONTINUE
 
      DO 237 JL=KIDIA,KFDIA
      KKCRITH(JL)=MIN0(KKCRITH(JL),KKNU(JL))
 237  CONTINUE
c
c     directional info for flow blocking ************************* 
c
      do 251 jk=ilevh,klev    
      DO 252 JL=kidia,kfdia
      IF(jk.ge.kkenvh(jl)) THEN
      LO=(PUm1(JL,jk).LT.GVSEC).AND.(PUm1(JL,jk).GE.-GVSEC)
      IF(LO) THEN
        ZU=PUm1(JL,jk)+2.*GVSEC
      ELSE
        ZU=PUm1(JL,jk)
      ENDIF
       Zphi=ATAN(PVm1(JL,jk)/ZU)
       ppsi(jl,jk)=ptheta(jl)*pi/180.-zphi
      end if
 252  continue
 251  CONTINUE
c      forms the vertical 'leakiness' **************************

      alpha=3.
      
      DO 254  JK=ilevh,klev
      DO 253  JL=kidia,kfdia
      IF(jk.ge.kkenvh(jl)) THEN
        zggeenv=AMAX1(1.,
     *          (pgeom1(jl,kkenvh(jl))+pgeom1(jl,kkenvh(jl)-1))/2.)      
        zggeom1=AMAX1(pgeom1(jl,jk),1.)
        zgvar=amax1(pvaror(jl)*g,1.)     
        pzdep(jl,jk)=SQRT((zggeenv-zggeom1)/(zggeom1+zgvar))      
      end if
 253  CONTINUE
 254  CONTINUE

 260  CONTINUE

      RETURN
      END
