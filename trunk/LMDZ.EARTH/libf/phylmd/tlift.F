!
! $Header$
!
        SUBROUTINE TLIFT(P,T,RR,RS,GZ,PLCL,ICB,NK,
     .                  TVP,TPK,CLW,ND,NL,
     .                  DTVPDT1,DTVPDQ1)
C
C     Argument NK ajoute (jyg) = Niveau de depart de la
C     convection
C
        PARAMETER (NA=60)
        REAL GZ(ND),TPK(ND),CLW(ND)
        REAL T(ND),RR(ND),RS(ND),TVP(ND),P(ND)
        REAL DTVPDT1(ND),DTVPDQ1(ND)   ! Derivatives of parcel virtual
C                                   temperature wrt T1 and Q1
C
        REAL CLW_NEW(NA),QI(NA)
        REAL DTPDT1(NA),DTPDQ1(NA)      ! Derivatives of parcel temperature
C                                   wrt T1 and Q1
 
C
        LOGICAL ICE_CONV
C
C   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***
C
c sb        CPD=1005.7
c sb      CPV=1870.0
c sb        CL=4190.0
c sb        CPVMCL=2320.0
c sb        RV=461.5
c sb        RD=287.04
c sb        EPS=RD/RV
c sb        ALV0=2.501E6
ccccccccccccccccccccccc
c constantes coherentes avec le modele du Centre Europeen
c sb      RD = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 28.9644
c sb      RV = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 18.0153
c sb      CPD = 3.5 * RD
c sb      CPV = 4.0 * RV
c sb      CL = 4218.0
c sb      CI=2090.0
c sb      CPVMCL=CL-CPV
c sb      CLMCI=CL-CI
c sb      EPS=RD/RV
c sb      ALV0=2.5008E+06
c sb      ALF0=3.34E+05
 
cccccccccccc
c on utilise les constantes thermo du Centre Europeen: (SB)
c
#include "YOMCST.h"
       GRAVITY = RG !sb: Pr que gravite ne devienne pas humidite!
c
       CPD = RCPD
       CPV = RCPV
       CL = RCW
       CI = RCS
       CPVMCL = CL-CPV
       CLMCI = CL-CI
       EPS = RD/RV
       ALV0 = RLVTT
       ALF0 = RLMLT ! (ALF0 = RLSTT-RLVTT)
c 
cccccccccccccccccccccc
C
C   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***
C
        ICB1=MAX(ICB,2)
        ICB1=MIN(ICB,NL)
C
Cjyg1
CC      CPP=CPD*(1.-RR(1))+RR(1)*CPV
      CPP=CPD*(1.-RR(NK))+RR(NK)*CPV
Cjyg2
      CPINV=1./CPP
Cjyg1
C         ICB may be below condensation level
CCC        DO 100 I=1,ICB1-1
CCC         TPK(I)=T(1)-GZ(I)*CPINV
CCC         TVP(I)=TPK(I)*(1.+RR(1)/EPS)
        DO 50 I=1,ICB1
         CLW(I)=0.0
50      CONTINUE
C
        DO 100 I=NK,ICB1
         TPK(I)=T(NK)-(GZ(I) - GZ(NK))*CPINV
Cjyg1
CCC         TVP(I)=TPK(I)*(1.+RR(NK)/EPS)
         TVP(I)=TPK(I)*(1.+RR(NK)/EPS-RR(NK))
Cjyg2
         DTVPDT1(I) = 1.+RR(NK)/EPS-RR(NK)
         DTVPDQ1(I) = TPK(I)*(1./EPS-1.)
C
Cjyg2
 
  100   CONTINUE
 
C
C    ***  FIND LIFTED PARCEL TEMPERATURE AND MIXING RATIO    ***
C
Cjyg1
CC        AH0=(CPD*(1.-RR(1))+CL*RR(1))*T(1)
CC     $     +RR(1)*(ALV0-CPVMCL*(T(1)-273.15))
        AH0=(CPD*(1.-RR(NK))+CL*RR(NK))*T(NK)
     $     +RR(NK)*(ALV0-CPVMCL*(T(NK)-273.15)) + GZ(NK)
Cjyg2
C
Cjyg1
        IMIN = ICB1
C         If ICB is below LCL, start loop at ICB+1
        IF (PLCL .LT. P(ICB1)) IMIN = MIN(IMIN+1,NL)
C
CCC        DO 300 I=ICB1,NL
        DO 300 I=IMIN,NL
Cjyg2
         ALV=ALV0-CPVMCL*(T(I)-273.15)
         ALF=ALF0+CLMCI*(T(I)-273.15)
 
        RG=RS(I)
        TG=T(I)
C       S=CPD+ALV*ALV*RG/(RV*T(I)*T(I))
Cjyg1
CC        S=CPD*(1.-RR(1))+CL*RR(1)+ALV*ALV*RG/(RV*T(I)*T(I))
        S=CPD*(1.-RR(NK))+CL*RR(NK)+ALV*ALV*RG/(RV*T(I)*T(I))
Cjyg2
        S=1./S
 
        DO 200 J=1,2
Cjyg1
CC         AHG=CPD*TG+(CL-CPD)*RR(1)*TG+ALV*RG+GZ(I)
         AHG=CPD*TG+(CL-CPD)*RR(NK)*TG+ALV*RG+GZ(I)
Cjyg2
        TG=TG+S*(AH0-AHG)
        TC=TG-273.15
        DENOM=243.5+TC
        DENOM=MAX(DENOM,1.0)
C
C       FORMULE DE BOLTON POUR PSAT
C
        ES=6.112*EXP(17.67*TC/DENOM)
        RG=EPS*ES/(P(I)-ES*(1.-EPS))
 
 
  200   CONTINUE
 
Cjyg1
CC        TPK(I)=(AH0-GZ(I)-ALV*RG)/(CPD+(CL-CPD)*RR(1))
        TPK(I)=(AH0-GZ(I)-ALV*RG)/(CPD+(CL-CPD)*RR(NK))
Cjyg2
C       TPK(I)=(AH0-GZ(I)-ALV*RG-(CL-CPD)*T(I)*RR(1))/CPD
 
Cjyg1
CC        CLW(I)=RR(1)-RG
        CLW(I)=RR(NK)-RG
Cjyg2
        CLW(I)=MAX(0.0,CLW(I))
Cjyg1
CCC        TVP(I)=TPK(I)*(1.+RG/EPS)
        TVP(I)=TPK(I)*(1.+RG/EPS-RR(NK))
Cjyg2
C
Cjyg1       Derivatives
C
        DTPDT1(I) = CPD*S
        DTPDQ1(I) = ALV*S
C
         DTVPDT1(I) = DTPDT1(I)*(1. + RG/EPS -
     .           RR(NK) + ALV*RG/(RD*TPK(I)) )
        DTVPDQ1(I) = DTPDQ1(I)*(1. + RG/EPS -
     .           RR(NK) + ALV*RG/(RD*TPK(I)) ) - TPK(I)
C
Cjyg2
 
  300   CONTINUE
C
      ICE_CONV = .FALSE.

      IF (ICE_CONV) THEN
C
CJAM
C       RAJOUT DE LA PROCEDURE ICEFRAC
C
c sb        CALL ICEFRAC(T,CLW,CLW_NEW,QI,ND,NL)
 
        DO 400 I=ICB1,NL
        IF (T(I).LT.263.15) THEN
        TG=TPK(I)
        TC=TPK(I)-273.15
        DENOM=243.5+TC
        ES=6.112*EXP(17.67*TC/DENOM)
        ALV=ALV0-CPVMCL*(T(I)-273.15)
        ALF=ALF0+CLMCI*(T(I)-273.15)
 
        DO J=1,4
        ESI=EXP(23.33086-(6111.72784/TPK(I))+0.15215*LOG(TPK(I)))
        QSAT_NEW=EPS*ESI/(P(I)-ESI*(1.-EPS))
CCC        SNEW= CPD*(1.-RR(1))+CL*RR(1)+ALV*ALV*QSAT_NEW/(RV*TPK(I)*TPK(I))
        SNEW= CPD*(1.-RR(NK))+CL*RR(NK)
     .        +ALV*ALV*QSAT_NEW/(RV*TPK(I)*TPK(I))
C
        SNEW=1./SNEW
        TPK(I)=TG+(ALF*QI(I)+ALV*RG*(1.-(ESI/ES)))*SNEW
c@$$        PRINT*,'################################'
c@$$        PRINT*,TPK(I)
c@$$        PRINT*,(ALF*QI(I)+ALV*RG*(1.-(ESI/ES)))*SNEW
        ENDDO
CCC        CLW(I)=RR(1)-QSAT_NEW
        CLW(I)=RR(NK)-QSAT_NEW
        CLW(I)=MAX(0.0,CLW(I))
Cjyg1
CCC        TVP(I)=TPK(I)*(1.+QSAT_NEW/EPS)
        TVP(I)=TPK(I)*(1.+QSAT_NEW/EPS-RR(NK))
Cjyg2
        ELSE
        CONTINUE
        ENDIF
 
  400   CONTINUE
C
      ENDIF
C
 
******************************************************
** BK :  RAJOUT DE LA TEMPERATURE DES ASCENDANCES
**   NON DILUES AU  NIVEAU KLEV = ND
**   POSONS LE ENVIRON EGAL A CELUI DE KLEV-1
********************************************************
 
      TPK(NL+1)=TPK(NL)
 
*******************************************************

      RG = GRAVITY  ! RG redevient la gravite de YOMCST (sb)
 
 
        RETURN
        END
 
 
 
 
 
 
 
 
