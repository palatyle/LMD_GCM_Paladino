      SUBROUTINE SWR_TOON ( KDLON, KFLEV, KNU
     S     ,  aerosol,QVISsQREF3d,omegaVIS3d,gVIS3d
     &     ,  albedo,PDSIG,PPSOL,PRMU,PSEC
     S     ,  PFD,PFU )

      use dimradmars_mod, only: sunfr, ndlo2, nsun, nflev, 
     &                          ndlon, naerkind
      use yomlw_h, only: nlaylte
      
      IMPLICIT NONE
C     
#include "callkeys.h"

C     
C   SWR - Continuum scattering computations
C     
C     PURPOSE.
C     --------
C     Computes the reflectivity and transmissivity in case oF
C     Continuum scattering
c     F. Forget (1999)
c    
c     Modified by Tran The Trung, using radiative transfer code 
c     of Toon 1981.  
C     
C     IMPLICIT ARGUMENTS :
C     --------------------
C     
C     ==== INPUTS ===
c
c    KDLON :  number of horizontal grid points
c    KFLEV :  number of vertical layers
c    KNU   :   Solar band # (1 or 2)
c   aerosol               aerosol extinction optical depth
c                         at reference wavelength "longrefvis" set
c                         in dimradmars_mod , in each layer, for one of
c                         the "naerkind" kind of aerosol optical properties.
c    albedo   hemispheric surface albedo
c                         albedo (i,1) : mean albedo for solar band#1
c                                        (see below)
c                         albedo (i,2) : mean albedo for solar band#2
c                                        (see below)
c    PDSIG      layer thickness in sigma coordinates
c    PPSOL       Surface pressure (Pa)
c    PRMU:  cos of solar zenith angle (=1 when sun at zenith)
c           (CORRECTED for high zenith angle (atmosphere), unlike mu0)
c    PSEC   =1./PRMU

C     ==== OUTPUTS ===
c
c    PFD : downward flux in spectral band #INU in a given mesh
c         (normalized to the total incident flux at the top of the atmosphere)
c    PFU : upward flux in specatral band #INU in a given mesh
c         (normalized to the total incident flux at the top of the atmosphere)
C
C     
C     METHOD.
C     -------
C     
C     Computes continuum fluxes corresponding to aerosoL
C     Or/and rayleigh scattering (no molecular gas absorption)
C     
C-----------------------------------------------------------------------
C     
C     
C-----------------------------------------------------------------------
C     
     
C     ARGUMENTS
C     ---------
      INTEGER KDLON, KFLEV, KNU
      REAL aerosol(NDLO2,KFLEV,naerkind), albedo(NDLO2,2), 
     S     PDSIG(NDLO2,KFLEV),PSEC(NDLO2)

      REAL QVISsQREF3d(NDLO2,KFLEV,nsun,naerkind)  
      REAL omegaVIS3d(NDLO2,KFLEV,nsun,naerkind)   
      REAL gVIS3d(NDLO2,KFLEV,nsun,naerkind)

      REAL PPSOL(NDLO2)
      REAL PFD(NDLO2,KFLEV+1),PFU(NDLO2,KFLEV+1)
      REAL PRMU(NDLO2)

C     LOCAL ARRAYS
C     ------------
 
      INTEGER jk,jl,jae
      REAL ZTRAY, ZRATIO,ZGAR, ZFF
      REAL ZCGAZ(NDLO2,NFLEV),ZPIZAZ(NDLO2,NFLEV),ZTAUAZ(NDLO2,NFLEV)
      REAL  ZRAYL(NDLON)

c     Part added by Tran The Trung
c     inputs to gfluxv 
      REAL*8 DTDEL(nlaylte), WDEL(nlaylte), CDEL(nlaylte)
c     outputs of gfluxv
      REAL*8 FP(nlaylte+1), FM(nlaylte+1)      
c     normalization of top downward flux
      REAL*8 norm
c     End part added by Tran The Trung

c     Function
c     --------
      real CVMGT

c Computing TOTAL single scattering parameters by adding
c  properties of all the NAERKIND kind of aerosols

      DO JK = 1 , nlaylte
         DO  JL = 1 , KDLON
            ZCGAZ(JL,JK) = 0.
            ZPIZAZ(JL,JK) =  0.
            ZTAUAZ(JL,JK) = 0.
         END DO 
         DO 106 JAE=1,naerkind
            DO 105 JL = 1 , KDLON
c              Mean Extinction optical depth in the spectral band
c              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ZTAUAZ(JL,JK)=ZTAUAZ(JL,JK)
     S              +aerosol(JL,JK,JAE)*QVISsQREF3d(JL,JK,KNU,JAE)
c              Single scattering albedo
c              ~~~~~~~~~~~~~~~~~~~~~~~~ 
               ZPIZAZ(JL,JK)=ZPIZAZ(JL,JK)+aerosol(JL,JK,JAE)*
     S           QVISsQREF3d(JL,JK,KNU,JAE)*
     &           omegaVIS3d(JL,JK,KNU,JAE)
c              Assymetry factor
c              ~~~~~~~~~~~~~~~~
               ZCGAZ(JL,JK) =  ZCGAZ(JL,JK) +aerosol(JL,JK,JAE)*
     S           QVISsQREF3d(JL,JK,KNU,JAE)*
     &           omegaVIS3d(JL,JK,KNU,JAE)*gVIS3d(JL,JK,KNU,JAE)
 105        CONTINUE
 106     CONTINUE
      END DO
C     
      DO JK = 1 , nlaylte
         DO JL = 1 , KDLON
            ZCGAZ(JL,JK) = CVMGT( 0., ZCGAZ(JL,JK) / ZPIZAZ(JL,JK),
     S            (ZPIZAZ(JL,JK).EQ.0) )
            ZPIZAZ(JL,JK) = CVMGT( 1., ZPIZAZ(JL,JK) / ZTAUAZ(JL,JK),
     S           (ZTAUAZ(JL,JK).EQ.0) )
         END DO
      END DO

C     --------------------------------
C     INCLUDING RAYLEIGH SCATERRING 
C     -------------------------------
      if (rayleigh) then 

        call swrayleigh(kdlon,knu,ppsol,prmu,ZRAYL)

c       Modifying mean aerosol parameters to account rayleigh scat by gas:

        DO JK = 1 , nlaylte
           DO JL = 1 , KDLON
c             Rayleigh opacity in each layer :
              ZTRAY = ZRAYL(JL) * PDSIG(JL,JK)
c             ratio Tau(rayleigh) / Tau (total)
              ZRATIO = ZTRAY / (ZTRAY + ZTAUAZ(JL,JK))
              ZGAR = ZCGAZ(JL,JK)
              ZFF = ZGAR * ZGAR
                ZTAUAZ(JL,JK)=ZTRAY+ZTAUAZ(JL,JK)*(1.-ZPIZAZ(JL,JK)*ZFF)
              ZCGAZ(JL,JK) = ZGAR * (1. - ZRATIO) / (1. + ZGAR)
              ZPIZAZ(JL,JK) =ZRATIO+(1.-ZRATIO)*ZPIZAZ(JL,JK)*(1.-ZFF)
     S           / (1. -ZPIZAZ(JL,JK) * ZFF)
           END DO
        END DO
      end if

c     Part added by Tran The Trung
c     new radiative transfer

      do JL = 1, KDLON

c     assign temporary inputs
        do JK = 1, nlaylte
           jae = nlaylte+1-JK
           DTDEL(JK) = real(ZTAUAZ(JL,jae),8)
           WDEL(JK) = real(ZPIZAZ(JL,jae),8)
           CDEL(JK) = real(ZCGAZ(JL,jae),8)
        end do

c     call gfluxv
        call gfluxv(nlaylte, DTDEL,WDEL,CDEL,
     S              real(PRMU(JL),8),
     S              real(albedo(JL,KNU),8),
     S              FP,FM)

c     assign output
        norm = FM(1)
c     here we can have a check of norm not being 0.0
c     however it would never happen in practice, 
c     so we can comment out
c        if (norm .gt. 0.0) then
          do JK = 1, nlaylte+1
             jae = nlaylte+2-JK
             PFU(JL,JK) = sunfr(KNU)*real(FP(jae)/norm,4)
             PFD(JL,JK) = sunfr(KNU)*real(FM(jae)/norm,4) 
          end do
c        elseif (norm .eq. 0.0) then
c          do JK = 1, nlaylte+1
c             PFU(JL,JK) = 0.0 
c             PFD(JL,JK) = 0.0
c          end do
c        else
c          stop "Error: top downward visible flux is negative!"
c        end if
 
      end do
c     End part added by Tran The Trung

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GFLUXV(NAYER,DTDEL,WDEL,CDEL,UBAR0,RSF
     &                 ,FP,FM)
      IMPLICIT NONE
C  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITONS
C  FOR THE VISIBLE  FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
C  THE LEVELS. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
C  MEASURED FROM THE TOP OF EACH LAYER.  (DTAU) TOP OF EACH LAYER HAS
C  OPTICAL DEPTH TAU(N). IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
C  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
C  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.
C  THIS SUBROUTINE DIFFERS FROM ITS IR CONTERPART IN THAT HERE WE SOLVE FOR
C  THE FLUXES DIRECTLY USING THE GENERALIZED NOTATION OF MEADOR AND WEAVOR
C  J.A.S., 37, 630-642, 1980.
C  THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY
C  VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER

C  THIS VERSION HAS BEEN MODIFIED BY TRAN THE TRUNG WITH:
C   1. Simplified input & output for swr.F subroutine in LMDZ.MARS gcm model
C   2. Use delta function to modify optical properties
C   3. Use delta-eddington G1, G2, G3 parameters
C   4. Optimized for speed

C  INPUTS: 
      INTEGER NAYER     
c   NAYER = number of layer
c            first layer is at top 
c            last layer is at bottom
      REAL*8 DTDEL(NAYER), WDEL(NAYER), CDEL(NAYER) 
c   DTDEL = optical thickness of layer
c   WDEL = single scattering of layer
c   CDEL = assymetry parameter
      REAL*8 UBAR0, RSF
c   UBAR0 = absolute value of cosine of solar zenith angle
c   RSF = surface albedo
C  OUTPUTS:
      REAL*8 FP(NAYER+1), FM(NAYER+1) 
c   FP = flux up
c   FM = flux down
C  PRIVATES:
      INTEGER J,NL,NLEV
!!!! AS+JBM 03/2010 BUG BUG si trop niveaux verticaux (LES)
!!!!                ET PAS BESOIN DE HARDWIRE SALE ICI  !   
!!!! CORRIGER CE BUG AMELIORE EFFICACITE ET FLEXIBILITE      
      !! PARAMETER (NL=201) 
      !! C THIS VALUE (201) MUST BE .GE. 2*NAYER
      REAL*8 BSURF,AP,AM,DENOM,EM,EP,G4
      !! REAL*8 W0(NL), COSBAR(NL), DTAU(NL), TAU(NL)
      !! REAL*8 LAMDA(NL),XK1(NL),XK2(NL)
      !! REAL*8 G1(NL),G2(NL),G3(NL)
      !! REAL*8 GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL)
      !! REAL*8 E1(NL),E2(NL),E3(NL),E4(NL)
      REAL*8 W0(2*NAYER), COSBAR(2*NAYER), DTAU(2*NAYER), TAU(2*NAYER)  
      REAL*8 LAMDA(2*NAYER),XK1(2*NAYER),XK2(2*NAYER)
      REAL*8 G1(2*NAYER),G2(2*NAYER),G3(2*NAYER)
      REAL*8 GAMA(2*NAYER),CP(2*NAYER),CM(2*NAYER),CPM1(2*NAYER)
      REAL*8 CMM1(2*NAYER)
      REAL*8 E1(2*NAYER),E2(2*NAYER),E3(2*NAYER),E4(2*NAYER)

      NL = 2*NAYER  !!! AS+JBM 03/2010 
      NLEV = NAYER+1
      
C  TURN ON THE DELTA-FUNCTION IF REQUIRED HERE
c      TAU(1) = 0.0
c      DO J=1,NAYER
c        W0(J)=WDEL(J)
c        COSBAR(J)=CDEL(J)
c        DTAU(J)=DTDEL(J)
c        TAU(J+1)=TAU(J)+DTAU(J)
c      END DO
C  FOR THE DELTA FUNCTION  HERE...
      TAU(1) = 0.0
      DO J=1,NAYER
        COSBAR(J)=CDEL(J)/(1.+CDEL(J))
        W0(J)=1.-WDEL(J)*CDEL(J)**2
        DTAU(J)=DTDEL(J)*W0(J)
        W0(J)=WDEL(J)*(1.-CDEL(J)**2)/W0(J)
        TAU(J+1)=TAU(J)+DTAU(J)
      END DO
      
c     Optimization, this is the major speed gain 
      TAU(1) = 1.0
      do J = 2, NAYER+1
        TAU(J) = EXP(-TAU(J)/UBAR0)
      end do
      BSURF = RSF*UBAR0*TAU(NLEV)
C  WE GO WITH THE HEMISPHERIC CONSTANT APPROACH
C  AS DEFINED BY M&W - THIS IS THE WAY THE IR IS DONE
      DO J=1,NAYER
c        Optimization: ALPHA not used 
c        ALPHA(J)=SQRT( (1.-W0(J))/(1.-W0(J)*COSBAR(J)) )
C  SET OF CONSTANTS DETERMINED BY DOM
c      G1(J)= (SQRT(3.)*0.5)*(2. - W0(J)*(1.+COSBAR(J)))
c      G2(J)= (SQRT(3.)*W0(J)*0.5)*(1.-COSBAR(J))
c      G3(J)=0.5*(1.-SQRT(3.)*COSBAR(J)*UBAR0)
c  We use delta-Eddington instead
        G1(J)=0.25*(7.-W0(j)*(4.+3*cosbar(j)))
        G2(J)=-0.25*(1.-W0(j)*(4.-3*cosbar(j)))
        G3(J)=0.5*(1.-SQRT(3.)*COSBAR(J)*UBAR0)
        LAMDA(J)=SQRT(G1(J)**2 - G2(J)**2)
        GAMA(J)=(G1(J)-LAMDA(J))/G2(J)
      END DO 

      DO J=1,NAYER
        G4=1.-G3(J)
        DENOM=LAMDA(J)**2 - 1./UBAR0**2
C  NOTE THAT THE ALGORITHM DONOT ACCEPT UBAR0 .eq. 0
C  THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C  THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN
C  THE SCATTERING GOES TO ZERO
C  PREVENT THIS WITH AN IF STATEMENT
        IF ( DENOM .EQ. 0.) THEN
          DENOM=1.E-6
        END IF
        DENOM = W0(J)/DENOM
        AM=DENOM*(G4   *(G1(J)+1./UBAR0) +G2(J)*G3(J) )
        AP=DENOM*(G3(J)*(G1(J)-1./UBAR0) +G2(J)*G4    )
C  CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
C  AT THE TOP OF THE LAYER, THAT IS LOWER OPTICAL DEPTH TAU(J)
        CPM1(J)=AP*TAU(J)
        CMM1(J)=AM*TAU(J)
C  CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C  BOTTOM OF THE LAYER. THAT IS AT HIGHER OPTICAL DEPTH TAU(J+1)
        CP(J)=AP*TAU(J+1)
        CM(J)=AM*TAU(J+1)
      END DO

C  NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
C  FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
C  WARNING IF DTAU(J) IS GREATER THAN ABOUT 35
C  WE CLIPP IT TO AVOID OVERFLOW.
C  EXP (TAU) - EXP(-TAU) WILL BE NONSENSE THIS IS
C  CORRECTED IN THE DSOLVER ROUTINE. ??FLAG?
      DO J=1,NAYER
c        EXPTRM(J) = MIN(35.,LAMDA(J)*DTAU(J))
        EP=EXP(MIN(35.0_8,LAMDA(J)*DTAU(J)))
        EM=1.0/EP
        AM = GAMA(J)*EM
        E1(J)=EP+AM
        E2(J)=EP-AM
        AP = GAMA(J)*EP
        E3(J)=AP+EM
        E4(J)=AP-EM
      END DO

      CALL DSOLVER(NAYER,GAMA,CP,CM,CPM1,CMM1
     &            ,E1,E2,E3,E4,0.0_8,BSURF,RSF,XK1,XK2)

C  EVALUATE THE NAYER FLUXES THROUGH THE NAYER LAYERS
C  USE THE TOP (TAU=0) OPTICAL DEPTH EXPRESSIONS TO EVALUATE FP AND FM
C  AT THE THE TOP OF EACH LAYER,J = LEVEL J
      DO J=1,NAYER
        FP(J)= XK1(J) + GAMA(J)*XK2(J) + CPM1(J)
        FM(J)= GAMA(J)*XK1(J) + XK2(J) + CMM1(J)
      END DO

C  USE EXPRESSION FOR BOTTOM FLUX TO GET THE FP AND FM AT NLEV
c     Optimization: no need this step since result of last
c     loop at about EP above give this
c      EP=EXP(EXPTRM(NAYER))
c      EM=1.0/EP
      FP(NLEV)=XK1(NAYER)*EP + XK2(NAYER)*AM + CP(NAYER)
      FM(NLEV)=XK1(NAYER)*AP + XK2(NAYER)*EM + CM(NAYER)

C  ADD THE DIRECT FLUX TERM TO THE DOWNWELLING RADIATION, LIOU 182
      DO J=1,NLEV
        FM(J)=FM(J)+UBAR0*TAU(J)
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DSOLVER(NL,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,
     *                   BSURF,RSF,XK1,XK2)

C DOUBLE PRECISION VERSION OF SOLVER

cc      PARAMETER (NMAX=201)
cc AS+JBM 03/2010
      IMPLICIT REAL*8  (A-H,O-Z)
      DIMENSION GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL),XK1(NL),
     *          XK2(NL),E1(NL),E2(NL),E3(NL),E4(NL)
cc AS+JBM 03/2010      
cc      DIMENSION AF(NMAX),BF(NMAX),CF(NMAX),DF(NMAX),XK(NMAX)
      DIMENSION AF(2*NL),BF(2*NL),CF(2*NL),DF(2*NL),XK(2*NL)

C*********************************************************
C* THIS SUBROUTINE SOLVES FOR THE COEFFICIENTS OF THE    *
C* TWO STREAM SOLUTION FOR GENERAL BOUNDARY CONDITIONS   *
C* NO ASSUMPTION OF THE DEPENDENCE ON OPTICAL DEPTH OF   *
C* C-PLUS OR C-MINUS HAS BEEN MADE.                      *
C* NL     = NUMBER OF LAYERS IN THE MODEL                *
C* CP     = C-PLUS EVALUATED AT TAO=0 (TOP)              *
C* CM     = C-MINUS EVALUATED AT TAO=0 (TOP)             *
C* CPM1   = C-PLUS  EVALUATED AT TAOSTAR (BOTTOM)        *
C* CMM1   = C-MINUS EVALUATED AT TAOSTAR (BOTTOM)        *
C* EP     = EXP(LAMDA*DTAU)                              *
C* EM     = 1/EP                                         *
C* E1     = EP + GAMA *EM                                *
C* E2     = EP - GAMA *EM                                *
C* E3     = GAMA*EP + EM                                 *
C* E4     = GAMA*EP - EM                                 *
C* BTOP   = THE DIFFUSE RADIATION INTO THE MODEL AT TOP  *
C* BSURF  = THE DIFFUSE RADIATION INTO THE MODEL AT      *
C*          THE BOTTOM: INCLUDES EMMISION AND REFLECTION *
C*          OF THE UNATTENUATED PORTION OF THE DIRECT    *
C*          BEAM. BSTAR+RSF*FO*EXP(-TAOSTAR/U0)          *
C* RSF    = REFLECTIVITY OF THE SURFACE                  *
C* XK1    = COEFFICIENT OF THE POSITIVE EXP TERM         *
C* XK2    = COEFFICIENT OF THE NEGATIVE EXP TERM         *
C*********************************************************
C THIS ROUTINE CALLS ROUTINE DTRIDGL TO SOLVE TRIDIAGONAL
C SYSTEMS
C======================================================================C

      L=2*NL
 
C     ************MIXED COEFFICENTS**********
C     THIS VERSION AVOIDS SINGULARITIES ASSOC.
C     WITH W0=0 BY SOLVING FOR XK1+XK2, AND XK1-XK2.

      AF(1) = 0.0
      BF(1) = GAMA(1)+1.
      CF(1) = GAMA(1)-1.
      DF(1) = BTOP-CMM1(1)
      N     = 0
      LM2   = L-2

C     EVEN TERMS
 
      DO I=2,LM2,2
        N     = N+1
        AF(I) = (E1(N)+E3(N))*(GAMA(N+1)-1.)       
        BF(I) = (E2(N)+E4(N))*(GAMA(N+1)-1.)
        CF(I) = 2.0*(1.-GAMA(N+1)**2)
        DF(I) = (GAMA(N+1)-1.) * (CPM1(N+1) - CP(N)) +
     *            (1.-GAMA(N+1))* (CM(N)-CMM1(N+1))
      END DO
 
      N   = 0
      LM1 = L-1
      DO I=3,LM1,2
        N     = N+1
        AF(I) = 2.0*(1.-GAMA(N)**2)
        BF(I) = (E1(N)-E3(N))*(1.+GAMA(N+1))
        CF(I) = (E1(N)+E3(N))*(GAMA(N+1)-1.)
        DF(I) = E3(N)*(CPM1(N+1) - CP(N)) + E1(N)*(CM(N) - CMM1(N+1))
      END DO
 
      AF(L) = E1(NL)-RSF*E3(NL)
      BF(L) = E2(NL)-RSF*E4(NL)
      CF(L) = 0.0
      DF(L) = BSURF-CP(NL)+RSF*CM(NL)
 
      CALL DTRIDGL(L,AF,BF,CF,DF,XK)
 
C     ***UNMIX THE COEFFICIENTS****

      DO 28 N=1,NL
        XK1(N) = XK(2*N-1)+XK(2*N)
        XK2(N) = XK(2*N-1)-XK(2*N)

C       NOW TEST TO SEE IF XK2 IS REALLY ZERO TO THE LIMIT OF THE
C       MACHINE ACCURACY  = 1 .E -30
C       XK2 IS THE COEFFICIENT OF THE GROWING EXPONENTIAL AND MUST
C       BE TREATED CAREFULLY

        IF(XK2(N) .EQ. 0.0) GO TO 28
        IF (ABS (XK2(N)/XK(2*N-1)) .LT. 1.E-30) XK2(N)=0.0

   28 CONTINUE
 
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DTRIDGL(L,AF,BF,CF,DF,XK)

C     DOUBLE PRECISION VERSION OF TRIDGL

cc AS+JBM 03/2010 : OBSOLETE MAINTENANT      
cc      PARAMETER (NMAX=201)
      IMPLICIT REAL*8  (A-H,O-Z)
      DIMENSION AF(L),BF(L),CF(L),DF(L),XK(L)
cc AS+JBM 03/2010 : OBSOLETE MAINTENANT
cc      DIMENSION AS(NMAX),DS(NMAX)
      DIMENSION AS(L),DS(L)

C*    THIS SUBROUTINE SOLVES A SYSTEM OF TRIDIAGIONAL MATRIX
C*    EQUATIONS. THE FORM OF THE EQUATIONS ARE:
C*    A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = D(I)
C*    WHERE I=1,L  LESS THAN 103.
C* ..............REVIEWED -CP........

C======================================================================C

      AS(L) = AF(L)/BF(L)
      DS(L) = DF(L)/BF(L)

      DO I=2,L
        X         = 1./(BF(L+1-I) - CF(L+1-I)*AS(L+2-I))
        AS(L+1-I) = AF(L+1-I)*X
        DS(L+1-I) = (DF(L+1-I)-CF(L+1-I)*DS(L+2-I))*X
      END DO
 
      XK(1)=DS(1)
      DO I=2,L
        XK(I) = DS(I)-AS(I)*XK(I-1)
      END DO

      RETURN
      END
      
