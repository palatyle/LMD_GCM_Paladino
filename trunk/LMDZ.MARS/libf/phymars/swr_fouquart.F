      SUBROUTINE SWR_FOUQUART ( KDLON, KFLEV, KNU
     S     ,  aerosol,QVISsQREF3d,omegaVIS3d,gVIS3d
     &     ,  albedo,PDSIG,PPSOL,PRMU,PSEC
     S     ,  PFD,PFU )

      use dimradmars_mod, only: sunfr, ndlo2, nsun, 
     &                          ndlon, nflev, naerkind
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
c      BASED ON MORCRETTE EARTH MODEL
C     (See radiation's part of the ecmwf research department
C     documentation, and Fouquart and BonneL (1980)
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
 
      INTEGER jk,ja,jl,jae, jkl,jklp1,jkm1,jaj
      REAL ZTRAY, ZRATIO,ZGAR, ZFF
      real zfacoa,zcorae
      real ZMUE, zgap,zbmu0, zww,zto,zden,zmu1,zbmu1,zden1,zre11

      REAL ZC1I(NDLON,NFLEV+1), ZGG(NDLON), ZREF(NDLON)
     S ,  ZRE1(NDLON), ZRE2(NDLON)
     S ,  ZRMUZ(NDLON), ZRNEB(NDLON), ZR21(NDLON)
     S ,  ZR23(NDLON),  ZSS1(NDLON), ZTO1(NDLON), ZTR(NDLON,2,NFLEV+1)
     S ,  ZTR1(NDLON), ZTR2(NDLON), ZW(NDLON)

      REAL ZRAY1(NDLO2,NFLEV+1), ZRAY2(NDLO2,NFLEV+1)
     s   ,  ZREFZ(NDLO2,2,NFLEV+1)
     S   ,  ZRMUE(NDLO2,NFLEV+1)
     S   ,  ZCGAZ(NDLO2,NFLEV),ZPIZAZ(NDLO2,NFLEV),ZTAUAZ(NDLO2,NFLEV)

      REAL  ZRAYL(NDLON)
     S     ,   ZRJ(NDLON,6,NFLEV+1)
     S     ,  ZRK(NDLON,6,NFLEV+1)
     S     ,  ZTRA1(NDLON,NFLEV+1), ZTRA2(NDLON,NFLEV+1)

c     Function
c     --------
      real CVMGT

C    --------------------------------
C     OPTICAL PARAMETERS FOR AEROSOLS 
C     -------------------------------
C     
      DO  JK = 1 , nlaylte+1
         DO  JA = 1 , 6
            DO JL = 1 , KDLON
               ZRJ(JL,JA,JK) = 0.
               ZRK(JL,JA,JK) = 0.
            END DO
         END DO
      END DO

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

     
C     ----------------------------------------------
C     TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
C     ----------------------------------------------
C     
 200  CONTINUE
      
      DO JL = 1 , KDLON
         ZR23(JL) = 0.
         ZC1I(JL,nlaylte+1) = 0.
      END DO
      
      DO JK = 1 , nlaylte
         JKL = nlaylte+1 - JK
         JKLP1 = JKL + 1
         DO JL = 1 , KDLON
            ZFACOA = 1.-ZPIZAZ(JL,JKL)*ZCGAZ(JL,JKL)*ZCGAZ(JL,JKL)
            ZCORAE = ZFACOA * ZTAUAZ(JL,JKL) * PSEC(JL)
            ZR21(JL) = EXP(-ZCORAE   )
            ZSS1(JL) =  1.0-ZR21(JL)
            ZC1I(JL,JKL) = 1.0-(1.0-ZSS1(JL))*(1.0-ZC1I(JL,JKLP1))
         END DO
      END DO

C     -----------------------------------------------
C     REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
C     -----------------------------------------------
C     
      DO JL = 1 , KDLON
         ZRAY1(JL,nlaylte+1) = 0.
         ZRAY2(JL,nlaylte+1) = 0.
         ZREFZ(JL,2,1) = albedo(JL,KNU)
         ZREFZ(JL,1,1) = albedo(JL,KNU)
         ZTRA1(JL,nlaylte+1) = 1.
         ZTRA2(JL,nlaylte+1) = 1.
      END DO

      DO JK = 2 , nlaylte+1
         JKM1 = JK-1
         DO 342 JL = 1 , KDLON
            ZRNEB(JL)= 1.e-5   ! used to be "cloudiness" (PCLDSW in Morcrette)

            ZRE1(JL)=0.
            ZTR1(JL)=0.
            ZRE2(JL)=0.
            ZTR2(JL)=0.
      
C           EQUIVALENT ZENITH ANGLE
c           ~~~~~~~~~~~~~~~~~~~~~~~
            ZMUE = (1.-ZC1I(JL,JK)) * PSEC(JL)
     S           + ZC1I(JL,JK) * 1.66
            ZRMUE(JL,JK) = 1./ZMUE

C     ------------------------------------------------------------------
C          REFLECT./TRANSMISSIVITY DUE TO AEROSOLS (and rayleigh ?)
C     ------------------------------------------------------------------

            ZGAP = ZCGAZ(JL,JKM1)
            ZBMU0 = 0.5 - 0.75 * ZGAP / ZMUE
            ZWW =ZPIZAZ(JL,JKM1)
            ZTO = ZTAUAZ(JL,JKM1)
            ZDEN = 1. + (1. - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE
     S           + (1-ZWW) * (1. - ZWW +2.*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
            ZRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
            ZTRA1(JL,JKM1) = 1. / ZDEN
C     
            ZMU1 = 0.5
            ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1
            ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1
     S         + (1-ZWW) * (1. - ZWW +2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
            ZRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
            ZTRA2(JL,JKM1) = 1. / ZDEN1

            ZGG(JL) =  ZCGAZ(JL,JKM1)
            ZW(JL) =ZPIZAZ(JL,JKM1)
            ZREF(JL) = ZREFZ(JL,1,JKM1)
            ZRMUZ(JL) = ZRMUE(JL,JK)
            ZTO1(JL) =  ZTAUAZ(JL,JKM1)/ZPIZAZ(JL,JKM1)

 342     CONTINUE

C     
         CALL DEDD ( KDLON
     S        , ZGG,ZREF,ZRMUZ,ZTO1,ZW
     S        , ZRE1,ZRE2,ZTR1,ZTR2     )
C     
         DO JL = 1 , KDLON
C     
            ZREFZ(JL,1,JK) = (1.-ZRNEB(JL)) * (ZRAY1(JL,JKM1)
     S           + ZREFZ(JL,1,JKM1) * ZTRA1(JL,JKM1)
     S           * ZTRA2(JL,JKM1)
     S           / (1.-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))
     S           + ZRNEB(JL) * ZRE2(JL)
C     
            ZTR(JL,1,JKM1) = ZRNEB(JL) * ZTR2(JL) + (ZTRA1(JL,JKM1)
     S           / (1.-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))
     S           * (1.-ZRNEB(JL))
C     
            ZREFZ(JL,2,JK) = (1.-ZRNEB(JL)) * (ZRAY1(JL,JKM1)
     S           + ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)
     S           * ZTRA2(JL,JKM1) )
     S           + ZRNEB(JL) * ZRE1(JL)
C     
            ZTR(JL,2,JKM1) = ZRNEB(JL) * ZTR1(JL)
     S           + ZTRA1(JL,JKM1) * (1.-ZRNEB(JL))
C     
         END DO
      END DO 
C     
C     
C     ------------------------------------------------------------------
C     
C     *         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C     -------------------------------------------------
C     
 350  CONTINUE
C     
      IF (KNU.EQ.1) THEN
         JAJ = 2
         DO 351 JL = 1 , KDLON
            ZRJ(JL,JAJ,nlaylte+1) = 1.
            ZRK(JL,JAJ,nlaylte+1) = ZREFZ(JL, 1,nlaylte+1)
 351     CONTINUE
C     
         DO 353 JK = 1 , nlaylte
            JKL = nlaylte+1 - JK
            JKLP1 = JKL + 1
            DO 352 JL = 1 , KDLON
               ZRE11= ZRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
               ZRJ(JL,JAJ,JKL) = ZRE11
               ZRK(JL,JAJ,JKL) = ZRE11 * ZREFZ(JL,  1,JKL)
 352        CONTINUE
 353     CONTINUE
 354     CONTINUE
C     
      ELSE
C     
         DO 358 JAJ = 1 , 2
            DO 355 JL = 1 , KDLON
               ZRJ(JL,JAJ,nlaylte+1) = 1.
               ZRK(JL,JAJ,nlaylte+1) = ZREFZ(JL,JAJ,nlaylte+1)
 355        CONTINUE
C     
            DO 357 JK = 1 , nlaylte
               JKL = nlaylte+1 - JK
               JKLP1 = JKL + 1
               DO 356 JL = 1 , KDLON
                  ZRE11= ZRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
                  ZRJ(JL,JAJ,JKL) = ZRE11
                  ZRK(JL,JAJ,JKL) = ZRE11 * ZREFZ(JL,JAJ,JKL)
 356           CONTINUE
 357        CONTINUE
 358     CONTINUE
      END IF

C     
C     
C     
C     ------------------------------------------------------------------
C     ---------------
C     DOWNWARD FLUXES
C     ---------------
C    
      JAJ = 2

      do JK = 1 , nlaylte+1
        JKL = nlaylte+1 - JK + 1
        DO  JL = 1 , KDLON
            PFD(JL,JKL) =   ZRJ(JL,JAJ,JKL) * sunfr(KNU)
        end do
      end do
C    
C  -------------
C  UPWARD FLUXES
C  -------------
      DO JK = 1 , nlaylte+1
         DO  JL = 1 , KDLON
c           ZRK = upward flux / incident top flux
            PFU(JL,JK) =    ZRK(JL,JAJ,JK) * sunfr(KNU)
         END DO
      END DO

C     
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DEDD (KDLON,PGG,PREF,PRMUZ,PTO1,PW
     S                ,      PRE1,PRE2,PTR1,PTR2         )
      use dimradmars_mod, only: ndlo2
      implicit none
C
C
C**** *DEDD* - DELTA-EDDINGTON IN A CLOUDY LAYER
C
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
C     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.
C
C**   INTERFACE.
C     ----------
C          *DEDD* IS CALLED BY *SW*.
C
C     SUBROUTINE DEDD (KDLON,PGG,PREF,PRMUZ,PTO1,PW
C    S                ,      PRE1,PRE2,PTR1,PTR2         )
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C PGG    : (NDLON)             ; ASSYMETRY FACTOR
C PREF   : (NDLON)             ; REFLECTIVITY OF THE UNDERLYING LAYER
C PRMUZ  : (NDLON)             ; COSINE OF SOLAR ZENITH ANGLE
C PTO1   : (NDLON)             ; OPTICAL THICKNESS
C PW     : (NDLON)             ; SINGLE SCATTERING ALBEDO
C     ==== OUTPUTS ===
C PRE1   : (NDLON)             ; LAYER REFLECTIVITY ASSUMING NO
C                              ; REFLECTION FROM UNDERLYING LAYER
C PTR1   : (NDLON)             ; LAYER TRANSMISSIVITY ASSUMING NO
C                              ; REFLECTION FROM UNDERLYING LAYER
C PRE2   : (NDLON)             ; LAYER REFLECTIVITY ASSUMING
C                              ; REFLECTION FROM UNDERLYING LAYER
C PTR2   : (NDLON)             ; LAYER TRANSMISSIVITY ASSUMING
C                              ; REFLECTION FROM UNDERLYING LAYER
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.
C
C     EXTERNALS.
C     ----------
C
C          NONE
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "IN CORE MODEL"
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C     ------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
      INTEGER KDLON
C
      REAL PGG(NDLO2),PREF(NDLO2),PRMUZ(NDLO2),PTO1(NDLO2),PW(NDLO2)
      REAL PRE1(NDLO2),PRE2(NDLO2),PTR1(NDLO2),PTR2(NDLO2)

c   local
      integer jl
      real*8 ZFF,ZGP,ZTOP,ZWCP,ZDT,ZX1,ZWM,ZRM2,ZRK,ZX2,ZRP,ZALPHA
      real*8 ZBETA,ZEXMU0,ZEXKP,ZEXKM,ZXP2P,ZXM2P,ZAP2B,ZAM2B
      real*8 ZA11,ZA12,ZA13,ZA22,ZA21,ZA23,ZDENA,ZC1A,ZC2A
      real*8 ZRI0A,ZRI1A,ZRI0B,ZRI1B
      real*8 ZB21,ZB22,ZB23,ZDENB,ZC1B,ZC2B
      real*8 ZRI0C,ZRI1C,ZRI0D,ZRI1D
C
C     ------------------------------------------------------------------
C
C*         1.      DELTA-EDDINGTON CALCULATIONS
C
 100  CONTINUE
C
      DO 131 JL   =   1 , KDLON
C
C*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS
C
 110  CONTINUE
C
      ZFF = PGG(JL)*PGG(JL)
      ZGP = PGG(JL)/(1.+PGG(JL))
      ZTOP = (1.- PW(JL) * ZFF) * PTO1(JL)
      ZWCP = (1-ZFF)* PW(JL) /(1.- PW(JL) * ZFF)
      ZDT = 2./3.
      ZX1 = 1.-ZWCP*ZGP
      ZWM = 1.-ZWCP
      ZRM2 =  PRMUZ(JL) * PRMUZ(JL)
      ZRK = SQRT(3.*ZWM*ZX1)
      ZX2 = 4.*(1.-ZRK*ZRK*ZRM2)
      ZRP = SQRT(3.*ZWM/ZX1)
      ZALPHA = 3.*ZWCP*ZRM2*(1.+ZGP*ZWM)/ZX2
      ZBETA = 3.*ZWCP* PRMUZ(JL) *(1.+3.*ZGP*ZRM2*ZWM)/ZX2
      ZEXMU0 = EXP(-ZTOP/ PRMUZ(JL) )
      ZEXKP = EXP(ZRK*ZTOP)
      ZEXKM = 1./ZEXKP
      ZXP2P = 1.+ZDT*ZRP
      ZXM2P = 1.-ZDT*ZRP
      ZAP2B = ZALPHA+ZDT*ZBETA
      ZAM2B = ZALPHA-ZDT*ZBETA
C
C*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER
C
 120  CONTINUE
C
      ZA11 = ZXP2P
      ZA12 = ZXM2P
      ZA13 = ZAP2B
      ZA22 = ZXP2P*ZEXKP
      ZA21 = ZXM2P*ZEXKM
      ZA23 = ZAM2B*ZEXMU0
      ZDENA = ZA11 * ZA22 - ZA21 * ZA12
      ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA
      ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA
      ZRI0A = ZC1A+ZC2A-ZALPHA
      ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA
      PRE1(JL) = (ZRI0A-ZDT*ZRI1A)/ PRMUZ(JL)
      ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMU0
      ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMU0
      PTR1(JL) = ZEXMU0+(ZRI0B+ZDT*ZRI1B)/ PRMUZ(JL)
C
C*         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER
C
 130  CONTINUE
C
      ZB21 = ZA21- PREF(JL) *ZXP2P*ZEXKM
      ZB22 = ZA22- PREF(JL) *ZXM2P*ZEXKP
      ZB23 = ZA23- PREF(JL) *ZEXMU0*(ZAP2B - PRMUZ(JL) )
      ZDENB = ZA11 * ZB22 - ZB21 * ZA12
      ZC1B = (ZB22*ZA13-ZA12*ZB23)/ZDENB
      ZC2B = (ZA11*ZB23-ZB21*ZA13)/ZDENB
      ZRI0C = ZC1B+ZC2B-ZALPHA
      ZRI1C = ZRP*(ZC1B-ZC2B)-ZBETA
      PRE2(JL) = (ZRI0C-ZDT*ZRI1C) / PRMUZ(JL)
      ZRI0D = ZC1B*ZEXKM + ZC2B*ZEXKP - ZALPHA*ZEXMU0
      ZRI1D = ZRP * (ZC1B*ZEXKM - ZC2B*ZEXKP) - ZBETA*ZEXMU0
      PTR2(JL) = ZEXMU0 + (ZRI0D + ZDT*ZRI1D) / PRMUZ(JL)
C
 131  CONTINUE
      RETURN
      END

