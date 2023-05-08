      SUBROUTINE SWR ( KDLON, KFLEV, KNU
     S     ,  aerosol,albedo,PDSIG,PPSOL,PRMU,PSEC
     S     ,  PFD,PFU )

      IMPLICIT NONE
C     
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "callkeys.h"

#include "yomaer.h"
#include "yomlw.h"

#include "fisice.h"

#include "aerice.h"
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
c                         in dimradmars.h , in each layer, for one of
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


      real ray,coefsizew

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
     S              +aerosol(JL,JK,JAE)*QVISsQREF(KNU,JAE)
c              Single scattering albedo
c              ~~~~~~~~~~~~~~~~~~~~~~~~ 
c TEST : to account for the varying w with varying crystal size
               if (activice.and.JAE.eq.naerkind.and.KNU.eq.2) then
                 ray=min( max(rice(JL,JK)*1.e+6, 1.),10.) 
	         coefsizew=(0.0001417*ray**2.-0.00328*ray+0.99667)
     &             /omegavis(KNU,JAE)
	       else
		 coefsizew=1.
               endif
               ZPIZAZ(JL,JK)=ZPIZAZ(JL,JK)+aerosol(JL,JK,JAE)
     S           * QVISsQREF(KNU,JAE)*omegavis(KNU,JAE)*coefsizew
c              Assymetry factor 
c              ~~~~~~~~~~~~~~~~
               ZCGAZ(JL,JK) =  ZCGAZ(JL,JK) +aerosol(JL,JK,JAE)
     S            * QVISsQREF(KNU,JAE)*omegavis(KNU,JAE)*gvis(KNU,JAE)
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
