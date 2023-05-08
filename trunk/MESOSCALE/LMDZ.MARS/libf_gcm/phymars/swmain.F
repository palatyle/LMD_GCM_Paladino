      SUBROUTINE SWMAIN ( KDLON, KFLEV,
     $                PCST, albedo,
     $                PRMU0, PDP, PPLEV, aerosol,PFRACT,
     $                PHEAT, PFLUXD,PFLUXU)

      IMPLICIT NONE
      
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"

#include "yomaer.h"
#include "yomlw.h"
#include "callkeys.h"
c     
c     PURPOSE.
c     --------
c     
c     This routine computes the shortwave (solar wavelength)
c     radiation fluxes in two spectral intervals 
c     and heating rate on the first "nlaylte" layers.
C     
c     Francois Forget (2000), adapted from 
C     Fouquart and Bonnel's ECMWF program
c 
C     IMPLICIT ARGUMENTS :
C     --------------------
C     
C     ==== INPUTS ===
c
c    KDLON :  number of horizontal grid points
c    PCST  :  Solar constant on Mars (W.m-2)
c    albedo   hemispheric surface albedo
c                         albedo (i,1) : mean albedo for solar band#1
c                                        (see below)
c                         albedo (i,2) : mean albedo for solar band#2
c                                        (see below)
c    PRMU0 :         cos of solar zenith angle (=1 when sun at zenith)
c    PDP :           Layer thickness (Pa)
c    PPLEV           pressure (Pa) at boundaries of each layer 
c   aerosol               aerosol extinction optical depth
c                         at reference wavelength "longrefvis" set
c                         in dimradmars.h , in each layer, for one of
c                         the "naerkind" kind of aerosol optical properties.
c    Pfract :        day fraction of the time interval
c                          =1 during the full day ; =0 during the night
C     ==== OUTPUTS ===
c    PHEAT :         Heating rate (K/s)
c    PFLUXD :        SW downward flux at boundaries of each layer (W.m-2)
c    PFLUXU :        SW upward flux at boundaries of each layer (W.m-2)
C     
C     ----------
C     
C-----------------------------------------------------------------------
C     
C     
C-----------------------------------------------------------------------
C     
C     ARGUMENTS
C     ---------
      
      INTEGER KDLON, KFLEV
      REAL ZPSOL(NDLO2), aerosol(NDLO2,KFLEV,naerkind),PRMU0(NDLO2)
      real PCST
      REAL albedo(NDLO2,2)
      REAL PDP(NDLO2,KFLEV)
      REAL PPLEV(NDLO2,KFLEV+1)
      REAL PHEAT(NDLO2,KFLEV)
      REAL PFRACT(NDLO2)
      real PFLUXD(NDLON,NFLEV+1,2)
      real PFLUXU(NDLON,NFLEV+1,2)
      
C     LOCAL ARRAYS
C     ------------
      
      REAL ZDSIG(NDLON,NFLEV), ZFACT(NDLON)
     S     ,  ZFD(NDLON,NFLEV+1)
     S     ,  ZFU(NDLON,NFLEV+1)
     S     ,  ZRMU(NDLON), ZSEC(NDLON)
     S     ,  ZUD(NDLON,3,NFLEV+1), ZUM(NDLON,NFLEV+1)
      REAL ZSIGN(NDLON),  ZSIGO(NDLON)

c following line has been changed, kflev--->nflev (to avoid error message 
c when compiling on NASA Ames Sun)
      REAL ZFDOWN(NDLO2,NFLEV+1),ZFUP(NDLO2,NFLEV+1)

      integer jl, jk, jkp1, jkl
      integer INU 
      real  zdfnet

C     ------------------------------------------------------------------
C     Initializations :
C     ------------------------------------------------------------------

c     Incident Solar flux and corrected angle in the atmosphere 
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO JL = 1 , KDLON
c        Incident Flux at the top of the atmosphere
         ZFACT(JL)= PRMU0(JL) * PCST * PFRACT(JL)

c        Cos of solar zenith angle CORRECTED for high zenith angle
       if (PRMU0(JL).GT.0) then
         ZRMU(JL)=SQRT(1224.* PRMU0(JL) * PRMU0(JL) + 1.) / 35.
       else
         ZRMU(JL)= 1. / 35.
       endif
         ZSEC(JL)=1./ZRMU(JL)
      END DO

c     Calcul of ZDSIG  (thickness of layers in sigma coordinates)
c     ~~~~~~~~~~~~~~~
      DO JL = 1 , KDLON
         ZSIGO(JL) = 1.0
         ZPSOL(JL) = PPLEV(JL,1)
      END DO
      DO JK = 1 , nlaylte
         JKP1 = JK + 1
         JKL = nlaylte+1 - JK
         DO  JL = 1 , KDLON
            ZSIGN(JL) = PPLEV(JL,JKP1) / PPLEV(JL,1)
            ZDSIG(JL,JK) = ZSIGO(JL) - ZSIGN(JL)
            ZSIGO(JL) = ZSIGN(JL)
         END DO
      END DO

C------------------------------------------------------------------
C          LOOP ON SPECTRAL INTERVAL in solar spectrum
C------------------------------------------------------------------
c  2 spectral interval in solar spectrum :
c  - INU=1: between wavelength "long1vis" and "long2vis" set in dimradmars.h
c  - INU=2: between wavelength "long2vis" and "long3vis" set in dimradmars.h

      DO INU = 1,2
         CALL SWR ( KDLON, kflev, INU
     S     ,  aerosol,albedo,ZDSIG,ZPSOL,ZRMU,ZSEC
     S     ,  ZFD,ZFU )

         DO JK = 1 , nlaylte+1
            DO  JL = 1 , KDLON
              PFLUXD(JL,JK,INU)=ZFD(JL,JK)*ZFACT(JL)
              PFLUXU(JL,JK,INU)=ZFU(JL,JK)*ZFACT(JL)
            END DO
         END DO
      END DO

C     ------------------------------------------------------
C     HEATING RATES
C     ------------------------------------------------------

      DO JK = 1 , nlaylte+1
         DO  JL = 1 , KDLON
c           wavelength integrated flux at every level:
            ZFUP(JL,JK)= (PFLUXU(JL,JK,1)+ PFLUXU(JL,JK,2))
            ZFDOWN(JL,JK)= (PFLUXD(JL,JK,1)+ PFLUXD(JL,JK,2))
         END DO
      END DO
     
      DO JK = 1 , nlaylte
         DO  JL = 1 , KDLON
            ZDFNET = ZFUP (JL,JK  ) - ZFDOWN(JL,JK  )
     S              -ZFUP (JL,JK+1) + ZFDOWN(JL,JK+1)
c           Heating rate
            PHEAT(JL,JK) = gcp * ZDFNET / PDP(JL,JK)
c TEST pour diminuer erreur eddington pour diffuseurs simples	    
c           if (activice) then   
c	    if(aerosol(JL,JK,2)/aerosol(JL,JK,1).gt.1) 
c     s 	    PHEAT(JL,JK) =PHEAT(JL,JK)/10.
c           endif 
         END DO
      END DO

      RETURN
      END
