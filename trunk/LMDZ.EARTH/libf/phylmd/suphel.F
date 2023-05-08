!
! $Header$
!
      SUBROUTINE suphel
C
#include "YOMCST.h"
#include "YOETHF.h"
cIM cf. JLD
       LOGICAL firstcall
       SAVE firstcall
c$OMP THREADPRIVATE(firstcall)
       DATA firstcall /.TRUE./
       
       IF (firstcall) THEN
         PRINT*, 'suphel initialise les constantes du GCM'
         firstcall = .FALSE.
       ELSE
         PRINT*, 'suphel DEJA APPELE '
         RETURN
       ENDIF
C      -----------------------------------------------------------------
C
C*       1.    DEFINE FUNDAMENTAL CONSTANTS.
C              -----------------------------
C
      WRITE(UNIT=6,FMT='(''0*** Constants of the ICM   ***'')')
      RPI=2.*ASIN(1.)
      RCLUM=299792458.
      RHPLA=6.6260755E-34
      RKBOL=1.380658E-23
      RNAVO=6.0221367E+23
      WRITE(UNIT=6,FMT='('' *** Fundamental constants ***'')')
      WRITE(UNIT=6,FMT='(''           PI = '',E13.7,'' -'')')RPI
      WRITE(UNIT=6,FMT='(''            c = '',E13.7,''m s-1'')')
     S RCLUM
      WRITE(UNIT=6,FMT='(''            h = '',E13.7,''J s'')')
     S RHPLA
      WRITE(UNIT=6,FMT='(''            K = '',E13.7,''J K-1'')')
     S RKBOL
      WRITE(UNIT=6,FMT='(''            N = '',E13.7,''mol-1'')')
     S RNAVO
C
C     ----------------------------------------------------------------
C
C*       2.    DEFINE ASTRONOMICAL CONSTANTS.
C              ------------------------------
C
      RDAY=86400.
      REA=149597870000.
      REPSM=0.409093
C
      RSIYEA=365.25*RDAY*2.*RPI/6.283076
      RSIDAY=RDAY/(1.+RDAY/RSIYEA)
      ROMEGA=2.*RPI/RSIDAY
c
c exp1      R_ecc = 0.05
c exp1      R_peri = 102.04
c exp1      R_incl = 22.5
c exp1      print*, 'Parametres orbitaux modifies'
c ref      R_ecc = 0.016724
c ref      R_peri = 102.04
c ref      R_incl = 23.5
c
cIM 161002 : pour avoir les ctes AMIP II
cIM 161002   R_ecc = 0.016724
cIM 161002   R_peri = 102.04
cIM 161002   R_incl = 23.5
cIM on mets R_ecc, R_peri, R_incl dans conf_phys.F90
c     R_ecc = 0.016715
c     R_peri = 102.7
c     R_incl = 23.441
c
      WRITE(UNIT=6,FMT='('' *** Astronomical constants ***'')')
      WRITE(UNIT=6,FMT='(''          day = '',E13.7,'' s'')')RDAY
      WRITE(UNIT=6,FMT='('' half g. axis = '',E13.7,'' m'')')REA
      WRITE(UNIT=6,FMT='('' mean anomaly = '',E13.7,'' -'')')REPSM
      WRITE(UNIT=6,FMT='('' sideral year = '',E13.7,'' s'')')RSIYEA
      WRITE(UNIT=6,FMT='(''  sideral day = '',E13.7,'' s'')')RSIDAY
      WRITE(UNIT=6,FMT='(''        omega = '',E13.7,'' s-1'')')
     S                  ROMEGA
c     write(unit=6,fmt='('' excentricite = '',e13.7,''-'')')R_ecc
c     write(unit=6,fmt='(''     equinoxe = '',e13.7,''-'')')R_peri
c     write(unit=6,fmt='(''  inclinaison = '',e13.7,''-'')')R_incl
C
C     ------------------------------------------------------------------
C
C*       3.    DEFINE GEOIDE.
C              --------------
C
      RG=9.80665
      RA=6371229.
      R1SA=SNGL(1.D0/DBLE(RA))
      WRITE(UNIT=6,FMT='('' ***         Geoide         ***'')')
      WRITE(UNIT=6,FMT='(''      Gravity = '',E13.7,'' m s-2'')')
     S      RG
      WRITE(UNIT=6,FMT='('' Earth radius = '',E13.7,'' m'')')RA
      WRITE(UNIT=6,FMT='('' Inverse E.R. = '',E13.7,'' m'')')R1SA
C
C     -----------------------------------------------------------------
C
C*       4.    DEFINE RADIATION CONSTANTS.
C              ---------------------------
C
c z.x.li      RSIGMA=2. * RPI**5 * RKBOL**4 /(15.* RCLUM**2 * RHPLA**3)
      rsigma = 2.*rpi**5 * (rkbol/rhpla)**3 * rkbol/rclum/rclum/15.
cIM init. dans conf_phys.F90   RI0=1365.
      WRITE(UNIT=6,FMT='('' ***        Radiation       ***'')')
      WRITE(UNIT=6,FMT='('' Stefan-Bol.  = '',E13.7,'' W m-2 K-4''
     S )')  RSIGMA
cIM init. dans conf_phys.F90   WRITE(UNIT=6,FMT='('' Solar const. = '',E13.7,'' W m-2'')')
cIM init. dans conf_phys.F90  S      RI0
C
C     -----------------------------------------------------------------
C
C*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
C              ------------------------------------------
C
      R=RNAVO*RKBOL
      RMD=28.9644
      RMO3=47.9942
      RMV=18.0153
      RD=1000.*R/RMD
      RV=1000.*R/RMV
      RCPD=3.5*RD
      RCVD=RCPD-RD
      RCPV=4. *RV
      RCVV=RCPV-RV
      RKAPPA=RD/RCPD
      RETV=RV/RD-1.
      WRITE(UNIT=6,FMT='('' *** Thermodynamic, gas     ***'')')
      WRITE(UNIT=6,FMT='('' Perfect gas  = '',e13.7)') R
      WRITE(UNIT=6,FMT='('' Dry air mass = '',e13.7)') RMD
      WRITE(UNIT=6,FMT='('' Ozone   mass = '',e13.7)') RMO3
      WRITE(UNIT=6,FMT='('' Vapour  mass = '',e13.7)') RMV
      WRITE(UNIT=6,FMT='('' Dry air cst. = '',e13.7)') RD
      WRITE(UNIT=6,FMT='('' Vapour  cst. = '',e13.7)') RV
      WRITE(UNIT=6,FMT='(''         Cpd  = '',e13.7)') RCPD
      WRITE(UNIT=6,FMT='(''         Cvd  = '',e13.7)') RCVD
      WRITE(UNIT=6,FMT='(''         Cpv  = '',e13.7)') RCPV
      WRITE(UNIT=6,FMT='(''         Cvv  = '',e13.7)') RCVV
      WRITE(UNIT=6,FMT='(''      Rd/Cpd  = '',e13.7)') RKAPPA
      WRITE(UNIT=6,FMT='(''     Rv/Rd-1  = '',e13.7)') RETV
C
C     ----------------------------------------------------------------
C
C*       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
C              ---------------------------------------------
C
      RCW=RCPV
      WRITE(UNIT=6,FMT='('' *** Thermodynamic, liquid  ***'')')
      WRITE(UNIT=6,FMT='(''         Cw   = '',E13.7)') RCW
C
C     ----------------------------------------------------------------
C
C*       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
C              --------------------------------------------
C
      RCS=RCPV
      WRITE(UNIT=6,FMT='('' *** thermodynamic, solid   ***'')')
      WRITE(UNIT=6,FMT='(''         Cs   = '',E13.7)') RCS
C
C     ----------------------------------------------------------------
C
C*       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
C              ----------------------------------------------------
C
      RTT=273.16
      RLVTT=2.5008E+6
      RLSTT=2.8345E+6
      RLMLT=RLSTT-RLVTT
      RATM=100000.
      WRITE(UNIT=6,FMT='('' *** Thermodynamic, trans.  ***'')')
      WRITE(UNIT=6,FMT='('' Fusion point  = '',E13.7)') RTT
      WRITE(UNIT=6,FMT='(''        RLvTt  = '',E13.7)') RLVTT
      WRITE(UNIT=6,FMT='(''        RLsTt  = '',E13.7)') RLSTT
      WRITE(UNIT=6,FMT='(''        RLMlt  = '',E13.7)') RLMLT
      WRITE(UNIT=6,FMT='('' Normal press. = '',E13.7)') RATM
      WRITE(UNIT=6,FMT='('' Latent heat :  '')')
C
C     ----------------------------------------------------------------
C
C*       9.    SATURATED VAPOUR PRESSURE.
C              --------------------------
C
      RESTT=611.14
      RGAMW=(RCW-RCPV)/RV
      RBETW=RLVTT/RV+RGAMW*RTT
      RALPW=LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT)
      RGAMS=(RCS-RCPV)/RV
      RBETS=RLSTT/RV+RGAMS*RTT
      RALPS=LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT)
      RGAMD=RGAMS-RGAMW
      RBETD=RBETS-RBETW
      RALPD=RALPS-RALPW
C
C     ------------------------------------------------------------------
c
c calculer les constantes pour les fonctions thermodynamiques
c
      RVTMP2=RCPV/RCPD-1.
      RHOH2O=RATM/100.
      R2ES=RESTT*RD/RV
      R3LES=17.269
      R3IES=21.875
      R4LES=35.86
      R4IES=7.66
      R5LES=R3LES*(RTT-R4LES)
      R5IES=R3IES*(RTT-R4IES)
C
      RETURN
      END
