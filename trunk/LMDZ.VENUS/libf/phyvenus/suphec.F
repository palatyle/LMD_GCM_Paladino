!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/suphec.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
      SUBROUTINE suphec(cpp_)
C
#include "YOMCST.h"
cIM cf. JLD
       REAL cpp_   ! from dynamics
       LOGICAL firstcall
       SAVE firstcall
       DATA firstcall /.TRUE./
       IF (firstcall) THEN
         PRINT*, 'suphec initialise les constantes du GCM'
         firstcall = .FALSE.
       ELSE
         PRINT*, 'suphec DEJA APPELE '
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
c TERRE
c     RDAY=86400.
c     REA=149597870000.
c     REPSM=0.409093
C
c     RSIYEA=365.25*RDAY*2.*RPI/6.283076
c 1/(duree du jour) = 1/(periode rotation) - 1/(periode revolution)
c     RSIDAY=RDAY/(1.+RDAY/RSIYEA)
c     ROMEGA=2.*RPI/RSIDAY

c VENUS
      RSIDAY=20.9961e6   ! 243.01 j
      RSIYEA=19.4141e6   ! 224.7 j
      ROMEGA=2.*RPI/RSIDAY
c 1/(duree du jour) = 1/(periode rotation) + 1/(periode revolution)
      RDAY=RSIDAY/(1.+RSIDAY/RSIYEA) ! 116.748 j
      REA=108.15e9 
      REPSM=0.  ! 0. veut dire qu'on commence au point vernal
c
cIM on mets R_ecc, R_peri, R_incl dans conf_phys.F90

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
c TERRE
c     RG=9.80665
c     RA=6371229.

c VENUS
      RG=8.87
      RA=6051300.

      R1SA=SNGL(1.D0/DBLE(RA))
      WRITE(UNIT=6,FMT='('' ***         Geoide         ***'')')
      WRITE(UNIT=6,FMT='(''       Gravity = '',E13.7,'' m s-2'')')
     S      RG
      WRITE(UNIT=6,FMT='('' Planet radius = '',E13.7,'' m'')')RA
      WRITE(UNIT=6,FMT='(''  Inverse P.R. = '',E13.7,'' m-1'')')R1SA
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
c TERRE
c     RMD=28.9644
      RMV=18.0153

c VENUS
      RMD=43.44

      RD=1000.*R/RMD
      RV=1000.*R/RMV
c TERRE
c     RCPD=3.5*RD
      RCPV=4. *RV
c VENUS
! ADAPTATION GCM POUR CP(T)
! VENUS: Cp(T) = RCPD*(T/T0)^nu (RCPD phys = cpp dyn)
! avec RCPD=1000., T0=460. et nu=0.35
      RCPD=cpp_

      RCVD=RCPD-RD
      RCVV=RCPV-RV
      RKAPPA=RD/RCPD
      RETV=RV/RD-1.
      WRITE(UNIT=6,FMT='('' *** Thermodynamic, gas     ***'')')
      WRITE(UNIT=6,FMT='('' Perfect gas  = '',e13.7)') R
      WRITE(UNIT=6,FMT='('' Dry air mass = '',e13.7)') RMD
      WRITE(UNIT=6,FMT='('' Vapour  mass = '',e13.7)') RMV
      WRITE(UNIT=6,FMT='('' Dry air cst. = '',e13.7)') RD
      WRITE(UNIT=6,FMT='('' Vapour  cst. = '',e13.7)') RV
      WRITE(UNIT=6,FMT='(''        Cpd0  = '',e13.7)') RCPD
      WRITE(UNIT=6,FMT='(''         Cvd  = '',e13.7)') RCVD
      WRITE(UNIT=6,FMT='(''         Cpv  = '',e13.7)') RCPV
      WRITE(UNIT=6,FMT='(''         Cvv  = '',e13.7)') RCVV
      WRITE(UNIT=6,FMT='(''     Rd/Cpd0  = '',e13.7)') RKAPPA
      WRITE(UNIT=6,FMT='(''     Rv/Rd-1  = '',e13.7)') RETV
C
      RETURN
      END
