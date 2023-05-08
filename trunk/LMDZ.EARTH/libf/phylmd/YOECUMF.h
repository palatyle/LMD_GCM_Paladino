!
! $Id: YOECUMF.h 1279 2009-12-10 09:02:56Z fairhead $
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez n'utiliser que des ! pour les commentaires
!                 et bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!     ----------------------------------------------------------------
!*    *COMMON* *YOECUMF* - PARAMETERS FOR CUMULUS MASSFLUX SCHEME
!     ----------------------------------------------------------------
!
      COMMON /YOECUMF/                                                  &
     &                 LMFPEN,LMFSCV,LMFMID,LMFDD,LMFDUDV,              &
     &                 ENTRPEN,ENTRSCV,ENTRMID,ENTRDD,CMFCTOP,          &
     &                 CMFCMAX,CMFCMIN,CMFDEPS,RHCDD,CPRCON

      LOGICAL          LMFPEN,LMFSCV,LMFMID,LMFDD,LMFDUDV
      REAL ENTRPEN, ENTRSCV, ENTRMID, ENTRDD
      REAL CMFCTOP, CMFCMAX, CMFCMIN, CMFDEPS, RHCDD, CPRCON
!$OMP THREADPRIVATE(/YOECUMF/)
!
!*if (DOC,declared) <> 'UNKNOWN'
!*    *COMMON* *YOECUMF* - PARAMETERS FOR CUMULUS MASSFLUX SCHEME
!
!     M.TIEDTKE       E. C. M. W. F.      18/1/89
!
!     NAME      TYPE      PURPOSE
!     ----      ----      -------
!
!     LMFPEN    LOGICAL  TRUE IF PENETRATIVE CONVECTION IS SWITCHED ON
!     LMFSCV    LOGICAL  TRUE IF SHALLOW     CONVECTION IS SWITCHED ON
!     LMFMID    LOGICAL  TRUE IF MIDLEVEL    CONVECTION IS SWITCHED ON
!     LMFDD     LOGICAL  TRUE IF CUMULUS DOWNDRAFT      IS SWITCHED ON
!     LMFDUDV   LOGICAL  TRUE IF CUMULUS FRICTION       IS SWITCHED ON
!     ENTRPEN   REAL     ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV   REAL     ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID   REAL     ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD    REAL     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     CMFCTOP   REAL     RELAT. CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANC
!     CMFCMAX   REAL     MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     CMFCMIN   REAL     MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     CMFDEPS   REAL     FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     RHCDD     REAL     RELATIVE SATURATION IN DOWNDRAFTS
!     CPRCON    REAL     COEFFICIENTS FOR DETERMINING CONVERSION
!                        FROM CLOUD WATER TO RAIN
!*ifend
!     ----------------------------------------------------------------
