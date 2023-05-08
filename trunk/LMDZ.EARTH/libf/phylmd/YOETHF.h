!
! $Header$
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!*    COMMON *YOETHF* DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
!
!     *R__ES*   *CONSTANTS USED FOR COMPUTATION OF SATURATION
!                MIXING RATIO OVER LIQUID WATER(*R_LES*) OR
!                ICE(*R_IES*).
!     *RVTMP2*  *RVTMP2=RCPV/RCPD-1.
!     *RHOH2O*  *DENSITY OF LIQUID WATER.   (RATM/100.)
!
      REAL R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES
      REAL RVTMP2, RHOH2O
      COMMON /YOETHF/R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES,    &
     &               RVTMP2, RHOH2O
!$OMP THREADPRIVATE(/YOETHF/)
