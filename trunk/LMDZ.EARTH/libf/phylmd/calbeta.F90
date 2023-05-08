!
! $Header$
!
SUBROUTINE calbeta(dtime,indice,knon,snow,qsol, &
     vbeta,vcal,vdif)

  USE dimphy
  IMPLICIT none
!======================================================================
! Auteur(s): Z.X. Li (LMD/CNRS) (adaptation du GCM au LMD)
! date: 19940414
!======================================================================
!
! Calculer quelques parametres pour appliquer la couche limite
! ------------------------------------------------------------
  INCLUDE "indicesol.h"
  
! Variables d'entrees
!****************************************************************************************
  REAL, INTENT(IN)                   :: dtime
  INTEGER, INTENT(IN)                :: indice
  INTEGER, INTENT(IN)                :: knon
  REAL, DIMENSION(klon), INTENT(IN)  :: snow
  REAL, DIMENSION(klon), INTENT(IN)  :: qsol

  
! Variables de sorties
!****************************************************************************************
  REAL, DIMENSION(klon), INTENT(OUT) :: vbeta
  REAL, DIMENSION(klon), INTENT(OUT) :: vcal
  REAL, DIMENSION(klon), INTENT(OUT) :: vdif

! Variables locales
!****************************************************************************************
  REAL, PARAMETER :: tau_gl=86400.0*5.0 ! temps de relaxation pour la glace de mer
!cc      PARAMETER (tau_gl=86400.0*30.0)
  REAL, PARAMETER :: mx_eau_sol=150.0
  REAL, PARAMETER :: calsol=1.0/(2.5578E+06*0.15)
  REAL, PARAMETER :: calsno=1.0/(2.3867E+06*0.15)
  REAL, PARAMETER :: calice=1.0/(5.1444E+06*0.15)
  
  INTEGER         :: i

!****************************************************************************************  
   
  vbeta(:) = 0.0
  vcal(:) = 0.0
  vdif(:) = 0.0
  
  IF (indice.EQ.is_oce) THEN
     DO i = 1, knon
        vcal(i)   = 0.0
        vbeta(i)  = 1.0
        vdif(i) = 0.0
     ENDDO
  ENDIF
  
  IF (indice.EQ.is_sic) THEN
     DO i = 1, knon
        vcal(i) = calice
        IF (snow(i) .GT. 0.0) vcal(i) = calsno
        vbeta(i)  = 1.0
        vdif(i) = 1.0/tau_gl
!          vdif(i) = calice/tau_gl ! c'etait une erreur
     ENDDO
  ENDIF
  
  IF (indice.EQ.is_ter) THEN
     DO i = 1, knon
        vcal(i) = calsol
        IF (snow(i) .GT. 0.0) vcal(i) = calsno
        vbeta(i)  = MIN(2.0*qsol(i)/mx_eau_sol, 1.0)
        vdif(i) = 0.0
     ENDDO
  ENDIF
  
  IF (indice.EQ.is_lic) THEN
     DO i = 1, knon
        vcal(i) = calice
        IF (snow(i) .GT. 0.0) vcal(i) = calsno
        vbeta(i)  = 1.0
        vdif(i) = 0.0
     ENDDO
  ENDIF
  
END SUBROUTINE calbeta

















