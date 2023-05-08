  SUBROUTINE stratocu_if(klon,klev,pctsrf,paprs, pplay,t &
,seuil_inversion,weak_inversion,dthmin)
implicit none

!======================================================================
! J'introduit un peu de diffusion sauf dans les endroits
! ou une forte inversion est presente
! On peut dire qu'il represente la convection peu profonde
!
! Arguments:
! klon-----input-I- nombre de points a traiter
! paprs----input-R- pression a chaque intercouche (en Pa)
! pplay----input-R- pression au milieu de chaque couche (en Pa)
! t--------input-R- temperature (K)
!
! weak_inversion-----logical
!======================================================================
!
! Arguments:
!
    INTEGER, INTENT(IN)                       :: klon,klev
    REAL, DIMENSION(klon, klev+1), INTENT(IN) ::  paprs
    REAL, DIMENSION(klon, klev), INTENT(IN)   ::  pplay
    REAL, DIMENSION(klon, 4), INTENT(IN)   ::  pctsrf
    REAL, DIMENSION(klon, klev), INTENT(IN)   :: t
    
    REAL, DIMENSION(klon), INTENT(OUT)  :: weak_inversion
!
! Quelques constantes et options:
!
    REAL seuil_inversion ! au-dela l'inversion est consideree trop faible
!    PARAMETER (seuil=-0.1)

!
! Variables locales:
!
    INTEGER i, k, invb(klon)
    REAL zl2(klon)
    REAL dthmin(klon), zdthdp

    INCLUDE "indicesol.h"
    INCLUDE "YOMCST.h"

!
! Chercher la zone d'inversion forte
!

    DO i = 1, klon
       invb(i) = klev
       dthmin(i)=0.0
    ENDDO
    DO k = 2, klev/2-1
       DO i = 1, klon
          zdthdp = (t(i,k)-t(i,k+1))/(pplay(i,k)-pplay(i,k+1)) &
               - RD * 0.5*(t(i,k)+t(i,k+1))/RCPD/paprs(i,k+1)
          zdthdp = zdthdp * 100.0
          IF (pplay(i,k).GT.0.8*paprs(i,1) .AND. &
               zdthdp.LT.dthmin(i) ) THEN
             dthmin(i) = zdthdp
             invb(i) = k
          ENDIF
       ENDDO
    ENDDO


!
! Introduire une diffusion:
!
    DO i = 1, klon
       IF ( (pctsrf(i,is_oce) < 0.5) .OR. &
          (invb(i) == klev) .OR. (dthmin(i) > seuil_inversion) ) THEN 
          weak_inversion(i)=1.
       ELSE
          weak_inversion(i)=0.
       ENDIF
    ENDDO

  END SUBROUTINE stratocu_if
