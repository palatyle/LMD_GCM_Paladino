!
! $Id $
!
SUBROUTINE radio_decay(radio,rnpb,dtime,tautr,tr,d_tr) 
!
! Caluclate radioactive decay for all tracers with radio(it)=true
!
  USE dimphy
  USE infotrac, ONLY : nbtr
  USE traclmdz_mod, ONLY : id_rn, id_pb
  IMPLICIT NONE
!-----------------------------------------------------------------------
! Auteur(s): AA + CG (LGGE/CNRS) Date 24-06-94
! Objet: Calcul de la tendance radioactive des traceurs type radioelements
!        Cas particulier pour le couple radon-plomb : Le radon decroit en plomb
!-----------------------------------------------------------------------
!
! Entrees
!
  LOGICAL,DIMENSION(nbtr),INTENT(IN)        :: radio ! .true. = traceur radioactif  
  LOGICAL,INTENT(IN)                        :: rnpb  ! .true. = decroissance RN = source PB
  REAL,INTENT(IN)                           :: dtime ! Pas de temps physique (secondes)
  REAL,DIMENSION(nbtr),INTENT(IN)           :: tautr ! Constante de decroissance radioactive
  REAL,DIMENSION(klon,klev,nbtr),INTENT(IN) :: tr    ! Concentrations traceurs U/kgA
!
! Sortie
!
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT) :: d_tr  ! Tendance de decroissance radioactive
!
! Locales
!
  INTEGER  :: i,k,it


  DO it = 1,nbtr
     d_tr(:,:,it) = 0.
     IF ( radio(it) ) THEN
        IF (tautr(it) .GT. 0.) THEN
           DO k = 1,klev
              DO i = 1,klon
                 d_tr(i,k,it) = - tr(i,k,it) * dtime / tautr(it)
              END DO
           END DO
        END IF
     END IF
  END DO

!-------------------------------------------------------
! Cas particulier radon (id_rn) => plomb (id_pb)
!-------------------------------------------------------
  IF ( rnpb ) THEN
     DO k = 1,klev
        DO i = 1,klon
           d_tr(i,k,id_pb) = d_tr(i,k,id_pb) - d_tr(i,k,id_rn)
        ENDDO
     ENDDO
  ENDIF

END SUBROUTINE radio_decay
