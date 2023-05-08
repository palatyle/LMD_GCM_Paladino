!
! $Id: initrrnpb.F90 1409 2010-07-08 11:46:29Z jghattas $
!
SUBROUTINE  initrrnpb(ftsol,pctsrf,masktr,fshtr,hsoltr,tautr,vdeptr,scavtr)
  USE dimphy
  USE infotrac, ONLY : nbtr
  USE traclmdz_mod, ONLY : id_rn, id_pb
  IMPLICIT NONE
!======================================================================
! Auteur(s): AA + CG (LGGE/CNRS) Date 24-06-94
! Objet: initialisation des constantes des traceurs
! id_rn : identificateur du traceur radon
! id_pb : identificateur du traceur plomb
!======================================================================
! Arguments:
! nbtr.............. nombre de vrais traceurs (sans l'eau)
! ftsol....input-R-  Temperature du sol (Kelvin)
! pctsrf...input-R-  Nature de sol (pourcentage de sol)
! masktr...output-R- Masque reservoir de sol traceur (1 = reservoir)
! fshtr....output-R- Flux surfacique de production dans le reservoir de sol
! hsoltr...output-R- Epaisseur equivalente du reservoir de sol
! tautr....output-R- Constante de decroissance radioactive du traceur
! vdeptr...output-R- Vitesse de depot sec dans la couche Brownienne
! scavtr...output-R- Coefficient de lessivage
!======================================================================
  INCLUDE "indicesol.h"
!======================================================================

  REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: pctsrf
  REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: ftsol
  REAL,DIMENSION(klon,nbtr),INTENT(OUT) :: masktr
  REAL,DIMENSION(klon,nbtr),INTENT(OUT) :: fshtr
  REAL,DIMENSION(nbtr),INTENT(OUT)      :: hsoltr
  REAL,DIMENSION(nbtr),INTENT(OUT)      :: tautr
  REAL,DIMENSION(nbtr),INTENT(OUT)      :: vdeptr
  REAL,DIMENSION(nbtr),INTENT(OUT)      :: scavtr
  INTEGER                               :: i, it
  REAL                                  :: s

  CHARACTER (LEN=20) :: modname='initrrnpb'
  CHARACTER (LEN=80) :: abort_message

!
! Radon it = id_rn
!----------------
  IF (id_rn /= 0) THEN
     it = id_rn
     s = 1.E4             ! Source: atome par m2
     hsoltr(it) = 0.1     ! Hauteur equivalente du reservoir : 
                          ! 1 m * porosite 0.1
     tautr(it) = 4.765E5  ! Decroissance du radon, secondes
     vdeptr(it) = 0.      ! Pas de depot sec pour le radon
     scavtr(it) = 0.      ! Pas de lessivage pour le radon
     
     WRITE(*,*)'-------------- SOURCE DU RADON ------------------------ '
     WRITE(*,*)'it = ',it
     WRITE(*,*)'Source : ', s
     WRITE(*,*)'Hauteur equivalente du reservoir de sol: ',hsoltr(it) 
     WRITE(*,*)'Decroissance (s): ', tautr(it)
     WRITE(*,*)'Vitesse de depot sec: ',vdeptr(it) 
     WRITE(*,*)'Facteur de lessivage: ',scavtr(it)

     DO i = 1,klon
        masktr(i,it) = 0.
        IF ( NINT(pctsrf(i,1)) .EQ. 1 ) masktr(i,it) = 1.
        fshtr(i,it) = s * masktr(i,it)
     END DO

  END IF ! id_rn /= 0

!
! 210Pb it = id_pb
!----------------
  IF (id_pb /= 0) THEN
     it = id_pb
     s = 0.                ! Pas de source 
     hsoltr(it) = 10.      ! Hauteur equivalente du reservoir 
                           ! a partir duquel le depot Brownien a lieu
     tautr(it) = 1.028E9   ! Decroissance du Pb210, secondes
     vdeptr(it) = 1.E-3    ! 1 mm/s pour le 210Pb
     scavtr(it) =  .5      ! Lessivage du Pb210
     DO i = 1,klon
        masktr(i,it) = 1.  ! Le depot sec peut avoir lieu partout
        fshtr(i,it) = s * masktr(i,it)
     END DO
     WRITE(*,*)'-------------- SOURCE DU PLOMB ------------------------ '
     WRITE(*,*)'it = ',it
     WRITE(*,*)'Source : ', s
     WRITE(*,*)'Hauteur equivalente du reservoir : ',hsoltr(it) 
     WRITE(*,*)'Decroissance (s): ', tautr(it)
     WRITE(*,*)'Vitesse de depot sec: ',vdeptr(it) 
     WRITE(*,*)'Facteur de lessivage: ',scavtr(it)
     
  END IF
     
END SUBROUTINE initrrnpb
