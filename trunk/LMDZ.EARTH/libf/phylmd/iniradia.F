      SUBROUTINE iniradia (klon,klev,pres)
  
      IMPLICIT none
c======================================================================
c
c Auteur(s) MP Lefebvre        date: 20080827
c
c Objet: initialise le rayonnement RRTM           
c======================================================================
c  Arguments:
c
c klon----input-I-nombre de points horizontaux
c klev----input-I-nombre de couches verticales
c pres----input-R-pression pour chaque inter-couche (en Pa)
c======================================================================
c
      INTEGER klon
      INTEGER klev
      REAL pres(klev+1)

!         CALL suphel     ! initialiser constantes et parametres phys.
!     print*,'Physiq: apres suphel '
!        CALL SUINIT(klon,klev)
!     print*,'iniradia: apres suinit '
! calcul des niveaux de pression de reference au bord des couches pour
! l'intialisation des aerosols. Momentannement, on passe un point de
! grille du profil de pression.
!        CALL SURAYOLMD(pres(klev+1))  ! initialiser le rayonnement RRTM
!     print*,'iniradia: apres surayolmd '

      RETURN
      END
