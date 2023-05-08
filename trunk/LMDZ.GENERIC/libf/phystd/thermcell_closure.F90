!
!
!
SUBROUTINE thermcell_closure(ngrid,nlay,ptimestep,rho,zlev,                   &
                             lmax,alim_star,zmin,zmax,wmax,f)
      
      
!===============================================================================
!  Purpose: fermeture, determination de f
!
! Modification 7 septembre 2009
! 1. On enleve alim_star_tot des arguments pour le recalculer et etre ainis
! coherent avec l'integrale au numerateur.
! 2. On ne garde qu'une version des couples wmax,zmax et wmax_sec,zmax_sec
! l'idee etant que le choix se fasse a l'appel de thermcell_closure
! 3. Vectorisation en mettant les boucles en l a l'exterieur avec des if
!===============================================================================
      
      USE thermcell_mod
      
      IMPLICIT NONE
      
      
!===============================================================================
! Declaration
!===============================================================================
      
!     Inputs:
!     -------
      
      INTEGER, INTENT(in) :: ngrid
      INTEGER, INTENT(in) :: nlay
      INTEGER, INTENT(in) :: lmax(ngrid)
      
      REAL, INTENT(in) :: ptimestep
      REAL, INTENT(in) :: rho(ngrid,nlay)
      REAL, INTENT(in) :: zlev(ngrid,nlay)
      REAL, INTENT(in) :: alim_star(ngrid,nlay)
      REAL, INTENT(in) :: zmin(ngrid)
      REAL, INTENT(in) :: zmax(ngrid)
      REAL, INTENT(in) :: wmax(ngrid)
      
!     Outputs:
!     --------
      
      REAL, INTENT(out) :: f(ngrid)
      
!     Local:
!     ------
      
      INTEGER ig, l 
      INTEGER llmax
      
      REAL alim_star_tot(ngrid)
      REAL alim_star2(ngrid)
      REAL plume_height(ngrid)
      
!===============================================================================
! Initialization
!===============================================================================
      
      alim_star2(:) = 0.
      alim_star_tot(:) = 0.
      
      f(:) = 0.
      
      llmax = 1
      
      DO ig=1,ngrid
         plume_height(ig) = zmax(ig) - zmin(ig)
      ENDDO
      
!===============================================================================
! Closure
!===============================================================================
      
!-------------------------------------------------------------------------------
! Indice vertical max atteint par les thermiques sur le domaine
!-------------------------------------------------------------------------------
      
      DO ig=1,ngrid
         IF (lmax(ig) > llmax) THEN
            llmax = lmax(ig)
         ENDIF
      ENDDO
      
!-------------------------------------------------------------------------------
! Calcul des integrales sur la verticale de alim_star et de alim_star^2/(rho dz)
!-------------------------------------------------------------------------------
      
      DO l=1,llmax-1
         DO ig=1,ngrid
               alim_star2(ig) = alim_star2(ig) + alim_star(ig,l)**2           &
               &              / (rho(ig,l) * (zlev(ig,l+1) - zlev(ig,l)))        ! => integration is ok because alim_star = a* dz
               alim_star_tot(ig) = alim_star_tot(ig) + alim_star(ig,l)
         ENDDO
      ENDDO
      
      DO ig=1,ngrid
         IF ((alim_star2(ig) > 0.).and.(plume_height(ig) > 0.)) THEN
            f(ig) = wmax(ig) * alim_star_tot(ig)                              &  ! => normalization is ok
            &     / (plume_height(ig) * r_aspect_thermals * alim_star2(ig))
         ELSE
            f(ig) = 0.
         ENDIF
      ENDDO
      
      
RETURN
END
