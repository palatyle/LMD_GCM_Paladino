!
! $Id: gr_dyn_fi_p.F 1615 2012-02-10 15:42:26Z emillour $
!
      SUBROUTINE gr_dyn_fi_p(nfield,im,jm,ngrid,pdyn,pfi)
#ifdef CPP_PARA
! Interface with parallel physics,
      USE mod_interface_dyn_phys
      USE dimphy
      USE parallel_lmdz
      IMPLICIT NONE
c=======================================================================
c   passage d'un champ de la grille scalaire a la grille physique
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

      INTEGER im,jm,ngrid,nfield
      REAL pdyn(im,jm,nfield)
      REAL pfi(ngrid,nfield)

      INTEGER i,j,ig,l

c-----------------------------------------------------------------------
c   calcul:
c   -------

c      IF(ngrid.NE.2+(jm-2)*(im-1)) STOP 'probleme de dim'
c   traitement des poles
c   traitement des point normaux
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,nfield    
       DO ig=1,klon
         i=index_i(ig)
         j=index_j(ig)
         pfi(ig,l)=pdyn(i,j,l)
       ENDDO
      ENDDO
c$OMP END DO NOWAIT
#endif
! of #ifdef CPP_PARA
      RETURN
      END
