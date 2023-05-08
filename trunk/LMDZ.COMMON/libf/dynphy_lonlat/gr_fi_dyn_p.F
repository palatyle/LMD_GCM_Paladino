!
! $Id: gr_fi_dyn_p.F 1615 2012-02-10 15:42:26Z emillour $
!
      SUBROUTINE gr_fi_dyn_p(nfield,ngrid,im,jm,pfi,pdyn)
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

      INTEGER i,j,ifield,ig

c-----------------------------------------------------------------------
c   calcul:
c   -------
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO ifield=1,nfield

        do ig=1,klon
          i=index_i(ig)
          j=index_j(ig)
          pdyn(i,j,ifield)=pfi(ig,ifield)
          if (i==1) pdyn(im,j,ifield)=pdyn(i,j,ifield)
	enddo

c   traitement des poles
      if (pole_nord) then
        do i=1,im
	  pdyn(i,1,ifield)=pdyn(1,1,ifield)
	enddo
      endif
       
      if (pole_sud) then
        do i=1,im
	  pdyn(i,jm,ifield)=pdyn(1,jm,ifield)
	enddo
      endif
      
      ENDDO
c$OMP END DO NOWAIT
#endif
! of #ifdef CPP_PARA
      RETURN
      END
