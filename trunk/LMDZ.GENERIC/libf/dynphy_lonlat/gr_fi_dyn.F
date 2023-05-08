      SUBROUTINE gr_fi_dyn(nfield,ngrid,im,jm,pfi,pdyn)
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
      EXTERNAL SCOPY

c-----------------------------------------------------------------------
c   calcul:
c   -------

      DO ifield=1,nfield
c   traitement des poles
         DO i=1,im
            pdyn(i,1,ifield)=pfi(1,ifield)
            pdyn(i,jm,ifield)=pfi(ngrid,ifield)
         ENDDO

c   traitement des point normaux
         DO j=2,jm-1
	    ig=2+(j-2)*(im-1)
            CALL SCOPY(im-1,pfi(ig,ifield),1,pdyn(1,j,ifield),1)
	    pdyn(im,j,ifield)=pdyn(1,j,ifield)
         ENDDO
      ENDDO

      RETURN
      END
