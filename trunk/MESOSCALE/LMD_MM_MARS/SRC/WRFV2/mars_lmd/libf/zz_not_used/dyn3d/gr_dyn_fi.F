      SUBROUTINE gr_dyn_fi(nfield,im,jm,ngrid,pdyn,pfi)

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

      INTEGER j,ifield,ig
      EXTERNAL SCOPY

c-----------------------------------------------------------------------
c   calcul:
c   -------

      IF(ngrid.NE.2+(jm-2)*(im-1)) STOP 'probleme de dim'
c   traitement des poles
      CALL SCOPY(nfield,pdyn,im*jm,pfi,ngrid)
      CALL SCOPY(nfield,pdyn(1,jm,1),im*jm,pfi(ngrid,1),ngrid)

c   traitement des point normaux
      DO ifield=1,nfield
         DO j=2,jm-1
            ig=2+(j-2)*(im-1)
            CALL SCOPY(im-1,pdyn(1,j,ifield),1,pfi(ig,ifield),1)
         ENDDO
      ENDDO

      RETURN
      END
