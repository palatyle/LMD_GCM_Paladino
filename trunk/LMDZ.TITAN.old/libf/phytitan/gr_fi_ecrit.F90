SUBROUTINE gr_fi_ecrit(nfield,nlon,iim,jjmp1,fi,ecrit)
  IMPLICIT none
  !
  ! Tranformer une variable de la grille physique a
  ! la grille d'ecriture
  !
  ! WARNING: This only works on the full global grid
  !          (ie for GCM in serial mode)
  INTEGER nfield,nlon,iim,jjmp1, jjm
  REAL fi(nlon,nfield), ecrit(iim*jjmp1,nfield)
  !
  INTEGER i, n, ig
  !
  jjm = jjmp1 - 1
  DO n = 1, nfield
     DO i=1,iim
        ecrit(i,n) = fi(1,n)
        ecrit(i+jjm*iim,n) = fi(nlon,n)
     ENDDO
     DO ig = 1, nlon - 2
        ecrit(iim+ig,n) = fi(1+ig,n)
     ENDDO
  ENDDO
END SUBROUTINE gr_fi_ecrit
