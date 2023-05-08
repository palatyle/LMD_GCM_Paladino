










!
! $Header$
!
	SUBROUTINE gr_ecrit_fi(nfield,nlon,iim,jjmp1,ecrit,fi)

	IMPLICIT none

c Transformer une variable de la grille d'ecriture a la grille physique
	
	INTEGER nfield,nlon,iim,jjmp1, jjm
      REAL fi(nlon,nfield), ecrit(iim,jjmp1,nfield)
c
      INTEGER i, j, n, ig
c
c	print*,'iim jjm ',iim,jjm

c modif par abd 21 02 01

        jjm = jjmp1 - 1
	do n = 1, nfield
	    fi(1,n) = ecrit(1,1,n)
            fi(nlon,n) = ecrit(1,jjm+1,n)
         DO j = 2, jjm
            ig = 2+(j-2)*iim
            DO i = 1, iim
	     fi(ig-1+i,n) = ecrit(i,j,n)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END

