










!
! $Header$
!
       SUBROUTINE minmax2(imax, jmax, lmax, xi, zmin, zmax )
c
       INTEGER lmax,jmax,imax
       REAL xi(imax*jmax*lmax) 
       REAL zmin,zmax
       INTEGER i
    
       zmin = xi(1)
       zmax = xi(1)

       DO i = 2, imax*jmax*lmax
         zmin = MIN( zmin,xi(i) )
         zmax = MAX( zmax,xi(i) )
       ENDDO

       RETURN
       END
