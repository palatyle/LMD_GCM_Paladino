       SUBROUTINE minmax(imax, xi, zmin, zmax )
c
c      P. Le Van

       INTEGER imax
       REAL    xi(imax)
       REAL    zmin,zmax
       INTEGER i

       zmin = xi(1)
       zmax = xi(1)

       DO i = 2, imax
         zmin = MIN( zmin,xi(i) )
         zmax = MAX( zmax,xi(i) )
       ENDDO

       RETURN
       END

