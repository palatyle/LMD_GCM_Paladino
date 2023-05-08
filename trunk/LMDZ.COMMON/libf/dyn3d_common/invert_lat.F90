
SUBROUTINE invert_lat(xsize,ysize,vsize,field)

    IMPLICIT NONE
 
! Input variables
    INTEGER, INTENT(IN) :: xsize,ysize,vsize
    REAL, DIMENSION (xsize,ysize,vsize), INTENT(INOUT) :: field
! Local variables
    REAL, DIMENSION (xsize,ysize,vsize)                :: f_aux
    INTEGER :: l,j
 
    DO l=1,vsize
        DO j=1,ysize
            f_aux(:,j,l)=field(:,ysize+1-j,l)
	END DO
    END DO
    
    field=f_aux

    END SUBROUTINE invert_lat
