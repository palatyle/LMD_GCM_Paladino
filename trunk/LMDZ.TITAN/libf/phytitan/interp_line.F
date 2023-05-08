      subroutine interp_line(x1,y1,len1,x2,y2,len2)
      implicit none
!-----------------------------------------------------------------------
!
!  Purpose: Do a series of linear interpolations
!  Data sets are organized as vectors (see below)
!  If x2(:), and abscissa at which interpolation is requiered, lies
!  outside of the interval covered by x1(:), instead of doing an
!  extrapolation, y2() is set to the value y1() corresponding to
!  the nearby x1(:) point
!  
c-----------------------------------------------------------------------
!  arguments
!  ---------
!  inputs:
      real x1(len1) ! ordered list of abscissas
      real y1(len1) ! values at x1(:)
      integer len1  ! length of x1(:) and y1(:)
      real x2(len2) !ordered list of abscissas at which interpolation is done
      integer len2  ! length of x2(:) and y2(:)
!  outputs:
      real y2(len2) ! interpolated values
!-----------------------------------------------------------------------

! local variables:
      integer i,j
      

      do i=1,len2
        ! check if x2(i) lies outside of the interval covered by x1()
        if(((x2(i).le.x1(1)).and.(x2(i).le.x1(len1))).or.
     &     ((x2(i).ge.x1(1)).and.(x2(i).ge.x1(len1)))) then
	  ! set y2(i) to y1(1) or y1(len1)
	  if (abs(x2(i)-x1(1)).lt.abs(x2(i)-x1(len1))) then
	    ! x2(i) lies closest to x1(1)
	    y2(i)=y1(1)
	  else
	    ! x2(i) lies closest to x1(len1)
	    y2(i)=y1(len1)
	  endif

	else
        ! find the nearest neigbours and do a linear interpolation
	 do j=1,len1-1
          if(((x2(i).ge.x1(j)).and.(x2(i).le.x1(j+1))).or.
     &       ((x2(i).le.x1(j)).and.(x2(i).ge.x1(j+1)))) then
	    y2(i)=((x2(i)-x1(j))/(x1(j+1)-x1(j)))*y1(j+1)+
     &            ((x2(i)-x1(j+1))/(x1(j)-x1(j+1)))*y1(j)
	  endif
	 enddo
	endif

      enddo

      end
