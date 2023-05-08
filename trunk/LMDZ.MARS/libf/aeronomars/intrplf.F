c******************************************************
      SUBROUTINE intrplf(x,y,xd,yd,nd)
c interpolation, give y = f(x) with array xd,yd known, size nd
 
c  Version with CONSTANT values outside limits
c**********************************************************
 
c Variable declaration
c --------------------
c  Arguments :
      real x,y
      real xd(nd),yd(nd)
      integer nd
c  internal
      integer i,j
      real y_undefined
 
c run
c ---
      y_undefined=1.e20
 
      y=0.
      if ((x.le.xd(1)).and.(x.le.xd(nd))) then
        if (xd(1).lt.xd(nd)) y = yd(1) ! yd(1)
        if (xd(1).ge.xd(nd)) y = yd(nd) ! yd(1)
      else if ((x.ge.xd(1)).and.(x.ge.xd(nd))) then
        if (xd(1).lt.xd(nd)) y = yd(nd) ! yd(1)
        if (xd(1).ge.xd(nd)) y = yd(1) ! yd(1)
c        y = yd (nd)
      else
        do i=1,nd-1
         if ( ( (x.ge.xd(i)).and.(x.lt.xd(i+1)) )
     &     .or. ( (x.le.xd(i)).and.(x.gt.xd(i+1)) ) ) then
           y=yd(i)+(x-xd(i))*(yd(i+1)-yd(i))/(xd(i+1)-xd(i))
           goto 99
         end if
        end do
      end if
 
 99   continue
      return
      end                    
