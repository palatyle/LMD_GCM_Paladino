      subroutine bilinearbig(nX,nY,x_arr,y_arr,f2d_arr,x_in,y_in,f,ind)

!     Necessary for interpolation of continuum data
!     optimized by A. Spiga 01/2013 

      implicit none

      integer nX,nY,i,j,ind,b

      real*8 x_in,y_in,x1,x2,y1,y2
      real*8 f,f11,f12,f21,f22,fA,fB
      real*8 x_arr(nX)
      real*8 y_arr(nY)
      real*8 f2d_arr(nX,nY)
      real*8,save :: x,y
!$OMP THREADPRIVATE(x,y)

      integer strlen
      character*100 label
      label='subroutine bilinear'


      x=x_in
      y=y_in

   !! AS: important to optimize here because the array is quite large
   !! ... and actually calculations only need to be done once
   !! IF ind=-9999 we have not calculated yet
   if ( ind == -9999) then
      !1st check we're within the wavenumber range
      if ((x.lt.x_arr(2)).or.(x.gt.x_arr(nX-2))) then
         ind=-1
      else
        i=1
        x2=x_arr(i)
        do while ( x2 .le. x )
          x1=x2
          i=i+1
          x2=x_arr(i)
          ind=i-1
        end do
      endif
   endif

   !! Either we already saw we are out of wavenumber range
   !! ... and we just have to set f=0 and exit
   if ( ind == -1) then 
      f=0.0D+0
      return
   !! Or we already determined ind -- so we just proceed
   else
      x1=x_arr(ind)
      x2=x_arr(ind+1)
   endif

!     ... and for y within the temperature range
      if ((y.le.y_arr(1)).or.(y.ge.y_arr(nY))) then
         !print*,y_arr(1),y_arr(nY)
         !write(*,*) 'Warning from bilinearbig routine:'
         !write(*,*) 'Outside continuum temperature range!'
         if(y.le.y_arr(1))then
            y=y_arr(1)+0.01
            b=1
            y1=y_arr(b)
            y2=y_arr(b+1)
         endif
         if(y.ge.y_arr(nY))then
            y=y_arr(nY)-0.01
            b=nY-1
            y1=y_arr(b)
            y2=y_arr(b+1)
         endif
      else
        j=1
        y2=y_arr(j)
        do while ( y2 .lt. y )
          y1=y2
          j=j+1
          y2=y_arr(j)
          b=j-1
        end do
      endif
      
      f11=f2d_arr(ind,b)
      f21=f2d_arr(ind+1,b)
      f12=f2d_arr(ind,b+1)
      f22=f2d_arr(ind+1,b+1)

      call bilinear(f,f11,f21,f12,f22,x,x1,x2,y,y1,y2)

      return
    end subroutine bilinearbig
