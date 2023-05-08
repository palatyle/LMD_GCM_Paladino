      subroutine spline(x,y,n,yp1,ypn,y2)
     
c
     
c     Routine to set up the interpolating function for a cubic spline
     
c     interpolation (see "Numerical Recipes" for details).
     
c
      implicit real (a-h,o-z)
      implicit integer (i-n)
     
      parameter(nllm=4096)
     
      dimension x(n),y(n),y2(n),u(nllm)
     
c
c	write(6,*)(x(i),i=1,n)
c	write(6,*)(y(i),i=1,n)
     
      if(yp1.gt.0.99E30) then
c the lower boundary condition is set
       y2(1)=0.
c either to be "natural"
       u(1)=0.
     
      else
c or else to have a specified first
       y2(1)=-0.5
c derivative
       u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
     
      end if
     
      do 11 i=2,n-1
c decomposition loop of the tridiagonal
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
c algorithm. Y2 and U are used
       p=sig*y2(i-1)+2.
c for temporary storage of the decompo-
       y2(i)=(sig-1.)/p
c sed factors
       u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     
     . /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
     
 11   continue
     
      if(ypn.gt.0.99E30) then
c the upper boundary condition is set
       qn=0.
c either to be "natural"
       un=0.
     
      else
c or else to have a specified first
       qn=0.5
c derivative
       un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
     
      end if
     
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
     
      do 12 k=n-1,1,-1
c this is the backsubstitution loop of
       y2(k)=y2(k)*y2(k+1)+u(k)
c the tridiagonal algorithm
 12   continue
     
c
     
      return
     
      end
     
