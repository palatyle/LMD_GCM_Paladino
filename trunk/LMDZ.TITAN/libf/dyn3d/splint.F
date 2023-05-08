     
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
     
c
c     Routine to compute a cubic-spline interpolated value Y given the
c     value of X, the arrays XA, YA and the 2nd derivative array Y2A
c     computed by SUBROUTINE SPLINE. See "Numerical Recipes" for details
c
     
      IMPLICIT REAL (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      DIMENSION xa(n),ya(n),y2a(n)
     
      kl0=1
     
      khi=n
c means of bisection
 1    IF(khi-kl0.gt.1) THEN
     
       k=(khi+kl0)/2
     
       IF(xa(k).gt.x) THEN
     
        khi=k
     
       ELSE
     
        kl0=k
     
       END IF
     
       GO TO 1
     
      END IF
c KL0 and KHI now bracket the X
      h=xa(khi)-xa(kl0)
     
      IF(h.eq.0.0) STOP
      a=(xa(khi)-x)/h
c evaluation of cubic spline polynomial
      b=(x-xa(kl0))/h
     
      y=a*ya(kl0)+b*ya(khi)+((a**3-a)*y2a(kl0)+(b**3-b)*y2a(khi))*(h**2)
     
     ./6.
     
c
     
      RETURN
     
      END
     

