!
! $Header$
!
      subroutine zilch(x,m)
c
c Zero the real array x dimensioned m.
c
      implicit none
c
      integer m,i
      real x(m)
      do 1 i=1,m
      x(i)= 0.0  
    1 continue
      return
      end
