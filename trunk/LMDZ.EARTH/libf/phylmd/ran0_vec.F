!
! $Header$
!
      subroutine ran0_vec(npoints,idum,ran0)

!     $Id: ran0_vec.F 524 2004-05-19 12:53:04Z lmdzadmin $
!     Platform independent random number generator from
!     Numerical Recipies
!     Mark Webb July 1999
      
      implicit none

      integer j,npoints,idum(npoints),IA,IM,IQ,IR,k(npoints)
      real ran0(npoints),AM

      parameter (IA=16807, IM=2147483647, AM=1.0/IM, IQ=127773, IR=2836)
      
c     do j=1,npoints
c       if (idum(j).eq.0) then
c     	  write(6,*) 'idum=',idum
c  write(6,*) 'ZERO seed not allowed'
c  stop
c       endif
c     enddo

      do j=1,npoints
        k(j)=idum(j)/IQ
        idum(j)=IA*(idum(j)-k(j)*IQ)-IR*k(j)
        if (idum(j).lt.0) idum(j)=idum(j)+IM
        ran0(j)=AM*idum(j)
      enddo

      end
      
