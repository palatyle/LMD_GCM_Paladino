      subroutine eofdump(ngrid,nlayer,u,v,t,rho,ps)

      implicit none
c
c     Dumps profiles for calculation of variability EOFs 
c     Modified to include rho, FF 09/2004
c     Corrected small bug in sampling rate/count, EM 11/2007
c
#include "dimensions.h"
#include "dimphys.h"
#include "eofdump.h"
c

      integer ngrid
      integer nlayer
      real*4 u(ngrid,nlayer)
      real*4 v(ngrid,nlayer)
      real*4 t(ngrid,nlayer)
      real*4 rho(ngrid,nlayer)
      real*4 ps(ngrid)
      integer count
      data count /0/
      save count
      integer i,j,l, ig

      LOGICAL firstcall
      SAVE firstcall
      DATA firstcall/.true./

c-------------------------------------------------------
c     Initialization at first call:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (firstcall) THEN
        write(*,*) 'CALL ineofdump'
        CALL ineofdump
        firstcall=.false.
      END IF

c-------------------------------------------------------
c     Dumps every ieofs physics timesteps
c
!      write(*,*)'eofdump:count=',count,' ps(1)=',ps(1)
!      if ((ieofs.gt.0).and.(mod(count,ieofs).eq.0)) then
      if (mod(count+1,ieofs).eq.0) then
!        write(*,*)'eofdump: dump --> ps(1)=',ps(1)
        do i=1,iim,eofskip
          do j=1+eofskip/2,jjm+1,eofskip
            ig = 1+ (j-2)*iim +i
            write(uedata) (u(ig,l),l=1,nlayer)
            write(uedata) (v(ig,l),l=1,nlayer)
            write(uedata) (t(ig,l),l=1,nlayer)
            write(uedata) (rho(ig,l),l=1,nlayer)
            write(uedata) ps(ig)
          enddo
        enddo
      endif
      count=count+1
      return
      end
