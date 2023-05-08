      subroutine diagadv(q,masse,descript)
      implicit none
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
      real q(iip1,jjp1,llm),masse(iip1,jjp1,llm)
      integer i,j,l
      character*10 descript
      real sss,smas,smin,smax
 
      external ismin,ismax
      integer ismin,ismax
 
      sss=0.
      smas=0.
      smin=1.e33
      smax=-1.e33
      do l=1,llm
        do j=1,jjp1
           do i=1,iim
              sss=sss+masse(i,j,l)*q(i,j,l)
              smas=smas+masse(i,j,l)
              smin=min(smin,q(i,j,l))
              smax=max(smax,q(i,j,l))
           enddo
        enddo
      enddo

      write(*,1000) descript,smin,smax,sss/smas,sss
 
      return
1000  format(a10,'   MIN:',e15.4,'   MAX:',e15.4,'   MOY:',e15.4
     s   ,'   TOT:',e15.4)
      end

      FUNCTION ismin(n,sx,incx)
c
      IMPLICIT NONE
c
      integer n,i,incx,ismin,ix
      real sx((n-1)*incx+1),sxmin
c
      ix=1
      ismin=1
      sxmin=sx(1)
      DO i=1,n-1
         ix=ix+incx
         if(sx(ix).lt.sxmin) then
             sxmin=sx(ix)
             ismin=i+1
         endif
      ENDDO
c
      return
      end
C
      function ismax(n,sx,incx)
c
      IMPLICIT NONE
c
      INTEGER n,i,incx,ismax,ix
      real sx((n-1)*incx+1),sxmax
c
      ix=1
      ismax=1
      sxmax=sx(1)
      do 10 i=1,n-1
       ix=ix+incx
       if(sx(ix).gt.sxmax) then
         sxmax=sx(ix)
         ismax=i+1
       endif
10    continue
c
      return
      end

