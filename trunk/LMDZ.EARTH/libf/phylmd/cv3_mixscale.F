        SUBROUTINE cv3_mixscale(nloc,ncum,na,ment,m)
***************************************************************
*                                                             *
* CV3_MIXSCALE                                                *
*                                                             *
*                                                             *
* written by   : Jean-Yves Grandpeix, 30/05/2003, 16.34.37    *
* modified by :                                               *
***************************************************************
*
      implicit none

#include "cv3param.h"

      integer nloc,ncum,na
      integer i,j,il
      real ment(nloc,na,na),m(nloc,na)
c
      do 100 j=1,nl
        do 101 i=1,nl
          do 102 il=1,ncum
             ment(il,i,j) = m(il,i)*ment(il,i,j)
102      continue
101    continue
100   continue

c
      return
      end
