      SUBROUTINE gr_v_scal(nx,x_v,x_scal)
c%W%    %G%
c=======================================================================
c
c   Author:    Frederic Hourdin      original: 11/11/92
c   -------
c
c   Subject:
c   ------
c
c   Method:
c   --------
c
c   Interface:
c   ----------
c
c      Input:
c      ------
c
c      Output:
c      -------
c
c=======================================================================
      IMPLICIT NONE
c-----------------------------------------------------------------------
c   Declararations:
c   ---------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

c   Arguments:
c   ----------

      INTEGER nx
      REAL x_v(ip1jm,nx),x_scal(ip1jmp1,nx)

c   Local:
c   ------

      INTEGER l,ij

c-----------------------------------------------------------------------

      DO l=1,nx
         DO ij=iip2,ip1jm
            x_scal(ij,l)=
     s      (airev(ij-iip1)*x_v(ij-iip1,l)+airev(ij)*x_v(ij,l))
     s      /(airev(ij-iip1)+airev(ij))
         ENDDO
         DO ij=1,iip1
            x_scal(ij,l)=0.
         ENDDO
         DO ij=ip1jm+1,ip1jmp1
            x_scal(ij,l)=0.
         ENDDO
      ENDDO

      RETURN
      END
