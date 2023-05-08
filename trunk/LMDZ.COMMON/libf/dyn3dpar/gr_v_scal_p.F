!
! $Header$
!
      SUBROUTINE gr_v_scal_p(nx,x_v,x_scal)
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
      USE parallel_lmdz
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
      INTEGER :: ijb,ije
c-----------------------------------------------------------------------
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
      DO l=1,nx
         DO ij=ijb,ije
            x_scal(ij,l)=
     s      (airev(ij-iip1)*x_v(ij-iip1,l)+airev(ij)*x_v(ij,l))
     s      /(airev(ij-iip1)+airev(ij))
         ENDDO
      ENDDO
      
      if (pole_nord) then
        DO l=1,nx
           DO ij=1,iip1
              x_scal(ij,l)=0.
           ENDDO
        ENDDO
      endif
    
      if (pole_sud) then
        DO l=1,nx
           DO ij=ip1jm+1,ip1jmp1
              x_scal(ij,l)=0.
           ENDDO
        ENDDO
      endif

      RETURN
      END
