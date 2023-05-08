!
! $Id: adaptdt.F 1403 2010-07-01 09:02:53Z fairhead $
!
      subroutine adaptdt(nadv,dtbon,n,pbaru,
     c                   masse)

      USE control_mod
      USE comconst_mod, ONLY: dtvr
      IMPLICIT NONE

#include "dimensions.h"
c#include "paramr2.h"
#include "paramet.h"
#include "comdissip.h"
#include "comgeom2.h"

c----------------------------------------------------------
c     Arguments
c----------------------------------------------------------
      INTEGER n,nadv
      REAL dtbon 
      REAL pbaru(iip1,jjp1,llm)
      REAL masse(iip1,jjp1,llm)
c----------------------------------------------------------    
c     Local
c----------------------------------------------------------
      INTEGER i,j,l
      REAL CFLmax,aaa,bbb
      
        CFLmax=0.
        do l=1,llm
         do j=2,jjm
          do i=1,iim
             aaa=pbaru(i,j,l)*dtvr/masse(i,j,l)
             CFLmax=max(CFLmax,aaa)
             bbb=-pbaru(i,j,l)*dtvr/masse(i+1,j,l)
             CFLmax=max(CFLmax,bbb)
          enddo
         enddo
        enddo              
        n=int(CFLmax)+1
c pour reproduire cas VL du code qui appele x,y,z,y,x
c        if (nadv.eq.30) n=n/2   ! Pour Prather
        dtbon=dtvr/n
        
       return
       end







