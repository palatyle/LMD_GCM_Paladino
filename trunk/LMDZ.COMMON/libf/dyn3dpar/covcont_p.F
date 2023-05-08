      SUBROUTINE covcont_p (klevel,ucov, vcov, ucont, vcont )
      USE parallel_lmdz, ONLY: ij_begin,ij_end,OMP_CHUNK,
     &                         pole_nord, pole_sud
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c  *********************************************************************
c    calcul des compos. contravariantes a partir des comp.covariantes
c  ********************************************************************
c
c=======================================================================

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

      INTEGER klevel
      REAL ucov( ip1jmp1,klevel ),  vcov( ip1jm,klevel )
      REAL ucont( ip1jmp1,klevel ), vcont( ip1jm,klevel )
      INTEGER   l,ij
      INTEGER ijb_u,ijb_v,ije_u,ije_v

      
      ijb_u=ij_begin-iip1
      ijb_v=ij_begin-iip1
      ije_u=ij_end+iip1
      ije_v=ij_end+iip1
      
      if (pole_nord) then 
        ijb_u=ij_begin+iip1
        ijb_v=ij_begin
      endif
      
      if (pole_sud) then
        ije_u=ij_end-iip1
        ije_v=ij_end-iip1
      endif

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
      DO 10 l = 1,klevel

      DO 2  ij = ijb_u,ije_u
      ucont( ij,l ) = ucov( ij,l ) * unscu2( ij )
   2  CONTINUE

      DO 4 ij = ijb_v,ije_v
      vcont( ij,l ) = vcov( ij,l ) * unscv2( ij )
   4  CONTINUE

  10  CONTINUE
c$OMP END DO NOWAIT
      RETURN
      END
