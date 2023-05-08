      SUBROUTINE nxgrad_gam_p( klevel, rot, x, y )
c
c  P. Le Van
c
c   ********************************************************************
c      calcul du gradient tourne de pi/2 du rotationnel du vect.v
c   ********************************************************************
c       rot          est un argument  d'entree pour le s-prog
c       x  et y    sont des arguments de sortie pour le s-prog
c
      USE parallel_lmdz
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
      INTEGER klevel
      REAL rot( ip1jm,klevel ),x( ip1jmp1,klevel ),y(ip1jm,klevel )
      INTEGER   l,ij
      integer ismin,ismax
      external ismin,ismax
      INTEGER :: ijb,ije
c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 10 l = 1,klevel
c
      ijb=ij_begin
      ije=ij_end
      if(pole_sud) ije=ij_end-iip1
      
      DO 1  ij = ijb+1, ije
      y( ij,l ) = (rot( ij,l ) - rot( ij-1,l )) * cvscuvgam( ij )
   1  CONTINUE
c
c    ..... correction pour  y ( 1,j,l )  ......
c
c    ....    y(1,j,l)= y(iip1,j,l) ....
CDIR$ IVDEP
      DO 2  ij = ijb, ije, iip1
      y( ij,l ) = y( ij +iim,l )
   2  CONTINUE
c
      ijb=ij_begin
      ije=ij_end+iip1
      if(pole_nord) ijb=ij_begin+iip1
      if(pole_sud) ije=ij_end-iip1
      
      DO 4  ij = ijb,ije
      x( ij,l ) = (rot( ij,l ) - rot( ij -iip1,l )) * cuscvugam( ij )
   4  CONTINUE
    
      if (pole_nord) then
        DO  ij = 1,iip1
         x(    ij    ,l ) = 0.
        ENDDO
      endif

      if (pole_sud) then
        DO  ij = 1,iip1
         x( ij +ip1jm,l ) = 0.
        ENDDO
      endif
c
  10  CONTINUE
c$OMP END DO NOWAIT
      RETURN
      END
