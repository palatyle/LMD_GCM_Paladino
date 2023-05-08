      SUBROUTINE massbarxy_p(  masse, massebxy )
      USE parallel_lmdz
      implicit none
c **********************************************************************
c
c  Calcule les moyennes en x et  y de la masse d'air dans chaque maille.
c **********************************************************************
c    Auteurs : P. Le Van , Fr. Hourdin  .
c   ..........
c
c  ..  masse          est  un  argum. d'entree  pour le s-pg ...
c  ..  massebxy       est  un  argum. de sortie pour le s-pg ...
c     
c
c     IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
c
       REAL  masse( ip1jmp1,llm ), massebxy( ip1jm,llm )
c
      INTEGER ij,l,ijb,ije

      
      ijb=ij_begin-iip1
      ije=ij_end
      
      if (pole_nord) ijb=ijb+iip1
      if (pole_sud)  ije=ije-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
      DO   100    l = 1 , llm
c
      DO 5 ij = ijb, ije - 1
      massebxy( ij,l ) = masse(    ij  ,l ) * alpha2(   ij    )   +
     +                   masse(   ij+1 ,l ) * alpha3(  ij+1   )   +
     +                   masse( ij+iip1,l ) * alpha1( ij+iip1 )   +
     +                   masse( ij+iip2,l ) * alpha4( ij+iip2 )
   5  CONTINUE

c    ....  correction pour     massebxy( iip1,j )  ........

CDIR$ IVDEP

      DO 7 ij = ijb+iip1-1, ije+iip1-1, iip1
      massebxy( ij,l ) = massebxy( ij - iim,l )
   7  CONTINUE

100   CONTINUE
c$OMP END DO NOWAIT
c
      RETURN
      END
