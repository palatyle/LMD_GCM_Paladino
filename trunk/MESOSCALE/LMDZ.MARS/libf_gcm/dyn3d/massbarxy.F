      SUBROUTINE massbarxy(  masse, massebxy )
c
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
#include "comconst.h"
#include "comgeom.h"
c
       REAL  masse( ip1jmp1,llm ), massebxy( ip1jm,llm )
c

      DO l = 1 , llm
c
         DO ij = 1, ip1jm - 1
            massebxy( ij,l ) = masse(    ij  ,l ) * alpha2(   ij    ) +
     +                      masse(   ij+1 ,l ) * alpha3(  ij+1   ) +
     +                      masse( ij+iip1,l ) * alpha1( ij+iip1 ) +
     +                      masse( ij+iip2,l ) * alpha4( ij+iip2 )
         ENDDO

c    ....  correction pour     massebxy( iip1,j )  ........

CDIR$ IVDEP

         DO ij = iip1, ip1jm, iip1
            massebxy( ij,l ) = massebxy( ij - iim,l )
         ENDDO

      ENDDO
c
      RETURN
      END
