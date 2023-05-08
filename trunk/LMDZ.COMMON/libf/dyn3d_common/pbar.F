!
! $Header$
!
      SUBROUTINE pbar ( pext, pbarx, pbary, pbarxy )
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c **********************************************************************
c calcul des moyennes en x et en y de (pression au sol*aire variable) ..
c *********************************************************************
c
c          pext               est  un argum. d'entree  pour le s-pg ..
c     pbarx,pbary et pbarxy  sont des argum. de sortie pour le s-pg ..
c
c   Methode:
c   --------
c
c    A chaque point scalaire P (i,j) est affecte 4 coefficients d'aires
c       alpha1(i,j)  calcule  au point ( i+1/4,j-1/4 )
c       alpha2(i,j)  calcule  au point ( i+1/4,j+1/4 )
c       alpha3(i,j)  calcule  au point ( i-1/4,j+1/4 )
c       alpha4(i,j)  calcule  au point ( i-1/4,j-1/4 )
c
c    Avec  alpha1(i,j) = aire(i+1/4,j-1/4)/ aire(i,j)        
c
c    N.B .  Pour plus de details, voir s-pg  ...  iniconst ...
c
c
c
c   alpha4 .         . alpha1    . alpha4
c    (i,j)             (i,j)       (i+1,j)
c
c             P .        U .          . P
c           (i,j)       (i,j)         (i+1,j)
c
c   alpha3 .         . alpha2    .alpha3 
c    (i,j)              (i,j)     (i+1,j)
c
c             V .        Z .          . V
c           (i,j)
c
c   alpha4 .         . alpha1    .alpha4
c   (i,j+1)            (i,j+1)   (i+1,j+1) 
c
c             P .        U .          . P
c          (i,j+1)                    (i+1,j+1)
c
c
c
c
c                       On  a :
c
c    pbarx(i,j) = Pext(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))      +
c                 Pext(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
c     localise  au point  ... U (i,j) ...
c
c    pbary(i,j) = Pext(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )     +
c                 Pext(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
c     localise  au point  ... V (i,j) ...
c
c  pbarxy(i,j)= Pext(i,j) *alpha2(i,j) + Pext(i+1,j) *alpha3(i+1,j) +
c               Pext(i,j+1)*alpha1(i,j+1)+ Pext(i+1,j+1)*alpha4(i+1,j+1)
c     localise  au point  ... Z (i,j) ...
c
c
c
c=======================================================================


#include "dimensions.h"
#include "paramet.h"

#include "comgeom.h"

      REAL pext( ip1jmp1 ),  pbarx ( ip1jmp1 )
      REAL pbary(  ip1jm  ),  pbarxy(  ip1jm  )

      INTEGER   ij



      DO 1 ij = 1, ip1jmp1 - 1
      pbarx( ij ) = pext(ij) * alpha1p2(ij) + pext(ij+1)*alpha3p4(ij+1)
   1  CONTINUE

c    .... correction pour pbarx( iip1,j) .....

c    ...    pbarx(iip1,j)= pbarx(1,j) ...
CDIR$ IVDEP
      DO 2 ij = iip1, ip1jmp1, iip1
      pbarx( ij ) = pbarx( ij - iim )
   2  CONTINUE


      DO 3 ij = 1,ip1jm
      pbary( ij ) = pext(   ij  )   * alpha2p3(   ij   )     +
     *              pext( ij+iip1 ) * alpha1p4( ij+iip1 )
   3  CONTINUE


      DO 5 ij = 1, ip1jm - 1
      pbarxy( ij ) = pext(ij)*alpha2(ij) + pext(ij+1)*alpha3(ij+1) +
     *   pext(ij+iip1)*alpha1(ij+iip1) + pext(ij+iip2)*alpha4(ij+iip2)
   5  CONTINUE


c    ....  correction pour     pbarxy( iip1,j )  ........

CDIR$ IVDEP

      DO 7 ij = iip1, ip1jm, iip1
      pbarxy( ij ) = pbarxy( ij - iim )
   7  CONTINUE


      RETURN
      END
