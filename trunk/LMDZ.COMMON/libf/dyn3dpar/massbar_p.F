      SUBROUTINE massbar_p(  masse, massebx, masseby )
     
c
c **********************************************************************
c
c  Calcule les moyennes en x et  y de la masse d'air dans chaque maille.
c **********************************************************************
c    Auteurs : P. Le Van , Fr. Hourdin  .
c   ..........
c
c  ..  masse                 est  un argum. d'entree  pour le s-pg ...
c  ..  massebx,masseby      sont des argum. de sortie pour le s-pg ...
c     
c
      USE parallel_lmdz
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
c
      REAL    masse( ip1jmp1,llm ), massebx( ip1jmp1,llm )  ,
     *      masseby(   ip1jm,llm )
      INTEGER ij,l,ijb,ije
c
c
c   Methode pour calculer massebx et masseby .
c   ----------------------------------------
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
c                       On  a :
c
c    massebx(i,j) = masse(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))   +
c                   masse(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
c     localise  au point  ... U (i,j) ...
c
c    masseby(i,j) = masse(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )  +
c                   masse(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
c     localise  au point  ... V (i,j) ...
c
c
c=======================================================================
      
      
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)  
      DO   100    l = 1 , llm
c
        ijb=ij_begin
        ije=ij_end+iip1
        if (pole_sud) ije=ije-iip1
        
        DO  ij = ijb, ije - 1
         massebx(ij,l) =  masse( ij, l) * alpha1p2( ij  )     + 
     *                   masse(ij+1, l) * alpha3p4(ij+1 )
        ENDDO

c    .... correction pour massebx( iip1,j) .....
c    ...    massebx(iip1,j)= massebx(1,j) ...
c
CDIR$ IVDEP

        

        DO  ij = ijb+iim, ije+iim, iip1
         massebx( ij,l ) = massebx( ij - iim,l )
        ENDDO


      
        ijb=ij_begin-iip1
        ije=ij_end+iip1
        if (pole_nord) ijb=ij_begin
        if (pole_sud) ije=ij_end-iip1

         DO  ij = ijb,ije
         masseby( ij,l ) = masse(  ij   , l ) * alpha2p3(   ij    )  +
     *                     masse(ij+iip1, l ) * alpha1p4( ij+iip1 )
         ENDDO

100   CONTINUE
c$OMP END DO NOWAIT
c
      RETURN
      END
