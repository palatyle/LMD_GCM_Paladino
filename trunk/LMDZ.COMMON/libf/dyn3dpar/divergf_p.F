      SUBROUTINE divergf_p(klevel,x,y,div)
c
c     P. Le Van
c
c  *********************************************************************
c  ... calcule la divergence a tous les niveaux d'1 vecteur de compos. 
c     x et y...
c              x et y  etant des composantes covariantes   ...
c  *********************************************************************
      USE parallel_lmdz
      IMPLICIT NONE
c
c      x  et  y  sont des arguments  d'entree pour le s-prog
c        div      est  un argument  de sortie pour le s-prog
c
c
c   ---------------------------------------------------------------------
c
c    ATTENTION : pendant ce s-pg , ne pas toucher au COMMON/scratch/  .
c
c   ---------------------------------------------------------------------
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
c
c    ..........          variables en arguments    ...................
c
      INTEGER klevel
      REAL x( ip1jmp1,klevel ),y( ip1jm,klevel ),div( ip1jmp1,klevel )
      INTEGER   l,ij
c
c    ...............     variables  locales   .........................

      REAL aiy1( iip1 ) , aiy2( iip1 )
      REAL sumypn,sumyps
c    ...................................................................
c
      EXTERNAL  SSUM
      REAL      SSUM
      INTEGER :: ijb,ije,jjb,jje
c
c
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ij_begin+iip1
      if(pole_sud)  ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO 10 l = 1,klevel
c
        DO  ij = ijb, ije - 1
         div( ij + 1, l )     =  
     *   cvusurcu( ij+1 ) * x( ij+1,l ) - cvusurcu( ij ) * x( ij , l) +
     *   cuvsurcv(ij-iim) * y(ij-iim,l) - cuvsurcv(ij+1) * y(ij+1,l) 
        ENDDO

c
c     ....  correction pour  div( 1,j,l)  ......
c     ....   div(1,j,l)= div(iip1,j,l) ....
c
CDIR$ IVDEP
        DO  ij = ijb,ije,iip1
         div( ij,l ) = div( ij + iim,l )
        ENDDO
c
c     ....  calcul  aux poles  .....
c
        if (pole_nord) then
        
          DO  ij  = 1,iim
           aiy1(ij) =    cuvsurcv(    ij       ) * y(     ij     , l )
          ENDDO
          sumypn = SSUM ( iim,aiy1,1 ) / apoln

c
          DO  ij = 1,iip1
           div(     ij    , l ) = - sumypn
          ENDDO
          
        endif
        
        if (pole_sud) then
        
          DO  ij  = 1,iim
           aiy2(ij) =    cuvsurcv( ij+ ip1jmi1 ) * y( ij+ ip1jmi1, l )
          ENDDO
          sumyps = SSUM ( iim,aiy2,1 ) / apols
c
          DO  ij = 1,iip1
           div( ij + ip1jm, l ) =   sumyps
          ENDDO
          
        endif
        
  10    CONTINUE
c$OMP END DO NOWAIT

c
        jjb=jj_begin
        jje=jj_end
        if (pole_sud) jje=jj_end-1
        
        CALL filtreg_p( div,jjb,jje, jjp1, klevel, 2, 2, .TRUE., 1 )
      
c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l = 1, klevel
           DO ij = ijb,ije
            div(ij,l) = div(ij,l) * unsaire(ij) 
          ENDDO
        ENDDO
c$OMP END DO NOWAIT
c
       RETURN
       END
