      SUBROUTINE convflu_p( xflu,yflu,nbniv,convfl )
c
c  P. Le Van
c
c
c    *******************************************************************
c  ... calcule la (convergence horiz. * aire locale)du flux ayant pour
c      composantes xflu et yflu ,variables extensives .  ......
c    *******************************************************************
c      xflu , yflu et nbniv sont des arguments d'entree pour le s-pg ..
c      convfl                est  un argument de sortie pour le s-pg .
c
c     njxflu  est le nombre de lignes de latitude de xflu, 
c     ( = jjm ou jjp1 )
c     nbniv   est le nombre de niveaux vert. de  xflu et de yflu .
c
      USE parallel_lmdz
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
      REAL       xflu,yflu,convfl,convpn,convps
      INTEGER    l,ij,nbniv
      DIMENSION  xflu( ip1jmp1,nbniv ),yflu( ip1jm,nbniv ) ,
     *         convfl( ip1jmp1,nbniv )
c
      INTEGER ijb,ije
      EXTERNAL   SSUM
      REAL       SSUM
c
c
#include "comgeom.h"
c
     
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)          
      DO 5 l = 1,nbniv
c
        ijb=ij_begin
        ije=ij_end+iip1
      
        IF (pole_nord) ijb=ij_begin+iip1
        IF (pole_sud)  ije=ij_end-iip1
        
        DO 2  ij = ijb , ije - 1
          convfl(ij+1,l) = xflu(ij,l) - xflu(ij+ 1,l)   +
     *                     yflu(ij +1,l ) - yflu( ij -iim,l )
   2    CONTINUE
c
c

c     ....  correction pour  convfl( 1,j,l)  ......
c     ....   convfl(1,j,l)= convfl(iip1,j,l) ...
c
CDIR$ IVDEP
        DO 3 ij = ijb,ije,iip1
          convfl( ij,l ) = convfl( ij + iim,l )
   3    CONTINUE
c
c     ......  calcul aux poles  .......
c
        IF (pole_nord) THEN
      
          convpn =   SSUM( iim, yflu(     1    ,l ),  1 )
        
          DO ij = 1,iip1
            convfl(ij,l) = convpn * aire(ij) / apoln
          ENDDO
        
        ENDIF
      
        IF (pole_sud) THEN
        
          convps = - SSUM( iim, yflu( ip1jm-iim,l ),  1 )
        
          DO ij = 1,iip1
            convfl(ij+ip1jm,l) = convps * aire(ij+ ip1jm) / apols
          ENDDO
        
        ENDIF
      
   5  CONTINUE
c$OMP END DO NOWAIT   
      RETURN
      END
