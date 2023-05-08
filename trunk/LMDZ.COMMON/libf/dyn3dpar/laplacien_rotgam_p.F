      SUBROUTINE laplacien_rotgam_p ( klevel, rotin, rotout )
c
c     P. Le Van
c
c   ************************************************************
c   ... calcul de  (rotat x nxgrad)_gam  du rotationnel rotin ..
c   ************************************************************
c     klevel et teta  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
      USE parallel_lmdz
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

c
c    .............   variables  en  arguments    ...........
c
      INTEGER klevel
      REAL rotin( ip1jm,klevel ), rotout( ip1jm,klevel )
c
c   ............     variables   locales     ...............
c
      INTEGER l, ij
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)
c   ........................................................
c
      INTEGER :: ijb,ije
      
c

      CALL   nxgrad_gam_p ( klevel, rotin,   ghx ,   ghy  )
      CALL   rotat_nfil_p ( klevel, ghx  ,   ghy , rotout )
c
      ijb=ij_begin
      ije=ij_end
      if(pole_sud) ije=ij_end-iip1
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1, klevel
        DO ij = ijb, ije
         rotout(ij,l) = rotout(ij,l) * unsairz_gam(ij)
        ENDDO
      ENDDO
c$OMP END DO NOWAIT
      RETURN
      END
