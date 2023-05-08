      SUBROUTINE laplacien_rot_p ( klevel, rotin, rotout,ghx,ghy )
c
c    P. Le Van
c
c   ************************************************************
c    ...  calcul de  ( rotat x nxgrad )  du rotationnel rotin  .
c   ************************************************************
c
c     klevel et rotin  sont des arguments  d'entree pour le s-prog
c      rotout           est  un argument  de sortie pour le s-prog
c
      USE parallel_lmdz
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

c 
c   ..........    variables  en  arguments     .............
c
      INTEGER klevel
      REAL rotin( ip1jm,klevel ), rotout( ip1jm,klevel )
c
c   ..........    variables   locales       ................
c
      REAL ghy(ip1jm,klevel), ghx(ip1jmp1,klevel)
c   ........................................................
c
c
      INTEGER :: ijb,ije,jjb,jje
      
      jjb=jj_begin-1
      jje=jj_end+1
      
      if (pole_nord) jjb=jj_begin
      if (pole_sud) jje=jj_end-1
      
      CALL  filtreg_p ( rotin ,jjb,jje,jjm, klevel,2, 1, .FALSE., 1)

      CALL   nxgrad_p ( klevel, rotin,   ghx ,  ghy               )
      CALL   rotatf_p  ( klevel, ghx  ,   ghy , rotout             )
c
      RETURN
      END
