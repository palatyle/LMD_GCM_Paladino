      SUBROUTINE laplacien_rot ( klevel, rotin, rotout,ghx,ghy )
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
      EXTERNAL  filtreg, nxgrad, rotatf
c
      CALL filtreg ( rotin ,   jjm, klevel,   2, 1, .FALSE., 1 )

      CALL nxgrad ( klevel, rotin,   ghx ,  ghy               )
      CALL rotatf  ( klevel, ghx  ,   ghy , rotout             )
c
      RETURN
      END
