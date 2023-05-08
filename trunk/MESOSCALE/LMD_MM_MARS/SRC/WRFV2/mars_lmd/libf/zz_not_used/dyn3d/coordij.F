      SUBROUTINE coordij(lon,lat,ilon,jlat)

c=======================================================================
c
c   calcul des coordonnees i et j de la maille scalaire dans
c   laquelle se trouve le point (lon,lat) en radian
c
c=======================================================================

      IMPLICIT NONE
      REAL lon,lat
      INTEGER ilon,jlat
      INTEGER i,j

#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comgeom.h"
#include "serre.h"

      real zlon,zlat

      zlon=lon
      zlat=lat

      DO i=1,iim+1
        print*,i,' ',rlonu(i),' ',zlon 
         IF (rlonu(i).GT.zlon) THEN
            ilon=i
            GOTO 10
         ENDIF
      ENDDO
10    CONTINUE

      j=0
      DO j=1,jjm
        print*,j,' ',rlatv(j),' ' ,zlat 
         IF(rlatv(j).LT.zlat) THEN
           
            jlat=j
            GOTO 20
         ENDIF
      ENDDO
20    CONTINUE
      IF(j.EQ.0) j=jjm+1

      RETURN
      END
