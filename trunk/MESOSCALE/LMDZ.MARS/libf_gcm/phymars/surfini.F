      SUBROUTINE surfini(ngrid,piceco2,qsurf,psolaralb)
      IMPLICIT NONE
c=======================================================================
c
c   creation des calottes pour l'etat initial
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------
#include "dimensions.h"
#include "dimphys.h"
#include "surfdat.h"
#include "callkeys.h"
#include "tracer.h"
c
      INTEGER ngrid,ig,icap
      REAL  piceco2(ngrid),psolaralb(ngrid,2)
      REAL qsurf(ngrid,nqmx) !tracer on surface (kg/m2)

      EXTERNAL ISMIN,ISMAX
      INTEGER ISMIN,ISMAX
c
c=======================================================================

c
c     calcul de piceco2 (kg/m2) a l'etat initial
c     ------------------------------------------

      DO 100 ig=1,ngrid
         psolaralb(ig,1)=albedodat(ig)
         psolaralb(ig,2)=albedodat(ig)
100   CONTINUE

      PRINT*,'minimum des donnees albedo',
     s     albedodat(ISMIN(ngrid,albedodat,1))
      PRINT*,'maximum des donnees albedo',
     s     albedodat(ISMAX(ngrid,albedodat,1))
c      calcul de psolaralb 
c      -------------------
      DO 115 ig=1,ngrid

c	 IF (water) THEN 
c	   if (qsurf(ig,nqmx).gt.0.005) then
c             psolaralb(ig,1) = 0.4
c             psolaralb(ig,2) = 0.4 
c           endif 
c         ENDIF
c IF there is more than 5 pr. um of h2o ice but no C02 ice, surface albedo is set to 0.4. 
         IF (piceco2(ig) .GT. 0.) THEN
             IF(ig.GT.ngrid/2+1) THEN
                icap=2
             ELSE
                icap=1
             ENDIF
             psolaralb(ig,1) = albedice(icap)
             psolaralb(ig,2) =  albedice(icap)
         END IF
115   CONTINUE      

      PRINT*,'minimum des donnees albedo',
     s     psolaralb(ISMIN(ngrid,psolaralb,1),1)
      PRINT*,'maximum des donnees albedo',
     s     psolaralb(ISMAX(ngrid,psolaralb,1),1)

      RETURN
      END
