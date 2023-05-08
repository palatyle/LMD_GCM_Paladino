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
#include "comgeomfi.h"
#include "comcstfi.h"
#include "watercap.h"

c
      INTEGER ngrid,ig,icap,iq
      REAL  piceco2(ngrid),psolaralb(ngrid,2)
      REAL qsurf(ngrid,nqmx) !tracer on surface (kg/m2)

      EXTERNAL ISMIN,ISMAX
      INTEGER ISMIN,ISMAX
c
c=======================================================================

      do ig=1,ngridmx

         !write(*,*) "all qsurf to zero. dirty."
         do iq=1,nqmx
         qsurf(ig,iq)=0.  !! on jette les inputs GCM
                          !! on regle juste watercaptag
                          !! il faudrait garder les inputs GCM 
                          !! si elles sont consequentes
         enddo

!! AS: my previous stuff
!         if ( ( lati(ig)*180./pi .gt. 70. ) .and.
!     .        ( albedodat(ig) .ge. 0.26   ) )  then
!                 write(*,*)"outlier ",ig,albedodat(ig),inertiedat(ig)
!                 watercaptag(ig)  = .true.
!                 dryness(ig)      = 1.
!                 albedodat(ig)    = 0.45 !albedo_h2o_ice  !! pour output 
!                 inertiedat(ig)   = 800.
!                 write(*,*)"new values ",ig,albedodat(ig),inertiedat(ig)
!         else
!                 watercaptag(ig)  = .false.
!                 dryness(ig)      = 1.
!         endif 
!         if ( inertiedat(ig) .ge. 800. ) then
!            write(*,*)"change inertie from ",inertiedat(ig)," to 800."
!            inertiedat(ig)   = 800.
!         endif

!! using Tyler and Barnes maps.
         dryness(ig)      = 1.
         if ( inertiedat(ig) .ge. 600. ) then
                 write(*,*)"ice ",ig,albedodat(ig),inertiedat(ig)
                 watercaptag(ig)  = .true.
         else
                 watercaptag(ig)  = .false.
         endif
         
      enddo


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
