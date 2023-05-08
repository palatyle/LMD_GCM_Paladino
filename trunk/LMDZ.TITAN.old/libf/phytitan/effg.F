      FUNCTION EFFG(Z)
#include "YOMCST.h"
! RA en m, Z en m...

! Quand on prendra atmosphere epaisse dans dynamique 
!    (et dans physique, attention a clmain et autres...)

      EFFG = RG * (RA/(RA + Z ) )**2

! Pour l'instant:
!     EFFG = RG 
      RETURN
      END

