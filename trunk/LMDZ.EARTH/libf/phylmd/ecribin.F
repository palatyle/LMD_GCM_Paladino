!
! $Header$
!
      SUBROUTINE ecribins(unit,pz)
      USE dimphy
      IMPLICIT none
c-----------------------------------------------------------------------
#include "dimensions.h"
cccc#include "dimphy.h"
#include "paramet.h"
#include "comgeom.h"
#include "comvert.h"
c
c   arguments:
c   ----------
      INTEGER unit
      REAL pz(klon)
c
c   local:
c   ------
      INTEGER i,j, ig
      REAL zz(iim +1,jjm+1)
c-----------------------------------------------------------------------
c   passage a la grille dynamique:
c   ------------------------------
         DO i=1,iim +1
            zz(i,1)=pz(1)
            zz(i,jjm+1)=pz(klon)
         ENDDO
c   traitement des point normaux
         DO j=2,jjm
            ig=2+(j-2)*iim
            CALL SCOPY(iim,pz(ig),1,zz(1,j),1)
            zz(iim+1,j)=zz(1,j)
         ENDDO
c-----------------------------------------------------------------------
#ifdef VPP
      CALL ecriture(unit,zz,(iim+1)*(jjm+1))
#else
      WRITE(unit) zz
#endif
c

      RETURN
      END
      SUBROUTINE ecribina(unit,pz)
      USE dimphy
      IMPLICIT none
c-----------------------------------------------------------------------
#include "dimensions.h"
cccc#include "dimphy.h"
#include "paramet.h"
#include "comgeom.h"
#include "comvert.h"
c
c   arguments:
c   ----------
      INTEGER unit
      REAL pz(klon,klev)
c
c   local:
c   ------
      INTEGER i,j,ilay,ig
      REAL zz(iim+1,jjm+1,llm)
c-----------------------------------------------------------------------
c   passage a la grille dynamique:
c   ------------------------------
      DO ilay=1,llm
c   traitement des poles
         DO i=1,iim +1
            zz(i,1,ilay)=pz(1,ilay)
            zz(i,jjm+1,ilay)=pz(klon,ilay)
         ENDDO
c   traitement des point normaux
         DO j=2,jjm
            ig=2+(j-2)*iim
            CALL SCOPY(iim,pz(ig,ilay),1,zz(1,j,ilay),1)
            zz(iim+1,j,ilay)=zz(1,j,ilay)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
      DO ilay = 1, llm
#ifdef VPP
         CALL ecriture(unit, zz(1,1,ilay), (iim+1)*(jjm+1))
#else
         WRITE(unit) ((zz(i,j,ilay),i=1,iim +1),j=1,jjm+1)
#endif
      ENDDO
c
      RETURN
      END
#ifdef VPP
@OPTIONS NODOUBLE
      SUBROUTINE ecriture(nunit, r8, n)
      INTEGER nunit, n, i
      REAL(KIND=8) r8(n)
      REAL r4(n)
      DO i = 1, n
         r4(i) = r8(i)
      ENDDO
      WRITE(nunit)r4
      RETURN
      END
#endif
