!
! $Header$
!
      SUBROUTINE ecriregs(unit,pz)
      use dimphy
      IMPLICIT none
c-----------------------------------------------------------------------
#include "dimensions.h"
cccc#include "dimphy.h"
#include "paramet.h"
#include "comgeom.h"
#include "comvert.h"
#include "regdim.h"
c
c   arguments:
c   ----------
      INTEGER unit
      REAL pz(klon)
c
c   local:
c   ------
      INTEGER i,j, ig
      REAL zz(iim,jjm+1)
      INTEGER nleng
      PARAMETER (nleng=(i2_fin-i2_deb+1+i1_fin-i1_deb+1)
     .                *(j_fin-j_deb+1))
      REAL zzz(nleng)
c
c-----------------------------------------------------------------------
c   passage a la grille dynamique:
c   ------------------------------
         DO i=1,iim
            zz(i,1)=pz(1)
            zz(i,jjm+1)=pz(klon)
         ENDDO
c
c   traitement des point normaux
         DO j=2,jjm
            ig=2+(j-2)*iim
            CALL SCOPY(iim,pz(ig),1,zz(1,j),1)
         ENDDO
c-----------------------------------------------------------------------
      ig = 0
      DO j = j_deb, j_fin
         DO i=i1_deb,i1_fin
            ig = ig + 1
            zzz(ig) = zz(i,j)
         ENDDO
         DO i=i2_deb,i2_fin
            ig = ig + 1
            zzz(ig) = zz(i,j)
         ENDDO
      ENDDO
#ifdef VPP
      CALL ecriture(unit,zzz,nleng)
#else
      WRITE(unit) zzz
#endif
      RETURN
      END
      SUBROUTINE ecrirega(unit,pz)
      USE dimphy
      IMPLICIT none
c-----------------------------------------------------------------------
#include "dimensions.h"
cccc#include "dimphy.h"
#include "paramet.h"
#include "comgeom.h"
#include "comvert.h"
#include "regdim.h"
c
c   arguments:
c   ----------
      INTEGER unit
      REAL pz(klon,klev)
c
c   local:
c   ------
      INTEGER i,j,ilay,ig
      REAL zz(iim,jjm+1,llm)
      INTEGER nleng
      PARAMETER (nleng=(i2_fin-i2_deb+1+i1_fin-i1_deb+1)
     .                *(j_fin-j_deb+1))
      REAL zzz(nleng)
c-----------------------------------------------------------------------
c   passage a la grille dynamique:
c   ------------------------------
      DO ilay=1,llm
c   traitement des poles
         DO i=1,iim
            zz(i,1,ilay)=pz(1,ilay)
            zz(i,jjm+1,ilay)=pz(klon,ilay)
         ENDDO
c   traitement des point normaux
         DO j=2,jjm
            ig=2+(j-2)*iim
            CALL SCOPY(iim,pz(ig,ilay),1,zz(1,j,ilay),1)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
      DO ilay = 1, llm
      ig = 0
      DO j = j_deb, j_fin
         DO i=i1_deb,i1_fin
            ig = ig + 1
            zzz(ig) = zz(i,j,ilay)
         ENDDO
         DO i=i2_deb,i2_fin
            ig = ig + 1
            zzz(ig) = zz(i,j,ilay)
         ENDDO
      ENDDO
#ifdef VPP
      CALL ecriture(unit,zzz,nleng)
#else
      WRITE(unit) zzz
#endif
      ENDDO

      RETURN
      END
