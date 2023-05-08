      SUBROUTINE groupeun_p(jjmax,llmax,jjb,jje,q)
      USE parallel_lmdz
      USE comconst_mod, ONLY: ngroup
      USE Write_Field_p
      IMPLICIT NONE

#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"

      INTEGER jjmax,llmax,jjb,jje
      REAL q(iip1,jjmax,llmax)

!      INTEGER ngroup
!      PARAMETER (ngroup=3)

      REAL airecn,qn
      REAL airecs,qs

      INTEGER i,j,l,ig,ig2,j1,j2,i0,jd

c--------------------------------------------------------------------c 
c Strategie d'optimisation                                           c
c stocker les valeurs systematiquement recalculees                   c
c et identiques d'un pas de temps sur l'autre. Il s'agit des         c
c aires des cellules qui sont sommees. S'il n'y a pas de changement  c
c de grille au cours de la simulation tout devrait bien se passer.   c
c Autre optimisation : determination des bornes entre lesquelles "j" c
c varie, au lieu de faire un test à chaque fois...
c--------------------------------------------------------------------c 

      INTEGER j_start, j_finish

      REAL, SAVE :: airen_tab(iip1,jjp1,0:1)
      REAL, SAVE :: aires_tab(iip1,jjp1,0:1)
!$OMP THREADPRIVATE(airen_tab, aires_tab)

      LOGICAL, SAVE :: first = .TRUE.
!$OMP THREADPRIVATE(first)
!      INTEGER,SAVE :: i_index(iim,ngroup)
      INTEGER      :: offset
!      REAL         :: qsum(iim/ngroup)

      IF (first) THEN
         CALL INIT_GROUPEUN_P(airen_tab, aires_tab)
         first = .FALSE.
      ENDIF

c Champs 3D
      jd=jjp1-jjmax
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         j1=1+jd
         j2=2
         DO ig=1,ngroup

c     Concerne le pole nord
            j_start  = MAX(jjb, j1-jd)
            j_finish = MIN(jje, j2-jd)
            DO ig2=1,ngroup-ig+1
              offset=2**(ig2-1)
              DO j=j_start, j_finish
!CDIR NODEP
!CDIR ON_ADB(q)
                 DO i0=1,iim,2**ig2
                   q(i0,j,l)=q(i0,j,l)+q(i0+offset,j,l) 
                 ENDDO
              ENDDO
            ENDDO
            
            DO j=j_start, j_finish
!CDIR NODEP
!CDIR ON_ADB(q)
               DO i=1,iim
                 q(i,j,l)=q(i-MOD(i-1,2**(ngroup-ig+1)),j,l)
               ENDDO
            ENDDO

            DO j=j_start, j_finish
!CDIR ON_ADB(airen_tab)
!CDIR ON_ADB(q)
               DO i=1,iim
                 q(i,j,l)=q(i,j,l)*airen_tab(i,j,jd)
               ENDDO
               q(iip1,j,l)=q(1,j,l)
            ENDDO
       
!c     Concerne le pole sud
            j_start  = MAX(1+jjp1-jje-jd, j1-jd)
            j_finish = MIN(1+jjp1-jjb-jd, j2-jd)
            DO ig2=1,ngroup-ig+1
              offset=2**(ig2-1)
              DO j=j_start, j_finish
!CDIR NODEP
!CDIR ON_ADB(q)
                 DO i0=1,iim,2**ig2
                   q(i0,jjp1-j+1-jd,l)= q(i0,jjp1-j+1-jd,l)
     &                                 +q(i0+offset,jjp1-j+1-jd,l) 
                 ENDDO
              ENDDO
            ENDDO


            DO j=j_start, j_finish
!CDIR NODEP
!CDIR ON_ADB(q)
               DO i=1,iim
                 q(i,jjp1-j+1-jd,l)=q(i-MOD(i-1,2**(ngroup-ig+1)),
     &                                jjp1-j+1-jd,l)
               ENDDO
            ENDDO

            DO j=j_start, j_finish
!CDIR ON_ADB(aires_tab)
!CDIR ON_ADB(q)
               DO i=1,iim
                 q(i,jjp1-j+1-jd,l)=q(i,jjp1-j+1-jd,l)*  
     &                              aires_tab(i,jjp1-j+1,jd)
               ENDDO
               q(iip1,jjp1-j+1-jd,l)=q(1,jjp1-j+1-jd,l)
            ENDDO

        
            j1=j2+1
            j2=j2+2**ig
         ENDDO
      ENDDO
!$OMP END DO NOWAIT

      RETURN
      END



      SUBROUTINE INIT_GROUPEUN_P(airen_tab, aires_tab)

      USE parallel_lmdz
      USE comconst_mod, ONLY: ngroup
      IMPLICIT NONE

#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"

!      INTEGER ngroup
!      PARAMETER (ngroup=3)

      REAL airen,airecn
      REAL aires,airecs

      INTEGER i,j,l,ig,j1,j2,i0,jd

      INTEGER j_start, j_finish

      REAL :: airen_tab(iip1,jjp1,0:1)
      REAL :: aires_tab(iip1,jjp1,0:1)

      DO jd=0, 1
         j1=1+jd
         j2=2
         DO ig=1,ngroup
            
!     c     Concerne le pole nord
            j_start = j1-jd
            j_finish = j2-jd
            DO j=j_start, j_finish
               DO i0=1,iim,2**(ngroup-ig+1)
                  airen=0.
                  DO i=i0,i0+2**(ngroup-ig+1)-1
                     airen = airen+aire(i,j)
                  ENDDO
                  DO i=i0,i0+2**(ngroup-ig+1)-1
                     airen_tab(i,j,jd) = 
     &                    aire(i,j) / airen
                  ENDDO
               ENDDO
            ENDDO
            
!     c     Concerne le pole sud
            j_start = j1-jd
            j_finish = j2-jd
            DO j=j_start, j_finish
               DO i0=1,iim,2**(ngroup-ig+1)
                  aires=0.
                  DO i=i0,i0+2**(ngroup-ig+1)-1
                     aires=aires+aire(i,jjp1-j+1)
                  ENDDO
                  DO i=i0,i0+2**(ngroup-ig+1)-1
                     aires_tab(i,jjp1-j+1,jd) = 
     &                    aire(i,jjp1-j+1) / aires
                  ENDDO
               ENDDO
            ENDDO
            
            j1=j2+1
            j2=j2+2**ig
         ENDDO
      ENDDO
      
      RETURN
      END
