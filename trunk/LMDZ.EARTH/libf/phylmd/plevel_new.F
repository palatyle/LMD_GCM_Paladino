!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/plevel.F,v 1.1.1.1.10.1 2006/08/17 15:41:51 fairhead Exp $
!
c================================================================
c================================================================
      SUBROUTINE plevel_new(ilon,ilev,klevSTD,lnew,pgcm,pres,Qgcm,Qpres)
c================================================================
c================================================================
      USE netcdf
      USE dimphy
      IMPLICIT none

cym#include "dimensions.h"
cy#include "dimphy.h"

c================================================================
c
c Interpoler des champs 3-D u, v et g du modele a un niveau de
c pression donnee (pres)
c
c INPUT:  ilon ----- nombre de points
c         ilev ----- nombre de couches
c         lnew ----- true si on doit reinitialiser les poids
c         pgcm ----- pressions modeles
c         pres ----- pression vers laquelle on interpolle
c         Qgcm ----- champ GCM
c         Qpres ---- champ interpolle au niveau pres
c
c================================================================
c
c   arguments :
c   -----------

      INTEGER ilon, ilev, klevSTD
      logical lnew
      
      REAL pgcm(ilon,ilev)
      REAL Qgcm(ilon,ilev)
      real pres(klevSTD)
      REAL Qpres(ilon, klevSTD)

c   local :
c   -------

cym      INTEGER lt(klon), lb(klon)
cym      REAL ptop, pbot, aist(klon), aisb(klon)

cym      save lt,lb,ptop,pbot,aist,aisb
      INTEGER,ALLOCATABLE,SAVE,DIMENSION(:,:) :: lt,lb
      REAL,ALLOCATABLE,SAVE,DIMENSION(:,:) :: aist,aisb
c$OMP THREADPRIVATE(lt,lb,aist,aisb)      
      REAL,SAVE :: ptop, pbot
c$OMP THREADPRIVATE(ptop, pbot)      
      LOGICAL,SAVE :: first = .true.
      INTEGER :: nlev
c$OMP THREADPRIVATE(first)
      INTEGER i, k
c
      REAL missing_val
c
      missing_val=nf90_fill_real
c
      if (first) then
         allocate(lt(klon,klevSTD),lb(klon,klevSTD))
         allocate(aist(klon,klevSTD),aisb(klon, klevSTD))
         first=.false.
      endif
      
c=====================================================================
      if (lnew) then
c   on reinitialise les reindicages et les poids
c=====================================================================


c Chercher les 2 couches les plus proches du niveau a obtenir
c
c Eventuellement, faire l'extrapolation a partir des deux couches
c les plus basses ou les deux couches les plus hautes:
c
c
         DO nlev = 1, klevSTD
            DO i = 1, klon
               IF ( ABS(pres(nlev)-pgcm(i,ilev) ) .LT.
     &              ABS(pres(nlev)-pgcm(i,1)) ) THEN
                  lt(i,nlev) = ilev  ! 2
                  lb(i,nlev) = ilev-1 ! 1
               ELSE
                  lt(i,nlev) = 2
                  lb(i,nlev) = 1
               ENDIF
            ENDDO
            DO k = 1, ilev-1
               DO i = 1, klon
                  pbot = pgcm(i,k)
                  ptop = pgcm(i,k+1)
                  IF (ptop.LE.pres(nlev) .AND. pbot.GE.pres(nlev)) THEN
                     lt(i,nlev) = k+1
                     lb(i,nlev) = k
                  ENDIF
               ENDDO
            ENDDO
            
c     Interpolation lineaire:
            DO i = 1, klon
c     interpolation en logarithme de pression:
c     
c     ...   Modif . P. Le Van    ( 20/01/98) ....
c     Modif Frederic Hourdin (3/01/02)

               aist(i,nlev) = LOG( pgcm(i,lb(i,nlev))/ pres(nlev) )
     &              / LOG( pgcm(i,lb(i,nlev))/ pgcm(i,lt(i,nlev)) )
               aisb(i,nlev) = LOG( pres(nlev) / pgcm(i,lt(i,nlev)) )
     &              / LOG( pgcm(i,lb(i,nlev))/ pgcm(i,lt(i,nlev)))
            ENDDO
         ENDDO

      ENDIF ! lnew

c======================================================================
c    inteprollation
c    ET je mets les vents a zero quand je rencontre une montagne
c======================================================================

      DO nlev = 1, klevSTD
         DO i=1,klon
            IF (pgcm(i,1).LT.pres(nlev)) THEN
               Qpres(i,nlev) = missing_val
            ELSE
               Qpres(i,nlev) = 
     &              Qgcm(i,lb(i,nlev))*aisb(i,nlev) +
     &              Qgcm(i,lt(i,nlev))*aist(i,nlev)
            ENDIF
         ENDDO
      ENDDO

c     
      RETURN
      END
