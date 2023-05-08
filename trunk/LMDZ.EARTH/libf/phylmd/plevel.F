!
! $Header$
!
c================================================================
c================================================================
      SUBROUTINE plevel(ilon,ilev,lnew,pgcm,pres,Qgcm,Qpres)
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

      INTEGER ilon, ilev
      logical lnew

      REAL pgcm(ilon,ilev)
      REAL Qgcm(ilon,ilev)
      real pres
      REAL Qpres(ilon)

c   local :
c   -------

cym      INTEGER lt(klon), lb(klon)
cym      REAL ptop, pbot, aist(klon), aisb(klon)

cym      save lt,lb,ptop,pbot,aist,aisb
      INTEGER,ALLOCATABLE,SAVE,DIMENSION(:) :: lt,lb
      REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: aist,aisb
c$OMP THREADPRIVATE(lt,lb,aist,aisb)      
      REAL,SAVE :: ptop, pbot
c$OMP THREADPRIVATE(ptop, pbot)      
      LOGICAL,SAVE :: first = .true.
c$OMP THREADPRIVATE(first)
      INTEGER i, k
c
      REAL missing_val
c
      missing_val=nf90_fill_real
c
      if (first) then
        allocate(lt(klon),lb(klon),aist(klon),aisb(klon))
	first=.false.
      endif
      
c=====================================================================
      if (lnew) then
c   on r�nitialise les r�ndicages et les poids
c=====================================================================


c Chercher les 2 couches les plus proches du niveau a obtenir
c
c Eventuellement, faire l'extrapolation a partir des deux couches
c les plus basses ou les deux couches les plus hautes:
      DO 130 i = 1, klon
         IF ( ABS(pres-pgcm(i,ilev) ) .LT.
     .        ABS(pres-pgcm(i,1)) ) THEN
            lt(i) = ilev     ! 2
            lb(i) = ilev-1   ! 1
         ELSE
            lt(i) = 2
            lb(i) = 1
         ENDIF
  130 CONTINUE
      DO 150 k = 1, ilev-1
         DO 140 i = 1, klon
            pbot = pgcm(i,k)
            ptop = pgcm(i,k+1)
            IF (ptop.LE.pres .AND. pbot.GE.pres) THEN
               lt(i) = k+1
               lb(i) = k
            ENDIF
  140    CONTINUE
  150 CONTINUE
c
c Interpolation lineaire:
c
      DO i = 1, klon
c interpolation en logarithme de pression:
c
c ...   Modif . P. Le Van    ( 20/01/98) ....
c       Modif Fr��ic Hourdin (3/01/02)

        aist(i) = LOG( pgcm(i,lb(i))/ pres )
     .       / LOG( pgcm(i,lb(i))/ pgcm(i,lt(i)) )
        aisb(i) = LOG( pres / pgcm(i,lt(i)) )
     .       / LOG( pgcm(i,lb(i))/ pgcm(i,lt(i)))
      enddo


      endif ! lnew

c======================================================================
c    inteprollation
c======================================================================

      do i=1,klon
         Qpres(i)= Qgcm(i,lb(i))*aisb(i)+Qgcm(i,lt(i))*aist(i)
      enddo
c
c Je mets les vents a zero quand je rencontre une montagne
      do i = 1, klon
         if (pgcm(i,1).LT.pres) THEN
            Qpres(i)=missing_val
         endif
      enddo

c
      RETURN
      END
