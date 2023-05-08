










c================================================================
c================================================================
      SUBROUTINE tetaleveli1j1(ilon,ilev,lnew,pgcm,pres,Qgcm,Qpres)
c================================================================
c================================================================

! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
!      USE dimphy
      IMPLICIT none

!-----------------------------------------------------------------------
!   INCLUDE 'dimensions.h'
!
!   dimensions.h contient les dimensions du modele
!   ndm est tel que iim=2**ndm
!-----------------------------------------------------------------------

      INTEGER iim,jjm,llm,ndm

      PARAMETER (iim= 128,jjm=96,llm=23,ndm=1)

!-----------------------------------------------------------------------
cccc#include "dimphy.h"

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

cIM 211004
c     INTEGER lt(klon), lb(klon)
c     REAL ptop, pbot, aist(klon), aisb(klon)
c
!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
!-----------------------------------------------------------------------
!   INCLUDE 'paramet.h'

      INTEGER  iip1,iip2,iip3,jjp1,llmp1,llmp2,llmm1
      INTEGER  kftd,ip1jm,ip1jmp1,ip1jmi1,ijp1llm
      INTEGER  ijmllm,mvar
      INTEGER jcfil,jcfllm

      PARAMETER( iip1= iim+1,iip2=iim+2,iip3=iim+3                       &
     &    ,jjp1=jjm+1-1/jjm)
      PARAMETER( llmp1 = llm+1,  llmp2 = llm+2, llmm1 = llm-1 )
      PARAMETER( kftd  = iim/2 -ndm )
      PARAMETER( ip1jm  = iip1*jjm,  ip1jmp1= iip1*jjp1 )
      PARAMETER( ip1jmi1= ip1jm - iip1 )
      PARAMETER( ijp1llm= ip1jmp1 * llm, ijmllm= ip1jm * llm )
      PARAMETER( mvar= ip1jmp1*( 2*llm+1) + ijmllm )
      PARAMETER( jcfil=jjm/2+5, jcfllm=jcfil*llm )

!-----------------------------------------------------------------------
c
      INTEGER lt(ip1jmp1), lb(ip1jmp1)
      REAL ptop, pbot, aist(ip1jmp1), aisb(ip1jmp1)
cMI 211004
      save lt,lb,ptop,pbot,aist,aisb

      INTEGER i, k
c
c     PRINT*,'tetalevel pres=',pres
c=====================================================================
      if (lnew) then
c   on réinitialise les réindicages et les poids
c=====================================================================


c Chercher les 2 couches les plus proches du niveau a obtenir
c
c Eventuellement, faire l'extrapolation a partir des deux couches
c les plus basses ou les deux couches les plus hautes:
      DO 130 i = 1, ilon
cIM      IF ( ABS(pres-pgcm(i,ilev) ) .LT.
         IF ( ABS(pres-pgcm(i,ilev) ) .GT.
     .        ABS(pres-pgcm(i,1)) ) THEN
            lt(i) = ilev     ! 2
            lb(i) = ilev-1   ! 1
         ELSE
            lt(i) = 2
            lb(i) = 1
         ENDIF
cIM   PRINT*,'i, ABS(pres-pgcm),ABS(pres-pgcm)',
cIM  .i, ABS(pres-pgcm(i,ilev)),ABS(pres-pgcm(i,1))
  130 CONTINUE
      DO 150 k = 1, ilev-1
         DO 140 i = 1, ilon
            pbot = pgcm(i,k)
            ptop = pgcm(i,k+1)
cIM         IF (ptop.LE.pres .AND. pbot.GE.pres) THEN
            IF (ptop.GE.pres .AND. pbot.LE.pres) THEN
               lt(i) = k+1
               lb(i) = k
            ENDIF
  140    CONTINUE
  150 CONTINUE
c
c Interpolation lineaire:
c
      DO i = 1, ilon
c interpolation en logarithme de pression:
c
c ...   Modif . P. Le Van    ( 20/01/98) ....
c       Modif Frédéric Hourdin (3/01/02)

        IF(pgcm(i,lb(i)).EQ.0.OR.
     $     pgcm(i,lt(i)).EQ.0.) THEN
c
        PRINT*,'i,lb,lt,2pgcm,pres',i,lb(i),
     .  lt(i),pgcm(i,lb(i)),pgcm(i,lt(i)),pres
c
        ENDIF 
c
        aist(i) = LOG( pgcm(i,lb(i))/ pres )
     .       / LOG( pgcm(i,lb(i))/ pgcm(i,lt(i)) )
        aisb(i) = LOG( pres / pgcm(i,lt(i)) )
     .       / LOG( pgcm(i,lb(i))/ pgcm(i,lt(i)))
      enddo


      endif ! lnew

c======================================================================
c    inteprollation
c======================================================================

      do i=1,ilon
         Qpres(i)= Qgcm(i,lb(i))*aisb(i)+Qgcm(i,lt(i))*aist(i)
cIM      PRINT*,'i,Qgcm,Qpres',i,Qgcm(i,lb(i)),aisb(i),
cIM  $   Qgcm(i,lt(i)),aist(i),Qpres(i)
      enddo
c
c Je mets les vents a zero quand je rencontre une montagne
      do i = 1, ilon
cIM      if (pgcm(i,1).LT.pres) THEN
         if (pgcm(i,1).GT.pres) THEN
c           Qpres(i)=1e33
            Qpres(i)=1e+20
cIM         PRINT*,'i,pgcm(i,1),pres =',i,pgcm(i,1),pres
         endif
      enddo

c
      RETURN
      END
