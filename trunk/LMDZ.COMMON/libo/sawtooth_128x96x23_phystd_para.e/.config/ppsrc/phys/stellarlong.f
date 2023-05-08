










      SUBROUTINE stellarlong(pday,pstellong)
      
      USE planete_mod, ONLY: year_day, peri_day, e_elips, timeperi
      use comcstfi_mod, only: pi
      IMPLICIT NONE

c=======================================================================
c
c   Objet:
c   ------
c
c      Calcul de la distance soleil-planete et de la declinaison
c   en fonction du jour de l'annee.
c
c
c   Methode:
c   --------
c
c      Calcul complet de l'elipse
c
c   Interface:
c   ----------
c
c      Uncommon comprenant les parametres orbitaux.
c
c   Arguments:
c   ----------
c
c   Input:
c   ------
c   pday          jour de l'annee (le jour 0 correspondant a l'equinoxe)
c
c   Output:
c   -------
c   pdist_star     distance entre le soleil et la planete
c                 ( en unite astronomique pour utiliser la constante 
c                  solaire terrestre 1370 Wm-2 )
c   pdecli        declinaison ( en radians )
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

c arguments:
c ----------

      REAL pday,pdist_star,pdecli,pstellong
      LOGICAL lwrite

c Local:
c ------

      REAL zanom,xref,zx0,zdx,zteta,zz
      INTEGER iter


c-----------------------------------------------------------------------
c calcul de l'angle polaire et de la distance au soleil :
c -------------------------------------------------------

c  calcul de l'zanomalie moyenne

      zz=(pday-peri_day)/year_day
      pi=2.*asin(1.)
      zanom=2.*pi*(zz-nint(zz))
      xref=abs(zanom)

c  resolution de l'equation horaire  zx0 - e * sin (zx0) = xref
c  methode de Newton

      zx0=xref+e_elips*sin(xref)
      DO 110 iter=1,10
         zdx=-(zx0-e_elips*sin(zx0)-xref)/(1.-e_elips*cos(zx0))
         if(abs(zdx).le.(1.e-7)) goto 120
         zx0=zx0+zdx
110   continue
120   continue
      zx0=zx0+zdx
      if(zanom.lt.0.) zx0=-zx0

c zteta est la longitude solaire

      zteta=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))

      pstellong=zteta-timeperi

      IF(pstellong.LT.0.) pstellong=pstellong+2.*pi
      IF(pstellong.GT.2.*pi) pstellong=pstellong-2.*pi
c-----------------------------------------------------------------------
c   sorties eventuelles:
c   ---------------------

c     IF (lwrite) THEN
c        PRINT*,'jour de l"annee   :',pday
c        PRINT*,'distance au soleil (en unite astronomique) :',pdist_star
c        PRINT*,'declinaison (en degres) :',pdecli*180./pi
c     ENDIF

      RETURN
      END
