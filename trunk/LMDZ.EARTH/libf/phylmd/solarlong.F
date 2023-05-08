      SUBROUTINE solarlong(pday,psollong,pdist_sol)

      USE ioipsl

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
c   lwrite        clef logique pour sorties de controle
c
c   Output:
c   -------
c   pdist_sol     distance entre le soleil et la planete
c                 ( en unite astronomique pour utiliser la constante 
c                  solaire terrestre 1370 Wm-2 )
c   pdecli        declinaison ( en radians )
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "planete.h"
#include "YOMCST.h"
      include 'iniprint.h'

c arguments:
c ----------

      REAL pday,pdist_sol,pdecli,psollong
      LOGICAL lwrite

c Local:
c ------

      REAL zanom,xref,zx0,zdx,zteta,zz,pi
      INTEGER iter
      REAL :: pyear_day,pperi_day
      REAL :: jD_eq, jD_peri
      LOGICAL,SAVE :: first=.TRUE.
c$OMP THREADPRIVATE(first)

c-----------------------------------------------------------------------
c calcul de l'angle polaire et de la distance au soleil :
c -------------------------------------------------------

c   Initialisation eventuelle:
      if(first) then
        call ioget_calendar(pyear_day)
        call ymds2ju(2000, 3, 21, 0., jD_eq)
        call ymds2ju(2001, 1, 4, 0., jD_peri)
        pperi_day = jD_peri - jD_eq
        pperi_day = R_peri + 180.
        write(lunout,*)' Number of days in a year = ',pyear_day
c         call iniorbit(249.22,206.66,669.,485.,25.2)
         call iniorbit(152.59,146.61,pyear_day,pperi_day,R_incl)
         first=.FALSE.
      endif

c  calcul de l'zanomalie moyenne

      zz=(pday-peri_day)/year_day
      pi=2.*asin(1.)
      zanom=2.*pi*(zz-nint(zz))
      xref=abs(zanom)

c  resolution de l'equation horaire  zx0 - e * sin (zx0) = xref
c  methode de Newton

!      zx0=xref+e_elips*sin(xref)
      zx0=xref+R_ecc*sin(xref)
      DO 110 iter=1,10
!         zdx=-(zx0-e_elips*sin(zx0)-xref)/(1.-e_elips*cos(zx0))
         zdx=-(zx0-R_ecc*sin(zx0)-xref)/(1.-R_ecc*cos(zx0))
         if(abs(zdx).le.(1.e-7)) goto 120
         zx0=zx0+zdx
110   continue
120   continue
      zx0=zx0+zdx
      if(zanom.lt.0.) zx0=-zx0

c zteta est la longitude solaire

!      zteta=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))
      zteta=2.*atan(sqrt((1.+R_ecc)/(1.-R_ecc))*tan(zx0/2.))

      psollong=zteta-timeperi

      IF(psollong.LT.0.) psollong=psollong+2.*pi
      IF(psollong.GT.2.*pi) psollong=psollong-2.*pi

      psollong = psollong * 180. / pi

c distance soleil

      pdist_sol = (1-R_ecc*R_ecc)
     &      /(1+R_ecc*COS(pi/180.*(psollong-(R_peri+180.0))))
!      pdist_sol = (1-e_elips*e_elips)
!     &      /(1+e_elips*COS(pi/180.*(psollong-(R_peri+180.0))))
c-----------------------------------------------------------------------
c   sorties eventuelles:
c   ---------------------

c     IF (lwrite) THEN
c        PRINT*,'jour de l"annee   :',pday
c        PRINT*,'distance au soleil (en unite astronomique) :',pdist_sol
c        PRINT*,'declinaison (en degres) :',pdecli*180./pi
c     ENDIF

      RETURN
      END
