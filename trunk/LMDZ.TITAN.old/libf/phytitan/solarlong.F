      SUBROUTINE solarlong(pday,psollong)
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
c      Calcul complet de l'ellipse
c
c   Interface:
c   ----------
c
c      Un common comprenant les parametres orbitaux.
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
c   psollong     Longitude solaire
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "comorbit.h"

c arguments:
c ----------

      REAL pday,psollong

c Local:
c ------

      REAL zanom,xref,zx0,zdx,zteta,zz
      INTEGER iter


c-----------------------------------------------------------------------
c calcul de l'angle polaire et de la distance au soleil :
c -------------------------------------------------------

c  calcul de l'zanomalie moyenne

      zz=(pday-peri_day)/year_day
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

      psollong=zteta-timeperi

      IF(psollong.LT.0.) psollong=psollong+2.*pi
      IF(psollong.GT.2.*pi) psollong=psollong-2.*pi
c-----------------------------------------------------------------------
c   sorties eventuelles:
c   ---------------------

      RETURN
      END
