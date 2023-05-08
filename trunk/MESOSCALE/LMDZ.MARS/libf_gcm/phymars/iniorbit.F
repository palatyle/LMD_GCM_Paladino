      SUBROUTINE iniorbit
     $     (paphelie,pperiheli,pyear_day,pperi_day,pobliq)
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:
c   -------
c     Frederic Hourdin      22 Fevrier 1991
c
c   Objet:
c   ------
c    Initialisation du sous programme orbite qui calcule
c    a une date donnee de l'annee de duree year_day commencant
c    a l'equinoxe de printemps et dont le perihelie se situe
c    a la date peri_day, la distance au soleil et la declinaison.
c
c   Interface:
c   ----------
c   - Doit etre appele avant d'utiliser orbite.
c   - initialise une partie du common planete.h
c
c   Arguments:
c   ----------
c
c   Input:
c   ------
c   aphelie       \   aphelie et perihelie de l'orbite
c   periheli      /   en millions de kilometres.
c
c=======================================================================

c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "planete.h"
#include "comcstfi.h"

c   Arguments:
c   ----------

      REAL paphelie,pperiheli,pyear_day,pperi_day,pobliq

c   Local:
c   ------

      REAL zxref,zanom,zz,zx0,zdx
      INTEGER iter

c-----------------------------------------------------------------------

      pi=2.*asin(1.)

      aphelie =paphelie
      periheli=pperiheli
      year_day=pyear_day
      obliquit=pobliq
      peri_day=pperi_day

      PRINT*,'Perihelie en Mkm  ',periheli
      PRINT*,'Aphelise  en Mkm  ',aphelie 
      PRINT*,'obliquite en degres  :',obliquit
      unitastr=149.597927
      e_elips=(aphelie-periheli)/(periheli+aphelie)
      p_elips=0.5*(periheli+aphelie)*(1-e_elips*e_elips)/unitastr

      print*,'e_elips',e_elips
      print*,'p_elips',p_elips

c-----------------------------------------------------------------------
c calcul de l'angle polaire et de la distance au soleil :
c -------------------------------------------------------

c  calcul de l'zanomalie moyenne

      zz=(year_day-pperi_day)/year_day
      zanom=2.*pi*(zz-nint(zz))
      zxref=abs(zanom)
      PRINT*,'zanom  ',zanom

c  resolution de l'equation horaire  zx0 - e * sin (zx0) = zxref
c  methode de Newton

      zx0=zxref+e_elips*sin(zxref)
      DO 110 iter=1,100
         zdx=-(zx0-e_elips*sin(zx0)-zxref)/(1.-e_elips*cos(zx0))
         if(abs(zdx).le.(1.e-12)) goto 120
         zx0=zx0+zdx
110   continue
120   continue
      zx0=zx0+zdx
      if(zanom.lt.0.) zx0=-zx0
      PRINT*,'zx0   ',zx0

c zteta est la longitude solaire

      timeperi=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))
      PRINT*,'longitude solaire du perihelie timeperi = ',timeperi

      RETURN
      END
