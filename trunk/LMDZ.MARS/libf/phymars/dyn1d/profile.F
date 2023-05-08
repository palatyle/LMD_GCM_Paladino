!      SUBROUTINE profile(unit,nlev,zkm,temp)
      SUBROUTINE profile(nlev,zkm,temp)
! to use  'getin'
      USE ioipsl_getincom 
      IMPLICIT NONE
c=======================================================================
c     Subroutine utilisee dans le modele 1-D  "testphys1d"
c     pour l'initialisation du profil atmospherique
c=======================================================================
c
c   differents profils d'atmospheres. T=f(z)
c   entree:
c     nlev    nombre de niveaux
c     zkm     alititudes en km
c     ichoice choix de l'atmosphere:
c             1 Temperature constante
c             2 profil Savidjari
c             3 Lindner (profil polaire)
c             4 Inversion pour Francois
c             5 Seiff (moyen)
c             6 T constante + perturbation gauss (level) (christophe 10/98)
c             7 T constante + perturbation gauss   (km)  (christophe 10/98)
c             8 Lecture du profile dans un fichier ASCII (profile)
c     tref    temperature de reference
c     isin    ajout d'une perturbation (isin=1)
c     pic     pic perturbation gauss pour ichoice = 6 ou 7
c     largeur largeur de la perturbation gauss pour ichoice = 6 ou 7
c     hauteur hauteur de la perturbation gauss pour ichoice = 6 ou 7
c
c   sortie:
c     temp    temperatures en K
c     
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------

c   arguments:
c   ----------

       INTEGER nlev
       REAL zkm(nlev),temp(nlev)

c   local:
c   ------

      INTEGER il,ichoice,nseiff,iseiff,isin,iter
      REAL pi
      PARAMETER(nseiff=37)
      REAL tref,t1,t2,t3,ww
      REAL tseiff(nseiff)
      DATA tseiff/214.,213.8,213.4,212.4,209.3,205.,201.4,197.8,
     $           194.6,191.4,188.2,185.2,182.5,180.,177.5,175,
     $           172.5,170.,167.5,164.8,162.4,160.,158.,156.,
     $           154.1,152.2,150.3,148.7,147.2,145.7,144.2,143.,
     $           142.,141.,140,139.5,139./
      REAL pic,largeur
      REAL hauteur,tmp

c-----------------------------------------------------------------------
c   read input profile type:
c--------------------------

      ichoice=1 ! default value for ichoice
      call getin("ichoice",ichoice)
      tref=200 ! default value for tref
      call getin("tref",tref)
      isin=0 ! default value for isin (=0 means no perturbation)
      call getin("isin",isin)
      pic=26.522 ! default value for pic
      call getin("pic",pic)
      largeur=10 ! default value for largeur
      call getin("largeur",largeur)
      hauteur=30 ! default value for hauteur
      call getin("hauteur",hauteur)

c-----------------------------------------------------------------------
c   ichoice=1 temperature constante:
c   --------------------------------

      IF(ichoice.EQ.1) THEN
         DO il=1,nlev
            temp(il)=tref
         ENDDO

c-----------------------------------------------------------------------
c   ichoice=2 savidjari:
c   --------------------

      ELSE IF(ichoice.EQ.2) THEN
         DO il=1,nlev
            temp(il)=AMAX1(219.-2.5*zkm(il),140.)
         ENDDO

c-----------------------------------------------------------------------
c   ichoice=3 Lindner:
c   ------------------

      ELSE IF(ichoice.EQ.3) THEN
         DO il=1,nlev
            IF(zkm(il).LT.2.5) THEN
               temp(il)=150.+30.*zkm(il)/2.5
            ELSE IF(zkm(il).LT.5.) THEN
               temp(il)=180.
            ELSE
               temp(il)=AMAX1(180.-2.*(zkm(il)-5.),130.)
            ENDIF
         ENDDO

c-----------------------------------------------------------------------
c   ichoice=4 Inversion pour Francois:
c   ------------------

      ELSE IF(ichoice.EQ.4) THEN
         DO il=1,nlev
            IF(zkm(il).LT.20.) THEN
               temp(il)=135.
            ELSE
               temp(il)=AMIN1(135.+5.*(zkm(il)-20.),200.)
            ENDIF
         ENDDO


c-----------------------------------------------------------------------
c   ichoice=5 Seiff:
c   ----------------

      ELSE IF(ichoice.EQ.5) THEN
	 DO il=1,nlev
	    iseiff=INT(zkm(il)/2.)+1
	    IF(iseiff.LT.nseiff-1) THEN
	       temp(il)=tseiff(iseiff)+(zkm(il)-2.*(iseiff-1))*
     $         (tseiff(iseiff+1)-tseiff(iseiff))/2.
	    ELSE
	       temp(il)=tseiff(nseiff)
	    ENDIF
	 ENDDO
c IF(ichoice.EQ.6) THEN
c	    DO iter=1,6
c	    t2=temp(1)
c	    t3=temp(2)
c	    DO il=2,nlev-1
c	       t1=t2
c	       t2=t3
c	       t3=temp(il+1)
c	       ww=tanh(zkm(il)/20.)
c	       ww=ww*ww*ww
c	       temp(il)=t2+ww*.5*(t1-2.*t2+t3)
c	    ENDDO
c	    ENDDO
c	 ENDIF

c-----------------------------------------------------------------------
c   ichoice=6 
c   ---------

      ELSE IF(ichoice.EQ.6) THEN
      DO il=1,nlev
        tmp=il-pic
        temp(il)=tref + hauteur*exp(-tmp*tmp/largeur/largeur)
      ENDDO


c-----------------------------------------------------------------------
c   ichoice=7
c   ---------

      ELSE IF(ichoice.EQ.7) THEN
      DO il=1,nlev
        tmp=zkm(il)-pic
        temp(il)=tref + hauteur*exp(-tmp*tmp*4/largeur/largeur)
      ENDDO

c-----------------------------------------------------------------------
c   ichoice=8
c   ---------

      ELSE IF(ichoice.GE.8) THEN
      OPEN(11,file='profile',status='old',form='formatted',err=101)
      DO il=1,nlev
        READ (11,*) temp(il)
      ENDDO

      GOTO 201
101   STOP'fichier profile inexistant'
201   CONTINUE
      CLOSE(10)

c-----------------------------------------------------------------------

      ENDIF

c-----------------------------------------------------------------------
c   rajout eventuel d'une perturbation:
c   -----------------------------------

      IF(isin.EQ.1) THEN
	 pi=2.*ASIN(1.)
	 DO il=1,nlev
c       if (nlev.EQ.501) then
c         if (zkm(il).LE.70.5) then
c       temp(il)=temp(il)+(1.-1000./(1000+zkm(il)*zkm(il)))*(
c    s      6.*SIN(zkm(il)*pi/6.)+9.*SIN(zkm(il)*pi/10.3) )
c         endif
c       else
        temp(il)=temp(il)+(1.-1000./(1000+zkm(il)*zkm(il)))*(
     s      6.*SIN(zkm(il)*pi/6.)+9.*SIN(zkm(il)*pi/10.3) )
c       endif
	 ENDDO
      ENDIF


c-----------------------------------------------------------------------
c   Ecriture du profil de temperature dans un fichier profile.out
c   -------------------------------------------------------------


      OPEN(12,file='profile.out',form='formatted')
	  DO il=1,nlev
	    write(12,*) temp(il)
	  ENDDO
      CLOSE(12)

      RETURN
      END
