      SUBROUTINE profile(unit,nlev,dzst,temp)
      IMPLICIT NONE
c=======================================================================
c     Subroutine utilisee dans le modele 1-D  "rcm1d"
c     pour l'initialisation du profil atmospherique
c=======================================================================
c
c   VERSION VENUS
c
c   differents profils d'atmospheres. T=f(z)
c   entree:
c     unit    unite de lecture de "rcm1d.def"
c     nlev    nombre de niveaux (nlev=llm+1, surf + couches 1 a llm)
c     dzst    dz/T (avec dz = epaisseur de la couche en m)
c     ichoice choix de l'atmosphere:
c             1 Temperature constante
c             2 profil VIRA lisse
c             3 
c             4 
c             5 
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

       INTEGER nlev, unit
       REAL dzst(nlev),temp(nlev)

c   local:
c   ------

      INTEGER il,ichoice,isin,iter
      REAL pi
      REAL tref,t1,t2,t3,ww
      REAL pic,largeur
      REAL hauteur,tmp
      REAL zkm(nlev)    ! altitude en km

      isin = 0

c-----------------------------------------------------------------------
c   choix du profil:
c   ----------------

c la lecture se fait dans le rcm1d.def, ouvert par rcm1d.F
      READ(unit,*) 
      READ(unit,*) 
      READ(unit,*) 
      READ(unit,*) ichoice
      READ(unit,*) tref
      READ(unit,*) isin
      READ(unit,*) pic
      READ(unit,*) largeur
      READ(unit,*) hauteur

c-----------------------------------------------------------------------
c   ichoice=1 temperature constante:
c   --------------------------------

      IF(ichoice.EQ.1) THEN
         temp(1) = tref
         zkm(1)  = 0.0
         DO il=2,nlev
            temp(il)= tref
            zkm(il) = zkm(il-1)+temp(il)*dzst(il)/1000.
         ENDDO

c-----------------------------------------------------------------------
c   ichoice=2 VIRA lisse:
c   ---------------------

      ELSE IF(ichoice.EQ.2) THEN
         temp(1) = 735.
         zkm(1)  = 0.0
         DO il=2,nlev
            zkm(il) = zkm(il-1)+temp(il-1)*dzst(il)/1000. ! approx avec T(l-1)
            if(zkm(il).lt.60.) then
              temp(il)=735.-7.95*zkm(il)
            else
              temp(il)=AMAX1(258.-3.*(zkm(il)-60.),168.)
            endif
            zkm(il) = zkm(il-1)+(temp(il-1)+temp(il))/2.*dzst(il)/1000.
         ENDDO

c-----------------------------------------------------------------------
c   ichoice=3 VIRA lisse - 135K:
c   ----------------------------

      ELSE IF(ichoice.EQ.3) THEN
         temp(1) = 600.
         zkm(1)  = 0.0
         DO il=2,nlev
            zkm(il) = zkm(il-1)+temp(il-1)*dzst(il)/1000. ! approx avec T(l-1)
            temp(il)=AMAX1(600.-7.95*zkm(il),168.)
            zkm(il) = zkm(il-1)+(temp(il-1)+temp(il))/2.*dzst(il)/1000.
         ENDDO

c-----------------------------------------------------------------------
c   ichoice=4 :
c   ------------------

      ELSE IF(ichoice.EQ.4) THEN
         print*,"Cas non defini..."
         print*,"Stop dans profile.F"
         STOP

c-----------------------------------------------------------------------
c   ichoice=5 :
c   ----------------

      ELSE IF(ichoice.EQ.5) THEN
         print*,"Cas non defini..."
         print*,"Stop dans profile.F"
         STOP

c-----------------------------------------------------------------------
c   ichoice=6 
c   ---------

      ELSE IF(ichoice.EQ.6) THEN
      temp(1) = tref
      zkm(1)  = 0.0
      DO il=2,nlev
        tmp=il-pic
        temp(il)= tref + hauteur*exp(-tmp*tmp/largeur/largeur)
        zkm(il) = zkm(il-1)+temp(il)*dzst(il)/1000.
      ENDDO


c-----------------------------------------------------------------------
c   ichoice=7
c   ---------

      ELSE IF(ichoice.EQ.7) THEN
      temp(1) = tref
      zkm(1)  = 0.0
      DO il=2,nlev
        zkm(il) = zkm(il-1)+tref*dzst(il)/1000. ! approx
        tmp=zkm(il)-pic
        temp(il)= tref + hauteur*exp(-tmp*tmp*4/largeur/largeur)
        zkm(il) = zkm(il-1)+(temp(il-1)+temp(il))/2.*dzst(il)/1000.
      ENDDO

c-----------------------------------------------------------------------
c   ichoice=8
c   ---------

      ELSE IF(ichoice.GE.8) THEN
      OPEN(11,file='profile',status='old',form='formatted',err=101)
      DO il=1,nlev
        READ (11,*) temp(il)
      ENDDO
      zkm(1) = 0.0
      DO il=2,nlev
        zkm(il) = zkm(il-1)+(temp(il-1)+temp(il))/2.*dzst(il)/1000.
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
        temp(il)=temp(il)+(1.-1000./(1000+zkm(il)*zkm(il)))*(
     s      6.*SIN(zkm(il)*pi/6.)+9.*SIN(zkm(il)*pi/10.3) )
	 ENDDO
      ENDIF


c-----------------------------------------------------------------------
c   Ecriture du profil de temperature dans un fichier profile.out
c   -------------------------------------------------------------


c     OPEN(12,file='profile.out')
c         DO il=1,nlev
c            write(12,*) temp(il)
c           write(12,*) temp(il),zkm(il)
c         ENDDO
c     CLOSE(12)

      RETURN
      END
