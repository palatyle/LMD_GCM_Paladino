      PROGRAM anldoppler2
      IMPLICIT NONE
c======================================================================
c
c
c Programme pour analyser vent a heure locale fixee
c a partir d'un fichier Netcdf de la MCD
c
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------


#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comdissip.h"
#include "comvert.h"
#include "comgeom2.h"
#include "logic.h"
#include "temps.h"
#include "control.h"
#include "ener.h"
#include "netcdf.inc"
#include "description.h"


      INTEGER point
      PARAMETER(point=11)     		! nbre de points d''observation
      real olon(point),olat(point)	! Coordonnees des points observes
      INTEGER iraie                	! pour le choix de la raie du CO

c--------------------------------------------------------
c coordonees lon/lat des points observes
c--------------------------------------------------------

c     ANNEE 2001

      DATA olon/0,-30,30,0,-90,90,-75,75,0,-70,70/
      DATA olat/-1,0,0,40,0,0,40,40,90,-40,40/

c     ANNEE 1999

c      DATA olon/0,0,-30,30,-85,85,0,-90,90,0,-65,65,0/
c      DATA olat/18,-10,0,0,0,0,55,40,40,90,-40,-40,-70/

c     ANNEE 1997

c      DATA olon/0,-85,85,0,0/
c      DATA olat/24,0,0,90,-70/

c     ANNEE 1992

c      DATA olon/0,-85,85,0,0,-75,75,-90,90,0,0/
c      DATA olat/8,0,0,-40,40,-40,-40,40,40,90,-90/

c     ANNEE 1990

c      DATA olon/0,-90,90,0,0/
c      DATA olat/-1,0,0,70,-90/

c     ANNEE 1988

c      DATA olon/0,-90,90,0,0,0,-35,35,0,0,-45,45,-60,60/
c      DATA olat/-21,0,0,-90,60,0,0,0,-35,25,-35,-35,25,25/

c----------------------------------------------------------

      INTEGER mcdreadnc
      external mcdreadnc
      INTEGER itau,nbpas,nbpasmx
      PARAMETER(nbpasmx=1000000)
      REAL temps(nbpasmx)
      INTEGER unitlec
      INTEGER i,j,l,itautot

c   Declarations DRS:
c   -----------------
      INTEGER ierr, dimtype
      CHARACTER*120 dimsou, dimnam*16, dimtit*80
      CHARACTER*40 dimunit
      CHARACTER*100  varname

c   variables meteo dans la grille verticale GCM
c   --------------------------------------------
      REAL v(iip1,jjp1,llm),u(iip1,jjp1,llm)
      REAL t(iip1,jjp1,llm)
      REAL ps(iip1,jjp1) , tsurf(iip1,jjp1)

c   variables meteo en coordonnee de pression
c   --------------------------------------------
      REAL vp(iip1,jjp1,llm),up(iip1,jjp1,llm)
      REAL tp(iip1,jjp1,llm)

c     Niveaux de pression (Pa)
      real pref(llm)

c     version 25 couches

c      data pref/1291.37,1005.720,783.255,600.,475.069,350.,
c     &          288.144,174.768,106.002,64.2935,38.9960,
c     &          23.6523,14.3458,8.70000,5.28000,3.20000,
c     &          1.94000,1.18000,0.714000,0.433000,0.263000,
c     &          0.118000,0.0614000,0.0333000,0.0163000/

c     version 32 couches

      data pref/1291.37,1005.720,783.255,600.,475.069,350.,
     &          288.144,174.768,106.002,64.2935,38.9960,
     &          23.6523,14.3458,8.70000,5.28000,3.20000,
     &          1.94000,1.18000,0.714000,0.433000,0.263000,
     &          0.118000,0.0614000,0.0333000,0.0163000, 
     &          8.e-3,4.e-3,2.e-3,1.e-3,5.e-4,2.e-4,1.e-4/

c   Local:
c   ------
      REAL phis(iip1,jjp1)
      real utime,latst,lonst,tamp,hlst
      real localtime(iip1),heure(iip1)
      common/temporaire/localtime

      INTEGER*4 day0
c      integer nmoy(jjp1,llm),tp1,tp2,pass
	

      CHARACTER*1 yes

      LOGICAL ldrs
      integer iformat,nid
       
c   declarations de l'interface avec mywrite:
c   -----------------------------------------

      CHARACTER file*80
      CHARACTER nomfich*60

c   externe:
c   --------

      EXTERNAL iniconst,inigeom,covcont,mywrite
      EXTERNAL inifilr,exner
      EXTERNAL solarlong,coordij,moy2
      EXTERNAL SSUM
      REAL SSUM
      EXTERNAL lnblnk
      INTEGER lnblnk
      logical icdf,idrs

c   Dust
c   ----

#include "fxyprim.h"

c-----------------------------------------------------------------------
c   initialisations:
c   ----------------

      ldrs=.true.
      unitlec=11
      itautot=0
      idrs=.false.
      icdf=.false.

c     Lecture du fichier a lire
c     -------------------------

      PRINT*,'entrer le nom du fichier NetCDF de la MCD'
      READ(5,'(a)') nomfich
      file=nomfich(1:lnblnk(nomfich))
      PRINT*,'file',file


c     Lecture des autres parametres du programme
c     ------------------------------------------

      write(*,*) 'Raie a simuler ? '
      write(*,*)'[1]12CO(2-1) [2]:13CO(2-1) ',
     &     ' [3]:12CO(1-0) [4]13CO(1-0) [5]H2O(183)'  
      read(*,*) iraie

      write(*,*) 'heure local du point subterrestre ? '
      read(*,*) hlst

      write(*,*) 'latitude du point subterrestre ? '
      read(*,*) latst

c---------------------------------------------------------
c   Lecture des fichiers stats , diagfi ou histmoy
c---------------------------------------------------------

800   continue ! LOOP SUR LES FICHIERS

c     Ouverture fichiers :
c     ------------------

c        Ouverture fichier NetCDF
            ierr = NF_OPEN (file(1:lnblnk(file))//'.nc', NF_NOWRITE,nid)
            IF (ierr.NE.NF_NOERR) THEN
              write(6,*)' Pb d''ouverture du fichier ',file
              write(6,*)' ierr = ', ierr
              CALL abort
            ENDIF
     
c----------------------------------------------------------------------
c   initialisation de la physique:
c   ------------------------------

c initialisations aux valeurs martiennes

      rad=3397200.              ! rayon de mars (m)  ~3397200 m
      daysec=88775.             ! duree du sol (s)  ~88775 s
      omeg=4.*asin(1.)/(daysec) ! vitesse de rotation (rad.s-1)
      g=3.72                    ! gravite (m.s-2) ~3.72
      mugaz=43.49               ! Masse molaire de l''atm (g.mol-1) ~43.49
      kappa=.256793             ! = r/cp  ~0.256793
      r = 191.18213
      pi=2.*asin(1.)
      ecritphy =1
      iphysiq=1
      day_ini=0.
      day_step=1
   
      WRITE (*,*) 'day0 = ' , day0

      CALL iniconst
      CALL inigeom
      CALL inifilr

c   sigma aux niveaux s:

      DO l=1,llm
         print*,'sig_s(',l,')  =  ',sig_s(l)
      ENDDO

c      Format Grads (Stats) : On suppose que l'on a 12 pas de temps :
c      --------------------

      nbpas=12
      print*,'nbpas=',nbpas

c======================================================================
c   debut de la boucle sur les etats dans un fichier histoire:
c======================================================================

       itau= 1

801    continue ! LOOP SUR LES PAS DE TEMPS START HERE

       varname='u'
       ierr =mcdreadnc(nid,varname,llm,u,itau)
       varname='v'
       ierr =mcdreadnc(nid,varname,llm,v,itau)
       varname='temp'
       ierr =mcdreadnc(nid,varname,llm,t,itau)
       varname='ps'
       ierr =mcdreadnc(nid,varname,1,ps,itau)
c      PRINT*,'ps',ps(iip1/2,jjp1/2)
       varname='tsurf'
       ierr =mcdreadnc(nid,varname,1,tsurf,itau)
c      PRINT*,'tsurf',tsurf(iip1/2,jjp1/2)

c**********************************************************
c  Traitement des donnees : (we are in LOOP on timestep)
c**********************************************************


c universal time (in stats file)
c-------------------------------
      utime = (itau-1)*2.

c Longitude du point subterrestre
c -------------------------------

      tamp=(hlst-utime)*(iip1/24.)
      if(tamp.GT.iip1) then
      tamp=tamp-iip1
      endif
      if(tamp.LT.0)  then
      tamp=iip1+tamp
      endif
      lonst=rlonv(nint(tamp))*180./pi
      write(*,*) hlst,latst,lonst

       lonst=(hlst-utime)*15
       if(lonst.gt.180) lonst=lonst-180.
       if(lonst.lt.-180) lonst=lonst+180.
       write(*,*) hlst,latst,lonst

c Local time
c ----------

c local time (in stats file)
            DO i = 1,iim
              localtime(i)=utime + 12.*rlonv(i)/pi
              if(localtime(i).gt.24) localtime(i)=localtime(i)-24.
              if(localtime(i).lt.0) localtime(i)=localtime(i)+24.
	      heure(i)=localtime(i)
c	      if(itau.EQ.1) heure(i)=0.375*(i-1)
c	      if(itau.EQ.1) print*,heure(i)
            ENDDO

c              if(itau.EQ.1) heure(iip1)=24.
	       heure(iip1)=localtime(1)

c ******************************************


c c Passage en niveaux de pression
c --------------------------------

          call xsig2p (ps,T,pref,llm,Tp)
          call xsig2p (ps,u,pref,llm,up)
          call xsig2p (ps,v,pref,llm,vp)

c         Test temporaire
c          do l=1, llm
c             write(77,*) pref(l),up(1,1,l)
c          end do

c Traitement des donnees -> creation de spectres synthetiques
c -----------------------------------------------------------

c calcul des vitesses projetees

      do l=1,llm
         do j=1,jjp1
            do i=1,iip1


               up(i,j,l)=(242.*cos(rlatu(j))+up(i,j,l))*
     .              sin(rlonv(i)-lonst*pi/180)*cos(latst*pi/180)


           vp(i,j,l)=vp(i,j,l)*(cos(rlonv(i)-lonst*pi/180)*sin(rlatu(j))
     .    *cos(latst*pi/180)-sin(latst*pi/180)*cos(rlatu(j)))


            enddo
         enddo
      enddo

c synthese des spectres

c	print*,'itau=',itau,nbpas

      call spectre2(itau,iraie,latst,lonst,olon,olat,
     .vp,up,tp,ps,tsurf)


	
c  Fin de la boucle sur les pas de temps ? : (we are in LOOP on timestep)
c  ****************************************

      itau=itau+1

      if(itau.gt.nbpas) then
         ierr = NF_CLOSE(nid)
         write(*,*) 'Si vous voulez entrer un nouveau fichier'
         print*,'entrez son nom. Sinon, RETURN'
            READ(5,'(a)') nomfich
         if (nomfich(1:lnblnk(nomfich)).eq."") then
            print*,'OK, on termine la'
         else
            file=nomfich(1:lnblnk(nomfich))
            PRINT*,'file',file
            goto 800
         endif
      else
         goto 801
      endif

 900   continue


9999  PRINT*,'Fin '
1000  format(a5,3x,i4,i3,x,a39)
7777  FORMAT ('latitude/longitude',4f7.1)

      END
