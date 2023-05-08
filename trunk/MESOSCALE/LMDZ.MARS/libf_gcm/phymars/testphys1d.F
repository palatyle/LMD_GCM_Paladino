
      PROGRAM testphys1d
      IMPLICIT NONE

c=======================================================================
c   subject:
c   --------
c   PROGRAM useful to run physical part of the martian GCM in a 1D column
c       
c Can be compiled with a command like (e.g. for 25 layers)
c  "makegcm -p mars -d 25 testphys1d"
c It requires the files "testphys1d.def" "callphys.def"
c      and a file describing the sigma layers (e.g. "z2sig.def")
c
c   author: Frederic Hourdin, R.Fournier,F.Forget
c   -------
c   
c   update: 12/06/2003 including chemistry (S. Lebonnois) 
c                            and water ice (F. Montmessin)
c 
c=======================================================================

#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "comgeomfi.h"
#include "surfdat.h"
#include "comdiurn.h"
#include "callkeys.h"
#include "comcstfi.h"
#include "planete.h"
#include "comsaison.h"
#include "yomaer.h"
#include "aerdust.h"
#include "control.h"
#include "comvert.h"
#include "netcdf.inc"
#include "comg1d.h"
#include "watercap.h"
#include "fisice.h"
#include "logic.h"

c --------------------------------------------------------------
c  Declarations
c --------------------------------------------------------------
c
      INTEGER unit           ! unite de lecture de "testphys1d.def"
      INTEGER unitstart      ! unite d'ecriture de "startfi"
      INTEGER nlayer,nlevel,nsoil,ndt
      INTEGER ilayer,ilevel,isoil,idt,iq
      LOGICAl firstcall,lastcall
c
      INTEGER day0          ! date initial (sol ; =0 a Ls=0)
      REAL day           ! date durant le run
      REAL time             ! time (0<time<1 ; time=0.5 a midi)
      REAL play(nlayermx)   ! Pressure at the middle of the layers (Pa)
      REAL plev(nlayermx+1) ! intermediate pressure levels (pa)
      REAL psurf,tsurf      
      REAL u(nlayermx),v(nlayermx)  ! zonal, meridional wind
      REAL gru,grv   ! prescribed "geostrophic" background wind
      REAL temp(nlayermx)   ! temperature at the middle of the layers
      REAL q(nlayermx,nqmx) ! tracer mixing ratio (e.g. kg/kg)
      REAL qsurf(nqmx)      ! tracer surface budget (e.g. kg.m-2)
      REAL tsoil(nsoilmx)   ! subsurface soik temperature (K)
      REAL co2ice           ! co2ice layer (kg.m-2)
      REAL emis             ! surface layer
      REAL q2(nlayermx+1)   ! Turbulent Kinetic Energy
      REAL zlay(nlayermx)   ! altitude estimee dans les couches (km)

c    Physical and dynamical tandencies (e.g.  m.s-2, K/s, Pa/s)
      REAL du(nlayermx),dv(nlayermx),dtemp(nlayermx)
      REAL dudyn(nlayermx),dvdyn(nlayermx),dtempdyn(nlayermx)
      REAL dpsurf    
      REAL dq(nlayermx,nqmx)
      REAL dqdyn(nlayermx,nqmx)

c   Various intermediate variables
      INTEGER thermo
      REAL zls
      REAL phi(nlayermx),h(nlayermx),s(nlayermx)
      REAL pks, ptif, w(nlayermx)
      REAL qtotinit, mqtot(nqmx),qtot
      INTEGER ierr, aslun
      REAL tmp1(0:nlayermx),tmp2(0:nlayermx)
      Logical  tracerdyn

      character*2 str2

c=======================================================================

c=======================================================================
c INITIALISATION
c=======================================================================

c ------------------------------------------------------
c  Constantes prescrites ICI
c ------------------------------------------------------

      pi=2.E+0*asin(1.E+0)

c     Constante de la Planete Mars
c     ----------------------------
      rad=3397200.               ! rayon de mars (m)  ~3397200 m
      daysec=88775.              ! duree du sol (s)  ~88775 s
      omeg=4.*asin(1.)/(daysec)  ! vitesse de rotation (rad.s-1)
      g=3.72                     ! gravite (m.s-2) ~3.72  
      mugaz=43.49                ! Masse molaire de l'atm (g.mol-1) ~43.49
      rcp=.256793                ! = r/cp  ~0.256793
      r= 8.314511E+0 *1000.E+0/mugaz
      cpp= r/rcp
      year_day = 669           ! duree de l'annee (sols) ~668.6
      periheli = 206.66          ! dist.min. soleil-mars (Mkm) ~206.66
      aphelie = 249.22           ! dist.max. soleil-mars (Mkm) ~249.22
      peri_day =  485.           ! date du perihelie (sols depuis printemps)
      obliquit = 25.2           ! Obliquite de la planete (deg) ~25.2         
 
c     Parametres Couche limite et Turbulence 
c     --------------------------------------
      z0 =  1.e-2                ! surface roughness (m) ~0.01 
      emin_turb = 1.e-6          ! energie minimale ~1.e-8
      lmixmin = 30               ! longueur de melange ~100
 
c     propriete optiques des calottes et emissivite du sol
c     ----------------------------------------------------
      emissiv= 0.95              ! Emissivite du sol martien ~.95
      emisice(1)=0.95            ! Emissivite calotte nord
      emisice(2)=0.95            ! Emissivite calotte sud  
      albedice(1)=0.5           ! Albedo calotte nord
      albedice(2)=0.5            ! Albedo calotte sud
      iceradius(1) = 100.e-6     ! mean scat radius of CO2 snow (north)
      iceradius(2) = 100.e-6     ! mean scat radius of CO2 snow (south)
      dtemisice(1) = 2.          ! time scale for snow metamorphism (north)
      dtemisice(2) = 2.          ! time scale for snow metamorphism (south
      hybrid=.false.

c ------------------------------------------------------
c  Lecture des parametres dans "testphys1d.def" 
c ------------------------------------------------------

c   Opening parameters file "testphys1d.def"
c   ---------------------------------------
      unit =97
      OPEN(unit,file='testphys1d.def',status='old',form='formatted'
     .     ,iostat=ierr)

      IF(ierr.ne.0) THEN
        write(*,*) 'Problem to open "testphys1d.def'
        write(*,*) 'Is it there ?'
        stop
      END IF

c  Date et heure locale du debut du run
c  ------------------------------------
c    Date (en sols depuis le solstice de printemps) du debut du run
      day0 = 0
      PRINT *,'date de depart ?'
      READ(unit,*) day0
      day=float(day0)
      PRINT *,day0
c  Heure de demarrage
      PRINT *,'heure de debut de simulation (entre 0 et 24) ?'
      READ(unit,*) time
      time=time/24.E+0

c  Discretisation (Definition de la grille et des pas de temps)
c  --------------
c
      nlayer=nlayermx
      nlevel=nlayer+1
      nsoil=nsoilmx
      PRINT *,'nombre de pas de temps par jour ?'
      READ(unit,*) day_step
      print*,day_step

      PRINT *,'nombre de jours simules ?'
      READ(unit,*) ndt
      print*,ndt

      ndt=ndt*day_step     
      dtphys=daysec/day_step  
c Pression de surface sur la planete
c ------------------------------------
c
      PRINT *,'pression au sol'
      READ(unit,*) psurf
      PRINT *,psurf
c Pression de reference
      pa=20.
      preff=610.      
 
c Proprietes des poussiere aerosol
c --------------------------------
      print *,'epaisseur optique dans le visible ?'
      READ(unit,*) tauvis
      print *,tauvis
 
c  latitude/longitude
c  ------------------
      PRINT *,'latitude en degres ?'
      READ(unit,*) lati(1)
      PRINT *,lati(1)
      lati(1)=lati(1)*pi/180.E+0
      long(1)=0.E+0
      long(1)=long(1)*pi/180.E+0

c  Initialisation albedo / inertie du sol
c  --------------------------------------
c
      PRINT *,'Albedo du sol nu ?'
      READ(unit,*) albedodat(1)
      PRINT *,albedodat(1)

      PRINT *,'Inertie thermique du sol ?'
      READ(unit,*) inertiedat(1)
      PRINT *,inertiedat(1)
c
c  pour le schema d'ondes de gravite
c  ---------------------------------
c
      zmea(1)=0.E+0
      zstd(1)=0.E+0
      zsig(1)=0.E+0
      zgam(1)=0.E+0
      zthe(1)=0.E+0



c   Initialisation speciales "physiq"
c   ---------------------------------
c   la surface de chaque maille est inutile en 1D --->
      area(1)=1.E+0

c   le geopotentiel au sol est inutile en 1D car tout est controle
c   par la pression de surface --->
      phisfi(1)=0.E+0

c  "inifis" reproduit un certain nombre d'initialisations deja faites
c  + lecture des clefs de callphys.def
c  + calcul de la frequence d'appel au rayonnement
c  + calcul des sinus et cosinus des longitude latitude

!Mars possible matter with dtphys in input and include!!!
      CALL inifis(1,llm,day0,daysec,dtphys,
     .            lati,long,area,rad,g,r,cpp)
c   Initialisation pour prendre en compte les vents en 1-D
c   ------------------------------------------------------
      ptif=2.E+0*omeg*sinlat(1)
 
c    vent geostrophique
      PRINT *,'composante vers l est du vent geostrophique (U) ?'
      READ(unit,*) gru
      PRINT *,'composante vers le nord du vent geostrophique (V) ?'
      READ(unit,*) grv

c     Initialisation des vents  au premier pas de temps
      DO ilayer=1,nlayer
         u(ilayer)=gru
         v(ilayer)=grv
      ENDDO

c     energie cinetique turbulente
      DO ilevel=1,nlevel
         q2(ilevel)=0.E+0
      ENDDO

c  glace de CO2 au sol
c  -------------------
      co2ice=0.E+0
      PRINT *,'co2ice (kg.m-2)'
      READ(unit,*) co2ice

c
c  emissivite
c  ----------
      emis=emissiv
      IF (co2ice.eq.1.E+0) THEN
         emis=emisice(1)
         IF(lati(1).LT.0) emis=emisice(2)
      ENDIF

 

c  calcul des pressions et altitudes en utilisant les niveaux sigma
c  ----------------------------------------------------------------

c    Vertical Coordinates
c    """"""""""""""""""""
      PRINT *,'Hybrid coordinates ?'
      READ(unit,*) hybrid
      PRINT *,'Hybrid =', hybrid

      CALL  disvert

      DO ilevel=1,nlevel
        plev(ilevel)=ap(ilevel)+psurf*bp(ilevel)
      ENDDO

      DO ilayer=1,nlayer
        play(ilayer)=aps(ilayer)+psurf*bps(ilayer)
      ENDDO

      DO ilayer=1,nlayer
        zlay(ilayer)=-200.E+0 *r*log(play(ilayer)/plev(1))
     &   /g
      ENDDO


c  profil de temperature au premier appel
c  --------------------------------------
      pks=psurf**rcp

c altitude en km dans profile: on divise zlay par 1000
      tmp1(0)=0.E+0
      DO ilayer=1,nlayer
        tmp1(ilayer)=zlay(ilayer)/1000.E+0
      ENDDO
      call profile(unit,nlayer+1,tmp1,tmp2)

      tsurf=tmp2(0)
      DO ilayer=1,nlayer
        temp(ilayer)=tmp2(ilayer)
      ENDDO
      


c     temperature du sous-sol
c     ~~~~~~~~~~~~~~~~~~~~~~~
      DO isoil=1,nsoil
         tsoil(isoil)=tsurf
      ENDDO

c    Initialisation des traceurs
c    ---------------------------

      DO iq=1,nqmx
        DO ilayer=1,nlayer
           q(ilayer,iq) = 0.
        ENDDO
      ENDDO

      if (photochem.or.callthermos) then
         write(*,*) 'Especes chimiques initialisees'
         ! thermo=0: initialisation pour toutes les couches 
         thermo=0
         call inichim_newstart(q,psurf,sig,nqmx,lati,long,area,
     $   thermo)
      endif
      watercaptag(ngridmx)=.false.
      
      DO iq=1,nqmx-1
        qsurf(iq) = 0.
      ENDDO


c    Initialisation pour les sorties GRADS dans "g1d.dat" et "g1d.ctl"
c    ----------------------------------------------------------------
c    (effectuee avec "writeg1d" a partir de "physiq.F"
c    ou tout autre subroutine physique, y compris celle ci).

        g1d_nlayer=nlayer
        g1d_nomfich='g1d.dat'
        g1d_unitfich=40
        g1d_nomctl='g1d.ctl'
        g1d_unitctl=41
        g1d_premier=.true.
        g2d_premier=.true.

c  Ecriture de "startfi"
c  --------------------
c  (Ce fichier sera aussitot relu au premier
c   appel de "physiq", mais il est necessaire pour passer
c   les variables purement physiques a "physiq"...

      call physdem1("startfi.nc",long,lati,nsoilmx,nqmx,
     .              dtphys,float(day0),
     .              time,tsurf,tsoil,co2ice,emis,q2,qsurf,
     .              area,albedodat,inertiedat,zmea,zstd,zsig,zgam,zthe)
c=======================================================================
c  BOUCLE TEMPORELLE DU MODELE 1D 
c=======================================================================
c
      firstcall=.true.
      lastcall=.false.

      DO idt=1,ndt
c        IF (idt.eq.ndt) lastcall=.true.
        IF (idt.eq.ndt-day_step-1) then       !test
         lastcall=.true.
         call solarlong(day*1.0,zls)
         write(103,*) 'Ls=',zls*180./pi
         write(103,*) 'Lat=', lati(1)*180./pi
         write(103,*) 'Tau=', tauvis/700*psurf
         write(103,*) 'RunEnd - Atmos. Temp. File'
         write(103,*) 'RunEnd - Atmos. Temp. File'
         write(104,*) 'Ls=',zls*180./pi
         write(104,*) 'Lat=', lati(1)
         write(104,*) 'Tau=', tauvis/700*psurf
         write(104,*) 'RunEnd - Atmos. Temp. File'
        ENDIF

c    calcul du geopotentiel 
c     ~~~~~~~~~~~~~~~~~~~~~
      DO ilayer=1,nlayer
        s(ilayer)=(aps(ilayer)/psurf+bps(ilayer))**rcp
        h(ilayer)=cpp*temp(ilayer)/(pks*s(ilayer))
      ENDDO
      phi(1)=pks*h(1)*(1.E+0-s(1))
      DO ilayer=2,nlayer
         phi(ilayer)=phi(ilayer-1)+
     &               pks*(h(ilayer-1)+h(ilayer))*.5E+0
     &                  *(s(ilayer-1)-s(ilayer))

      ENDDO

c       appel de la physique
c       --------------------

      CALL physiq (1,llm,nqmx,
     ,     firstcall,lastcall,
     ,     day,time,dtphys,
     ,     plev,play,phi,
     ,     u, v,temp, q,  
     ,     w,
C - sorties
     s     du, dv, dtemp, dq,dpsurf,tracerdyn)

c       evolution du vent : modele 1D
c       -----------------------------
 
c       la physique calcule les derivees temporelles de u et v.
c       on y rajoute betement un effet Coriolis.
c
c       DO ilayer=1,nlayer
c          du(ilayer)=du(ilayer)+ptif*(v(ilayer)-grv)
c          dv(ilayer)=dv(ilayer)+ptif*(-u(ilayer)+gru)
c       ENDDO

c       Pour certain test : pas de coriolis a l'equateur
c       if(lati(1).eq.0.) then
          DO ilayer=1,nlayer
             du(ilayer)=du(ilayer)+ (gru-u(ilayer))/1.e4
             dv(ilayer)=dv(ilayer)+ (grv-v(ilayer))/1.e4
          ENDDO
c       end if
c      
c
c       Calcul du temps au pas de temps suivant
c       ---------------------------------------
        firstcall=.false.
        time=time+dtphys/daysec
        IF (time.gt.1.E+0) then
            time=time-1.E+0
            day=day+1
        ENDIF

c       calcul des vitesses et temperature au pas de temps suivant
c       ----------------------------------------------------------

        DO ilayer=1,nlayer
           u(ilayer)=u(ilayer)+dtphys*du(ilayer)
           v(ilayer)=v(ilayer)+dtphys*dv(ilayer)
           temp(ilayer)=temp(ilayer)+dtphys*dtemp(ilayer)
        ENDDO

c       calcul des pressions au pas de temps suivant
c       ----------------------------------------------------------

           psurf=psurf+dtphys*dpsurf   ! evolution de la pression de surface
           DO ilevel=1,nlevel
             plev(ilevel)=ap(ilevel)+psurf*bp(ilevel)
           ENDDO
           DO ilayer=1,nlayer
             play(ilayer)=aps(ilayer)+psurf*bps(ilayer)
           ENDDO

      ENDDO   ! fin de la boucle temporelle

c    ========================================================
c    GESTION DES SORTIE
c    ========================================================

c    fermeture pour conclure les sorties format grads dans "g1d.dat"
c    et "g1d.ctl" (effectuee avec "writeg1d" a partir de "physiq.F"
c    ou tout autre subroutine physique, y compris celle ci).

c        CALL endg1d(1,nlayer,zphi/(g*1000.),ndt)
        CALL endg1d(1,nlayer,zlay/1000.,ndt)

c    ========================================================
      END
 
c***********************************************************************
c***********************************************************************
c     Subroutines Bidons utilise seulement en 3D, mais
c     necessaire a la compilation de testphys1d en 1-D

      subroutine gr_fi_dyn
      RETURN
      END
 
      subroutine iniwrite
      RETURN
      END
 
c***********************************************************************
c***********************************************************************

#include "../dyn3d/disvert.F"
