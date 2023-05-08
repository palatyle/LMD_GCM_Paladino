      PROGRAM rcm1d
      
      USE infotrac
      use control_mod, only: planet_type, day_step
      USE phys_state_var_mod
      use chemparam_mod
      USE comconst_mod, ONLY: cpp,t0_venus,nu_venus
      use cpdet_mod, only: ini_cpdet
      use moyzon_mod, only: tmoy
      USE comvert_mod, ONLY: ap,bp,presnivs,pa,preff,nivsigs,nivsig,
     .                       aps,bps,scaleheight,pseudoalt,
     .                       disvert_type,pressure_exner
      use conc, only: rho
      USE iniphysiq_mod, ONLY: iniphysiq
      USE mod_const_mpi, ONLY: comm_lmdz
      USE physiq_mod, ONLY: physiq
      USE logic_mod, ONLY: iflag_trac
      IMPLICIT NONE

c=======================================================================
c   subject:
c   --------
c   PROGRAM useful to run physical part of the venusian GCM in a 1D column
c       
c Can be compiled with a command like (e.g. for 50 layers)
c  "makelmdz -p venus -d 50 rcm1d"

c It requires the files "rcm1d.def" "physiq.def"
c      and a file describing the sigma layers (e.g. "z2sig.def")
c
c   author: Frederic Hourdin, R.Fournier,F.Forget (original Mars version)
c   ------- Sebastien Lebonnois (Venus version)
c   
c=======================================================================

#include "dimensions.h"
#include "dimsoil.h"
#include "comcstfi.h"
#include "netcdf.inc"
#include "clesphys.h"
#include "iniprint.h"
#include "tabcontrol.h"

c --------------------------------------------------------------
c  Declarations
c --------------------------------------------------------------
c
      INTEGER unit           ! unite de lecture de "rcm1d.def"
      INTEGER unitstart      ! unite d'ecriture de "startphy.nc"
      INTEGER nlayer,nlevel,nsoil,ndt
      INTEGER ilayer,ilevel,isoil,idt,iq
      LOGICAl firstcall,lastcall
c
      INTEGER day0          ! date initial (sol ; =0 a Ls=0)
      REAL day              ! date durant le run
      REAL time             ! time (0<time<1 ; time=0.5 a midi)
      REAL play(llm)   ! Pressure at the middle of the layers (Pa)
      REAL plev(llm+1) ! intermediate pressure levels (pa)
      REAL psurf      
      REAL u(llm),v(llm)  ! zonal, meridional wind
      REAL gru,grv   ! prescribed "geostrophic" background wind
      REAL temp(llm)   ! temperature at the middle of the layers
      REAL,allocatable :: q(:,:) ! tracer mixing ratio (e.g. kg/kg)
      REAL zlay(llm)   ! altitude estimee dans les couches (km)
      REAL long(1),lati(1),area(1)
      REAL cufi(1),cvfi(1)
      REAL phisfi(1)

c    Physical and dynamical tandencies (e.g.  m.s-2, K/s, Pa/s)
      REAL du(llm),dv(llm),dtemp(llm)
      REAL dudyn(llm),dvdyn(llm),dtempdyn(llm)
      REAL dpsurf(1)    
      REAL,allocatable :: dq(:,:)

c   Various intermediate variables
      REAL zls
      REAL phi(llm),s(llm)
      REAL pk(llm),pks, w(llm)
      INTEGER l, ierr, aslun
      REAL tmp1(0:llm),tmp2(0:llm)                        

      character*2 str2

      real pi

c=======================================================================

c=======================================================================
c INITIALISATION
c=======================================================================

      lunout = 6

c ------------------------------------------------------
c  Constantes prescrites ICI
c ------------------------------------------------------

      pi=2.E+0*asin(1.E+0)

c     Constante de la Planete Venus
c     -----------------------------
      planet_type = "venus"
      rad=6051300.               ! rayon de Venus (m)  ~6051300 m
      daysec= 1.0087e7           ! duree du sol (s)  ~1.e7 s
      omeg=4.*asin(1.)/19.4141e6 ! vitesse de rotation (rad.s-1)
      g= 8.87                    ! gravite (m.s-2) ~8.87
      mugaz=43.44                ! Masse molaire de l'atm (g.mol-1) ~43.44
! ADAPTATION GCM POUR CP(T)
! VENUS: Cp(T) = cpp*(T/T0)^nu 
! avec T0=460. et nu=0.35
      cpp=1.0e3
!     cpp=9.0e2      ! version constante
      r= 8.314511E+0 *1000.E+0/mugaz
      rcp= r/cpp

c-----------------------------------------------------------------------
c   Initialisation des traceurs
c   ---------------------------
c  Choix du nombre de traceurs et du schema pour l'advection
c  dans fichier traceur.def
      call infotrac_init
      if (nqtot.gt.1) iflag_trac=1

c Allocation de la tableau q : champs advectes   
      allocate(q(llm,nqtot))
      allocate(dq(llm,nqtot))

c ------------------------------------------------------
c  Lecture des parametres dans "rcm1d.def" 
c ------------------------------------------------------

c   Opening parameters file "rcm1d.def"
c   ---------------------------------------
      unit =97
      OPEN(unit,file='rcm1d.def',status='old',form='formatted'
     .     ,iostat=ierr)

      IF(ierr.ne.0) THEN
        write(*,*) 'Problem to open "rcm1d.def'
        write(*,*) 'Is it there ?'
        stop
      END IF

c  Date et heure locale du debut du run
c  ------------------------------------
c    Date (en sols depuis le solstice de printemps) du debut du run
      day0 = 0
      PRINT *,'date de depart ?'
      READ(unit,*) day0
      day=REAL(day0)
      PRINT *,day0
c  Heure de demarrage
      PRINT *,'heure de debut de simulation (entre 0 et 24) ?'
      READ(unit,*) time
      time=time/24.E+0

c  Discretisation (Definition de la grille et des pas de temps)
c  --------------
c
      nlayer=llm
      nlevel=nlayer+1
      nsoil=nsoilmx
      PRINT *,'nombre de pas de temps par jour ?'
      READ(unit,*) day_step
      print*,day_step

c     PRINT *,'nombre d appel au rayonnement par jour ?'
c     READ(unit,*) nbapp_rad
c     print*,nbapp_rad
c LU DANS PHYSIQ.DEF...
      nbapp_rad = 1000.

      PRINT *,'nombre de jours simules ?'
      READ(unit,*) ndt
      print*,ndt

      ndt=ndt*day_step     
      dtphys=daysec/day_step  
      dtime=dtphys

c Pression de surface sur la planete
c ------------------------------------
c
      PRINT *,'pression au sol'
      READ(unit,*) psurf
      PRINT *,psurf
c Pression de reference  ! voir dyn3d/etat0_venus
c     pa     =  5.e4 
      pa     =  1.e6 
      preff  = 9.2e6 ! 92 bars
c     preff  = psurf
 
c  latitude/longitude 
c  -------------------
      PRINT *,'latitude en degres ?'
      READ(unit,*) lati(1)
      PRINT *,lati(1)
      lati(1)=lati(1)*pi/180.  ! must be in radians.
      long(1)=0.E+0

c   Initialisation speciales "physiq"
c   ---------------------------------

!      CALL init_phys_lmdz(iim,jjm,llm,1,(/1/))

c   la surface de chaque maille est inutile en 1D --->
      area(1)=1.E+0
c de meme ?
      cufi(1)=1.E+0
      cvfi(1)=1.E+0

      call ini_cpdet

c Ehouarn: iniphysiq requires arrays related to (3D) dynamics grid,
c e.g. for cell boundaries, which are meaningless in 1D; so pad these 
c with '0.' when necessary
      CALL iniphysiq(1,1,llm,
     &            1,comm_lmdz,
     &            daysec,day0,dtphys,
     &            (/lati(1),0./),(/0./),
     &            (/0.,0./),(/long(1),0./),
     &            (/ (/area,0./),(/0.,0./) /),
     &            (/cufi,0.,0.,0./),
     &            (/cvfi,0./),
     &            rad,g,r,cpp,1)

c   le geopotentiel au sol est inutile en 1D car tout est controle
c   par la pression de surface --->
      phisfi(1)=0.E+0

c   Initialisation pour prendre en compte les vents en 1-D
c   ------------------------------------------------------
 
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

c  calcul des pressions et altitudes en utilisant les niveaux sigma
c  ----------------------------------------------------------------

c    Vertical Coordinates  (hybrids)
c    """"""""""""""""""""
      CALL  disvert_noterre
      
c     Calcul au milieu des couches : Vient de la version Mars
c     WARNING : le choix de placer le milieu des couches au niveau de
c     pression intermédiaire est arbitraire et pourrait etre modifié.
c     C'est fait de la meme facon dans disvert

      DO l = 1, llm
       aps(l) =  0.5 *( ap(l) +ap(l+1))
       bps(l) =  0.5 *( bp(l) +bp(l+1))
      ENDDO

      DO ilevel=1,nlevel
        plev(ilevel)=ap(ilevel)+psurf*bp(ilevel)
      ENDDO

      DO ilayer=1,nlayer
        play(ilayer)=aps(ilayer)+psurf*bps(ilayer)
        pk(ilayer)  =cpp*(play(ilayer)/preff)**rcp
c       write(120,*) ilayer,plev(ilayer),play(ilayer)
      ENDDO
c     write(120,*) nlevel,plev(nlevel)
c     stop
      
      pks=cpp*(psurf/preff)**rcp

c  init des variables pour phyredem
c  --------------------------------
      call phys_state_var_init

c  profil de temperature et altitude au premier appel
c  --------------------------------------------------

c modif par rapport a Mars: 
c   on envoie dz/T=-log(play/psurf)*r/g dans profile
      tmp1(0)=0.0
      tmp1(1)= -log(play(1)/psurf)*r/g
      DO ilayer=2,nlayer
        tmp1(ilayer)=-log(play(ilayer)/play(ilayer-1))*r/g
      ENDDO
      call profile(unit,nlayer+1,tmp1,tmp2)
      CLOSE(unit)

      print*,"               Pression        Altitude     Temperature"
      ilayer=1
      ftsol(1)=tmp2(0)
       temp(1)=tmp2(1)
       zlay(1)=tmp2(1)*tmp1(1)
      print*,"           0",ftsol(1)
      print*,ilayer,play(ilayer),zlay(ilayer),temp(ilayer)
      DO ilayer=2,nlayer
        temp(ilayer)=tmp2(ilayer)
        zlay(ilayer)=zlay(ilayer-1)+tmp2(ilayer)*tmp1(ilayer)
        print*,ilayer,play(ilayer),zlay(ilayer),temp(ilayer)
      ENDDO

      allocate(tmoy(llm))
      tmoy(:)=temp(:)
      
c     temperature du sous-sol
c     ~~~~~~~~~~~~~~~~~~~~~~~
      DO isoil=1,nsoil
         ftsoil(1,isoil)=ftsol(1)
      ENDDO

c    Initialisation des traceurs
c    ---------------------------

      DO iq=1,nqtot
        DO ilayer=1,nlayer
           q(ilayer,iq) = 0.
        ENDDO
      ENDDO

c FULL CHEMISTRY !! AJOUTER INIT AURELIEN...
C Faudrait lire les cles avant pour mettre ca en option....
c ou alors mettre ca dans physiq

c    Initialisation des parametres d'oro
c    -----------------------------------

      zmea(1) = 0.
      zstd(1) = 0.
      zsig(1) = 0.
      zgam(1) = 0.
      zthe(1) = 0.
      zpic(1) = 0.
      zval(1) = 0.

c  Initialisation albedo 
c  ----------------------

      falbe(1)=0.1

c  Ecriture de "startphy.nc"
c  -------------------------
c  (Ce fichier sera aussitot relu au premier
c   appel de "physiq", mais il est necessaire pour passer
c   les variables purement physiques a "physiq"...

      solsw(1)    = 0.
      sollw(1)    = 0.
      fder(1)     = 0.
      radsol(1)   = 0.
     
      radpas      = NINT(1.*day_step/nbapp_rad)
      soil_model  = .true.

      call phyredem("startphy.nc")

c  deallocation des variables phyredem
c  -----------------------------------
      call phys_state_var_end

c=======================================================================
c  BOUCLE TEMPORELLE DU MODELE 1D 
c=======================================================================
c
      firstcall=.true.
      lastcall=.false.

      DO idt=1,ndt
        IF (idt.eq.ndt) then
         lastcall=.true.
c toujours nulle dans le cas de Venus, pour l'instant...
         zls = 0.0
c        write(103,*) 'Ls=',zls*180./pi
c        write(103,*) 'Lat=', lati(1)
c        write(103,*) 'RunEnd - Atmos. Temp. File'
c        write(103,*) 'RunEnd - Atmos. Temp. File'
c        write(104,*) 'Ls=',zls*180./pi
c        write(104,*) 'Lat=', lati(1)
c        write(104,*) 'RunEnd - Atmos. Temp. File'
        ENDIF

c    calcul du geopotentiel 
c     ~~~~~~~~~~~~~~~~~~~~~
! ADAPTATION GCM POUR CP(T)
      DO ilayer=1,nlayer
        s(ilayer)=(play(ilayer)/psurf)**rcp
      ENDDO
      phi(1)=cpp*temp(1)*(1.E+0-s(1))
      DO ilayer=2,nlayer
         phi(ilayer)=phi(ilayer-1)+
     &     cpp*(temp(ilayer-1)/s(ilayer-1)+temp(ilayer)/s(ilayer))*0.5
     &        *(s(ilayer-1)-s(ilayer))

      ENDDO

c       appel de la physique
c       --------------------

      CALL physiq (1,llm,nqtot,
     ,     firstcall,lastcall,
     ,     day,time,dtphys,
     ,     plev,play,pk,phi,phisfi,
     ,     presnivs,
     ,     u,v,temp,q,
     ,     plev,temp, ! "planetary mean" plev and temperature
     ,     w,
C - sorties
     s     du,dv,dtemp,dq,dpsurf)

c     calcul de rho
       rho = 0.
c     print*,rho


c     print*,"DT APRES PHYSIQ=",day,time,dtime
c     print*,dtemp
c     print*,temp
c     print*," "
c     stop

c       evolution du vent : modele 1D
c       -----------------------------
 
c       la physique calcule les derivees temporelles de u et v.
c       Pas de coriolis 
          DO ilayer=1,nlayer
             du(ilayer)=du(ilayer)+ (gru-u(ilayer))/1.e4
             dv(ilayer)=dv(ilayer)+ (grv-v(ilayer))/1.e4
          ENDDO
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

c       calcul des traceurs au pas de temps suivant
c       -------------------------------------------
        if (iflag_trac.eq.1) then
         DO iq=1,nqtot
          DO ilayer=1,nlayer
           q(ilayer,iq)=q(ilayer,iq)+dq(ilayer,iq)*dtphys
          ENDDO
         ENDDO
        endif

c       calcul des pressions au pas de temps suivant
c       --------------------------------------------

           psurf=psurf+dtphys*dpsurf(1)   ! evolution de la pression de surface
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

        print*,"Temperature finale:"
        print*,temp
        
c stabilite
      DO ilayer=1,nlayer
        zlay(ilayer) = phi(ilayer)/g/1000.  !en km
      ENDDO
      DO ilayer=2,nlayer
        tmp1(ilayer) =
     .    (temp(ilayer)-temp(ilayer-1))/(zlay(ilayer)-zlay(ilayer-1)) 
     .   + 1000.*g/cpp
      ENDDO

      OPEN(11,file='profile.new')
      DO ilayer=1,nlayer
        write (11,*) zlay(ilayer),temp(ilayer),tmp1(ilayer)
      ENDDO

c    ========================================================
      END
 
c***********************************************************************
c***********************************************************************

!#include "../dyn3d_common/disvert_noterre.F"
!#include "../dyn3d/abort_gcm.F"

