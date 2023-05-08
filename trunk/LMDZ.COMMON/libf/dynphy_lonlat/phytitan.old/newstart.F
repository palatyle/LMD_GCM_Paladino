C======================================================================
      PROGRAM newstart
c=======================================================================
c
c
c   Auteur:   S. Lebonnois, 
c    a partir des newstart/start_archive/lect_start_archive martiens
c
c             Derniere modif : 02/09 (ecriture des q*)
c                              01/12 (inclusion dans svn dyn3d)
c
c   Objet:  Modify the grid for the initial state (LMD GCM VENUS/TITAN)
c   -----           (from file NetCDF start_archive.nc)
c
c
c=======================================================================

      use IOIPSL
      USE filtreg_mod
      USE startvar
      USE control_mod
      USE infotrac
      use cpdet_mod, only: ini_cpdet,t2tpot
      use exner_hyb_m, only: exner_hyb
      use exner_milieu_m, only: exner_milieu
      USE comconst_mod
      USE comvert_mod, ONLY: ap,bp,presnivs,pa,preff,nivsigs,nivsig,
     .                       aps,bps,scaleheight,pseudoalt,
     .                       disvert_type,pressure_exner
      USE serre_mod, ONLY: clon,clat,grossismx,grossismy,
     &                     dzoomx,dzoomy,taux,tauy
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0
      USE logic_mod, ONLY: iflag_trac,fxyhypb,ysinus
      USE temps_mod, ONLY: day_ref,annee_ref
      implicit none

#include "dimensions.h"
#include "paramet.h"
#include "comdissnew.h"
#include "comgeom2.h"
#include "description.h"
#include "dimsoil.h"
#include "netcdf.inc"

c-----------------------------------------------------------------------
c   Declarations
c-----------------------------------------------------------------------

c Variables pour fichier "ini"
c------------------------------------
      INTEGER   imold,jmold,lmold,nqold,ip1jmp1old
      INTEGER   length
      parameter (length = 100)
      real      tab_cntrl(2*length) 
      INTEGER isoil,iq,iqmax
      CHARACTER*2   str2

c Variable histoire 
c------------------
      REAL vcov(iip1,jjm,llm),ucov(iip1,jjp1,llm) ! vents covariants
      REAL teta(iip1,jjp1,llm),ps(iip1,jjp1)
      REAL phis(iip1,jjp1)                     ! geopotentiel au sol
      REAL masse(ip1jmp1,llm)                ! masse de l'atmosphere
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:):: q! champs advectes
      REAL tab_cntrl_dyn(length) ! tableau des parametres de start

c variable physique
c------------------
      integer    ngridmx
      parameter (ngridmx=(2+(jjm-1)*iim - 1/jjm))
      REAL tab_cntrl_fi(length) ! tableau des parametres de startfi
      real rlat(ngridmx),rlon(ngridmx)
      REAL tsurf(ngridmx),tsoil(ngridmx,nsoilmx)
      REAL albe(ngridmx),radsol(ngridmx),sollw(ngridmx)
      real solsw(ngridmx),dlw(ngridmx)
      REAL zmea(ngridmx), zstd(ngridmx)
      REAL zsig(ngridmx), zgam(ngridmx), zthe(ngridmx)
      REAL zpic(ngridmx), zval(ngridmx)
      real t_fi(ngridmx,llm)

c Variable nouvelle grille naturelle au point scalaire
c------------------------------------------------------
      real us(iip1,jjp1,llm),vs(iip1,jjp1,llm)
      REAL p3d(iip1,jjp1,llm+1)            ! pression aux interfaces
      REAL phisold_newgrid(iip1,jjp1)
      REAL T(iip1,jjp1,llm)
      real rlonS(iip1,jjp1),rlatS(iip1,jjp1)
      real tsurfS(iip1,jjp1),tsoilS(iip1,jjp1,nsoilmx)
      real albeS(ip1jmp1),radsolS(ip1jmp1),sollwS(ip1jmp1)
      real solswS(ip1jmp1),dlwS(ip1jmp1)
      real zmeaS(ip1jmp1),zstdS(ip1jmp1),zsigS(ip1jmp1)
      real zgamS(ip1jmp1),ztheS(ip1jmp1),zpicS(ip1jmp1)
      real zvalS(ip1jmp1)

      real ptotal

c Var intermediaires : vent naturel, mais pas coord scalaire
c-----------------------------------------------------------
      real vnat(iip1,jjm,llm),unat(iip1,jjp1,llm)

      REAL pks(iip1,jjp1)                      ! exner (f pour filtre)
      REAL pk(iip1,jjp1,llm)
      REAL pkf(iip1,jjp1,llm)
      REAL alpha(iip1,jjp1,llm),beta(iip1,jjp1,llm)


c Variable de l'ancienne grille 
c---------------------------------------------------------

      real, dimension(:),     allocatable :: rlonuold, rlatvold
      real, dimension(:),     allocatable :: rlonvold, rlatuold
      real, dimension(:),     allocatable :: nivsigsold,nivsigold
      real, dimension(:),     allocatable :: apold,bpold
      real, dimension(:),     allocatable :: presnivsold
      real, dimension(:,:,:), allocatable :: uold,vold,told
      real, dimension(:,:,:,:), allocatable :: qold
      real, dimension(:,:,:), allocatable :: tsoilold
      real, dimension(:,:),   allocatable :: psold,phisold
      real, dimension(:,:),   allocatable :: tsurfold
      real, dimension(:,:),   allocatable :: albeold,radsolold
      real, dimension(:,:),   allocatable :: sollwold,solswold
      real, dimension(:,:),   allocatable :: dlwold
      real, dimension(:,:),   allocatable :: zmeaold,zstdold,zsigold
      real, dimension(:,:),   allocatable :: zgamold,ztheold,zpicold
      real, dimension(:,:),   allocatable :: zvalold

      real ptotalold

c Variable intermediaires iutilise pour l'extrapolation verticale 
c----------------------------------------------------------------
      real, dimension(:,:,:), allocatable :: var,varp1 

c divers local
c-----------------

      integer ierr,nid,nvarid
      INTEGER ij, l,i,j
      character*80      fichnom      
      integer, dimension(4) :: start,counter
      REAL phisinverse(iip1,jjp1)  ! geopotentiel au sol avant inversion
      logical topoflag,albedoflag,razvitu,razvitv 
      real    albedo
      
c=======================================================================
c  INITIALISATIONS DIVERSES
c=======================================================================

c VENUS/TITAN

        iflag_trac = 1
c-----------------------------------------------------------------------
c   Initialisation des traceurs
c   ---------------------------
c  Choix du nombre de traceurs et du schema pour l'advection
c  dans fichier traceur.def, par default ou via INCA
      call infotrac_init

c Allocation de la tableau q : champs advectes   
      allocate(q(iip1,jjp1,llm,nqtot))

c-----------------------------------------------------------------------
c   Ouverture du fichier a modifier (start_archive.nc)
c-----------------------------------------------------------------------

        write(*,*) 'Creation d un etat initial a partir de'
        write(*,*) './start_archive.nc'
        write(*,*)
        fichnom = 'start_archive.nc'
        ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
        IF (ierr.NE.NF_NOERR) THEN
          write(6,*)' Pb d''ouverture du fichier ',fichnom
          write(6,*)' ierr = ', ierr
          CALL ABORT
        ENDIF
 
c-----------------------------------------------------------------------
c Lecture du tableau des parametres du run (pour la dynamique)
c-----------------------------------------------------------------------

        write(*,*) 'lecture tab_cntrl START_ARCHIVE'
c
        ierr = NF_INQ_VARID (nid, "controle", nvarid)
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
        ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
c
      write(*,*) 'Impression de tab_cntrl'
      do i=1,200
        write(*,*) i,tab_cntrl(i)
      enddo
      
c-----------------------------------------------------------------------
c		Initialisation des constantes
c-----------------------------------------------------------------------

      imold      = tab_cntrl(1)
      jmold      = tab_cntrl(2)
      lmold      = tab_cntrl(3)
      day_ref    = tab_cntrl(4)
      annee_ref  = tab_cntrl(5)
      rad        = tab_cntrl(6)
      omeg       = tab_cntrl(7)
      g          = tab_cntrl(8)
      cpp        = tab_cntrl(9)
      kappa      = tab_cntrl(10)
      daysec     = tab_cntrl(11)
      dtvr       = tab_cntrl(12)
      etot0      = tab_cntrl(13)
      ptot0      = tab_cntrl(14)
      ztot0      = tab_cntrl(15)
      stot0      = tab_cntrl(16)
      ang0       = tab_cntrl(17)
      pa         = tab_cntrl(18)
      preff      = tab_cntrl(19)
c
      clon       = tab_cntrl(20)
      clat       = tab_cntrl(21)
      grossismx  = tab_cntrl(22)
      grossismy  = tab_cntrl(23)

      IF ( tab_cntrl(24).EQ.1. )  THEN
        fxyhypb  = . TRUE .
        dzoomx   = tab_cntrl(25)
        dzoomy   = tab_cntrl(26)
        taux     = tab_cntrl(28)
        tauy     = tab_cntrl(29)
      ELSE
        fxyhypb = . FALSE .
        ysinus  = . FALSE .
        IF( tab_cntrl(27).EQ.1. ) ysinus = . TRUE.
      ENDIF

      ptotalold  = tab_cntrl(2*length)

      write(*,*) "Old dimensions:"
      write(*,*) "longitude: ",imold
      write(*,*) "latitude: ",jmold
      write(*,*) "altitude: ",lmold
      write(*,*) 

      ip1jmp1old = (imold+1-1/imold)*(jmold+1-1/jmold)
      
c dans run.def
      CALL getin('planet_type',planet_type)
      call ini_cpdet

c=======================================================================
c   CHANGEMENT DE CONSTANTES CONTENUES DANS tab_cntrl
c=======================================================================
c  changement de la resolution dans le fichier de controle
      tab_cntrl(1)  = REAL(iim)
      tab_cntrl(2)  = REAL(jjm)
      tab_cntrl(3)  = REAL(llm)

c--------------------------------
c We use a specific run.def to read new parameters that need to be changed
c--------------------------------
      
c Changement de cpp:
      CALL getin('cpp',cpp)
      kappa = (8.314511/43.44E-3)/cpp
      tab_cntrl(9)  = cpp
      tab_cntrl(10) = kappa

c CHANGING THE ZOOM PARAMETERS TO CHANGE THE GRID
      CALL getin('clon',clon)
      CALL getin('clat',clat)
      tab_cntrl(20) = clon
      tab_cntrl(21) = clat
      CALL getin('grossismx',grossismx)
      CALL getin('grossismy',grossismy)
      tab_cntrl(22) = grossismx
      tab_cntrl(23) = grossismy
      CALL getin('fxyhypb',fxyhypb)
      IF ( fxyhypb )  THEN
        CALL getin('dzoomx',dzoomx)
        CALL getin('dzoomy',dzoomy)
        tab_cntrl(25) = dzoomx
        tab_cntrl(26) = dzoomy
        CALL getin('taux',taux)
        CALL getin('tauy',tauy)
        tab_cntrl(28) = taux
        tab_cntrl(29) = tauy
      ELSE
        CALL getin('ysinus',ysinus)
        IF (ysinus) THEN
          tab_cntrl(27) = 1
        ELSE
          tab_cntrl(27) = 0
        ENDIF
      ENDIF

c changes must be done BEFORE these lines...
      DO l=1,length
         tab_cntrl_dyn(l)= tab_cntrl(l)
         tab_cntrl_fi(l) = tab_cntrl(l+length)
      ENDDO

c-----------------------------------------------------------------------
c	Autres initialisations 
c-----------------------------------------------------------------------

      CALL iniconst 
      CALL inigeom
      call inifilr

c-----------------------------------------------------------------------
c		Allocations des anciennes variables
c-----------------------------------------------------------------------

      allocate(rlonuold(imold+1), rlatvold(jmold))
      allocate(rlonvold(imold+1), rlatuold(jmold+1))
      allocate(nivsigsold(lmold+1),nivsigold(lmold))
      allocate(apold(lmold),bpold(lmold))
      allocate(presnivsold(lmold))
      allocate(uold(imold+1,jmold+1,lmold))
      allocate(vold(imold+1,jmold+1,lmold))
      allocate(told(imold+1,jmold+1,lmold))
      allocate(qold(imold+1,jmold+1,lmold,nqtot))
      allocate(psold(imold+1,jmold+1))
      allocate(phisold(imold+1,jmold+1))
      allocate(tsurfold(imold+1,jmold+1))
      allocate(tsoilold(imold+1,jmold+1,nsoilmx))
      allocate(albeold(imold+1,jmold+1),radsolold(imold+1,jmold+1))
      allocate(sollwold(imold+1,jmold+1),solswold(imold+1,jmold+1))
      allocate(dlwold(imold+1,jmold+1))
      allocate(zmeaold(imold+1,jmold+1),zstdold(imold+1,jmold+1))
      allocate(zsigold(imold+1,jmold+1),zgamold(imold+1,jmold+1))
      allocate(ztheold(imold+1,jmold+1),zpicold(imold+1,jmold+1))
      allocate(zvalold(imold+1,jmold+1))

      allocate(var (imold+1,jmold+1,llm))
      allocate(varp1 (imold+1,jmold+1,llm+1))

      print*,"Initialisations OK"

c-----------------------------------------------------------------------
c Lecture des longitudes et latitudes
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "rlonv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Le champ <rlonv> est absent de start.nc"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonvold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlonvold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Lecture echouee pour <rlonv>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlatu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Le champ <rlatu> est absent de start.nc"
         CALL abort
      ENDIF 
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatuold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlatuold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Lecture echouee pour <rlatu>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlonu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Le champ <rlonu> est absent de start.nc"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonuold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlonuold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Lecture echouee pour <rlonu>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlatv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Le champ <rlatv> est absent de start.nc"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatvold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlatvold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Lecture echouee pour <rlatv>"
         CALL abort
      ENDIF
c

      print*,"Lecture lon/lat OK"

c-----------------------------------------------------------------------
c Lecture des niveaux verticaux
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "nivsigs", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <nivsigs> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, nivsigsold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, nivsigsold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <nivsigs>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "nivsig", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <nivsig> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, nivsigold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, nivsigold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <nivsig>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "ap", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Le champ <ap> est absent de start.nc"
         CALL abort
      ELSE
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid, nvarid, apold)
#else
         ierr = NF_GET_VAR_REAL(nid, nvarid, apold)
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "new_grid: Lecture echouee pour <ap>"
         ENDIF
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "bp", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Le champ <bp> est absent de start.nc"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, bpold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, bpold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Lecture echouee pour <bp>"
         CALL abort
      END IF

      ierr = NF_INQ_VARID (nid, "presnivs", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <presnivs> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, presnivsold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, presnivsold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <presnivs>"
         CALL abort
      ENDIF

c-----------------------------------------------------------------------
c Lecture geopotentiel au sol
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "phisinit", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Le champ <phisinit> est absent de start.nc"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, phisold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, phisold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "new_grid: Lecture echouee pour <phisinit>"
         CALL abort
      ENDIF

      print*,"Lecture niveaux et geopot OK"

c-----------------------------------------------------------------------
c Lecture des champs 2D ()
c-----------------------------------------------------------------------

      start=(/1,1,1,0/)
      counter=(/imold+1,jmold+1,1,0/)

      ierr = NF_INQ_VARID (nid, "ps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <ps> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,psold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,psold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <ps>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "tsurf", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <tsurf> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,tsurfold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,tsurfold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <tsurf>"
         CALL abort
      ENDIF
c
      do isoil=1,nsoilmx
         write(str2,'(i2.2)') isoil
c
         ierr = NF_INQ_VARID (nid, "Tsoil"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .              Le champ <","Tsoil"//str2,"> est absent"
            CALL abort
         ENDIF
#ifdef NC_DOUBLE
         ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,
     .          tsoilold(1,1,isoil))
#else
         ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,
     .          tsoilold(1,1,isoil))
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .            Lecture echouee pour <","Tsoil"//str2,">"
            CALL abort
         ENDIF
      end do
c
      ierr = NF_INQ_VARID (nid, "albe", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <albe> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,albeold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,albeold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <albe>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "radsol", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <radsol> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,radsolold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,radsolold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <radsol>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "sollw", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <sollw> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,sollwold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,sollwold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <sollw>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "solsw", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <solsw> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,solswold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,solswold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <solsw>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "dlw", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <dlw> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,dlwold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,dlwold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <dlw>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "zmea", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <zmea> est absent"
         PRINT*, "          Il est donc initialise a zero"
         zmeaold=0.
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,zmeaold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,zmeaold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <zmea>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "zstd", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <zstd> est absent"
         PRINT*, "          Il est donc initialise a zero"
         zstdold=0.
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,zstdold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,zstdold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <zstd>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "zsig", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <zsig> est absent"
         PRINT*, "          Il est donc initialise a zero"
         zsigold=0.
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,zsigold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,zsigold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <zsig>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "zgam", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <zgam> est absent"
         PRINT*, "          Il est donc initialise a zero"
         zgamold=0.
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,zgamold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,zgamold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <zgam>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "zthe", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <zthe> est absent"
         PRINT*, "          Il est donc initialise a zero"
         ztheold=0.
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,ztheold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,ztheold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <zthe>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "zpic", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <zpic> est absent"
         PRINT*, "          Il est donc initialise a zero"
         zpicold=0.
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,zpicold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,zpicold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <zpic>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "zval", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <zval> est absent"
         PRINT*, "          Il est donc initialise a zero"
         zvalold=0.
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,zvalold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,zvalold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <zval>"
         CALL abort
      ENDIF
c

      print*,"Lecture champs 2D OK"

c-----------------------------------------------------------------------
c	Lecture des champs 3D ()
c-----------------------------------------------------------------------

      start=(/1,1,1,1/)
      counter=(/imold+1,jmold+1,lmold,1/)
c
      ierr = NF_INQ_VARID (nid,"temp", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <temp> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid, start, counter, told)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid, start, counter, told)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <temp>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid,"u", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <u> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,uold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,uold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <u>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid,"v", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <v> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,counter,vold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,vold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <v>"
         CALL abort
      ENDIF
c

c TNAME: IL EST LU A PARTIR DE traceur.def (mettre l'ancien si
c                changement du nombre de traceurs)

      IF(iflag_trac.eq.1) THEN
      DO iq=1,nqtot
        ierr =  NF_INQ_VARID (nid, tname(iq), nvarid)
        IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "dynetat0: Le champ <"//tname(iq)//"> est absent"
            PRINT*, "          Il est donc initialise a zero"
            qold(:,:,:,iq)=0.
        ELSE
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start,counter,qold(1,1,1,iq))
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,counter,qold(1,1,1,iq))
#endif
          IF (ierr .NE. NF_NOERR) THEN
             PRINT*, "dynetat0: Lecture echouee pour "//tname(iq)
             CALL abort
          ENDIF
        ENDIF
      ENDDO
      ENDIF


      print*,"Lecture champs 3D OK"

c=======================================================================
c   INTERPOLATION DANS LA NOUVELLE GRILLE et initialisation des variables
c=======================================================================
c  Interpolation horizontale puis passage dans la grille physique pour 
c  les variables physique 
c  Interpolation verticale puis horizontale pour chaque variable 3D
c=======================================================================

c-----------------------------------------------------------------------
c	Variable 2d :
c-----------------------------------------------------------------------

c Relief 
! topoflag = F: we keep the existing topography
!          = T: we read it from the Relief.nc file
! topoflag need to be in the specific run.def for newstart
       topoflag = . FALSE .
       CALL getin('topoflag',topoflag)

       IF ( topoflag ) then
         print*,"Topography (phis) read in file Relief.nc"
         print*,"For Venus, map directly inverted in Relief.nc"
         CALL startget_phys2d('surfgeo',iip1,jjp1,rlonv,rlatu,phis,0.0,
     .                jjm ,rlonu,rlatv,.true.)
c needed ?
          phis(iip1,:) = phis(1,:)

         CALL startget_phys1d('zmea',iip1,jjp1,rlonv,rlatu,ngridmx,zmea,
     .            0.0,jjm,rlonu,rlatv,.true.)
         CALL startget_phys1d('zstd',iip1,jjp1,rlonv,rlatu,ngridmx,zstd,
     .            0.0,jjm,rlonu,rlatv,.true.)
         CALL startget_phys1d('zsig',iip1,jjp1,rlonv,rlatu,ngridmx,zsig,
     .            0.0,jjm,rlonu,rlatv,.true.)
         CALL startget_phys1d('zgam',iip1,jjp1,rlonv,rlatu,ngridmx,zgam,
     .            0.0,jjm,rlonu,rlatv,.true.)
         CALL startget_phys1d('zthe',iip1,jjp1,rlonv,rlatu,ngridmx,zthe,
     .            0.0,jjm,rlonu,rlatv,.true.)
         CALL startget_phys1d('zpic',iip1,jjp1,rlonv,rlatu,ngridmx,zpic,
     .            0.0,jjm,rlonu,rlatv,.true.)
         CALL startget_phys1d('zval',iip1,jjp1,rlonv,rlatu,ngridmx,zval,
     .            0.0,jjm,rlonu,rlatv,.true.)

       ELSE
          print*,'Using existing topography'
          call interp_horiz (phisold,phis,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)

          call interp_horiz (zmeaold,zmeaS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
          call gr_dyn_fi (1,iip1,jjp1,ngridmx,zmeaS,zmea)
          call interp_horiz (zstdold,zstdS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
          call gr_dyn_fi (1,iip1,jjp1,ngridmx,zstdS,zstd)
          call interp_horiz (zsigold,zsigS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
          call gr_dyn_fi (1,iip1,jjp1,ngridmx,zsigS,zsig)
          call interp_horiz (zgamold,zgamS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
          call gr_dyn_fi (1,iip1,jjp1,ngridmx,zgamS,zgam)
          call interp_horiz (ztheold,ztheS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
          call gr_dyn_fi (1,iip1,jjp1,ngridmx,ztheS,zthe)
          call interp_horiz (zpicold,zpicS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
          call gr_dyn_fi (1,iip1,jjp1,ngridmx,zpicS,zpic)
          call interp_horiz (zvalold,zvalS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
          call gr_dyn_fi (1,iip1,jjp1,ngridmx,zvalS,zval)
       ENDIF

       print*,"New surface geopotential OK"

c Lat et lon pour physique
      do i=1,iip1
        rlatS(i,:)=rlatu(:)*180./pi
      enddo
      call gr_dyn_fi (1,iip1,jjp1,ngridmx,rlatS,rlat)

      do j=2,jjm
        rlonS(:,j)=rlonv(:)*180./pi
      enddo
      rlonS(:,1)=0.
      rlonS(:,jjp1)=0.
      call gr_dyn_fi (1,iip1,jjp1,ngridmx,rlonS,rlon)

c Temperature de surface
      call interp_horiz (tsurfold,tsurfS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iip1,jjp1,ngridmx,tsurfS,tsurf)
c     write(44,*) 'tsurf', tsurf

c Temperature du sous-sol
      call interp_horiz(tsoilold,tsoilS,
     &                  imold,jmold,iim,jjm,nsoilmx,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (nsoilmx,iip1,jjp1,ngridmx,tsoilS,tsoil)
c     write(45,*) 'tsoil',tsoil

! CHANGING ALBEDO: may be done through run.def
       albedoflag = . FALSE .
       CALL getin('albedoflag',albedoflag)

       IF ( albedoflag ) then
         print*,"Albedo is reinitialized to the albedo value in run.def"
         print*,"We may want to consider a map later on..."
         albedo=0.1
         CALL getin('albedo',albedo)
         albe=albedo
       ELSE
         call interp_horiz (albeold,albeS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
         call gr_dyn_fi (1,iip1,jjp1,ngridmx,albeS,albe)
       ENDIF

      call interp_horiz (radsolold,radsolS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iip1,jjp1,ngridmx,radsolS,radsol)

      call interp_horiz (sollwold,sollwS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iip1,jjp1,ngridmx,sollwS,sollw)

      call interp_horiz (solswold,solswS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iip1,jjp1,ngridmx,solswS,solsw)

      call interp_horiz (dlwold,dlwS,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iip1,jjp1,ngridmx,dlwS,dlw)

      print*,"Nouvelles var physiques OK"

c-----------------------------------------------------------------------
c	Traitement special de la pression au sol :
c-----------------------------------------------------------------------

c  Extrapolation la pression dans la nouvelle grille
      call interp_horiz(psold,ps,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)

c	On assure la conservation de la masse de l'atmosphere 
c--------------------------------------------------------------

!!! ATTENTION TEMPORAIRE
c     ps(:,:)=146700.

      ptotal =  0.
      DO j=1,jjp1
         DO i=1,iim
            ptotal=ptotal+ps(i,j)*aire(i,j)/g
         END DO
      END DO

      write(*,*)
      write(*,*)'Ancienne grille: masse de l atm :',ptotalold
      write(*,*)'Nouvelle grille: masse de l atm :',ptotal
      write (*,*) 'Ratio new atm./ old atm =', ptotal/ptotalold 
      write(*,*)


      DO j=1,jjp1
         DO i=1,iip1
            ps(i,j)=ps(i,j) * ptotalold/ptotal
         END DO
      END DO

c la pression de surface et les temperatures ne sont pas reequilibrees en fonction
c de la nouvelle topographie... 
c Si l'ajustement inevitable du debut pose des problemes, voir le newstart martien.

      print*,"Nouvelle ps OK"
      print*,"If changes done on topography, beware !"
      print*,"Some time may be needed for adjustments at the beginning"
      print*,"so if unstable, relax the timestep and/or dissipation."

c-----------------------------------------------------------------------
c	Variable 3d :
c-----------------------------------------------------------------------

      CALL pression(ip1jmp1, ap, bp, ps, p3d)
         if (disvert_type==1) then
           CALL exner_hyb(  ip1jmp1, ps, p3d, pks, pk, pkf )
         else ! we assume that we are in the disvert_type==2 case
           CALL exner_milieu( ip1jmp1, ps, p3d, pks, pk, pkf )
         endif
      
c temperatures atmospheriques

c ATTENTION: peut servir, mais bon...
c     do l=1,lmold
c      do j=1,jmold+1
c       do i=1,imold+1
c modif: profil uniforme
c          told(i,j,l)=told(1,jmold/2,l)
c mean T profile:
c       told(i,j,l) = 142.1*exp(-((p3d(i,j,l)/100.+21.45)/40.11)**2.) 
c    .              + 106.3*exp(-((p3d(i,j,l)/100.-3183.)/4737.)**2.)
c       enddo
c      enddo
c     enddo

      write (*,*) 'told ', told (1,jmold+1,1)  ! INFO
      call interp_vert
     &    (told,var,lmold,llm,apold,bpold,ap,bp,
     &     psold,(imold+1)*(jmold+1))
      write (*,*) 'var ', var (1,jmold+1,1)  ! INFO
      call interp_horiz(var,t,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'T ', T(1,jjp1,1)  ! INFO
! pour info:
! Si extension verticale, la T est extrapolee constante au-dessus de lmold

c passage grille physique pour restartphy.nc 
      call gr_dyn_fi (llm,iip1,jjp1,ngridmx,T,t_fi)

! ADAPTATION GCM POUR CP(T)
c passage en temperature potentielle
      call t2tpot(ip1jmp1*llm,T,teta,pk)
c on assure la periodicite
      teta(iip1,:,:) =  teta(1,:,:)

! RESETING U TO 0: may be done through run.def
       razvitu = . FALSE .
       CALL getin('razvitu',razvitu)
       razvitv = . FALSE .
       CALL getin('razvitv',razvitv)

c calcul des champ de vent; passage en vent covariant
      write (*,*) 'uold ', uold (1,2,1)  ! INFO
      call interp_vert
     & (uold,var,lmold,llm,apold,bpold,ap,bp,
     &  psold,(imold+1)*(jmold+1))
      write (*,*) 'var ', var (1,2,1)  ! INFO
      call interp_horiz(var,us,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'us ', us (1,2,1)   ! INFO

      call interp_vert
     & (vold,var,lmold,llm,
     &  apold,bpold,ap,bp,psold,(imold+1)*(jmold+1))
      call interp_horiz(var,vs,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call scal_wind(us,vs,unat,vnat)
! Reseting u=0
      if (razvitu) then
           unat(:,:,:) = 0.
      endif
      write (*,*) 'unat ', unat (1,2,1)    ! INFO
      do l=1,llm
        do j = 1, jjp1
          do i=1,iip1
            ucov( i,j,l ) = unat( i,j,l ) * cu(i,j)
! pour info:
! Si extension verticale, on impose u=0 au-dessus de lmold
            if (l.gt.lmold) ucov( i,j,l ) = 0
          end do
        end do
      end do 
      write (*,*) 'ucov ', ucov (1,2,1)  ! INFO
c     write(48,*) 'ucov',ucov
! Reseting v=0
      if (razvitv) then
           vnat(:,:,:) = 0.
      endif
      write (*,*) 'vnat ', vnat (1,2,1)    ! INFO
      do l=1,llm
        do j = 1, jjm
          do i=1,iim
            vcov( i,j,l ) = vnat( i,j,l ) * cv(i,j)
! pour info:
! Si extension verticale, on impose v=0 au-dessus de lmold
            if (l.gt.lmold) vcov( i,j,l ) = 0
          end do
          vcov( iip1,j,l ) = vcov( 1,j,l )
        end do
      end do
c     write(49,*) 'ucov',vcov

c masse: on la recalcule (ps a été ajustée pour conserver la masse totale)
      call massdair(p3d,masse)
      
c traceurs 3D
      do  iq = 1, nqtot
            call interp_vert(qold(1,1,1,iq),var,lmold,llm,
     &        apold,bpold,ap,bp,psold,(imold+1)*(jmold+1))
            call interp_horiz(var,q(1,1,1,iq),imold,jmold,iim,jjm,llm,
     &                  rlonuold,rlatvold,rlonu,rlatv)
      enddo

      print*,"Nouvelles var dynamiques OK"

c=======================================================================
c    Ecriture des restart.nc et restartphy.nc :
c=======================================================================

      call writerestart('restart.nc',tab_cntrl_dyn,
     .                    phis,vcov,ucov,teta,masse,q,ps)

      print*,"restart OK"

      call writerestartphy('restartphy.nc',tab_cntrl_fi,ngridmx,llm,
     .           rlat,rlon,tsurf,tsoil,albe,solsw, sollw,dlw,
     .           radsol,
     .           zmea,zstd,zsig,zgam,zthe,zpic,zval,
     .           t_fi)

      print*,"restartphy OK"

      end
