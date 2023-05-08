      SUBROUTINE lect_start_archive(ngrid,nlayer,nqtot,
     &     date,tsurf,tsoil,emis,q2,
     &     t,ucov,vcov,ps,co2ice,h,phisold_newgrid,
     &     q,qsurf,tauscaling,totcloudfrac,surfith,nid)
c=======================================================================
c
c
c   Auteur:    05/1997 , 12/2003 : coord hybride  FF
c   ------
c
c
c   Objet:     Lecture des variables d'un fichier "start_archive"
c              Plus besoin de régler ancienne valeurs grace
c              a l'allocation dynamique de memoire (Yann Wanherdrick)
c
c
c
c=======================================================================
      use infotrac, only: tname
      use comsoil_h, only: nsoilmx, layer, mlayer, volcapa, inertiedat
      use planete_h
      USE comvert_mod, ONLY: ap,bp,aps,bps,preff
      USE comconst_mod, ONLY: kappa,g,pi
      implicit none

#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
#include "netcdf.inc"
c=======================================================================
c   Declarations
c=======================================================================

! routine arguments
! -----------------
      integer,intent(in) :: ngrid ! number of atmosphferic columns
                                  ! on new physics grid
      integer,intent(in) :: nlayer ! number of atmospheric layers
                                   ! on new grid
      integer,intent(in) :: nqtot ! number of advected tracers
      REAL,INTENT(OUT) :: date
      REAL,INTENT(OUT) :: vcov(iip1,jjm,llm),ucov(iip1,jjp1,llm) ! vents covariants
      REAL,INTENT(OUT) :: h(iip1,jjp1,llm),ps(iip1,jjp1)
      REAL,INTENT(OUT) :: q(iip1,jjp1,llm,nqtot)
      REAL,INTENT(OUT) :: tsurf(ngrid) ! surface temperature
      REAL,INTENT(OUT) :: tsoil(ngrid,nsoilmx) ! soil temperature
      REAL,INTENT(OUT) :: co2ice(ngrid) ! CO2 ice layer
      REAL,INTENT(OUT) :: emis(ngrid)
      REAL,INTENT(OUT) :: q2(ngrid,nlayer+1),qsurf(ngrid,nqtot)
      REAL,INTENT(OUT) :: tauscaling(ngrid) ! dust conversion factor
      REAL,INTENT(OUT) :: totcloudfrac(ngrid) ! sub grid cloud fraction
      REAL,INTENT(OUT) :: phisold_newgrid(iip1,jjp1)
      REAL,INTENT(OUT) :: t(iip1,jjp1,llm)
      real,intent(in) :: surfith(iip1,jjp1) ! surface thermal inertia
      INTEGER,INTENT(IN) :: nid ! input NetCDF file identifier



c Old variables dimensions (from file)
c------------------------------------
      INTEGER   imold,jmold,lmold,nsoilold,nqold

c Variables pour les lectures des fichiers "ini" 
c--------------------------------------------------
!      INTEGER sizei,
      integer timelen,dimid
!      INTEGER length
!      parameter (length = 100)
      INTEGER tab0
      INTEGER isoil,iq,iqmax
      CHARACTER*2   str2

!      REAL dimfirst(4) ! tableau contenant les 1ers elements des dimensions

!      REAL dimlast(4) ! tableau contenant les derniers elements des dimensions

!      REAL dimcycl(4) ! tableau contenant les periodes des dimensions
!      CHARACTER*120 dimsource
!      CHARACTER*16 dimname
!      CHARACTER*80 dimtitle
!      CHARACTER*40 dimunits
!      INTEGER   dimtype

!      INTEGER dimord(4)  ! tableau contenant l''ordre
!      data dimord /1,2,3,4/ ! de sortie des dimensions

!      INTEGER vardim(4)
      INTEGER   memo
!      character (len=50) :: tmpname

c Variable histoire 
c------------------
      REAL ::qtot(iip1,jjp1,llm)

c autre variables dynamique nouvelle grille
c------------------------------------------

c!-*-
!      integer klatdat,klongdat
!      PARAMETER (klatdat=180,klongdat=360)

c Physique sur grille scalaire 
c----------------------------

c variable physique
c------------------
c     REAL phisfi(ngrid)

      INTEGER i,j,l
      INTEGER nvarid
c     REAL year_day,periheli,aphelie,peri_day
c     REAL obliquit,z0,emin_turb,lmixmin
c     REAL emissiv,emisice(2),albedice(2),tauvis
c     REAL iceradius(2) , dtemisice(2)

!      EXTERNAL RAN1
!      REAL RAN1
!      EXTERNAL geopot,inigeom
      integer ierr
!      integer ismin
!      external ismin
!      CHARACTER*80 datapath
      integer, dimension(4) :: start,count

c Variable nouvelle grille naturelle au point scalaire
c------------------------------------------------------
      real us(iip1,jjp1,llm),vs(iip1,jjp1,llm)
      real tsurfS(iip1,jjp1),tsoilS(iip1,jjp1,nsoilmx)
      real inertiedatS(iip1,jjp1,nsoilmx)
      real co2iceS(iip1,jjp1),emisS(iip1,jjp1)
      REAL q2S(iip1,jjp1,llm+1),qsurfS(iip1,jjp1,nqtot)
      real tauscalingS(iip1,jjp1)
      real totcloudfracS(iip1,jjp1)

      real ptotal, co2icetotal

c Var intermediaires : vent naturel, mais pas coord scalaire
c-----------------------------------------------------------
      real vnat(iip1,jjm,llm),unat(iip1,jjp1,llm)


c Variable de l'ancienne grille 
c---------------------------------------------------------

      real, dimension(:), allocatable :: timelist
      real, dimension(:), allocatable :: rlonuold, rlatvold
      real, dimension(:), allocatable :: rlonvold, rlatuold
      real, dimension(:), allocatable :: apsold,bpsold
      real, dimension(:), allocatable :: mlayerold
      real, dimension(:,:,:), allocatable :: uold,vold,told,q2old
      real, dimension(:,:,:), allocatable :: tsoilold,qsurfold
      real, dimension(:,:,:),allocatable :: tsoiloldnew
! tsoiloldnew: old soil values, but along new subterranean grid
      real, dimension(:,:,:), allocatable :: inertiedatold
! inertiedatoldnew: old inertia values, but along new subterranean grid
      real, dimension(:,:,:), allocatable :: inertiedatoldnew
      real, dimension(:,:), allocatable :: psold,phisold
      real, dimension(:,:), allocatable :: co2iceold,tsurfold
      real, dimension(:,:), allocatable :: emisold
      real, dimension(:,:,:,:), allocatable :: qold
      real, dimension(:,:), allocatable :: tauscalingold
      real, dimension(:,:), allocatable :: totcloudfracold

      real tab_cntrl(100)

      real ptotalold, co2icetotalold

      logical :: olddepthdef=.false. ! flag
! olddepthdef=.true. if soil depths are in 'old' (unspecified) format
      logical :: depthinterpol=.false. ! flag
! depthinterpol=.true. if interpolation will be requiered
      logical :: therminertia_3D=.true. ! flag
! therminertia_3D=.true. if thermal inertia is 3D and read from datafile
c Variable intermediaires iutilise pour l'extrapolation verticale 
c----------------------------------------------------------------
      real, dimension(:,:,:), allocatable :: var,varp1 
      real, dimension(:), allocatable :: oldgrid, oldval
      real, dimension(:), allocatable :: newval
!      real, dimension(:), allocatable :: oldmlayer

!      real surfithfi(ngrid)
      ! surface thermal inertia at old horizontal grid resolution
      real, dimension(:,:), allocatable :: surfithold 

! flag which identifies if archive file is using old tracer names (qsurf01,...)
      logical :: oldtracernames=.false.
      integer :: counter
      character(len=30) :: txt ! to store some text
      real :: tmpval ! to store a temporary variable/value

c=======================================================================

! 0. Preliminary stuff

! check if tracers follow old naming convention (q01, q02, q03, ...)
      counter=0
      do iq=1,nqtot
        txt= " "
        write(txt,'(a1,i2.2)')'q',iq
        ierr=NF_INQ_VARID(nid,txt,nvarid)
        if (ierr.ne.NF_NOERR) then
          ! did not find old tracer name
          exit ! might as well stop here
        else
          ! found old tracer name
          counter=counter+1
        endif
      enddo
      if (counter.eq.nqtot) then
        write(*,*) "lect_start_archive: tracers seem to follow old ",
     &             "naming convention (q01, q02,...)"
        oldtracernames=.true.
      endif


!-----------------------------------------------------------------------
! 1. Read data dimensions (i.e. size and length)
!-----------------------------------------------------------------------

! 1.2 Read the various dimension lengths of data in file 

      ierr= NF_INQ_DIMID(nid,"Time",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"temps",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,timelen)
      if (ierr.ne.NF_NOERR) then
        write(*,*) 'lect_start_archive error: cannot find Time length'
        stop
      else
        write(*,*) "lect_start_archive: timelen=",timelen
      endif

      ierr= NF_INQ_DIMID(nid,"latitude",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"rlatu",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,jmold)
      if (ierr.ne.NF_NOERR) then
        write(*,*) 'lect_start_archive error: cannot find lat length'
        stop
      else
        write(*,*) "lect_start_archive: jmold=",jmold
      endif
      jmold=jmold-1

      ierr= NF_INQ_DIMID(nid,"longitude",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"rlonv",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,imold)
      if (ierr.ne.NF_NOERR) then
        write(*,*) 'lect_start_archive error: cannot find lon length'
        stop
      else
        write(*,*) "lect_start_archive: imold=",imold
      endif
      imold=imold-1

      ierr= NF_INQ_DIMID(nid,"altitude",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"sig_s",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,lmold)
      if (ierr.ne.NF_NOERR) then
        write(*,*) 'lect_start_archive error: cannot find alt length'
        stop
      else
        write(*,*) "lect_start_archive: lmold=",lmold
      endif

      nqold=0
      do
         write(str2,'(i2.2)') nqold+1
         ierr= NF_INQ_VARID(nid,'q'//str2,dimid)
!        write(*,*) 'q'//str2
         if (ierr.eq.NF_NOERR) then
            nqold=nqold+1
         else
            exit
         endif
      enddo

! 1.2.1 find out the # of subsurface_layers
      nsoilold=0 !dummy initialisation
      ierr=NF_INQ_DIMID(nid,"subsurface_layers",dimid)
      if (ierr.eq.NF_NOERR) then
        ierr=NF_INQ_DIMLEN(nid,dimid,nsoilold)
	if (ierr.ne.NF_NOERR) then
	 write(*,*)'lec_start_archive: ',
     &              'Failed reading subsurface_layers length'
	endif
      else
        write(*,*)"lec_start_archive: did not find subsurface_layers"
      endif

      if (nsoilold.eq.0) then ! 'old' archive format;
      ! must use Tg//str2 fields to compute nsoilold
      write(*,*)"lec_start_archive: building nsoilold from Tg fields"
        do
	 write(str2,'(i2.2)') nsoilold+1
	 ierr=NF_INQ_VARID(nid,'Tg'//str2,dimid)
	 if (ierr.eq.NF_NOERR) then
	  nsoilold=nsoilold+1
	 else
	  exit
	 endif
	enddo
      endif


      if (nsoilold.ne.nsoilmx) then ! interpolation will be required
        depthinterpol=.true.
      endif

! 1.3 Report dimensions
      
      write(*,*) "lect_start_archive: Start_archive dimensions:"
      write(*,*) "longitude: ",imold
      write(*,*) "latitude: ",jmold
      write(*,*) "altitude: ",lmold
      if (nqold.gt.0) then
        write(*,*) "old tracers q*: ",nqold
      endif
      write(*,*) "subsurface_layers: ",nsoilold
      if (depthinterpol) then
      write(*,*) " => Warning, nsoilmx= ",nsoilmx
      write(*,*) '    which implies that you want subterranean interpola
     &tion.'
      write(*,*) '  Otherwise, set nsoilmx -in comsoil_h- to: ',nsoilold
      endif
      write(*,*) "time lenght: ",timelen
      write(*,*) 

!-----------------------------------------------------------------------
! 2. Allocate arrays to store datasets
!-----------------------------------------------------------------------

      allocate(timelist(timelen))
      allocate(rlonuold(imold+1), rlatvold(jmold))
      allocate(rlonvold(imold+1), rlatuold(jmold+1))
      allocate (apsold(lmold),bpsold(lmold))
      allocate(uold(imold+1,jmold+1,lmold))
      allocate(vold(imold+1,jmold+1,lmold))
      allocate(told(imold+1,jmold+1,lmold))
      allocate(psold(imold+1,jmold+1))
      allocate(phisold(imold+1,jmold+1))
      allocate(qold(imold+1,jmold+1,lmold,nqtot))
      allocate(co2iceold(imold+1,jmold+1))
      allocate(tsurfold(imold+1,jmold+1))
      allocate(emisold(imold+1,jmold+1))
      allocate(q2old(imold+1,jmold+1,lmold+1))
!      allocate(tsoilold(imold+1,jmold+1,nsoilmx))
      allocate(tsoilold(imold+1,jmold+1,nsoilold))
      allocate(tsoiloldnew(imold+1,jmold+1,nsoilmx))
      allocate(inertiedatold(imold+1,jmold+1,nsoilold)) ! soil thermal inertia
      allocate(inertiedatoldnew(imold+1,jmold+1,nsoilmx))
      ! surface thermal inertia at old horizontal grid resolution
      allocate(surfithold(imold+1,jmold+1))
      allocate(mlayerold(nsoilold))
      allocate(qsurfold(imold+1,jmold+1,nqtot))
      allocate(tauscalingold(imold+1,jmold+1))
      allocate(totcloudfracold(imold+1,jmold+1))

      allocate(var (imold+1,jmold+1,llm))
      allocate(varp1 (imold+1,jmold+1,llm+1))

      write(*,*) 'q2',ngrid,nlayer+1
      write(*,*) 'q2S',iip1,jjp1,llm+1
      write(*,*) 'q2old',imold+1,jmold+1,lmold+1

!-----------------------------------------------------------------------
! 3. Read time-independent data
!-----------------------------------------------------------------------

C-----------------------------------------------------------------------
c 3.1. Lecture du tableau des parametres du run 
c     (pour  la lecture ulterieure de "ptotalold" et "co2icetotalold")
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "controle", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Lect_start_archive: <controle> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <controle>"
         CALL abort
      ENDIF
c
      tab0 = 50

c-----------------------------------------------------------------------
c 3.2 Lecture des longitudes et latitudes
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "rlonv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <rlonv> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonvold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlonvold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <rlonv>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlatu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <rlatu> is missing"
         CALL abort
      ENDIF 
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatuold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlatuold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <rlatu>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlonu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <rlonu> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonuold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlonuold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <rlonu>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlatv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <rlatv> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatvold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlatvold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <rlatv>"
         CALL abort
      ENDIF
c

c-----------------------------------------------------------------------
c 3.3. Lecture des niveaux verticaux
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "aps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <aps> is missing"
         apsold=0
         PRINT*, "<aps> set to 0"
      ELSE
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid, nvarid, apsold)
#else
         ierr = NF_GET_VAR_REAL(nid, nvarid, apsold)
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: Failed loading <aps>"
         ENDIF
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "bps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <bps> is missing"
         PRINT*, "It must be an old start_archive, lets look for sig_s"
         ierr = NF_INQ_VARID (nid, "sig_s", nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "Nothing to do..."
            CALL abort
         ENDIF
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, bpsold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, bpsold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <bps>"
         CALL abort
      END IF

c-----------------------------------------------------------------------
c 3.4 Read Soil layers depths
c-----------------------------------------------------------------------
     
      ierr=NF_INQ_VARID(nid,"soildepth",nvarid)
      if (ierr.ne.NF_NOERR) then
       write(*,*)'lect_start_archive: Could not find <soildepth>'
       write(*,*)' => Assuming this is an archive in old format'
       olddepthdef=.true.
       depthinterpol=.true.
       ! this is how soil depth was defined in ye old days
	do isoil=1,nsoilold
	  mlayerold(isoil)=sqrt(887.75/3.14)*((2.**(isoil-0.5))-1.)
	enddo
      else
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid,nvarid,mlayerold)
#else
        ierr = NF_GET_VAR_REAL(nid,nvarid,mlayerold)
#endif
       if (ierr .NE. NF_NOERR) then
         PRINT*, "lect_start_archive: Failed reading <soildepth>"
         CALL abort
       endif

      endif !of if(ierr.ne.NF_NOERR)

      ! Read (or build) mlayer()
      if (depthinterpol) then
       ! Build (default) new soil depths (mlayer(:) is in comsoil.h),
       ! as in soil_settings.F
       write(*,*)' => Building default soil depths'
       do isoil=0,nsoilmx-1
         mlayer(isoil)=2.e-4*(2.**(isoil-0.5))
       enddo
       write(*,*)' => mlayer: ',mlayer
       ! Also build (default) new soil interlayer depth layer(:)
       do isoil=1,nsoilmx
         layer(isoil)=sqrt(mlayer(0)*mlayer(1))*
     &                      ((mlayer(1)/mlayer(0))**(isoil-1))
       enddo
       write(*,*)' =>  layer: ',layer
      else ! read mlayer() from file
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid,nvarid,mlayer)
#else
        ierr = NF_GET_VAR_REAL(nid,nvarid,mlayer)
#endif
       if (ierr .NE. NF_NOERR) then
         PRINT*, "lect_start_archive: Failed reading <soildepth>"
         CALL abort
       endif
       ! Also build (default) soil interlayer depth layer(:)
       do isoil=1,nsoilmx
         layer(isoil)=sqrt(mlayer(0)*mlayer(1))*
     &                      ((mlayer(1)/mlayer(0))**(isoil-1))
       enddo
      endif ! of if (depthinterpol)

c-----------------------------------------------------------------------
c 3.5 Read Soil thermal inertia
c-----------------------------------------------------------------------

      ierr=NF_INQ_VARID(nid,"inertiedat",nvarid)
      if (ierr.ne.NF_NOERR) then
       write(*,*)'lect_start_archive: Could not find <inertiedat>'
       write(*,*)' => Assuming this is an archive in old format'
       therminertia_3D=.false.
       write(*,*)' => Thermal inertia will be read from reference file'
       volcapa=1.e6
       write(*,*)'    and soil volumetric heat capacity is set to ',
     &           volcapa
      else
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid,nvarid,inertiedatold)
#else
        ierr = NF_GET_VAR_REAL(nid,nvarid,inertiedatold)
#endif
       if (ierr .NE. NF_NOERR) then
         PRINT*, "lect_start_archive: Failed reading <inertiedat>"
         CALL abort
       endif
      endif

c-----------------------------------------------------------------------
c 3.6 Lecture geopotentiel au sol
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "phisinit", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <phisinit> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, phisold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, phisold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <phisinit>"
         CALL abort
      ENDIF

C-----------------------------------------------------------------------
c   lecture de "ptotalold" et "co2icetotalold"
c-----------------------------------------------------------------------
      ptotalold = tab_cntrl(tab0+49)
      co2icetotalold = tab_cntrl(tab0+50)
 
c-----------------------------------------------------------------------
c 4. Lecture du temps et choix
c-----------------------------------------------------------------------
 
c  lecture du temps
c
      ierr = NF_INQ_DIMID (nid, "Time", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         ierr = NF_INQ_DIMID (nid, "temps", nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: <Time> is missing"
            CALL abort
         endif
      ENDIF

      ierr = NF_INQ_VARID (nid, "Time", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         ierr = NF_INQ_VARID (nid, "temps", nvarid)
      endif 
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, timelist)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, timelist)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <Time>"
         CALL abort
      ENDIF
c
      write(*,*)
      write(*,*)
      write(*,*) 'Dates of the stored initial states:'
      write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      pi=2.*ASIN(1.)
      do i=1,timelen
c       call solarlong(timelist(i),sollong(i))
c       sollong(i) = sollong(i)*180./pi
c        write(*,*) 'initial state at martian day: ',int(timelist(i))
        write(*,*) 'initial state at martian day: ',timelist(i)
c       write(*,6) nint(timelist(i)),nint(mod(timelist(i),669)),
c    .    sollong(i)
      end do

   6  FORMAT(i7,i7,f9.3)
 
      write(*,*)
      write(*,*) 'Choose the martian day to use'
 123  read(*,*,iostat=ierr) date
      if(ierr.ne.0) goto 123
      memo = 0
      do i=1,timelen
c        if (date.eq.int(timelist(i))) then
        if (abs(date-timelist(i)).lt.0.01) then
            memo = i
        endif
      end do
 
      if (memo.eq.0) then
        write(*,*)
        write(*,*)
        write(*,*) 'Wrong value for day number !!'
        write(*,*)
        write(*,*) 'Dates of the stored initial states:'
        write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        do i=1,timelen
          write(*,*) 'initial state at martian day: ',timelist(i)
c         write(*,6) nint(timelist(i)),nint(mod(timelist(i),669))
        end do
        goto 123
      endif
      
!-----------------------------------------------------------------------
! 5. Read (time-dependent) data from datafile
!-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c 5.1 Lecture des champs 2D (co2ice, emis,ps,tsurf,Tg[10], q2surf, tauscaling,totcloudfrac)
c-----------------------------------------------------------------------
 
      start=(/1,1,memo,0/)
      count=(/imold+1,jmold+1,1,0/)
       
      ierr = NF_INQ_VARID (nid, "co2ice", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <co2ice> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,co2iceold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,co2iceold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <co2ice>"
         PRINT*, NF_STRERROR(ierr)
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "emis", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <emis> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,emisold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,emisold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <emis>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "ps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <ps> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,psold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,psold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <ps>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "tsurf", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <tsurf> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,tsurfold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,tsurfold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <tsurf>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "q2surf", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <q2surf> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,q2old)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,q2old)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <q2surf>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "tauscaling", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <tauscaling> not in file"
         tauscalingold(:,:) = -1
      ELSE
#ifdef NC_DOUBLE
        ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,tauscalingold)
#else
        ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,tauscalingold)
#endif
        IF (ierr .NE. NF_NOERR) THEN
           PRINT*, "lect_start_archive: Failed loading <tauscaling>"
           PRINT*, NF_STRERROR(ierr)
           CALL abort
        ENDIF
      ENDIF

      ierr = NF_INQ_VARID (nid, "totcloudfrac", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <totcloudfrac> not in file"
         totcloudfracold(:,:) = 1
      ELSE
#ifdef NC_DOUBLE
        ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,
     &                            totcloudfracold)
#else
        ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,
     &                          totcloudfracold)
#endif
        IF (ierr .NE. NF_NOERR) THEN
           PRINT*, "lect_start_archive: Failed loading <totcloudfrac>"
           PRINT*, NF_STRERROR(ierr)
           CALL abort
        ENDIF
      ENDIF
c
      write(*,*)"lect_start_archive: rlonuold:"
     &           ,rlonuold," rlatvold:",rlatvold
      write(*,*)
c

c tracers: the 2 last ones are kept the 2 last one. 
c the others keep their rank. ! No longer true.
c -------------------------------------------
! Surface tracers:      
      qsurfold(1:imold+1,1:jmold+1,1:nqtot)=0

      DO iq=1,nqtot
        IF (oldtracernames) THEN
          txt=" "
          write(txt,'(a5,i2.2)')'qsurf',iq
        ELSE
          txt=trim(tname(iq))//"_surf"
          if (txt.eq."h2o_vap") then
            ! There is no surface tracer for h2o_vap;
            ! "h2o_ice" should be loaded instead
            txt="h2o_ice_surf"
            write(*,*) 'lect_start_archive: loading surface tracer',
     &                     ' h2o_ice instead of h2o_vap'
          endif
        ENDIF ! of IF (oldtracernames)
        write(*,*) "lect_start_archive: loading tracer ",trim(txt)
        ierr = NF_INQ_VARID (nid,txt,nvarid)
        IF (ierr .NE. NF_NOERR) THEN
          PRINT*, "lect_start_archive: ",
     &              " Tracer <",trim(txt),"> not found"
          print*, "which (constant) value should it be initialized to?"
          read(*,*) tmpval
          qsurfold(1:imold+1,1:jmold+1,iq)=tmpval
        ELSE ! tracer exists in file, load it
#ifdef NC_DOUBLE
          ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,
     &          qsurfold(1,1,iq))
#else
          ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,
     &          qsurfold(1,1,iq))
#endif
          IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: ",
     &             " Failed loading <",trim(txt),">"
            stop
          ENDIF
        ENDIF

      ENDDO ! of DO iq=1,nqtot

!-----------------------------------------------------------------------
! 5.2 Read 3D subterranean fields
!-----------------------------------------------------------------------

      start=(/1,1,1,memo/)
      count=(/imold+1,jmold+1,nsoilold,1/)
!
! Read soil temperatures
!
      if (olddepthdef) then ! tsoil stored using the 'old format'
         start=(/1,1,memo,0/)
         count=(/imold+1,jmold+1,1,0/) ! because the "Tg" are 2D datasets
       do isoil=1,nsoilold
!        write(*,*)'isoil',isoil
         write(str2,'(i2.2)') isoil
c
         ierr = NF_INQ_VARID (nid, "Tg"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: ",
     &              "Field <","Tg"//str2,"> not found"
            CALL abort
         ENDIF
#ifdef NC_DOUBLE
         ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,
     &          tsoilold(1,1,isoil))
#else
         ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,
     &          tsoilold(1,1,isoil))
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: ",
     &            "Failed reading <","Tg"//str2,">"
            CALL abort
         ENDIF
c
       enddo ! of do isoil=1,nsoilold
      
      ! reset 'start' and 'count' to "3D" behaviour
      start=(/1,1,1,memo/)
      count=(/imold+1,jmold+1,nsoilold,1/)
      
      else
       write(*,*) "lect_start_archive: loading tsoil "
       ierr=NF_INQ_VARID(nid,"tsoil",nvarid)
       if (ierr.ne.NF_NOERR) then
        write(*,*)"lect_start_archive: Cannot find <tsoil>"
	call abort
       else
#ifdef NC_DOUBLE
      ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start,count,tsoilold)
#else
      ierr=NF_GET_VARA_REAL(nid,nvarid,start,count,tsoilold)
#endif
       endif ! of if (ierr.ne.NF_NOERR)
       
      endif ! of if (olddepthdef)

!
! Read soil thermal inertias
!
!      if (.not.olddepthdef) then ! no thermal inertia data in "old" archives
!       ierr=NF_INQ_VARID(nid,"inertiedat",nvarid)
!       if (ierr.ne.NF_NOERR) then
!        write(*,*)"lect_start_archive: Cannot find <inertiedat>"
!	call abort
!       else
!#ifdef NC_DOUBLE
!      ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start,count,inertiedatold)
!#else
!      ierr=NF_GET_VARA_REAL(nid,nvarid,start,count,inertiedatold)
!#endif
!       endif ! of if (ierr.ne.NF_NOERR)
!      endif

c-----------------------------------------------------------------------
c 5.3	Lecture des champs 3D (t,u,v, q2atm,q)
c-----------------------------------------------------------------------

      start=(/1,1,1,memo/)
      count=(/imold+1,jmold+1,lmold,1/)

c
      ierr = NF_INQ_VARID (nid,"temp", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <temp> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid, start, count, told)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid, start, count, told)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <temp>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid,"u", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <u> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,uold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,uold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <u>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid,"v", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <v> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,vold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,vold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <v>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid,"q2atm", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: <q2atm> is missing"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,q2old(1,1,2))
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,q2old(1,1,2))
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Failed loading <q2atm>"
         CALL abort
      ENDIF
c

c tracers: the 2 last ones are kept the 2 last one. 
c the others keep their rank. ! No longer true.
c -------------------------------------------
! Tracers:
      qold(1:imold+1,1:jmold+1,1:lmold,1:nqtot)=0

      DO iq=1,nqtot
        IF (oldtracernames) THEN
          txt=" "
          write(txt,'(a1,i2.2)')'q',iq
        ELSE
          txt=tname(iq)
        ENDIF
        write(*,*)"lect_start_archive: loading tracer ",trim(txt)
        ierr = NF_INQ_VARID (nid,txt,nvarid)
        IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: ",
     &              " Tracer <",trim(txt),"> not found"
          print*, "which (constant) value should it be initialized to?"
          read(*,*) tmpval
          qold(1:imold+1,1:jmold+1,1:lmold,iq)=tmpval
        ELSE ! tracer exists in file, load it
#ifdef NC_DOUBLE
         ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start,count,qold(1,1,1,iq))
#else
         ierr=NF_GET_VARA_REAL(nid,nvarid,start,count,qold(1,1,1,iq))
#endif
          IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: ",
     &             "  Failed loading <",trim(txt),">"
            stop
          ENDIF
        ENDIF

      ENDDO ! of DO iq=1,nqtot


c Chemin pour trouver les donnees de surface (albedo, relief, th.inertia...)
c -------------------------------------------------------------------------

!      datapath = '/users/forget/gcm/data_mars_gcm'


!=======================================================================
! 6. Interpolation from old grid to new grid
!=======================================================================

c=======================================================================
c   INTERPOLATION DANS LA NOUVELLE GRILLE et initialisation des variables
c=======================================================================
c  Interpolation horizontale puis passage dans la grille physique pour 
c  les variables physique 
c  Interpolation verticale puis horizontale pour chaque variable 3D
c=======================================================================

c-----------------------------------------------------------------------
c 6.1	Variable 2d :
c-----------------------------------------------------------------------
c Relief 
      call interp_horiz (phisold,phisold_newgrid,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)

c Glace CO2
      call interp_horiz (co2iceold,co2ices,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)

c Temperature de surface
      call interp_horiz (tsurfold,tsurfs,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,tsurfs,tsurf)
c     write(44,*) 'tsurf', tsurf

c Temperature du sous-sol
!      call interp_horiz(tsoilold,tsoils,
!     &                  imold,jmold,iim,jjm,nsoilmx,
!     &                   rlonuold,rlatvold,rlonu,rlatv)
!      call gr_dyn_fi (nsoilmx,iim+1,jjm+1,ngrid,tsoils,tsoil)
c     write(45,*) 'tsoil',tsoil

c Emissivite de la surface
      call interp_horiz (emisold,emiss,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,emiss,emis)

c Dust conversion factor
      call interp_horiz (tauscalingold,tauscalings,imold,jmold,iim,jjm,
     &                   1,rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,tauscalings,tauscaling)
      
c Sub grid cloud fraction
      call interp_horiz (totcloudfracold,totcloudfracs,imold,jmold,iim,
     &                   jjm,1,rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,totcloudfracs,totcloudfrac)
      
c     write(46,*) 'emis',emis
c-----------------------------------------------------------------------
c 6.1.2	Traitement special de la pression au sol :
c-----------------------------------------------------------------------

c  Extrapolation la pression dans la nouvelle grille
      call interp_horiz(psold,ps,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)

c-----------------------------------------------------------------------
c	On assure la conservation de la masse de l'atmosphere + calottes
c-----------------------------------------------------------------------

      ptotal =  0.
      co2icetotal = 0.
      DO j=1,jjp1
         DO i=1,iim
            ptotal=ptotal+ps(i,j)*aire(i,j)/g
            co2icetotal = co2icetotal + co2iceS(i,j)*aire(i,j)
         END DO
      END DO

      write(*,*)
      write(*,*)'Old grid: mass of the atmosphere :',ptotalold
      write(*,*)'New grid: mass of the atmosphere :',ptotal
      write (*,*) 'Ratio new atm / old atm =', ptotal/ptotalold 
      write(*,*)
      write(*,*)'Old grid: mass of CO2 ice:',co2icetotalold
      write(*,*)'New grid: mass of CO2 ice:',co2icetotal
      if (co2icetotalold.ne.0.) then
       write(*,*)'Ratio new ice / old ice =',co2icetotal/co2icetotalold
      endif
      write(*,*)


      DO j=1,jjp1
         DO i=1,iip1
            ps(i,j)=ps(i,j) * ptotalold/ptotal
         END DO
      END DO

      if ( co2icetotalold.gt.0.) then 
         DO j=1,jjp1
            DO i=1,iip1
               co2iceS(i,j)=co2iceS(i,j) * co2icetotalold/co2icetotal
            END DO
         END DO
      end if

c-----------------------------------------------------------------------
c 6.2 Subterranean 3d variables:
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c 6.2.1 Thermal Inertia
c       Note: recall that inertiedat is a common in "comsoil.h"
c-----------------------------------------------------------------------

      ! depth-wise interpolation, if required
      if (depthinterpol.and.(.not.olddepthdef)) then
        allocate(oldval(nsoilold))
	allocate(newval(nsoilmx))
        write(*,*)'lect_start_archive: WARNING: vertical interpolation o
     &f soil thermal inertia; might be wiser to reset it.'
        write(*,*)
       
        do i=1,imold+1
         do j=1,jmold+1
	   !copy old values
	   oldval(1:nsoilold)=inertiedatold(i,j,1:nsoilold)
	   !interpolate
	   call interp_line(mlayerold,oldval,nsoilold,
     &                     mlayer,newval,nsoilmx)
           !copy interpolated values
           inertiedatoldnew(i,j,1:nsoilmx)=newval(1:nsoilmx)
	 enddo
        enddo
        ! cleanup
	deallocate(oldval)
	deallocate(newval)
      endif !of if (depthinterpol)

      if (therminertia_3D) then
        ! We have inertiedatold
       if((imold.ne.iim).or.(jmold.ne.jjm)) then
       write(*,*)'lect_start_archive: WARNING: horizontal interpolation 
     &of thermal inertia; might be better to reset it.'
       write(*,*)
       endif
       
        ! Do horizontal interpolation
	if (depthinterpol) then
	  call interp_horiz(inertiedatoldnew,inertiedatS,
     &                  imold,jmold,iim,jjm,nsoilmx,
     &                   rlonuold,rlatvold,rlonu,rlatv)
	else
          call interp_horiz(inertiedatold,inertiedatS,
     &                  imold,jmold,iim,jjm,nsoilold,
     &                   rlonuold,rlatvold,rlonu,rlatv)
        endif ! of if (depthinterpol)

      else ! no 3D thermal inertia data
       write(*,*)'lect_start_archive: using reference surface inertia'
        ! Use surface inertia (and extend it to all depths)
        do i=1,nsoilmx
         inertiedatS(1:iip1,1:jjp1,i)=surfith(1:iip1,1:jjp1)
        enddo
	! Build an old resolution surface thermal inertia
	! (will be needed for tsoil interpolation)
	call interp_horiz(surfith,surfithold,
     &                    iim,jjm,imold,jmold,1,
     &                    rlonu,rlatv,rlonuold,rlatvold)
      endif


      ! Reshape inertiedatS to scalar grid as inertiedat
      call gr_dyn_fi (nsoilmx,iim+1,jjm+1,ngrid,
     &                  inertiedatS,inertiedat)
      
c-----------------------------------------------------------------------
c 6.2.2 Soil temperature
c-----------------------------------------------------------------------
!      write(*,*) 'Soil'
      ! Recast temperatures along soil depth, if necessary
      if (olddepthdef) then
        allocate(oldgrid(nsoilold+1))
        allocate(oldval(nsoilold+1))
	allocate(newval(nsoilmx))
        do i=1,imold+1
	 do j=1,jmold+1
	   ! copy values
	   oldval(1)=tsurfold(i,j)
	   oldval(2:nsoilold+1)=tsoilold(i,j,1:nsoilold)
	   ! build vertical coordinate
	   oldgrid(1)=0. ! ground
	   oldgrid(2:nsoilold+1)=mlayerold(1:nsoilold)*
     &                (surfithold(i,j)/1.e6)
          ! Note; at this stage, we impose volcapa=1.e6 above
	  ! since volcapa isn't set in old soil definitions

	  ! interpolate
	  call interp_line(oldgrid,oldval,nsoilold+1,
     &                     mlayer,newval,nsoilmx)
	 ! copy result in tsoilold
	 tsoiloldnew(i,j,1:nsoilmx)=newval(1:nsoilmx)
	 enddo
	enddo
        ! cleanup
	deallocate(oldgrid)
	deallocate(oldval)
	deallocate(newval)

      else
       if (depthinterpol) then ! if vertical interpolation is required
        allocate(oldgrid(nsoilold+1))
        allocate(oldval(nsoilold+1))
	allocate(newval(nsoilmx))
        ! build vertical coordinate
	oldgrid(1)=0. ! ground
	oldgrid(2:nsoilold+1)=mlayerold(1:nsoilold)
        do i=1,imold+1
	 do j=1,jmold+1
	   ! copy values
	   oldval(1)=tsurfold(i,j)
	   oldval(2:nsoilold+1)=tsoilold(i,j,1:nsoilold)
	  ! interpolate
	  call interp_line(oldgrid,oldval,nsoilold+1,
     &                     mlayer,newval,nsoilmx)
	 ! copy result in tsoilold
	 tsoiloldnew(i,j,1:nsoilmx)=newval(1:nsoilmx)
	 enddo
	enddo
!	write(*,*)'tsoiloldnew(1,1,1):',tsoiloldnew(1,1,1)
        ! cleanup
	deallocate(oldgrid)
	deallocate(oldval)
	deallocate(newval)
       
       else
        tsoiloldnew(:,:,:)=tsoilold(:,:,:)
       endif ! of if (depthinterpol)
      endif ! of if (olddepthdef)

      ! Do the horizontal interpolation
       call interp_horiz(tsoiloldnew,tsoilS,
     &                  imold,jmold,iim,jjm,nsoilmx,
     &                   rlonuold,rlatvold,rlonu,rlatv)

      ! Reshape tsoilS to scalar grid as tsoil
       call gr_dyn_fi (nsoilmx,iim+1,jjm+1,ngrid,tsoilS,tsoil)



c-----------------------------------------------------------------------
c 6.3 Variable 3d :
c-----------------------------------------------------------------------
      
c temperatures atmospheriques
      write (*,*) 'lect_start_archive: told ', told (1,jmold+1,1)  ! INFO
      call interp_vert
     &    (told,var,lmold,llm,apsold,bpsold,aps,bps,
     &     psold,(imold+1)*(jmold+1))
      write (*,*) 'lect_start_archive: var ', var (1,jmold+1,1)  ! INFO
      call interp_horiz(var,t,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'lect_start_archive: t ', t(1,jjp1,1)  ! INFO

c q2 : pbl wind variance
      write (*,*) 'lect_start_archive: q2old ', q2old (1,2,1)  ! INFO
      call interp_vert (q2old,varp1,lmold+1,llm+1,
     &     apsold,bpsold,ap,bp,psold,(imold+1)*(jmold+1))
      write (*,*) 'lect_start_archive: varp1 ', varp1 (1,2,1)  ! INFO
      call interp_horiz(varp1,q2s,imold,jmold,iim,jjm,llm+1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'lect_start_archive: q2s ', q2s (1,2,1)  ! INFO
      call gr_dyn_fi (llm+1,iim+1,jjm+1,ngrid,q2s,q2)
      write (*,*) 'lect_start_archive: q2 ', q2 (1,2)  ! INFO
c     write(47,*) 'q2',q2

c calcul des champ de vent; passage en vent covariant
      write (*,*) 'lect_start_archive: uold ', uold (1,2,1)  ! INFO
      call interp_vert
     & (uold,var,lmold,llm,apsold,bpsold,aps,bps,
     &  psold,(imold+1)*(jmold+1))
      write (*,*) 'lect_start_archive: var ', var (1,2,1)  ! INFO
      call interp_horiz(var,us,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'lect_start_archive: us ', us (1,2,1)   ! INFO

      call interp_vert
     & (vold,var,lmold,llm,
     &  apsold,bpsold,aps,bps,psold,(imold+1)*(jmold+1))
      call interp_horiz(var,vs,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call scal_wind(us,vs,unat,vnat)
      write (*,*) 'lect_start_archive: unat ', unat (1,1,1)    ! INFO
      do l=1,llm
        do j = 1, jjp1
          do i=1,iip1
            ucov( i,j,l ) = unat( i,j,l ) * cu(i,j)
c           ucov( i,j,l ) = 0
          end do
        end do
      end do 
      write (*,*) 'lect_start_archive: ucov ', ucov (1,1,1)  ! INFO
c     write(48,*) 'ucov',ucov
      do l=1,llm
        do j = 1, jjm
          do i=1,iim
            vcov( i,j,l ) = vnat( i,j,l ) * cv(i,j)
c           vcov( i,j,l ) = 0
          end do
          vcov( iip1,j,l ) = vcov( 1,j,l )
        end do
      end do
c     write(49,*) 'ucov',vcov

c traceurs surface
      do iq = 1, nqtot
            call interp_horiz(qsurfold(1,1,iq) ,qsurfs(1,1,iq),
     &                  imold,jmold,iim,jjm,1,
     &                  rlonuold,rlatvold,rlonu,rlatv)
      enddo

      call gr_dyn_fi (nqtot,iim+1,jjm+1,ngrid,qsurfs,qsurf)

c traceurs 3D
      do  iq = 1, nqtot
            call interp_vert(qold(1,1,1,iq),var,lmold,llm,
     &        apsold,bpsold,aps,bps,psold,(imold+1)*(jmold+1))
            call interp_horiz(var,q(1,1,1,iq),imold,jmold,iim,jjm,llm,
     &                  rlonuold,rlatvold,rlonu,rlatv)
      enddo
cccccccccccccccccccccccccccccc      
c  make sure that sum of q = 1      
c dominent species is = 1 - sum(all other species)      
cccccccccccccccccccccccccccccc      
c      iqmax=1
c      
c      if (nqold.gt.10) then
c       do l=1,llm
c        do j=1,jjp1
c          do i=1,iip1
c           do iq=1,nqold
c            if (q(i,j,l,iq).gt.q(i,j,l,iqmax)) then
c              iqmax=iq
c            endif
c           enddo
c           q(i,j,l,iqmax)=1.
c           qtot(i,j,l)=0
c           do iq=1,nqold
c            if (iq.ne.iqmax) then        
c              q(i,j,l,iqmax)=q(i,j,l,iqmax)-q(i,j,l,iq)        
c            endif
c           enddo !iq
c           do iq=1,nqold
c            qtot(i,j,l)=qtot(i,j,l)+q(i,j,l,iq)
c            if (i.eq.1.and.j.eq.1.and.l.Eq.1) write(*,*)' qtot(i,j,l)',
c     $    qtot(i,j,l)
c           enddo !iq
c          enddo !i   
c         enddo !j   
c       enddo !l  
c      endif
ccccccccccccccccccccccccccccccc

c     Periodicite :
      do  iq = 1, nqtot
         do l=1, llm
            do j = 1, jjp1
               q(iip1,j,l,iq) = q(1,j,l,iq)
            end do
         end do
      enddo
      
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,co2ices,co2ice)

c-----------------------------------------------------------------------
c   Initialisation  h:	(passage de t -> h)
c-----------------------------------------------------------------------

      DO l=1,llm
         DO j=1,jjp1
            DO i=1,iim
               h(i,j,l) = t(i,j,l)*((ps(i,j)/preff)**kappa)
            END DO
            h(iip1,j,l) =  h(1,j,l)
         END DO
      END DO


c***********************************************************************
c***********************************************************************
c     Fin subroutine lecture ini
c***********************************************************************
c***********************************************************************

      deallocate(timelist)
      deallocate(rlonuold, rlatvold)
      deallocate(rlonvold, rlatuold)
      deallocate(apsold,bpsold)
      deallocate(uold)
      deallocate(vold)
      deallocate(told)
      deallocate(psold)
      deallocate(phisold)
      deallocate(qold)
      deallocate(co2iceold)
      deallocate(tsurfold)
      deallocate(emisold)
      deallocate(q2old)
      deallocate(tsoilold)
      deallocate(tsoiloldnew)
      deallocate(inertiedatold)
      deallocate(inertiedatoldnew)
      deallocate(surfithold)
      deallocate(mlayerold)
      deallocate(qsurfold)
      deallocate(tauscalingold)
      deallocate(totcloudfracold)
      deallocate(var,varp1)

!      write(*,*)'lect_start_archive: END'
      return
      end
