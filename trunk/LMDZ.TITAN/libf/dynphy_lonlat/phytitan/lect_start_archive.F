      SUBROUTINE lect_start_archive(ngrid,nlayer,
     &     date,tsurf,tsoil,emis,q2,
     &     t,ucov,vcov,ps,h,phisold_newgrid,
     &     q,qsurf,surfith,nid)

!      USE surfdat_h
      USE comchem_h, only : cnames, nkim
      USE comchem_newstart_h
      USE comsoil_h, ONLY: nsoilmx, layer, mlayer, volcapa, inertiedat
      USE infotrac, ONLY: tname, nqtot
!      USE control_mod
! to use  'getin'
      USE callkeys_mod, only: callchim
      USE comvert_mod, ONLY: ap,bp,aps,bps,preff
      USE comconst_mod, ONLY: kappa,g,pi

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

      implicit none

#include "dimensions.h"
!#include "dimphys.h"
!#include "planete.h"
#include "paramet.h"
#include "comgeom2.h"
!#include "control.h"
#include "netcdf.inc"
!#include"advtrac.h"
c=======================================================================
c   Declarations
c=======================================================================

      INTEGER,INTENT(IN) :: ngrid, nlayer

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
      REAL date
      INTEGER   memo
!      character (len=50) :: tmpname

c Variable histoire 
c------------------
      REAL vcov(iip1,jjm,llm),ucov(iip1,jjp1,llm) ! vents covariants
      REAL h(iip1,jjp1,llm),ps(iip1,jjp1)
      REAL q(iip1,jjp1,llm,nqtot),qtot(iip1,jjp1,llm)

c autre variables dynamique nouvelle grille
c------------------------------------------

c!-*-
!      integer klatdat,klongdat
!      PARAMETER (klatdat=180,klongdat=360)

c Physique sur grille scalaire 
c----------------------------

c variable physique
c------------------
      REAL tsurf(ngrid) ! surface temperature
      REAL tsoil(ngrid,nsoilmx) ! soil temperature
      REAL emis(ngrid)
      REAL q2(ngrid,llm+1),qsurf(ngrid,nqtot)
c     REAL phisfi(ngrid)

      INTEGER i,j,l
      INTEGER nid,nvarid
c     REAL year_day,periheli,aphelie,peri_day
c     REAL obliquit,z0,emin_turb,lmixmin
c     REAL emissiv,emisice(2),albedice(2)
c     REAL iceradius(2) , dtemisice(2)

      integer ierr
      integer, dimension(4) :: start,count

c Variable nouvelle grille naturelle au point scalaire
c------------------------------------------------------
      real us(iip1,jjp1,llm),vs(iip1,jjp1,llm)
      REAL phisold_newgrid(iip1,jjp1)
      REAL t(iip1,jjp1,llm)
      real tsurfS(iip1,jjp1),tsoilS(iip1,jjp1,nsoilmx)
      real inertiedatS(iip1,jjp1,nsoilmx)
      real emisS(iip1,jjp1)
      REAL q2S(iip1,jjp1,llm+1),qsurfS(iip1,jjp1,nqtot)
      
      real ptotal

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
      real, dimension(:,:), allocatable :: tsurfold
      real, dimension(:,:), allocatable :: emisold
      real, dimension(:,:,:,:), allocatable :: qold
      
      real tab_cntrl(100)

      real ptotalold

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

      real surfith(iip1,jjp1) ! surface thermal inertia
      ! surface thermal inertia at old horizontal grid resolution
      real, dimension(:,:), allocatable :: surfithold 

      character(len=30) :: txt ! to store some text

c=======================================================================

! 0. Preliminary stuff



!-----------------------------------------------------------------------
! 1. Read data dimensions (i.e. size and length)
!-----------------------------------------------------------------------

! 1.2 Read the various dimension lengths of data in file 

      ierr= NF_INQ_DIMID(nid,"Time",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"temps",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,timelen)

      ierr= NF_INQ_DIMID(nid,"latitude",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"rlatu",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,jmold)
      jmold=jmold-1

      ierr= NF_INQ_DIMID(nid,"longitude",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"rlonv",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,imold)
      imold=imold-1

      ierr= NF_INQ_DIMID(nid,"altitude",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"sig_s",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,lmold)

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

! 1.2.2 find out the # of upper chemistry layers

      ierr= NF_INQ_DIMID(nid,"upper_chemistry_layers",dimid)
      ierr= NF_INQ_DIMLEN(nid,dimid,nlaykimold)

! 1.3 Report dimensions
      
      write(*,*) "Start_archive dimensions:"
      write(*,*) "longitude: ",imold
      write(*,*) "latitude: ",jmold
      write(*,*) "altitude: ",lmold
      write(*,*) "tracers: ",nqold
      write(*,*) "subsurface_layers: ",nsoilold
      if (depthinterpol) then
      write(*,*) " => Warning, nsoilmx= ",nsoilmx
      write(*,*) '    which implies that you want subterranean interpola
     &tion.'
      write(*,*) '  Otherwise, set nsoilmx -in dimphys.h- to: ',nsoilold
      endif
      write(*,*) "upper_chemistry_layers: ",nlaykimold
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
      
      allocate(preskimold(nlaykimold))
      allocate(ykim_upS(nkim,iip1,jjp1,nlaykimold))
      allocate(ykim_upoldS(nkim,imold+1,jmold+1,nlaykimold))
      allocate(ykim_up_oldv(nkim,ngrid,nlaykimold))

      allocate(var (imold+1,jmold+1,llm))
      allocate(varp1 (imold+1,jmold+1,llm+1))

      write(*,*) 'q2',ngrid,llm+1
      write(*,*) 'q2S',iip1,jjp1,llm+1
      write(*,*) 'q2old',imold+1,jmold+1,lmold+1

!-----------------------------------------------------------------------
! 3. Read time-independent data
!-----------------------------------------------------------------------

C-----------------------------------------------------------------------
c 3.1. Lecture du tableau des parametres du run 
c     (pour  la lecture ulterieure de "ptotalold")
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "controle", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Lect_start_archive: champ <controle> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echoue pour <controle>"
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
         PRINT*, "lect_start_archive: Le champ <rlonv> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonvold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlonvold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <rlonv>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlatu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <rlatu> est absent"
         CALL abort
      ENDIF 
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatuold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlatuold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <rlatu>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlonu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <rlonu> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonuold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlonuold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <rlonu>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlatv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <rlatv> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatvold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlatvold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <rlatv>"
         CALL abort
      ENDIF
c

c-----------------------------------------------------------------------
c 3.3. Lecture des niveaux verticaux
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "aps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <aps> est absent"
         apsold=0
         PRINT*, "<aps> set to 0"
      ELSE
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid, nvarid, apsold)
#else
         ierr = NF_GET_VAR_REAL(nid, nvarid, apsold)
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: Lecture echouee pour <aps>"
         ENDIF
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "bps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <bps> est absent"
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
         PRINT*, "lect_start_archive: Lecture echouee pour <bps>"
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
c 3.6 Read upper chemistry mid-layer pressure
c-----------------------------------------------------------------------  
    
      ierr=NF_INQ_VARID(nid,"preskim",nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <preskim> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, preskimold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, preskimold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <preskim>"
         CALL abort
      ENDIF

c-----------------------------------------------------------------------
c 3.7 Lecture geopotentiel au sol
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "phisinit", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <phisinit> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, phisold)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, phisold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <phisinit>"
         CALL abort
      ENDIF

C-----------------------------------------------------------------------
c   lecture de "ptotalold"
c-----------------------------------------------------------------------
      ptotalold = tab_cntrl(tab0+49)
 
c-----------------------------------------------------------------------
c 4. Lecture du temps et choix
c-----------------------------------------------------------------------
 
c  lecture du temps
c
      ierr = NF_INQ_DIMID (nid, "Time", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         ierr = NF_INQ_DIMID (nid, "temps", nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: Le champ <Time> est absent"
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
         PRINT*, "lect_start_archive: Lecture echouee pour <Time>"
         CALL abort
      ENDIF
c
      write(*,*)
      write(*,*)
      write(*,*) 'Differentes dates des etats initiaux stockes:'
      write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      pi=2.*ASIN(1.)
      do i=1,timelen
c       call solarlong(timelist(i),sollong(i))
c       sollong(i) = sollong(i)*180./pi
        write(*,*) 'etat initial au jour martien' ,int(timelist(i))
c       write(*,6) nint(timelist(i)),nint(mod(timelist(i),669)),
c    .    sollong(i)
      end do

   6  FORMAT(i7,i7,f9.3)
 
      write(*,*)
      write(*,*) 'Choix de la date'
 123  read(*,*,iostat=ierr) date
      if(ierr.ne.0) goto 123
      memo = 0
      do i=1,timelen
        if (date.eq.int(timelist(i))) then
            memo = i
        endif
      end do
 
      if (memo.eq.0) then
        write(*,*)
        write(*,*)
        write(*,*) 'He alors... Y sait pas lire !?!'
        write(*,*)
        write(*,*) 'Differentes dates des etats initiaux stockes:'
        write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        do i=1,timelen
          write(*,*) 'etat initial au jour martien' ,nint(timelist(i))
c         write(*,6) nint(timelist(i)),nint(mod(timelist(i),669))
        end do
        goto 123
      endif

!-----------------------------------------------------------------------
! 5. Read (time-dependent) data from datafile
!-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c 5.1 Lecture des champs 2D (emis,ps,tsurf,Tg[10], qsurf)
c-----------------------------------------------------------------------
 
      start=(/1,1,memo,0/)
      count=(/imold+1,jmold+1,1,0/)
       

      ierr = NF_INQ_VARID (nid, "emis", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <emis> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,emisold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,emisold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <emis>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "ps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <ps> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,psold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,psold)
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
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,tsurfold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,tsurfold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <tsurf>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "q2surf", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <q2surf> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,q2old)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,q2old)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <q2surf>"
         CALL abort
      ENDIF

c
      write(*,*)"lect_start_archive: rlonuold:"
     &           ,rlonuold," rlatvold:",rlatvold
      write(*,*)

! Surface tracers:      
      do iq=1,nqtot
        ! initialize all surface tracers to zero
        call initial0((jmold+1)*(imold+1), qsurfold(1,1,iq))
      enddo


!      print*,'tname=',tname
!      print*,'nid',nid
!      print*,'nvarid',nvarid
!      stop

      DO iq=1,nqtot
        txt=trim(tname(iq))//"_surf"
        
        write(*,*) "lect_start_archive: loading tracer ",trim(txt)
        ierr = NF_INQ_VARID (nid,txt,nvarid)
        IF (ierr .NE. NF_NOERR) THEN
          PRINT*, "lect_start_archive: ",
     &              " Tracer <",trim(txt),"> not found"

!          print*,'RDW has added hack to let me continue...'
!          CALL abort
        ENDIF
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
          write (*,*) trim(txt),'    is set to 0'
!          call initial0((jmold+1)*(imold+1), qsurfold(1,1,iq))
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
c 5.3 Read 3D upper chemistry fields, if needed
c-----------------------------------------------------------------------  
       
      start=(/1,1,1,memo/)
      count=(/imold+1,jmold+1,nlaykimold,1/)
    
c NB : The sanity check on callchim is on H_up but could be on any chem. field
c as we assume we can't do incomplete chemistry - JVO 18
    
      PRINT*, "lect_start_archive: loading upper chemistry fields..."
      
      ierr=NF_INQ_VARID(nid,"H_up",nvarid)
      
      IF (ierr .NE. NF_NOERR) THEN ! H_up not found
      
         PRINT*, "lect_start_archive: Le champ <H_up> est absent..."        
         IF (callchim) THEN
           PRINT*, "... mais callchim=.TRUE. dans callphys.def !"
           PRINT*, "Verifiez start_archive.nc ou desactivez callchim !"
           CALL abort
         ELSE
           PRINT*, '... je suppose que les autres champs aussi et je
     &passerai donc mon chemin pour tout ce qui concerne la chimie !'
           WRITE (*,*)
         ENDIF
         
      ELSE ! H_up found
      
        IF (.not.callchim) THEN
          PRINT*, "lect_start_archive: Le champ <H_up> est present..."
          PRINT*, "... mais callchim=.FALSE. dans callphys.def !"
          PRINT*, "Si vous voulez gerer la chimie activez callchim !"
          ! CALL abort ! This is too violent to abort here we can start from an archive with chemistry and don't want to use it - JVO !
        ELSE
#ifdef NC_DOUBLE
          ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start,count,
     &                             ykim_upoldS(1,:,:,:))
#else
          ierr = NF_GET_VARA_REAL(nid,nvarid,start,count,
     &                            ykim_upoldS(1,:,:,:))
#endif
          IF (ierr .NE. NF_NOERR) THEN
             PRINT*, "lect_start_archive: Lecture echouee pour <H_up>"
             CALL abort
          ENDIF
                   
          ! Then read all the others by their name if needed
          DO iq=2,nkim
            ierr=NF_INQ_VARID(nid,trim(cnames(iq))//"_up",nvarid)
            IF (ierr .NE. NF_NOERR) THEN
              PRINT*, "lect_start_archive: Le champ <"
     &//trim(cnames(iq))//"_up> est absent..."
              CALL abort
            ENDIF
#ifdef NC_DOUBLE
            ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start,count,
     &                                ykim_upoldS(iq,:,:,:))
#else
            ierr = NF_GET_VARA_REAL(nid,nvarid,start,count,
     &                              ykim_upoldS(iq,:,:,:))
#endif
            IF (ierr .NE. NF_NOERR) THEN
               PRINT*, "lect_start_archive: Lecture echouee pour <"
     &//trim(cnames(iq))//"_up>"
               CALL abort
            ENDIF 
          ENDDO
      
         WRITE(*,*)     
         
         ENDIF ! if callchim
         
      ENDIF  ! if ierr.ne.nf_no_err

c-----------------------------------------------------------------------
c 5.4	Lecture des champs 3D (t,u,v, q2atm,q)
c-----------------------------------------------------------------------

      start=(/1,1,1,memo/)
      count=(/imold+1,jmold+1,lmold,1/)

c
      ierr = NF_INQ_VARID (nid,"temp", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <temp> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid, start, count, told)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid, start, count, told)
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
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,uold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,uold)
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
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,vold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,vold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <v>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid,"q2atm", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <q2atm> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,q2old(1,1,2))
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,q2old(1,1,2))
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <q2atm>"
         CALL abort
      ENDIF
c

! Tracers:      
      do iq=1,nqtot
         call initial0((jmold+1)*(imold+1)*lmold,qold(1,1,1,iq) )
      enddo

      DO iq=1,nqtot
        txt=tname(iq)
        write(*,*)"lect_start_archive: loading tracer ",trim(txt)
        ierr = NF_INQ_VARID (nid,txt,nvarid)
        IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: ",
     &              " Tracer <",trim(txt),"> not found"
!            CALL abort
        ENDIF
#ifdef NC_DOUBLE
        ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start,count,qold(1,1,1,iq))
#else
        ierr=NF_GET_VARA_REAL(nid,nvarid,start,count,qold(1,1,1,iq))
#endif
        IF (ierr .NE. NF_NOERR) THEN
          PRINT*, "lect_start_archive: ",
     &             "  Failed loading <",trim(txt),">"
          write (*,*) trim(txt),'      set to 1.E-30'
          do l=1,lmold
            do j=1,jmold+1
              do i=1,imold+1
                 qold(i,j,l,iq)=1.e-30
              end do
            end do
          end do
        ENDIF

      ENDDO ! of DO iq=1,nqtot


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

      !print*,'Problem in lect_start_archive interpolating'
      !print*,'to new resolution!!'

      ! Recast temperatures along soil depth, if necessary
      if (olddepthdef) then
        allocate(oldgrid(nsoilold+1))
        allocate(oldval(nsoilold+1))
	allocate(newval(nsoilmx))
        do i=1,imold+1
	 do j=1,jmold+1

            !if(i.gt.iip1 .or. j.gt.jjp1)then
               !print*,'Problem in lect_start_archive interpolating'
               !print*,'to new resolution!!'
               !call abort
            !endif

	   ! copy values
	   oldval(1)=tsurfold(i,j)
!	   oldval(1)=tsurfS(i,j)
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
!	   oldval(1)=tsurfS(i,j)
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
c 6.3 Upper chemistry 3d fields:
c
c NB : The vertical regriding, if needed cannot be done here since the
c new pressure grid is only computed at the end of newstart.
c Here we will just do the horizontal interpolation on scalar grid.
c
c-----------------------------------------------------------------------

      ! Do the horizontal interpolation
      DO i=1,nkim
        call interp_horiz(ykim_upoldS(i,:,:,:),ykim_upS(i,:,:,:),
     &                  imold,jmold,iim,jjm,nlaykimold,
     &                   rlonuold,rlatvold,rlonu,rlatv)

        ! Reshape ykim_upS to scalar grid as ykim_up_oldv
        call gr_dyn_fi(nlaykimold,iim+1,jjm+1,ngrid,
     &                 ykim_upS(i,:,:,:),ykim_up_oldv(i,:,:))
      ENDDO


c-----------------------------------------------------------------------
c 6.4 Variable 3d :
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
      write (*,*) 'lect_start_archive: unat ', unat (1,2,1)    ! INFO
      do l=1,llm
        do j = 1, jjp1
          do i=1,iip1
            ucov( i,j,l ) = unat( i,j,l ) * cu(i,j)
c           ucov( i,j,l ) = 0
          end do
        end do
      end do 
      write (*,*) 'lect_start_archive: ucov ', ucov (1,2,1)  ! INFO
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

      if (nqtot .gt. 0) then
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

      endif !! if nqtot .ne. 0

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
      deallocate(var,varp1)
      
      deallocate(ykim_upS)
      deallocate(ykim_upoldS)

!      write(*,*)'lect_start_archive: END'
      return
      end
