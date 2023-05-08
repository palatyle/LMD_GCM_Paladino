
      SUBROUTINE lect_start_archive(date,tsurf,tsoil,emis,q2,
     .     t,ucov,vcov,ps,co2ice,h,phisold_newgrid,q,qsurf,nid)
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
#include "dimphys.h"
#include "surfdat.h"
#include "dimradmars.h"
#include "yomaer.h"
#include "planete.h"
#include "paramet.h"
#include "comconst.h"
#include "comvert.h"
#include "comgeom2.h"
#include "control.h"
#include "logic.h"
#include "description.h"
#include "ener.h"
#include "temps.h"
#include "lmdstd.h"
#include "netcdf.inc"

c=======================================================================
c   Declarations
c=======================================================================

c Variables dimension du fichier "ini"
c------------------------------------
      INTEGER   imold,jmold,lmold,nqold

c et autres:
c----------
      INTEGER lnblnk
      EXTERNAL lnblnk

c Variables pour les lectures des fichiers "ini" 
c--------------------------------------------------
      INTEGER sizei,timelen,dimid
      INTEGER length
      parameter (length = 100)
      INTEGER tab0
      INTEGER isoil,iq,iqmax
      CHARACTER*2   str2

      REAL dimfirst(4) ! tableau contenant les 1ers elements des dimensions

      REAL dimlast(4) ! tableau contenant les derniers elements des dimensions

      REAL dimcycl(4) ! tableau contenant les periodes des dimensions
      CHARACTER*120 dimsource
      CHARACTER*16 dimname
      CHARACTER*80 dimtitle
      CHARACTER*40 dimunits
      INTEGER   dimtype

      INTEGER dimord(4)  ! tableau contenant l''ordre
      data dimord /1,2,3,4/ ! de sortie des dimensions

      INTEGER vardim(4)
      REAL date
      INTEGER   memo
      character (len=50) :: tmpname

c Variable histoire 
c------------------
      REAL vcov(iip1,jjm,llm),ucov(iip1,jjp1,llm) ! vents covariants
      REAL h(iip1,jjp1,llm),ps(iip1,jjp1)
      REAL q(iip1,jjp1,llm,nqmx),qtot(iip1,jjp1,llm)

c autre variables dynamique nouvelle grille
c------------------------------------------

c!-*-
      integer klatdat,klongdat
      PARAMETER (klatdat=180,klongdat=360)

c Physique sur grille scalaire 
c----------------------------

c variable physique
c------------------
      REAL tsurf(ngridmx),tsoil(ngridmx,nsoilmx),co2ice(ngridmx)
      REAL emis(ngridmx)
      REAL q2(ngridmx,nlayermx+1),qsurf(ngridmx,nqmx)
c     REAL phisfi(ngridmx)

      INTEGER i,j,l
      INTEGER nid,nvarid
c     REAL year_day,periheli,aphelie,peri_day
c     REAL obliquit,z0,emin_turb,lmixmin
c     REAL emissiv,emisice(2),albedice(2),tauvis
c     REAL iceradius(2) , dtemisice(2)

      EXTERNAL RAN1
      REAL RAN1
      EXTERNAL geopot,inigeom
      integer ierr
      integer ismin
      external ismin
      CHARACTER*80 datapath
      integer, dimension(4) :: start,count

c Variable nouvelle grille naturelle au point scalaire
c------------------------------------------------------
      real us(iip1,jjp1,llm),vs(iip1,jjp1,llm)
      REAL phisold_newgrid(iip1,jjp1)
      REAL t(iip1,jjp1,llm)
      real tsurfS(iip1,jjp1),tsoilS(iip1,jjp1,nsoilmx)
      real co2iceS(iip1,jjp1),emisS(iip1,jjp1)
      REAL q2S(iip1,jjp1,llm+1),qsurfS(iip1,jjp1,nqmx)

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
      real, dimension(:,:,:), allocatable :: uold,vold,told,q2old
      real, dimension(:,:,:), allocatable :: tsoilold,qsurfold
      real, dimension(:,:), allocatable :: psold,phisold
      real, dimension(:,:), allocatable :: co2iceold,tsurfold
      real, dimension(:,:), allocatable :: emisold
      real, dimension(:,:,:,:), allocatable :: qold

      real tab_cntrl(100)

      real ptotalold, co2icetotalold

c Variable intermediaires iutilise pour l'extrapolation verticale 
c----------------------------------------------------------------
      real, dimension(:,:,:), allocatable :: var,varp1 

c=======================================================================

c Catching the axis lenghts for dynamic memory allocation 

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

      write(*,*) "Start_archive dimensions:"
      write(*,*) "longitude: ",imold
      write(*,*) "latitude: ",jmold
      write(*,*) "altitude: ",lmold
      write(*,*) "tracers: ",nqold
      write(*,*) "time lenght: ",timelen
      write(*,*) 

      allocate(timelist(timelen))
      allocate(rlonuold(imold+1), rlatvold(jmold))
      allocate(rlonvold(imold+1), rlatuold(jmold+1))
      allocate (apsold(lmold),bpsold(lmold))
      allocate(uold(imold+1,jmold+1,lmold))
      allocate(vold(imold+1,jmold+1,lmold))
      allocate(told(imold+1,jmold+1,lmold))
      allocate(psold(imold+1,jmold+1))
      allocate(phisold(imold+1,jmold+1))
      allocate(qold(imold+1,jmold+1,lmold,nqmx))
      allocate(co2iceold(imold+1,jmold+1))
      allocate(tsurfold(imold+1,jmold+1))
      allocate(emisold(imold+1,jmold+1))
      allocate(q2old(imold+1,jmold+1,lmold+1))
      allocate(tsoilold(imold+1,jmold+1,nsoilmx))
      allocate(qsurfold(imold+1,jmold+1,nqmx))

      allocate(var (imold+1,jmold+1,llm))
      allocate(varp1 (imold+1,jmold+1,llm+1))

      write(*,*) 'q2',ngridmx,nlayermx+1
      write(*,*) 'q2S',iip1,jjp1,llm+1
      write(*,*) 'q2old',imold+1,jmold+1,lmold+1
 
C-----------------------------------------------------------------------
c Lecture du tableau des parametres du run 
c     (pour  la lecture ulterieure de "ptotalold" et "co2icetotalold")
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
c Lecture des longitudes et latitudes
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
c Lecture des niveaux verticaux
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
c Lecture geopotentiel au sol
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
c   lecture de "ptotalold" et "co2icetotalold"
c-----------------------------------------------------------------------
      ptotalold = tab_cntrl(tab0+49)
      co2icetotalold = tab_cntrl(tab0+50)
 
c-----------------------------------------------------------------------
c   Lecture du temps et choix
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
 
c-----------------------------------------------------------------------
c Lecture des champs 2D (co2ice, emis,ps,tsurf,Tg[10], q2surf)
c-----------------------------------------------------------------------
 
      start=(/1,1,memo,0/)
      count=(/imold+1,jmold+1,1,0/)
       
      ierr = NF_INQ_VARID (nid, "co2ice", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <co2ice> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,co2iceold)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,co2iceold)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <co2ice>"
         PRINT*, NF_STRERROR(ierr)
         CALL abort
      ENDIF
c
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
      do isoil=1,nsoilmx
         write(str2,'(i2.2)') isoil
c
         ierr = NF_INQ_VARID (nid, "Tg"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .              Le champ <","Tg"//str2,"> est absent"
            CALL abort
         ENDIF
#ifdef NC_DOUBLE
         ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,
     .          tsoilold(1,1,isoil))
#else
         ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,
     .          tsoilold(1,1,isoil))
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .            Lecture echouee pour <","Tg"//str2,">"
            CALL abort
         ENDIF
c
      end do

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
      write (*,*) "rlonuold,rlatvold"
      write (*,*) rlonuold
      write (*,*) rlatvold
      write (*,*)
c

c tracers: the 2 last ones are kept the 2 last one. 
c the others keep their rank.
c -------------------------------------------
      
      do iq=1,nqmx
        call initial0((jmold+1)*(imold+1), qsurfold(1,1,iq))
      enddo

      iq=nqold
        write(str2,'(i2.2)') iq
         ierr = NF_INQ_VARID (nid, "qsurf"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .               Le champ <","qsurf"//str2,"> est absent"
            CALL abort
         ENDIF
#ifdef NC_DOUBLE
         ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,
     .          qsurfold(1,1,nqmx))
#else
         ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,
     .          qsurfold(1,1,nqmx))
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .               Lecture echouee pour <","qsurf"//str2,">"
            write (*,*) 'qsurf'//str2,' set to 0'
            call initial0((jmold+1)*(imold+1), qsurfold(1,1,nqmx))
         ENDIF

      if ((nqold.gt.1).and.(nqmx.gt.1)) then
        iq=nqold-1
        write(str2,'(i2.2)') iq
         ierr = NF_INQ_VARID (nid, "qsurf"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .               Le champ <","qsurf"//str2,"> est absent"
            CALL abort
         ENDIF
#ifdef NC_DOUBLE
         ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,
     .          qsurfold(1,1,nqmx-1))
#else
         ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,
     .          qsurfold(1,1,nqmx-1))
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .               Lecture echouee pour <","qsurf"//str2,">"
            write (*,*) 'qsurf'//str2,' set to 0'
            call initial0((jmold+1)*(imold+1), qsurfold(1,1,nqmx-1))
         ENDIF
      endif

      if (nqold.gt.2) then
       do  iq = 1, nqold-2
       if (iq.lt.nqmx-1) then
         write(str2,'(i2.2)') iq
         ierr = NF_INQ_VARID (nid, "qsurf"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .               Le champ <","qsurf"//str2,"> est absent"
            CALL abort
         ENDIF
#ifdef NC_DOUBLE
         ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,
     .          qsurfold(1,1,iq))
#else
         ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,
     .          qsurfold(1,1,iq))
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .               Lecture echouee pour <","qsurf"//str2,">"
            write (*,*) 'qsurf'//str2,' set to 0'
            call initial0((jmold+1)*(imold+1), qsurfold(1,1,iq))
         ENDIF
       end if
       end do
      end if

c-----------------------------------------------------------------------
c	Lecture des champs 3D (t,u,v, q2atm,q)
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

c tracers: the 2 last ones are kept the 2 last one. 
c the others keep their rank.
c -------------------------------------------
      
      do iq=1,nqmx
         call initial0((jmold+1)*(imold+1)*lmold,qold(1,1,1,iq) )
      enddo

      iq=nqold
        write(str2,'(i2.2)') iq
         ierr = NF_INQ_VARID (nid, "q"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive:
     .               Le champ <","q"//str2,"> est absent"
            CALL abort
         ENDIF
#ifdef NC_DOUBLE
         ierr= NF_GET_VARA_DOUBLE(nid,nvarid,start,count,qold(1,1,1,nqmx))
#else
         ierr= NF_GET_VARA_REAL(nid,nvarid,start,count,qold(1,1,1,nqmx))
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .               Lecture echouee pour <","q"//str2,">"
            write (*,*) 'q'//str2,' set to 1.E-30'
            do l=1,lmold
              do j=1,jmold+1
                do i=1,imold+1
                   qold(1,1,1,nqmx)=1.e-30
                end do
              end do
            end do

         ENDIF

      if ((nqold.gt.1).and.(nqmx.gt.1)) then
        iq=nqold-1
        write(str2,'(i2.2)') iq
         ierr = NF_INQ_VARID (nid, "q"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive:
     .               Le champ <","q"//str2,"> est absent"
            CALL abort
         ENDIF
#ifdef NC_DOUBLE
         ierr= NF_GET_VARA_DOUBLE(nid,nvarid,start,count,
     .                            qold(1,1,1,nqmx-1))
#else
         ierr= NF_GET_VARA_REAL(nid,nvarid,start,count,
     .                            qold(1,1,1,nqmx-1))
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .               Lecture echouee pour <","q"//str2,">"
            write (*,*) 'q'//str2,' set to 1.E-30'
            do l=1,lmold
              do j=1,jmold+1
                do i=1,imold+1
                   qold(1,1,1,nqmx-1)=1.e-30
                end do
              end do
            end do

         ENDIF
      endif

      if (nqold.gt.2) then
       do  iq = 1, nqold-2
       if (iq.lt.nqmx-1) then
         write(str2,'(i2.2)') iq
         ierr = NF_INQ_VARID (nid, "q"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive:
     .               Le champ <","q"//str2,"> est absent"
            CALL abort
         ENDIF
#ifdef NC_DOUBLE
         ierr= NF_GET_VARA_DOUBLE(nid,nvarid,start,count,qold(1,1,1,iq))
#else
         ierr= NF_GET_VARA_REAL(nid,nvarid,start,count,qold(1,1,1,iq))
#endif
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: 
     .               Lecture echouee pour <","q"//str2,">"
            write (*,*) 'q'//str2,' set to 1.E-30 '
            do l=1,lmold
              do j=1,jmold+1
                do i=1,imold+1
                   qold(1,1,1,iq)=1.e-30
                end do
              end do
            end do

         ENDIF
       end if
       end do
      end if

c Chemin pour trouver les donnees de surface (albedo, relief, th.inertia...)
c -------------------------------------------------------------------------

      datapath = '/users/forget/gcm/data_mars_gcm'


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
      call interp_horiz (phisold,phisold_newgrid,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)

c Glace CO2
      call interp_horiz (co2iceold,co2ices,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)

c Temperature de surface
      call interp_horiz (tsurfold,tsurfs,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngridmx,tsurfs,tsurf)
c     write(44,*) 'tsurf', tsurf

c Temperature du sous-sol
      call interp_horiz(tsoilold,tsoils,
     &                  imold,jmold,iim,jjm,nsoilmx,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (nsoilmx,iim+1,jjm+1,ngridmx,tsoils,tsoil)
c     write(45,*) 'tsoil',tsoil

c Emissivite de la surface
      call interp_horiz (emisold,emiss,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngridmx,emiss,emis)
c     write(46,*) 'emis',emis
c-----------------------------------------------------------------------
c	Traitement special de la pression au sol :
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
      write(*,*)'Ancienne grille: masse de l atm :',ptotalold
      write(*,*)'Nouvelle grille: masse de l atm :',ptotal
      write (*,*) 'Ratio new atm./ old atm =', ptotal/ptotalold 
      write(*,*)
      write(*,*)'Ancienne grille: masse de la glace CO2:',co2icetotalold
      write(*,*)'Nouvelle grille: masse de la glace CO2:',co2icetotal
      write(*,*)'Ratio new ice./old ice =',co2icetotal/co2icetotalold
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
c	Variable 3d :
c-----------------------------------------------------------------------
      
c temperatures atmospheriques
      write (*,*) 'told ', told (1,jmold+1,1)  ! INFO
      call interp_vert
     &    (told,var,lmold,llm,apsold,bpsold,aps,bps,
     &     psold,(imold+1)*(jmold+1))
      write (*,*) 'var ', var (1,jmold+1,1)  ! INFO
      call interp_horiz(var,t,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 't ', t(1,jjp1,1)  ! INFO

c q2 : pbl wind variance
      write (*,*) 'q2old ', q2old (1,2,1)  ! INFO
      call interp_vert (q2old,varp1,lmold+1,llm+1,
     &     apsold,bpsold,aps,bps,psold,(imold+1)*(jmold+1))
      write (*,*) 'varp1 ', varp1 (1,2,1)  ! INFO
      call interp_horiz(varp1,q2s,imold,jmold,iim,jjm,llm+1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'q2s ', q2s (1,2,1)  ! INFO
      call gr_dyn_fi (llm+1,iim+1,jjm+1,ngridmx,q2s,q2)
      write (*,*) 'q2 ', q2 (1,2)  ! INFO
c     write(47,*) 'q2',q2

c calcul des champ de vent; passage en vent covariant
      write (*,*) 'uold ', uold (1,2,1)  ! INFO
      call interp_vert
     & (uold,var,lmold,llm,apsold,bpsold,aps,bps,
     &  psold,(imold+1)*(jmold+1))
      write (*,*) 'var ', var (1,2,1)  ! INFO
      call interp_horiz(var,us,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'us ', us (1,2,1)   ! INFO

      call interp_vert
     & (vold,var,lmold,llm,
     &  apsold,bpsold,aps,bps,psold,(imold+1)*(jmold+1))
      call interp_horiz(var,vs,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call scal_wind(us,vs,unat,vnat)
      write (*,*) 'unat ', unat (1,2,1)    ! INFO
      do l=1,llm
        do j = 1, jjp1
          do i=1,iip1
            ucov( i,j,l ) = unat( i,j,l ) * cu(i,j)
c           ucov( i,j,l ) = 0
          end do
        end do
      end do 
      write (*,*) 'ucov ', ucov (1,2,1)  ! INFO
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
      do iq = 1, nqmx
            call interp_horiz(qsurfold(1,1,iq) ,qsurfs(1,1,iq),
     &                  imold,jmold,iim,jjm,1,
     &                  rlonuold,rlatvold,rlonu,rlatv)
      enddo

      call gr_dyn_fi (nqmx,iim+1,jjm+1,ngridmx,qsurfs,qsurf)

c traceurs 3D
      do  iq = 1, nqmx
            call interp_vert(qold(1,1,1,iq),var,lmold,llm,
     &        apsold,bpsold,aps,bps,psold,(imold+1)*(jmold+1))
            call interp_horiz(var,q(1,1,1,iq),imold,jmold,iim,jjm,llm,
     &                  rlonuold,rlatvold,rlonu,rlatv)
      enddo
cccccccccccccccccccccccccccccc      
c  make sure that sum of q = 1      
c dominent species is = 1 - sum(all other species)      
cccccccccccccccccccccccccccccc      
c     iqmax=1
c     
c     if (nqold.gt.10) then
c      do l=1,llm
c       do j=1,jjp1
c         do i=1,iip1
c          do iq=1,nqold
c           if (q(i,j,l,iq).gt.q(i,j,l,iqmax)) then
c             iqmax=iq
c            endif
c          enddo
c          q(i,j,l,iqmax)=1.
c          qtot(i,j,l)=0
c          do iq=1,nqold
c           if (iq.ne.iqmax) then        
c             q(i,j,l,iqmax)=q(i,j,l,iqmax)-q(i,j,l,iq)        
c           endif
c          enddo !iq
c          do iq=1,nqold
c           qtot(i,j,l)=qtot(i,j,l)+q(i,j,l,iq)
c            if (i.eq.1.and.j.eq.1.and.l.Eq.1) write(*,*)' qtot(i,j,l)',
c     $    qtot(i,j,l)
c          enddo !iq
c         enddo !i   
c        enddo !j   
c      enddo !l  
c     endif
ccccccccccccccccccccccccccccccc

c     Periodicite :
      do  iq = 1, nqmx
         do l=1, llm
            do j = 1, jjp1
               q(iip1,j,l,iq) = q(1,j,l,iq)
            end do
         end do
      enddo
      
      call gr_dyn_fi (1,iim+1,jjm+1,ngridmx,co2ices,co2ice)

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
      deallocate(qsurfold)
      deallocate(var,varp1)

      return
      end
