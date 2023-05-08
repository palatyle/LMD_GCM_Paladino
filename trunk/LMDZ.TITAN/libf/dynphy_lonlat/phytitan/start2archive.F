c=======================================================================
      PROGRAM start2archive
c=======================================================================
c
c
c   Date:    01/1997
c   ----
c
c
c   Objet:   Passage des  fichiers netcdf d'etat initial "start" et
c   -----    "startfi" a un fichier netcdf unique "start_archive" 
c
c  "start_archive" est une banque d'etats initiaux:
c  On peut stocker plusieurs etats initiaux dans un meme fichier "start_archive"
c    (Veiller dans ce cas avoir un day_ini different pour chacun des start)
c 
c
c
c=======================================================================

      use infotrac, only: infotrac_init, nqtot, tname
      USE comsoil_h
      USE comchem_h, only : cnames, nkim, nlaykim_up, preskim, ykim_up
!      USE comgeomfi_h, ONLY: lati, long, area
!      use control_mod
!      use comgeomphy, only: initcomgeomphy
! to use  'getin'
      USE ioipsl_getincom
      USE planete_mod, only: year_day
      USE mod_const_mpi, ONLY: COMM_LMDZ
      USE control_mod, only: planet_type
      use filtreg_mod, only: inifilr
      USE comvert_mod, ONLY: ap,bp
      USE comconst_mod, ONLY: daysec,dtphys,rad,g,r,cpp
      USE temps_mod, ONLY: day_ini
      USE iniphysiq_mod, ONLY: iniphysiq
      use phyetat0_mod, only: phyetat0
      use tracer_h
      implicit none

      include "dimensions.h"
      integer, parameter :: ngridmx = (2+(jjm-1)*iim - 1/jjm) 
      include "paramet.h"
      include "comdissip.h"
      include "comgeom.h"
!#include "control.h"

!#include "dimphys.h"
!#include "planete.h"
!#include"advtrac.h"
      include "netcdf.inc"
c-----------------------------------------------------------------------
c   Declarations
c-----------------------------------------------------------------------

c variables dynamiques du GCM
c -----------------------------
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) ! vents covariants
      REAL teta(ip1jmp1,llm)                    ! temperature potentielle 
      REAL,ALLOCATABLE :: q(:,:,:)   ! champs advectes
      REAL pks(ip1jmp1)                      ! exner (f pour filtre)
      REAL pk(ip1jmp1,llm)
      REAL pkf(ip1jmp1,llm)
      REAL beta(iip1,jjp1,llm)
      REAL phis(ip1jmp1)                     ! geopotentiel au sol
      REAL masse(ip1jmp1,llm)                ! masse de l'atmosphere
      REAL ps(ip1jmp1)                       ! pression au sol
      REAL p3d(iip1, jjp1, llm+1)            ! pression aux interfaces
      
c Variable Physiques (grille physique)
c ------------------------------------
      REAL tsurf(ngridmx)        ! Surface temperature
      REAL,ALLOCATABLE :: tsoil(:,:) ! Soil temperature
      REAL q2(ngridmx,llm+1)
      REAL,ALLOCATABLE :: qsurf(:,:)
      REAL emis(ngridmx)
      INTEGER start,length
      PARAMETER (length = 100)
      REAL tab_cntrl_fi(length) ! tableau des parametres de startfi
      REAL tab_cntrl_dyn(length) ! tableau des parametres de start
      INTEGER*4 day_ini_fi

c     Added by JVO for Titan specifities
      REAL tankCH4(ngridmx) ! Depth of surface methane tank

c Variable naturelle / grille scalaire
c ------------------------------------
      REAL T(ip1jmp1,llm),us(ip1jmp1,llm),vs(ip1jmp1,llm)
      REAL tsurfS(ip1jmp1)
      REAL,ALLOCATABLE :: tsoilS(:,:)
      REAL,ALLOCATABLE :: ithS(:,:) ! Soil Thermal Inertia
      REAL q2S(ip1jmp1,llm+1)
      REAL,ALLOCATABLE :: qsurfS(:,:)
      REAL emisS(ip1jmp1)

c     Added by JVO for Titan specifities
      REAL tankCH4S(ip1jmp1)  ! Depth of surface methane tank

      REAL,ALLOCATABLE :: ykim_upS(:,:,:)  ! Upper chemistry fields
 
c Variables intermediaires : vent naturel, mais pas coord scalaire
c----------------------------------------------------------------
      REAL vn(ip1jm,llm),un(ip1jmp1,llm)

c Autres  variables
c -----------------
      LOGICAL startdrs
      INTEGER Lmodif

      LOGICAL kim

      REAL ptotal
      REAL timedyn,timefi !fraction du jour dans start, startfi
      REAL date

      CHARACTER*2 str2
      CHARACTER*80 fichier 
      data  fichier /'startfi'/

      INTEGER ij, l,i,j,isoil,iq
      character*80      fichnom
      integer :: ierr,ntime
      integer :: nq,numvanle
      character(len=30) :: txt ! to store some text

c Netcdf
c-------
      integer varid,dimid,timelen 
      INTEGER nid,nid1

c-----------------------------------------------------------------------
c   Initialisations 
c-----------------------------------------------------------------------

      CALL defrun_new(99, .TRUE. )

      planet_type="titan"

c=======================================================================
c Lecture des donnees
c=======================================================================
! Load tracer number and names:
      call infotrac_init

! allocate arrays:
      allocate(q(ip1jmp1,llm,nqtot))
      allocate(qsurf(ngridmx,nqtot))
      allocate(qsurfS(ip1jmp1,nqtot))
! other array allocations:
!      call ini_comsoil_h(ngridmx) ! done via iniphysiq

      fichnom = 'start.nc'
      CALL dynetat0(fichnom,vcov,ucov,teta,q,masse,
     .       ps,phis,timedyn)

! load 'controle' array from dynamics start file

       ierr = NF_OPEN (fichnom, NF_NOWRITE,nid1)
       IF (ierr.NE.NF_NOERR) THEN
         write(6,*)' Pb d''ouverture du fichier'//trim(fichnom)
        CALL ABORT
       ENDIF
                                                
      ierr = NF_INQ_VARID (nid1, "controle", varid)
      IF (ierr .NE. NF_NOERR) THEN
       PRINT*, "start2archive: Le champ <controle> est absent"
       CALL abort
      ENDIF
#ifdef NC_DOUBLE
       ierr = NF_GET_VAR_DOUBLE(nid1, varid, tab_cntrl_dyn)
#else
      ierr = NF_GET_VAR_REAL(nid1, varid, tab_cntrl_dyn)
#endif
       IF (ierr .NE. NF_NOERR) THEN
          PRINT*, "start2archive: Lecture echoue pour <controle>"
          CALL abort
       ENDIF

      ierr = NF_CLOSE(nid1)

! Get value of the "subsurface_layers" dimension from physics start file
      fichnom = 'startfi.nc'
      ierr = NF_OPEN (fichnom, NF_NOWRITE,nid1)
       IF (ierr.NE.NF_NOERR) THEN
         write(6,*)' Pb d''ouverture du fichier'//trim(fichnom)
        CALL ABORT
       ENDIF
      ierr = NF_INQ_DIMID(nid1,"subsurface_layers",varid)
      IF (ierr .NE. NF_NOERR) THEN
       PRINT*, "start2archive: No subsurface_layers dimension!!"
       CALL abort
      ENDIF
      ierr = NF_INQ_DIMLEN(nid1,varid,nsoilmx)
      IF (ierr .NE. NF_NOERR) THEN
       PRINT*, "start2archive: Failed reading subsurface_layers value!!"
       CALL abort
      ENDIF
      
      ! allocate arrays of nsoilmx size
      allocate(tsoil(ngridmx,nsoilmx))
      allocate(tsoilS(ip1jmp1,nsoilmx))
      allocate(ithS(ip1jmp1,nsoilmx))

! Get value of the "upper_chemistry_layers" dimension from physics start file

      ! NB : nlaykim_up is recalculated in phyetat0 ... this is redundant for newstart
      ! If you want to do it properly put ykim_upS in a specific separate header
      ! Anyway it should be safe. -- JVO 18

      ierr = NF_INQ_DIMID(nid1,"upper_chemistry_layers",varid)
      IF (ierr .NE. NF_NOERR) THEN
       PRINT*, "start2archive: No upper_chemistry_layers dimension!!"
       CALL abort
      ENDIF
      ierr = NF_INQ_DIMLEN(nid1,varid,nlaykim_up)
      IF (ierr .NE. NF_NOERR) THEN
       PRINT*, "start2archive: Failed reading 
     . upper_chemistry_layers value!!"
       CALL abort
      ENDIF

      ALLOCATE(preskim(nlaykim_up))
      
      ! Allocate other arrays of nlaykim_up size, only if they're present
      ! The test is on H but could be on any as we assume we can't do incomplete chemistry

      ierr = NF_INQ_VARID(nid1,'H_up',varid)
      IF (ierr .NE. NF_NOERR) THEN
        PRINT*, "start2archive: Missing field(s) for upper chemistry ...
     . I presume they're all absent !"
        kim=.FALSE.
      ELSE
        PRINT*,"start2archive: I found a field for upper chemistry ...
     . I presume they're all here as you can't do uncomplete chemistry!"
     
        ALLOCATE(ykim_up(nkim,ngridmx,nlaykim_up))
        ALLOCATE(ykim_upS(nkim,ip1jmp1,nlaykim_up))
        
        kim=.TRUE.
      ENDIF

      ierr = NF_CLOSE(nid1)

c-----------------------------------------------------------------------
c   Initialisations 
c-----------------------------------------------------------------------

      CALL defrun_new(99, .FALSE. )
      call iniconst
      call inigeom
      call inifilr

! Initialize the physics
         CALL iniphysiq(iim,jjm,llm,
     &                  (jjm-1)*iim+2,comm_lmdz,
     &                  daysec,day_ini,dtphys,
     &                  rlatu,rlatv,rlonu,rlonv,
     &                  aire,cu,cv,rad,g,r,cpp,
     &                  1)

      fichnom = 'startfi.nc'
      Lmodif=0

! Initialize tracer names, indexes and properties
      CALL initracer2(nqtot,tname)

      CALL phyetat0(.true.,ngridmx,llm,fichnom,0,Lmodif,nsoilmx,nqtot,
     .      day_ini_fi,timefi,
     .      tsurf,tsoil,emis,q2,qsurf,tankCH4)


! load 'controle' array from physics start file

       ierr = NF_OPEN (fichnom, NF_NOWRITE,nid1)
       IF (ierr.NE.NF_NOERR) THEN
         write(6,*)' Pb d''ouverture du fichier'//trim(fichnom)
        CALL ABORT
       ENDIF
                                                
      ierr = NF_INQ_VARID (nid1, "controle", varid)
      IF (ierr .NE. NF_NOERR) THEN
       PRINT*, "start2archive: Le champ <controle> est absent"
       CALL abort
      ENDIF
#ifdef NC_DOUBLE
       ierr = NF_GET_VAR_DOUBLE(nid1, varid, tab_cntrl_fi)
#else
      ierr = NF_GET_VAR_REAL(nid1, varid, tab_cntrl_fi)
#endif
       IF (ierr .NE. NF_NOERR) THEN
          PRINT*, "start2archive: Lecture echoue pour <controle>"
          CALL abort
       ENDIF

      ierr = NF_CLOSE(nid1)


c-----------------------------------------------------------------------
c Controle de la synchro
c-----------------------------------------------------------------------
!mars a voir      if ((day_ini_fi.ne.day_ini).or.(abs(timefi-timedyn).gt.1.e-10)) 
      if ((day_ini_fi.ne.day_ini)) 
     &  stop ' Probleme de Synchro entre start et startfi !!!'


c *****************************************************************
c    Option : Reinitialisation des dates dans la premieres annees :
       do while (day_ini.ge.year_day)
          day_ini=day_ini-year_day
       enddo
c *****************************************************************

      CALL pression(ip1jmp1, ap, bp, ps, p3d)
      call exner_hyb(ip1jmp1, ps, p3d, beta, pks, pk, pkf)

c=======================================================================
c Transformation EN VARIABLE NATURELLE / GRILLE SCALAIRE si necessaire
c=======================================================================
c  Les variables modeles dependent de la resolution. Il faut donc
c  eliminer les facteurs responsables de cette dependance
c  (pour utiliser newstart)
c=======================================================================

c-----------------------------------------------------------------------
c Vent   (depend de la resolution horizontale) 
c-----------------------------------------------------------------------
c
c ucov --> un  et  vcov --> vn
c un --> us  et   vn --> vs
c
c-----------------------------------------------------------------------

      call covnat(llm,ucov, vcov, un, vn) 
      call wind_scal(un,vn,us,vs) 

c-----------------------------------------------------------------------
c Temperature  (depend de la resolution verticale => de "sigma.def")
c-----------------------------------------------------------------------
c
c h --> T
c
c-----------------------------------------------------------------------

      DO l=1,llm
         DO ij=1,ip1jmp1
            T(ij,l)=teta(ij,l)*pk(ij,l)/cpp !mars deduit de l'equation dans newstart
         ENDDO
      ENDDO

c-----------------------------------------------------------------------
c Variable physique 
c-----------------------------------------------------------------------
c
c tsurf --> tsurfS
c tsoil --> tsoilS
c emis --> emisS
c q2 --> q2S
c qsurf --> qsurfS
c tankCH4 --> tankCH4S
c ykim_up --> ykim_upS ( nkim=44 upper chemistry fields )
c
c-----------------------------------------------------------------------

      call gr_fi_dyn(1,ngridmx,iip1,jjp1,tsurf,tsurfS)
      call gr_fi_dyn(nsoilmx,ngridmx,iip1,jjp1,tsoil,tsoilS)
      ! Note: thermal inertia "inertiedat" is in comsoil.h
      call gr_fi_dyn(nsoilmx,ngridmx,iip1,jjp1,inertiedat,ithS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,emis,emisS)
      call gr_fi_dyn(llm+1,ngridmx,iip1,jjp1,q2,q2S)
      call gr_fi_dyn(nqtot,ngridmx,iip1,jjp1,qsurf,qsurfS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,tankCH4,tankCH4S)

      IF (kim) THEN ! NB : fields are in comchem_startarch_h
         DO iq=1,nkim
           call gr_fi_dyn(nlaykim_up,ngridmx,iip1,jjp1,ykim_up(iq,:,:),
     &                                                 ykim_upS(iq,:,:))
         ENDDO
      ENDIF

c=======================================================================
c Info pour controler
c=======================================================================

      ptotal =  0.
      DO j=1,jjp1
         DO i=1,iim
           ptotal=ptotal+aire(i+(iim+1)*(j-1))*ps(i+(iim+1)*(j-1))/g
         ENDDO
      ENDDO
      write(*,*)'Ancienne grille : masse de l''atm :',ptotal

c-----------------------------------------------------------------------
c Passage de "ptotal"  par tab_cntrl_fi
c-----------------------------------------------------------------------

      tab_cntrl_fi(49) = ptotal
      tab_cntrl_fi(50) = 0.

c=======================================================================
c Ecriture dans le fichier  "start_archive"
c=======================================================================

c-----------------------------------------------------------------------
c Ouverture de "start_archive" 
c-----------------------------------------------------------------------

      ierr = NF_OPEN ('start_archive.nc', NF_WRITE,nid)
 
c-----------------------------------------------------------------------
c  si "start_archive" n'existe pas:
c    1_ ouverture
c    2_ creation de l'entete dynamique ("ini_archive")
c-----------------------------------------------------------------------
c ini_archive:
c On met dans l'entete le tab_cntrl dynamique (1 a 16) 
c  On y ajoute les valeurs du tab_cntrl_fi (a partir de 51)
c  En plus les deux valeurs ptotal et co2icetotal (99 et 100)
c-----------------------------------------------------------------------

      if (ierr.ne.NF_NOERR) then
         write(*,*)'OK, Could not open file "start_archive.nc"'
         write(*,*)'So let s create a new "start_archive"'
         ierr = NF_CREATE('start_archive.nc', NF_CLOBBER, nid)
         call ini_archive(nid,day_ini,phis,ithS,tab_cntrl_fi,
     &                                          tab_cntrl_dyn)
      endif

c-----------------------------------------------------------------------
c Ecriture de la coordonnee temps (date en jours)
c-----------------------------------------------------------------------

      date = day_ini
      ierr= NF_INQ_VARID(nid,"Time",varid)
      ierr= NF_INQ_DIMID(nid,"Time",dimid)
      ierr= NF_INQ_DIMLEN(nid,dimid,timelen)
      ntime=timelen+1

      write(*,*) "******************"
      write(*,*) "ntime",ntime
      write(*,*) "******************"
#ifdef NC_DOUBLE
      ierr= NF_PUT_VARA_DOUBLE(nid,varid,ntime,1,date)
#else
      ierr= NF_PUT_VARA_REAL(nid,varid,ntime,1,date)
#endif
      if (ierr.ne.NF_NOERR) then
         write(*,*) "time matter ",NF_STRERROR(ierr)
         stop
      endif

c-----------------------------------------------------------------------
c Ecriture des champs  (emis,ps,Tsurf,T,u,v,q2,q,qsurf,tankCH4)
c-----------------------------------------------------------------------
c ATTENTION: q2 a une couche de plus!!!!
c    Pour creer un fichier netcdf lisible par grads,
c    On passe donc une des couches de q2 a part
c    comme une variable 2D (la couche au sol: "q2surf")
c    Les lmm autres couches sont nommees "q2atm" (3D) 
c-----------------------------------------------------------------------

      call write_archive(nid,ntime,'emis','grd emis',' ',2,emisS)
      call write_archive(nid,ntime,'ps','Psurf','Pa',2,ps)
      call write_archive(nid,ntime,'tsurf','surf T','K',2,tsurfS)
      call write_archive(nid,ntime,'temp','temperature','K',3,t)
      call write_archive(nid,ntime,'u','Vent zonal','m.s-1',3,us)
      call write_archive(nid,ntime,'v','Vent merid','m.s-1',3,vs)
      call write_archive(nid,ntime,'q2surf','wind variance','m2.s-2',2,
     .              q2S)
      call write_archive(nid,ntime,'q2atm','wind variance','m2.s-2',3,
     .              q2S(1,2))

c-----------------------------------------------------------------------
c Ecriture du champs  q  ( q[1,nqtot] )
c-----------------------------------------------------------------------
      do iq=1,nqtot
        call write_archive(nid,ntime,tname(iq),'tracer','kg/kg',
     &         3,q(1,1,iq))
      end do
c-----------------------------------------------------------------------
c Ecriture du champs  qsurf  ( qsurf[1,nqtot] )
c-----------------------------------------------------------------------
      do iq=1,nqtot
        txt=trim(tname(iq))//"_surf"
        call write_archive(nid,ntime,txt,'Tracer on surface',
     &  'kg.m-2',2,qsurfS(1,iq))
      enddo


c-----------------------------------------------------------------------
c Ecriture du champs  tsoil  ( Tg[1,10] )
c-----------------------------------------------------------------------
c "tsoil" Temperature au sol definie dans 10 couches dans le sol
c   Les 10 couches sont lues comme 10 champs 
c  nommees Tg[1,10]

c      do isoil=1,nsoilmx
c       write(str2,'(i2.2)') isoil
c       call write_archive(nid,ntime,'Tg'//str2,'Ground Temperature ',
c     .   'K',2,tsoilS(1,isoil))
c      enddo

! Write soil temperatures tsoil
      call write_archive(nid,ntime,'tsoil','Soil temperature',
     &     'K',-3,tsoilS)

! Write soil thermal inertia
      call write_archive(nid,ntime,'inertiedat',
     &     'Soil thermal inertia',
     &     'J.s-1/2.m-2.K-1',-3,ithS)

! Write (0D) volumetric heat capacity (stored in comsoil.h)
!      call write_archive(nid,ntime,'volcapa',
!     &     'Soil volumetric heat capacity',
!     &     'J.m-3.K-1',0,volcapa)
! Note: no need to write volcapa, it is stored in "controle" table

c-----------------------------------------------------------------
c Ecriture du champs  tankCH4
c-----------------------------------------------------------------
      call write_archive(nid,ntime,'tankCH4',
     &         'Depth of surface methane tank','m',2,tankCH4S)

c-----------------------------------------------------------------
c Ecriture des champs upper_chemistry
c NB : We assume a given order of the nkim=44 chemistry tracers !
c ( H=1, H2=2 ..., C4N2=nkim=44) -> hardcoded in comchem_h.
c-----------------------------------------------------------------

      IF (kim) THEN
        DO iq=1,nkim
          call write_archive(nid,ntime,trim(cnames(iq))//'_up',
     .                trim(cnames(iq))//' in upper atmosphere',
     .                'mol/mol',4,ykim_upS(iq,:,:))
        ENDDO
      ENDIF

c Fin 
c-----------------------------------------------------------------------
      ierr=NF_CLOSE(nid)

      write(*,*) "start2archive: All is well that ends well."

      end 
