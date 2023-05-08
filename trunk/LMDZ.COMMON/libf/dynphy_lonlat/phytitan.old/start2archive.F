c=======================================================================
      PROGRAM start2archive
c=======================================================================
c
c
c   Date:    01/1997
c   ----
c
c   Version Venus: 09/2007
c           Titan: 02/2009
c
c   Objet:   Passage des  fichiers netcdf d etat initial "start" et
c   -----    "startphy" a un fichier netcdf unique "start_archive" 
c
c=======================================================================

      USE filtreg_mod
      USE infotrac
      USE control_mod
      use cpdet_mod, only: tpot2t,ini_cpdet
      use exner_hyb_m, only: exner_hyb
      use exner_milieu_m, only: exner_milieu
      USE comconst_mod
      USE comvert_mod, ONLY: ap,bp,presnivs,pa,preff,nivsigs,nivsig,
     .                       aps,bps,scaleheight,pseudoalt,
     .                       disvert_type,pressure_exner
      USE logic_mod, ONLY: iflag_trac

      implicit none

#include "dimensions.h"
#include "paramet.h"
#include "comdissnew.h"
#include "comgeom.h"
#include "description.h"
#include "dimsoil.h"
#include "netcdf.inc"

c-----------------------------------------------------------------------
c   Declarations
c-----------------------------------------------------------------------

c variables dynamiques du GCM
c -----------------------------
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) ! vents covariants
      REAL teta(ip1jmp1,llm)                    ! temperature potentielle 
      REAL, ALLOCATABLE, DIMENSION(:,:,:):: q! champs advectes
      REAL pks(ip1jmp1)                      ! exner (f pour filtre)
      REAL pk(ip1jmp1,llm)
      REAL pkf(ip1jmp1,llm)
      REAL alpha(iip1,jjp1,llm),beta(iip1,jjp1,llm)
      REAL phis(ip1jmp1)                     ! geopotentiel au sol
      REAL masse(ip1jmp1,llm)                ! masse de l'atmosphere
      REAL ps(ip1jmp1)                       ! pression au sol
      REAL p3d(iip1, jjp1, llm+1)            ! pression aux interfaces
      
c Variable Physiques (grille physique)
c ------------------------------------
      integer ngridmx,nlayermx
      parameter (ngridmx=(2+(jjm-1)*iim - 1/jjm))
      parameter (nlayermx=llm)

      real rlat(ngridmx),rlon(ngridmx)
      REAL tsurf(ngridmx),tsoil(ngridmx,nsoilmx)
      REAL albe(ngridmx),radsol(ngridmx),sollw(ngridmx)
      real solsw(ngridmx),dlw(ngridmx)
      REAL zmea(ngridmx), zstd(ngridmx)
      REAL zsig(ngridmx), zgam(ngridmx), zthe(ngridmx)
      REAL zpic(ngridmx), zval(ngridmx)
      
      INTEGER start,length
      PARAMETER (length = 100)
      REAL tab_cntrl_fi(length) ! tableau des parametres de startfi
      REAL tab_cntrl_dyn(length) ! tableau des parametres de start
      INTEGER*4 day_ini_fi

c Variable naturelle / grille scalaire
c ------------------------------------
      REAL T(ip1jmp1,llm),us(ip1jmp1,llm),vs(ip1jmp1,llm)
      REAL tsurfS(ip1jmp1),tsoilS(ip1jmp1,nsoilmx)
      real rlatS(ip1jmp1),rlonS(ip1jmp1)
      real albeS(ip1jmp1),radsolS(ip1jmp1),sollwS(ip1jmp1)
      real solswS(ip1jmp1),dlwS(ip1jmp1)
      real zmeaS(ip1jmp1),zstdS(ip1jmp1),zsigS(ip1jmp1)
      real zgamS(ip1jmp1),ztheS(ip1jmp1),zpicS(ip1jmp1)
      real zvalS(ip1jmp1)
      
c Variables intermediaires : vent naturel, mais pas coord scalaire
c----------------------------------------------------------------
      REAL vn(ip1jm,llm),un(ip1jmp1,llm)

c Autres  variables
c -----------------
      REAL ptotal

      CHARACTER*2 str2

      INTEGER ij, l,i,j,isoil,iq
      character*80      fichnom
      integer :: ierr

c Netcdf
c-------
      integer varid,dimid
      INTEGER nid

c-----------------------------------------------------------------------
c   Initialisations 
c-----------------------------------------------------------------------

c VENUS/TITAN

        iflag_trac = 1
c-----------------------------------------------------------------------
c   Initialisation des traceurs
c   ---------------------------
c  Choix du nombre de traceurs et du schema pour l advection
c  dans fichier traceur.def, par default ou via INCA
      call infotrac_init

c Allocation de la tableau q : champs advectes   
      allocate(q(ip1jmp1,llm,nqtot))

c=======================================================================
c Lecture des donnees
c=======================================================================

      fichnom = 'start.nc'
      CALL readstart(fichnom,nqtot,vcov,ucov,teta,q,masse,
     .       ps,phis,tab_cntrl_dyn)

      fichnom = 'startphy.nc'
      CALL readstartphy(fichnom,
     .       rlat,rlon,tsurf,tsoil,
     .       albe, solsw, sollw,
     .       dlw,radsol,
     .       zmea,zstd,zsig,zgam,zthe,zpic,zval,
     .       tab_cntrl_fi)

c-----------------------------------------------------------------------
c   Initialisations 
c-----------------------------------------------------------------------

      CALL conf_gcm( 99, .TRUE. )
      call iniconst
      call inigeom
      call inifilr
      call ini_cpdet

      CALL pression(ip1jmp1, ap, bp, ps, p3d)
         if (disvert_type==1) then
           CALL exner_hyb(  ip1jmp1, ps, p3d, pks, pk, pkf )
         else ! we assume that we are in the disvert_type==2 case
           CALL exner_milieu( ip1jmp1, ps, p3d, pks, pk, pkf )
         endif

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
! ADAPTATION GCM POUR CP(T)

      call tpot2t(ip1jmp1*llm,teta,T,pk)

c-----------------------------------------------------------------------
c Variable physique 
c-----------------------------------------------------------------------
c
c tsurf --> tsurfS
c et autres...
c
c-----------------------------------------------------------------------

      call gr_fi_dyn(1,ngridmx,iip1,jjp1,tsurf,tsurfS)
      call gr_fi_dyn(nsoilmx,ngridmx,iip1,jjp1,tsoil,tsoilS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,rlat,rlatS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,rlon,rlonS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,albe,albeS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,radsol,radsolS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,sollw,sollwS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,solsw,solswS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,dlw,dlwS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,zmea,zmeaS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,zstd,zstdS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,zsig,zsigS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,zgam,zgamS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,zthe,ztheS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,zpic,zpicS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,zval,zvalS)

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
c Passage de "ptotal" par tab_cntrl_fi
c-----------------------------------------------------------------------

      tab_cntrl_fi(length) = ptotal

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
c On met dans l'entete le tab_cntrl_dyn (1 a length) 
c  On y ajoute les valeurs du tab_cntrl_fi (length+1 a 2*length)
c-----------------------------------------------------------------------

      if (ierr.ne.NF_NOERR) then
         write(*,*)'OK, Could not open file "start_archive.nc"'
         write(*,*)'So let s create a new "start_archive"'
         ierr = NF_CREATE('start_archive.nc', NF_CLOBBER, nid)
         call ini_archive(nid,phis,tab_cntrl_dyn,tab_cntrl_fi)
      else
         write(*,*)'Attention, start_archive.nc existe deja...'
	 call abort
      endif

c-----------------------------------------------------------------------
c Ecriture des champs 
c-----------------------------------------------------------------------

      call write_archive(nid,'u','Vent zonal','m.s-1',3,us)
      call write_archive(nid,'v','Vent merid','m.s-1',3,vs)
      call write_archive(nid,'temp','temperature','K',3,T)
c-----------------------------------------------------------------------
c Ecriture du champs  q  ( q[1,nqtot] )
c-----------------------------------------------------------------------
       do iq=1,nqtot
        write(str2,'(i2.2)') iq
         call write_archive(nid,tname(iq),'tracer','kg/kg',
     .         3,q(1,1,iq))
       end do
c-----------------------------------------------------------------------
      call write_archive(nid,'masse','Masse','kg',3,masse)
      call write_archive(nid,'ps','Psurf','Pa',2,ps)
      call write_archive(nid,'tsurf','surf T','K',2,tsurfS)
c-----------------------------------------------------------------------
c Ecriture du champs  tsoil  ( Tsoil[1,nsoilmx] )
c-----------------------------------------------------------------------
c "tsoil" Temperature au sol definie dans nsoilmx couches dans le sol
c   Les nsoilmx couches sont lues comme nsoilmx champs 
c  nommees Tsoil[1,nsoilmx]
      do isoil=1,nsoilmx
       write(str2,'(i2.2)') isoil
       call write_archive(nid,'Tsoil'//str2,'Ground Temperature ',
     .   'K',2,tsoilS(1,isoil))
      enddo
c-----------------------------------------------------------------------
      call write_archive(nid,'rlat','Latitude','rad',2,rlatS)
      call write_archive(nid,'rlon','Longitude','rad',2,rlonS)
      call write_archive(nid,'albe','Albedo','',2,albeS)
      call write_archive(nid,'radsol',
     .             'Net flux at surface','W m-2',2,radsolS)
      call write_archive(nid,'sollw',
     .             'LW flux at surface','W m-2',2,sollwS)
      call write_archive(nid,'solsw',
     .             'SW flux at surface','W m-2',2,solswS)
      call write_archive(nid,'dlw','LW derive','?',2,dlwS)
      call write_archive(nid,'zmea','param oro sous-maille','m',2,zmeaS)
      call write_archive(nid,'zstd','param oro sous-maille','m',2,zstdS)
      call write_archive(nid,'zsig','param oro sous-maille','m',2,zsigS)
      call write_archive(nid,'zgam','param oro sous-maille','m',2,zgamS)
      call write_archive(nid,'zthe','param oro sous-maille','m',2,ztheS)
      call write_archive(nid,'zpic','param oro sous-maille','m',2,zpicS)
      call write_archive(nid,'zval','param oro sous-maille','m',2,zvalS)

      ierr=NF_CLOSE(nid)

c-----------------------------------------------------------------------
c Fin 
c-----------------------------------------------------------------------

      end 
