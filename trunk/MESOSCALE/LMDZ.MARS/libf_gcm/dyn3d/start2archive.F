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

      implicit none

#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comdissip.h"
#include "comvert.h"
#include "comgeom.h"
#include "logic.h"
#include "temps.h"
#include "control.h"
#include "ener.h"
#include "description.h"

#include "dimphys.h"
#include "netcdf.inc"

c-----------------------------------------------------------------------
c   Declarations
c-----------------------------------------------------------------------

c variables dynamiques du GCM
c -----------------------------
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) ! vents covariants
      REAL teta(ip1jmp1,llm)                    ! temperature potentielle 
      REAL q(ip1jmp1,llm,nqmx)               ! champs advectes
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
      REAL tsurf(ngridmx),tsoil(ngridmx,nsoilmx),co2ice(ngridmx)
      REAL q2(ngridmx,nlayermx+1),qsurf(ngridmx,nqmx)
      REAL emis(ngridmx)
      INTEGER start,length
      PARAMETER (length = 100)
      REAL tab_cntrl_fi(length) ! tableau des parametres de startfi
      INTEGER*4 day_ini_fi

c Variable naturelle / grille scalaire
c ------------------------------------
      REAL T(ip1jmp1,llm),us(ip1jmp1,llm),vs(ip1jmp1,llm)
      REAL tsurfS(ip1jmp1),tsoilS(ip1jmp1,nsoilmx),co2iceS(ip1jmp1)
      REAL q2S(ip1jmp1,llm+1),qsurfS(ip1jmp1,nqmx)
      REAL emisS(ip1jmp1)

c Variables intermediaires : vent naturel, mais pas coord scalaire
c----------------------------------------------------------------
      REAL vn(ip1jm,llm),un(ip1jmp1,llm)

c Autres  variables
c -----------------
      LOGICAL startdrs
      INTEGER Lmodif

      REAL ptotal, co2icetotal
      REAL timedyn,timefi !fraction du jour dans start, startfi
      REAL date

      CHARACTER*2 str2
      CHARACTER*80 fichier 
      data  fichier /'startfi'/

      INTEGER ij, l,i,j,isoil,iq
      character*80      fichnom
      integer :: ierr,ntime

c Netcdf
c-------
      integer varid,dimid,timelen 
      INTEGER nid,nid1

c-----------------------------------------------------------------------
c   Initialisations 
c-----------------------------------------------------------------------

      grireg   = .TRUE.

c=======================================================================
c Lecture des donnees
c=======================================================================

      fichnom = 'start.nc'
      CALL dynetat0(fichnom,nqmx,vcov,ucov,teta,q,masse,
     .       ps,phis,timedyn)

      fichnom = 'startfi.nc'
      Lmodif=0

      CALL phyetat0 (fichnom,0,Lmodif,nsoilmx,nqmx,day_ini_fi,timefi,
     .      tsurf,tsoil,emis,q2,qsurf,co2ice)

       ierr = NF_OPEN (fichnom, NF_NOWRITE,nid1)
       IF (ierr.NE.NF_NOERR) THEN
         write(6,*)' Pb d''ouverture du fichier'//fichnom
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
       do while (day_ini.ge.669)
          day_ini=day_ini-669
       enddo
c *****************************************************************

c-----------------------------------------------------------------------
c   Initialisations 
c-----------------------------------------------------------------------

      CALL defrun_new(99, .FALSE. )
      call iniconst
      call inigeom
      call inifilr
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
c co2ice --> co2iceS
c tsoil --> tsoilS
c emis --> emisS
c q2 --> q2S
c qsurf --> qsurfS
c
c-----------------------------------------------------------------------

      call gr_fi_dyn(1,ngridmx,iip1,jjp1,tsurf,tsurfS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,co2ice,co2iceS)
      call gr_fi_dyn(nsoilmx,ngridmx,iip1,jjp1,tsoil,tsoilS)
      call gr_fi_dyn(1,ngridmx,iip1,jjp1,emis,emisS)
      call gr_fi_dyn(llm+1,ngridmx,iip1,jjp1,q2,q2S)
      call gr_fi_dyn(nqmx,ngridmx,iip1,jjp1,qsurf,qsurfS)

c=======================================================================
c Info pour controler
c=======================================================================

      ptotal =  0.
      co2icetotal = 0.
      DO j=1,jjp1
         DO i=1,iim
           ptotal=ptotal+aire(i+(iim+1)*(j-1))*ps(i+(iim+1)*(j-1))/g
           co2icetotal = co2icetotal + 
     &            co2iceS(i+(iim+1)*(j-1))*aire(i+(iim+1)*(j-1))
         ENDDO
      ENDDO
      write(*,*)'Ancienne grille : masse de l''atm :',ptotal
      write(*,*)'Ancienne grille : masse de la glace CO2 :',co2icetotal

c-----------------------------------------------------------------------
c Passage de "ptotal" et "co2icetotal" par tab_cntrl_fi
c-----------------------------------------------------------------------

      tab_cntrl_fi(49) = ptotal
      tab_cntrl_fi(50) = co2icetotal

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
         call ini_archive(nid,day_ini,phis,tab_cntrl_fi)
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
c Ecriture des champs  (co2ice,emis,ps,Tsurf,T,u,v,q2,q,qsurf)
c-----------------------------------------------------------------------
c ATTENTION: q2 a une couche de plus!!!!
c    Pour creer un fichier netcdf lisible par grads,
c    On passe donc une des couches de q2 a part
c    comme une variable 2D (la couche au sol: "q2surf")
c    Les lmm autres couches sont nommees "q2atm" (3D) 
c-----------------------------------------------------------------------

      call write_archive(nid,ntime,'co2ice','couche de glace co2',
     &  'kg/m2',2,co2iceS)
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
c Ecriture du champs  tsoil  ( Tg[1,10] )
c-----------------------------------------------------------------------
c "tsoil" Temperature au sol definie dans 10 couches dans le sol
c   Les 10 couches sont lues comme 10 champs 
c  nommees Tg[1,10]

      do isoil=1,nsoilmx
       write(str2,'(i2.2)') isoil
       call write_archive(nid,ntime,'Tg'//str2,'Ground Temperature ',
     .   'K',2,tsoilS(1,isoil))
      enddo
c-----------------------------------------------------------------------
c Ecriture du champs  q  ( q[1,nqmx] )
c-----------------------------------------------------------------------
      do iq=1,nqmx
       write(str2,'(i2.2)') iq
        call write_archive(nid,ntime,'q'//str2,'tracer','kg/kg',
     .         3,q(1,1,iq))
      end do
c-----------------------------------------------------------------------
c Ecriture du champs  qsurf  ( qsurf[1,nqmx] )
c-----------------------------------------------------------------------
      do iq=1,nqmx
       write(str2,'(i2.2)') iq
       call write_archive(nid,ntime,'qsurf'//str2,'Tracer on surface',
     $  'kg.m-2',2,qsurfS(1,iq))
      enddo

      ierr=NF_CLOSE(nid)

c-----------------------------------------------------------------------
c Fin 
c-----------------------------------------------------------------------

      end 
