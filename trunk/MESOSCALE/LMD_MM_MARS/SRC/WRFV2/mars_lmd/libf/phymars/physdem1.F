      subroutine physdem1(filename,lonfi,latfi,nsoil,nq,
     .                   phystep,day_ini,
     .                   time,tsurf,tsoil,co2ice,emis,q2,qsurf,
     .                   airefi,alb,ith,pzmea,pzstd,pzsig,pzgam,pzthe)

      IMPLICIT none
c-------------------------------------------------------------
C Author : L. Fairhead
C Date   : 01/10/1999
C Objet  : Ecriture des etats initiaux physiques
c-------------------------------------------------------------
c
c
c
c
c
#include "dimensions.h"
#include "paramet.h"
c-----------------------------------------------------------------------
#include "comvert.h"
#include "comgeom2.h"
#include "control.h"
#include "comdissnew.h"
#include "logic.h"
#include "ener.h"
#include "netcdf.inc"
#include "dimphys.h"
c
      INTEGER nid,iq
      INTEGER, parameter :: ivap=1
      REAL, parameter :: qsolmax= 150.0
      character (len=*) :: filename
      character (len=7) :: str7

      REAL day_ini
      INTEGER nsoil,nq
      integer ierr,idim1,idim2,idim3,idim4,idim5,nvarid

c
      REAL phystep,time
      REAL latfi(ngridmx), lonfi(ngridmx)
      REAL champhys(ngridmx)
      REAL tsurf(ngridmx)
      INTEGER length
      PARAMETER (length=100)
      REAL tab_cntrl(length)

c

ccccccccccccccccccccc g95 errors
ccc probably looking for functions, 
ccc but disappointed when finding subroutines 
ccc
ccc obsolete
ccc
ccc      
c      EXTERNAL defrun_new,iniconst,geopot,inigeom,massdair,pression
c      EXTERNAL exner_hyb , SSUM
c
ccccccccccccccccccccc g95 errors

#include "serre.h"
#include "clesph0.h"
#include "fxyprim.h"
#include "comgeomfi.h"
#include "surfdat.h"
#include "planete.h"
#include "dimradmars.h"
#include "yomaer.h"
#include "comcstfi.h"

      real co2ice(ngridmx),tsoil(ngridmx,nsoil),emis(ngridmx)
      real q2(ngridmx, llm+1),qsurf(ngridmx,nq)
      real airefi(ngridmx)
      real alb(ngridmx),ith(ngridmx)
      real pzmea(ngridmx),pzstd(ngridmx)
      real pzsig(ngridmx),pzgam(ngridmx),pzthe(ngridmx)
      integer ig
      INTEGER lnblnk
      EXTERNAL lnblnk

c-----------------------------------------------------------------------

      CALL SCOPY(ngridmx,airefi,1,area,1)
      DO ig=1,ngridmx
         albedodat(ig)=alb(ig)
         inertiedat(ig)=ith(ig)
         zmea(ig)=pzmea(ig)
         zstd(ig)=pzstd(ig)
         zsig(ig)=pzsig(ig)
         zgam(ig)=pzgam(ig)
         zthe(ig)=pzthe(ig)
      ENDDO
c
c  stockage sur le fichier Physique:
c
      ierr = NF_CREATE(trim((filename)),NF_CLOBBER, nid)
      IF (ierr.NE.NF_NOERR) THEN
        WRITE(6,*)' Problem creating restartfi.nc'
        WRITE(6,*)' ierr = ', ierr
        CALL ABORT
      ENDIF
c
      ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 26,
     .                       "Fichier demarrage physique")
c
      ierr = NF_DEF_DIM (nid,"index",length,idim1)
      ierr = NF_DEF_DIM (nid,"physical_points",ngridmx,idim2)
      ierr = NF_DEF_DIM (nid,"subsurface_layers",nsoil,idim3)
      ierr = NF_DEF_DIM (nid,"nlayer+1",llm+1,idim4)
      ierr = NF_DEF_DIM (nid,"number_of_advected_fields",nq,idim5)
c
      ierr = NF_ENDDEF(nid)
c
      DO ierr = 1, length
         tab_cntrl(ierr) = 0.0
      ENDDO

      write(*,*) "ngridmx: ",ngridmx

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c pour la DOCUMENTATION    (fichier io/maj/fi_cntl)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Info sur la grille physique
      tab_cntrl(1) = float(ngridmx)  ! nombre de points de la grille physique
      tab_cntrl(2) = float(nlayermx) ! nombre de couches
      tab_cntrl(3) = day_ini + int(time)         ! jour initial
      tab_cntrl(4) = time -int(time)            ! heure initiale 0

c Info sur la Planete Mars pour la dynamique et la physique
      tab_cntrl(5) = rad      ! rayon de mars(m) ~3397200
      tab_cntrl(6) = omeg     ! vitesse de rotation (rad.s-1)
      tab_cntrl(7) = g        ! gravite (m.s-2) ~3.72
      tab_cntrl(8) = mugaz    ! Masse molaire de l''atm (g.mol-1) ~43.49
      tab_cntrl(9) = rcp      !  = r/cp  ~0.256793 (=kappa dans dynamique)
      tab_cntrl(10) = daysec  ! duree du sol (s)  ~88775

      tab_cntrl(11) = phystep  ! pas de temps de la physique 
      tab_cntrl(12) = 0.
      tab_cntrl(13) = 0.

c Info sur la Planete Mars pour la physique uniquement
      tab_cntrl(14) = year_day  ! duree de l''annee (sols) ~668.6
      tab_cntrl(15) = periheli  ! dist.min. soleil-mars (Mkm) ~206.66
      tab_cntrl(16) = aphelie   ! dist.max. soleil-mars (Mkm) ~249.22
      tab_cntrl(17) = peri_day  ! date du perihelie (sols depuis printemps)
      tab_cntrl(18) = obliquit  ! Obliquite de la planete (deg) ~23.98

c Couche limite et Turbulence
      tab_cntrl(19) = z0        ! surface roughness (m) ~0.01
      tab_cntrl(20) = lmixmin   ! longueur de melange ~100
      tab_cntrl(21) = emin_turb ! energie minimale ~1.e-8

c propriete optiques des calottes et emissivite du sol
      tab_cntrl(22) = albedice(1)  ! Albedo calotte nord ~0.5
      tab_cntrl(23) = albedice(2)  ! Albedo calotte sud ~0.5
      tab_cntrl(24) = emisice(1)   ! Emissivite calotte nord ~0.95
      tab_cntrl(25) = emisice(2)   ! Emissivite calotte sud ~0.95
      tab_cntrl(26) = emissiv      ! Emissivite du sol martien ~.95
      tab_cntrl(31) = iceradius(1) ! mean scat radius of CO2 snow (north)
      tab_cntrl(32) = iceradius(2) ! mean scat radius of CO2 snow (south)
      tab_cntrl(33) = dtemisice(1) ! time scale for snow metamorphism (north)
      tab_cntrl(34) = dtemisice(2) ! time scale for snow metamorphism (south)

c Proprietes des poussiere aerosol
      tab_cntrl(27) = tauvis      ! profondeur optique visible moyenne

      tab_cntrl(28) = 0. 
      tab_cntrl(29) = 0.
      tab_cntrl(30) = 0.

cc   ***    new_oliq   (  commentaires de L. LI dans routine physique ) 
cc   ***  ok_orodr  et ok_orolf   si on appelle l'orographie      ****
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "controle", NF_DOUBLE, 1, idim1,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "controle", NF_FLOAT, 1, idim1,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 22,
     .                        "Parametres de controle")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "longitude", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "longitude", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 32,
     .                        "Longitudes de la grille physique")
      ierr = NF_ENDDEF(nid)

#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,lonfi)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,lonfi)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "latitude", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "latitude", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 31,
     .                        "Latitudes de la grille physique")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,latfi)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,latfi)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "area", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "area", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 16,
     .                        "Aire des mailles")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,area)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,area)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "phisfi", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "phisfi", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                        "Geopotentiel au sol")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,phisfi)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,phisfi)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "albedodat", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "albedodat", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 16,
     .                        "Albedo du sol nu")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,albedodat)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,albedodat)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "inertiedat", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "inertiedat", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 24,
     .                        "Inertie thermique du sol")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,inertiedat)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,inertiedat)
#endif
c
c   fichier pour les programmes de Francois Lott

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZMEA", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZMEA", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12,
     .                        "Relief moyen")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zmea)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zmea)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZSTD", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZSTD", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 18,
     .                        "Ecartype du relief")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zstd)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zstd)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZSIG", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZSIG", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                        "Relief: parametre sigma")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zsig)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zsig)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZGAM", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZGAM", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                        "Relief: parametre gamma")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zgam)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zgam)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZTHE", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZTHE", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                        "Relief: parametre theta")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zthe)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zthe)
#endif

c Ecriture des champs physiques
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "co2ice", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "co2ice", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 13,
     .                        "CO2 ice cover")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,co2ice)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,co2ice)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "tsurf", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "tsurf", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                        "Surface temperature")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tsurf)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tsurf)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"tsoil",NF_DOUBLE,2,(/idim2,idim3/),nvarid)
#else
      ierr = NF_DEF_VAR (nid,"tsoil",NF_FLOAT,2,(/idim2,idim3/),nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 16,
     .                        "Soil temperature")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tsoil)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tsoil)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "emis", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "emis", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 18,
     .                        "Surface emissivity")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,emis)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,emis)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "q2", NF_DOUBLE, 2, 
     .  (/idim2,idim4/),nvarid)
c****WRF: ligne trop longue ....
#else
      ierr = NF_DEF_VAR (nid, "q2", NF_FLOAT, 2,(/idim2,idim4/),nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 17,
     .                        "pbl wind variance")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,q2)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,q2)
#endif
c
      IF(nq.GE.1) THEN
         DO iq=1,nq
            str7(1:5)='qsurf'
            WRITE(str7(6:7),'(i2.2)') iq
            ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
            ierr = NF_DEF_VAR (nid,str7,NF_DOUBLE,1,idim2,nvarid)
#else
            ierr = NF_DEF_VAR (nid,str7,NF_FLOAT,1,idim2,nvarid)
#endif
            ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 17,
     .                        "tracer on surface")
            ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
            ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,qsurf(1,iq))
#else
            ierr = NF_PUT_VAR_REAL (nid,nvarid,qsurf(1,iq))
#endif
         ENDDO
      ENDIF
c
      ierr = NF_CLOSE(nid)

      RETURN

      END
