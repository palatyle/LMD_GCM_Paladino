c=======================================================================
      subroutine ini_archive(nid,idayref,phis,ith,tab_cntrl_fi)
c=======================================================================
c
c
c   Date:    01/1997
c   ----
c
c   Objet:  ecriture de l'entete du fichier "start_archive"
c   -----
c
c	 Proche de iniwrite.F
c
c	 On ajoute dans le tableau "tab_cntrl" (dynamique), a partir de 51, 
c	 les valeurs de tab_cntrl_fi (les 38 parametres de controle physiques
c	 du RUN + ptotal et cotoicetotal)
c
c			tab_cntrl(50+l)=tab_cntrl_fi(l)
c
c   Arguments:
c   ---------
c
c	Inputs:
c   ------
c
c       nid            unite logique du fichier "start_archive"
c       idayref        Valeur du jour initial a mettre dans
c                      l'entete du fichier "start_archive"
c       phis           geopotentiel au sol
c       ith            soil thermal inertia
c       tab_cntrl_fi   tableau des param physiques
c

c=======================================================================
 
      use comsoil_h, only: nsoilmx, mlayer
      USE comvert_mod, ONLY: ap,bp,aps,bps,pa,preff,presnivs,pseudoalt
      USE comconst_mod, ONLY: daysec,dtvr,rad,omeg,g,cpp,kappa,pi
      USE logic_mod, ONLY: fxyhypb,ysinus
      USE serre_mod, ONLY: clon,clat,grossismx,grossismy,dzoomx,dzoomy
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0
      implicit none

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "netcdf.inc"

c-----------------------------------------------------------------------
c   Declarations
c-----------------------------------------------------------------------

c   Local:
c   ------
      INTEGER	length,l
      parameter (length = 100)
      REAL		tab_cntrl(length) ! tableau des parametres du run
      INTEGER	loop
      INTEGER	ierr, setvdim, putvdim, putdat, setname,cluvdb
      INTEGER	setdim
      INTEGER	ind1,indlast

c   Arguments:
c   ----------
      INTEGER*4	idayref
      REAL		phis(ip1jmp1)
      real ith(ip1jmp1,nsoilmx)
      REAL		tab_cntrl_fi(length)

!Mars --------Ajouts-----------
c   Variables locales pour NetCDF:
c
      INTEGER dims2(2), dims3(3) !, dims4(4)
      INTEGER idim_index
      INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
      INTEGER idim_llmp1,idim_llm
      INTEGER idim_tim
      INTEGER idim_nsoilmx ! "subsurface_layers" dimension ID #
      INTEGER nid,nvarid
      real sig_s(llm),s(llm)

      pi  = 2. * ASIN(1.)


c-----------------------------------------------------------------------
c   Remplissage du tableau des parametres de controle du RUN  (dynamique)
c-----------------------------------------------------------------------

      DO l=1,length
         tab_cntrl(l)=0.
      ENDDO

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      tab_cntrl(1)  = REAL(iim) ! nombre de points en longitude
      tab_cntrl(2)  = REAL(jjm) ! nombre de points en latitude
      tab_cntrl(3)  = REAL(llm) ! nombre de couches
      tab_cntrl(4)  = REAL(idayref) ! jour 0
      tab_cntrl(5)  = rad ! rayon de mars(m) ~3397200
      tab_cntrl(6)  = omeg ! vitesse de rotation (rad.s-1)
      tab_cntrl(7)  = g   ! gravite (m.s-2) ~3.72
      tab_cntrl(8)  = cpp 
      tab_cntrl(8)  = 43.49 !mars temporaire Masse molaire de l''atm (g.mol-1) ~43.49
      tab_cntrl(9)  = kappa ! = r/cp  ~0.256793 (=rcp dans physique)
      tab_cntrl(10) = daysec ! duree du sol (s)  ~88775
      tab_cntrl(11) = dtvr ! pas de temps de la dynamique (s)
      tab_cntrl(12) = etot0 ! energie totale    !
      tab_cntrl(13) = ptot0 ! pression totalei   !    variables
      tab_cntrl(14) = ztot0 ! enstrophie totale   !  de controle
      tab_cntrl(15) = stot0 ! enthalpie totale   !    globales
      tab_cntrl(16) = ang0 ! moment cinetique  !
      tab_cntrl(17) = pa
      tab_cntrl(18) = preff

c    .....    parametres  pour le zoom      ......   

      tab_cntrl(19)  = clon ! longitude en degres du centre du zoom
      tab_cntrl(20)  = clat ! latitude en degres du centre du zoom
      tab_cntrl(21)  = grossismx ! facteur de grossissement du zoom,selon longitude
      tab_cntrl(22)  = grossismy ! facteur de grossissement du zoom ,selon latitude

      IF ( fxyhypb )   THEN
       tab_cntrl(23) = 1.
       tab_cntrl(24) = dzoomx ! extension en longitude  de la zone du zoom
       tab_cntrl(25) = dzoomy ! extension en latitude  de la zone du zoom
      ELSE
       tab_cntrl(23) = 0.
       tab_cntrl(24) = dzoomx ! extension en longitude  de la zone du zoom
       tab_cntrl(25) = dzoomy ! extension en latitude  de la zone du zoom
       tab_cntrl(26) = 0.
       IF ( ysinus)  tab_cntrl(26) = 1.
      ENDIF

c-----------------------------------------------------------------------
c   Copie du tableau des parametres de controle du RUN  (physique)
c		dans le tableau dynamique
c-----------------------------------------------------------------------

      DO l=1,50
         tab_cntrl(50+l)=tab_cntrl_fi(l)
      ENDDO

c=======================================================================
c	Ecriture NetCDF de l''entete du fichier "start_archive"
c=======================================================================

c
c Preciser quelques attributs globaux:
c
      ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 21,
     &                       "Fichier start_archive")
c
c Definir les dimensions du fichiers:
c
c     CHAMPS AJOUTES POUR LA VISUALISATION T,ps, etc... avec Grads ou ferret:
      ierr = NF_DEF_DIM (nid, "latitude", jjp1, idim_rlatu)
      ierr = NF_DEF_DIM (nid, "longitude", iip1, idim_rlonv)
      ierr = NF_DEF_DIM (nid, "altitude", llm, idim_llm)
      ierr = NF_DEF_DIM (nid,"subsurface_layers",nsoilmx,idim_nsoilmx)

      ierr = NF_DEF_DIM (nid,"index", length, idim_index)
      ierr = NF_DEF_DIM (nid,"rlonu", iip1, idim_rlonu)
      ierr = NF_DEF_DIM (nid,"rlatv", jjm, idim_rlatv)
      ierr = NF_DEF_DIM (nid,"interlayer", llmp1, idim_llmp1)
      ierr = NF_DEF_DIM (nid,"Time", NF_UNLIMITED, idim_tim)

c
      ierr = NF_ENDDEF(nid) ! sortir du mode de definition

c-----------------------------------------------------------------------
c  Ecriture du tableau des parametres du run
c-----------------------------------------------------------------------

      call def_var(nid,"Time","Time","days since 00:00:00",1,
     .            idim_tim,nvarid,ierr)

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"controle",NF_DOUBLE,1,idim_index,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"controle",NF_FLOAT,1,idim_index,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Parametres de controle")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif

c-----------------------------------------------------------------------
c  Ecriture des longitudes et latitudes
c-----------------------------------------------------------------------

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlonu",NF_DOUBLE,1,idim_rlonu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlonu",NF_FLOAT,1,idim_rlonu,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                       "Longitudes des points U")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonu)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlatu",NF_DOUBLE,1,idim_rlatu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlatu",NF_FLOAT,1,idim_rlatu,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Latitudes des points U")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlonv",NF_DOUBLE,1,idim_rlonv,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlonv",NF_FLOAT,1,idim_rlonv,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                       "Longitudes des points V")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlatv",NF_DOUBLE,1,idim_rlatv,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlatv",NF_FLOAT,1,idim_rlatv,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Latitudes des points V")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatv)
#endif

c-----------------------------------------------------------------------
c  Ecriture des niveaux verticaux
c-----------------------------------------------------------------------

c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ap",NF_DOUBLE,1,idim_llmp1,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ap",NF_FLOAT,1,idim_llmp1,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 32,
     .                       "Coef A: niveaux pression hybride")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ap)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"bp",NF_DOUBLE,1,idim_llmp1,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"bp",NF_FLOAT,1,idim_llmp1,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 35,
     .                       "Coefficient B niveaux sigma hybride")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bp)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bp)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"aps",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"aps",NF_FLOAT,1,idim_llm,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 36,
     .      "Coef AS: hybrid pressure in midlayers")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aps)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"bps",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"bps",NF_FLOAT,1,idim_llm,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 30,
     .      "Coef BS: hybrid sigma midlayers")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bps)
#endif
c
c ----------------------

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"presnivs",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"presnivs",NF_FLOAT,1,idim_llm,nvarid)
#endif
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,presnivs)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,presnivs)
#endif
c ------------------------------------------------------------------
c  Variable uniquement pour visualisation avec Grads ou Ferret
c ------------------------------------------------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"latitude",NF_DOUBLE,1,idim_rlatu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"latitude",NF_FLOAT,1,idim_rlatu,nvarid)
#endif
      ierr =NF_PUT_ATT_TEXT(nid,nvarid,'units',13,"degrees_north")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 14,
     .      "North latitude")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu/pi*180)
#endif
c----------------------
       ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr =NF_DEF_VAR(nid,"longitude", NF_DOUBLE, 1, idim_rlonv,nvarid)
#else
      ierr = NF_DEF_VAR(nid,"longitude", NF_FLOAT, 1, idim_rlonv,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 14,
     .      "East longitude")
      ierr = NF_PUT_ATT_TEXT(nid,nvarid,'units',12,"degrees_east")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv/pi*180)
#endif
c--------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "altitude", NF_DOUBLE, 1,
     .       idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "altitude", NF_FLOAT, 1,
     .       idim_llm,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name",10,"pseudo-alt")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,'units',2,"km")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,'positive',2,"up")

      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,pseudoalt)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,pseudoalt)
#endif

!-------------------------------
! (soil) depth variable mlayer() (known from comsoil.h)
!-------------------------------
      ierr=NF_REDEF (nid) ! Enter NetCDF (re-)define mode
      ! define variable
#ifdef NC_DOUBLE
      ierr=NF_DEF_VAR(nid,"soildepth",NF_DOUBLE,1,idim_nsoilmx,nvarid)
#else
      ierr=NF_DEF_VAR(nid,"soildepth",NF_FLOAT,1,idim_nsoilmx,nvarid)
#endif
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 20,
     .                        "Soil mid-layer depth")
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"units",1,"m")
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"positive",4,"down")
      ierr=NF_ENDDEF(nid) ! Leave NetCDF define mode
      ! write variable
#ifdef NC_DOUBLE
      ierr=NF_PUT_VAR_DOUBLE (nid,nvarid,mlayer)
#else
      ierr=NF_PUT_VAR_REAL (nid,nvarid,mlayer)
#endif

!---------------------
! soil thermal inertia
!---------------------
      ierr=NF_REDEF (nid) ! Enter NetCDF (re-)define mode
      dims3(1)=idim_rlonv
      dims3(2)=idim_rlatu
      dims3(3)=idim_nsoilmx
      ! define variable
#ifdef NC_DOUBLE
      ierr=NF_DEF_VAR(nid,"inertiedat",NF_DOUBLE,3,dims3,nvarid)
#else
      ierr=NF_DEF_VAR(nid,"inertiedat",NF_FLOAT,3,dims3,nvarid)
#endif
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 20,
     &                        "Soil thermal inertia")
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"units",15,
     &                        "J.s-1/2.m-2.K-1")
      ierr=NF_ENDDEF(nid) ! Leave NetCDF define mode
      ! write variable
#ifdef NC_DOUBLE
      ierr=NF_PUT_VAR_DOUBLE (nid,nvarid,ith)
#else
      ierr=NF_PUT_VAR_REAL (nid,nvarid,ith)
#endif

c-----------------------------------------------------------------------
c  Ecriture aire et coefficients de passage cov. <-> contra. <--> naturel
c-----------------------------------------------------------------------

      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonu
      dims2(2) = idim_rlatu
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"cu",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"cu",NF_FLOAT,2,dims2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29,
     .                       "Coefficient de passage pour U")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cu)
#endif
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatv
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"cv",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"cv",NF_FLOAT,2,dims2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29,
     .                       "Coefficient de passage pour V")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cv)
#endif
c
c Aire de chaque maille:
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatu
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"aire",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"aire",NF_FLOAT,2,dims2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Aires de chaque maille")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aire)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aire)
#endif

c-----------------------------------------------------------------------
c  Ecriture du geopentiel au sol
c-----------------------------------------------------------------------

      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatu
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"phisinit",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"phisinit",NF_FLOAT,2,dims2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                       "Geopotentiel au sol")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,phis)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,phis)
#endif

      PRINT*,'iim,jjm,llm,idayref',iim,jjm,llm,idayref
      PRINT*,'rad,omeg,g,mugaz,kappa',
     s rad,omeg,g,43.49,kappa !mars temporaire (ecrire mugaz ensuite)
      PRINT*,'daysec,dtvr',daysec,dtvr

      RETURN
      END
