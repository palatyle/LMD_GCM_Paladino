c=======================================================================
      subroutine ini_archive(nid,phis,tab_cntrl_dyn,tab_cntrl_fi)
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
c   Arguments:
c   ---------
c
c	Inputs:
c   ------
c
c       nid            unite logique du fichier "start_archive"
c       phis           geopotentiel au sol
c       tab_cntrl_dyn  tableau des param dynamiques
c       tab_cntrl_fi   tableau des param physiques
c

c=======================================================================
 
      USE control_mod
      USE comconst_mod
      USE comvert_mod, ONLY: ap,bp,presnivs,pa,preff,nivsigs,nivsig,
     .                       aps,bps,scaleheight,pseudoalt,
     .                       disvert_type,pressure_exner

      implicit none

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "description.h"
#include "netcdf.inc"

c-----------------------------------------------------------------------
c   Declarations
c-----------------------------------------------------------------------

c   Local:
c   ------
      INTEGER	length,l
      parameter (length = 100)
      REAL	tab_cntrl(2*length) ! tableau des parametres du run
      INTEGER	loop
      INTEGER	ierr, setvdim, putvdim, putdat, setname,cluvdb
      INTEGER	setdim
      INTEGER	ind1,indlast

c   Arguments:
c   ----------
      REAL	phis(ip1jmp1)
      REAL	tab_cntrl_dyn(length)
      REAL	tab_cntrl_fi(length)

!Mars --------Ajouts-----------
c   Variables locales pour NetCDF:
c
      INTEGER dims2(2), dims3(3), dims4(4)
      INTEGER idim_index
      INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
      INTEGER idim_llmp1,idim_llm
      INTEGER idim_tim
      INTEGER nid,nvarid
      real sig_s(llm),s(llm)

      pi  = 2. * ASIN(1.)


c-----------------------------------------------------------------------
c   Remplissage du tableau des parametres de controle du RUN  
c-----------------------------------------------------------------------

      DO l=1,length
         tab_cntrl(l)       = tab_cntrl_dyn(l)
         tab_cntrl(length+l)= tab_cntrl_fi(l)
      ENDDO

c=======================================================================
c	Ecriture NetCDF de l''entete du fichier "start_archive"
c=======================================================================

c
c Preciser quelques attributs globaux:
c
      ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 27,
     .                       "Fichier start_archive")
c
c Definir les dimensions du fichiers:
c
c     CHAMPS AJOUTES POUR LA VISUALISATION T,ps, etc... avec Grads ou ferret:
      ierr = NF_DEF_DIM (nid,"index", 2*length, idim_index)
      ierr = NF_DEF_DIM (nid,"rlonu", iip1, idim_rlonu)
      ierr = NF_DEF_DIM (nid,"rlatu", jjp1, idim_rlatu)
      ierr = NF_DEF_DIM (nid,"rlonv", iip1, idim_rlonv)
      ierr = NF_DEF_DIM (nid,"rlatv", jjm, idim_rlatv)
      ierr = NF_DEF_DIM (nid,"sigs",  llm, idim_llm)
      ierr = NF_DEF_DIM (nid,"sig", llmp1, idim_llmp1)
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

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"nivsigs",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"nivsigs",NF_FLOAT,1,idim_llm,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 28,
     .                       "Numero naturel des couches s")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,nivsigs)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsigs)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"nivsig",NF_DOUBLE,1,idim_llmp1,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"nivsig",NF_FLOAT,1,idim_llmp1,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 28,
     .                       "Numero naturel des couches sigma")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,nivsig)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsig)
#endif
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

      RETURN
      END
