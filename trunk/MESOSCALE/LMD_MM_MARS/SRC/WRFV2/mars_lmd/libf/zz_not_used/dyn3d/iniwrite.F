      SUBROUTINE iniwrite(nid,idayref,phis)
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  L. Fairhead  ,  P. Le Van, Y. Wanherdrick, F. Forget
c   -------
c
c   Objet:
c   ------
c
c   ecriture du fichier histoire
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
!#include "comconst.h"
#include "comcstfi.h"
#include "comvert.h"
#include "comgeom.h"
#include "temps.h"
#include "ener.h"
#include "logic.h"
#include "netcdf.inc"
#include "description.h"
#include "serre.h"

c   Arguments:
c   ----------

      INTEGER*4 idayref
      REAL phis(ip1jmp1)

c   Local:
c   ------
      INTEGER length,l
      parameter (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      INTEGER loop
      INTEGER ierr
      INTEGER ind1,indlast

      integer :: nid,nvarid,idim_index,idim_rlonu,idim_rlonv
      integer :: idim_rlatu,idim_rlatv,idim_llmp1,idim_llm
      integer, dimension(2) :: id  
c-----------------------------------------------------------------------

      DO l=1,length
         tab_cntrl(l)=0.
      ENDDO
      tab_cntrl(1)  = FLOAT(iim)
      tab_cntrl(2)  = FLOAT(jjm)
      tab_cntrl(3)  = FLOAT(llm)
      tab_cntrl(4)  = FLOAT(idayref)
      tab_cntrl(5)  = rad
      tab_cntrl(6)  = omeg
      tab_cntrl(7)  = g
      tab_cntrl(8)  = mugaz
      tab_cntrl(9)  = rcp
      tab_cntrl(10) = daysec
      tab_cntrl(11) = dtphys
      tab_cntrl(12) = etot0
      tab_cntrl(13) = ptot0
      tab_cntrl(14) = ztot0
      tab_cntrl(15) = stot0
      tab_cntrl(16) = ang0
c
c    ..........    P.Le Van  ( ajout le 8/04/96 )    .........
c         .....        parametres  pour le zoom          ......   
      tab_cntrl(17)  = clon
      tab_cntrl(18)  = clat
      tab_cntrl(19)  = grossismx
      tab_cntrl(20)  = grossismy
c
c     .....   ajout  le 6/05/97 et le 15/10/97  .......
c
      IF ( fxyhypb )   THEN
        tab_cntrl(21) = 1.
        tab_cntrl(22) = dzoomx
        tab_cntrl(23) = dzoomy
      ELSE
        tab_cntrl(21) = 0.
        tab_cntrl(22) = dzoomx
        tab_cntrl(23) = dzoomy
        tab_cntrl(24) = 0.
        IF( ysinus )  tab_cntrl(24) = 1.
      ENDIF

c    .........................................................

! Definition des dimensions
    
      ierr = NF_REDEF (nid)

      ierr = NF_DEF_DIM (nid, "index", length, idim_index)
      ierr = NF_DEF_DIM (nid, "rlonu", iip1, idim_rlonu)
      ierr = NF_DEF_DIM (nid, "latitude", jjp1, idim_rlatu)
      ierr = NF_DEF_DIM (nid, "longitude", iip1, idim_rlonv)
      ierr = NF_DEF_DIM (nid, "rlatv", jjm, idim_rlatv)
      ierr = NF_DEF_DIM (nid, "interlayer", (llm+1), idim_llmp1)
      ierr = NF_DEF_DIM (nid, "altitude", llm, idim_llm)
c
      ierr = NF_ENDDEF(nid)

c  Ecriture DRS du tableau des parametres du run
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "controle", NF_DOUBLE, 1, 
     .       idim_index,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "controle", NF_FLOAT, 1, 
     .       idim_index,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 22,
     .                        "Parametres de controle")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif

c --------------------------
c  Ecriture des longitudes et latitudes
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "rlonu", NF_DOUBLE, 1, idim_rlonu,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "rlonu", NF_FLOAT, 1, idim_rlonu,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 23,
     .                        "Longitudes aux points u")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonu/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonu/pi*180)
#endif
c
c --------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr =NF_DEF_VAR(nid, "latitude", NF_DOUBLE, 1, idim_rlatu,nvarid)
#else
      ierr =NF_DEF_VAR (nid, "latitude", NF_FLOAT, 1, idim_rlatu,nvarid)
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
c
c --------------------------
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
c
c --------------------------
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
c
c --------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "rlatv", NF_DOUBLE, 1, idim_rlatv,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "rlatv", NF_FLOAT, 1, idim_rlatv,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 22,
     .                        "Latitudes aux points v")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatv/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatv/pi*180)
#endif
c
c --------------------------
c  Ecriture des niveaux verticaux
      call def_var(nid,"aps","hybrid pressure at midlayers "," ",
     .            1,idim_llm,nvarid,ierr)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aps)
#endif

      call def_var(nid,"bps","hybrid sigma at midlayers"," ",
     .            1,idim_llm,nvarid,ierr)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bps)
#endif

      call def_var(nid,"ap","hybrid pressure at interlayers"," ",
     .            1,idim_llmp1,nvarid,ierr)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ap)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
#endif

      call def_var(nid,"bp","hybrid sigma at interlayers"," ",
     .            1,idim_llmp1,nvarid,ierr)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bp)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bp)
#endif

c
c --------------------------
c  Ecriture aire et coefficients de passage cov. <-> contra. <--> naturel

      id(1)=idim_rlonu
      id(2)=idim_rlatu
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "cu", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "cu", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 39,
     .                      "Coefficients de passage cov <--> naturel")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cu)
#endif
c
      id(1)=idim_rlonv
      id(2)=idim_rlatv
c
c --------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "cv", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "cv", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 39,
     .                      "Coefficients de passage cov <--> naturel")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cv)
#endif
c
      id(1)=idim_rlonv
      id(2)=idim_rlatu
c
c --------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "aire", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "aire", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 17,
     .                        "Aires des mailles")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aire)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aire)
#endif
c
c  Ecriture du geopentiel au sol
      id(1)=idim_rlonv
      id(2)=idim_rlatu
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "phisinit", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "phisinit", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 24,
     .                        "Geopotentiel au sol init")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,phis)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,phis)
#endif
c
c  Ecriture aire et coefficients de passage cov. <-> contra. <--> naturel

      id(1)=idim_rlonu
      id(2)=idim_rlatu
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "cu", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "cu", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 39,
     .                      "Coefficients de passage cov <--> naturel")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cu)
#endif
c
      id(1)=idim_rlonv
      id(2)=idim_rlatv
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "cv", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "cv", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 39,
     .                      "Coefficients de passage cov <--> naturel")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cv)
#endif
c
      id(1)=idim_rlonv
      id(2)=idim_rlatu
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "aire", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "aire", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 17,
     .                        "Aires des mailles")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aire)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aire)
#endif
c
c  Ecriture du geopentiel au sol
      id(1)=idim_rlonv
      id(2)=idim_rlatu
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "phisinit", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "phisinit", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 24,
     .                        "Geopotentiel au sol init")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,phis)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,phis)
#endif
c
      PRINT*,'iim,jjm,llm,idayref',iim,jjm,llm,idayref
      PRINT*,'rad,omeg,g,mugaz,rcp',
     s rad,omeg,g,mugaz,rcp
      PRINT*,'daysec,dtphys',daysec,dtphys

      END
