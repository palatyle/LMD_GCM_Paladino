      SUBROUTINE writerestartphy(fichnom,tab_cntrl,klon,klev,
     .           rlat,rlon, tsol,tsoil,
     .           albedo, 
     .           solsw, sollw,fder,dlw,
     .           sollwdown,radsol,
     .    zmea, zstd, zsig, zgam, zthe, zpic, zval,
     .           t_ancien)


c======================================================================
c Ecriture de l'etat de redemarrage pour la physique (pour newstart)
c======================================================================
      IMPLICIT none
#include "dimensions.h"
#include "netcdf.inc"
#include "dimsoil.h"
c======================================================================
      CHARACTER*13 fichnom
      INTEGER length,klon,klev
      PARAMETER (length=100)
      REAL tab_cntrl(length)
      REAL rlat(klon), rlon(klon)
      REAL tsol(klon)
      REAL tsoil(klon,nsoilmx)
      REAL albedo(klon)
      real solsw(klon)
      real sollw(klon)
      real fder(klon)
      real dlw(klon)
      real sollwdown(klon)
      REAL radsol(klon)
      REAL zmea(klon), zstd(klon)
      REAL zsig(klon), zgam(klon), zthe(klon)
      REAL zpic(klon), zval(klon)
      REAL t_ancien(klon,klev)
c
      INTEGER nid, nvarid, idim1, idim2, idim3
      INTEGER ierr
c
      INTEGER isoil
      CHARACTER*2 str2
c
c-----------------------------------------------------------------------
c Creation du fichier:
c-----------------------------------------------------------------------
c
c     print*,'fichnom',fichnom
      ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
      IF (ierr.NE.NF_NOERR) THEN
        write(6,*)' Pb d''ouverture du fichier '//fichnom
        write(6,*)' ierr = ', ierr
        CALL ABORT
      ENDIF
c
      ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 28,
     .                       "Fichier redemarrage physique")
c
      ierr = NF_DEF_DIM (nid, "index", length, idim1)
      ierr = NF_DEF_DIM (nid, "points_physiques", klon, idim2)
      ierr = NF_DEF_DIM (nid, "horizon_vertical", klon*klev, idim3)
c

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
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlon)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlon)
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
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlat)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlat)
#endif
c
c-----------------------------------------------------------------------
c Ecriture des variables:
c-----------------------------------------------------------------------
c
        ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
        ierr = NF_DEF_VAR (nid, "TS", NF_DOUBLE, 1, idim2,nvarid)
#else
        ierr = NF_DEF_VAR (nid, "TS", NF_FLOAT, 1, idim2,nvarid)
#endif
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 22,
     .                        "Temperature de surface")
        ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
        ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tsol)
#else
        ierr = NF_PUT_VAR_REAL (nid,nvarid,tsol)
#endif
c
      DO isoil=1, nsoilmx
        IF (isoil.LE.99) THEN
        WRITE(str2,'(i2.2)') isoil
        ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
        ierr = NF_DEF_VAR (nid, "Tsoil"//str2,NF_DOUBLE,1,idim2,nvarid)
#else
        ierr = NF_DEF_VAR (nid, "Tsoil"//str2,NF_FLOAT,1,idim2,nvarid)
#endif
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 24,
     .                        "Temperature du sol No."//str2)
        ierr = NF_ENDDEF(nid)
        ELSE
        PRINT*, "Trop de couches"
        CALL abort
        ENDIF
#ifdef NC_DOUBLE
        ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tsoil(1,isoil))
#else
        ierr = NF_PUT_VAR_REAL (nid,nvarid,tsoil(1,isoil))
#endif
      ENDDO
c
        ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
        ierr = NF_DEF_VAR (nid,"ALBE",NF_DOUBLE,1,idim2,nvarid)
#else
        ierr = NF_DEF_VAR (nid,"ALBE",NF_FLOAT,1,idim2,nvarid)
#endif
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 18,
     .                        "albedo de surface")
        ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,albedo)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,albedo)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "solsw", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "solsw", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 32,
     .                        "Rayonnement solaire a la surface")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,solsw)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,solsw)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "sollw", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "sollw", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 27,
     .                        "Rayonnement IR a la surface")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,sollw)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,sollw)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "fder", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "fder", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 14,
     .                        "Derive de flux")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,fder)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,fder)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "dlw", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "dlw", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 14,
     .                        "Derivee flux IR")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,dlw)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,dlw)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "sollwdown", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "sollwdown", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 14,
     .                        "Flux IR vers le bas a la surface")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,sollwdown)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,sollwdown)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "RADS", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "RADS", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "Rayonnement net a la surface")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,radsol)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,radsol)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZMEA", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZMEA", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "zmea Orographie sous-maille")
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
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "zstd Orographie sous-maille")
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
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "zsig Orographie sous-maille")
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
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "zgam Orographie sous-maille")
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
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "zthe Orographie sous-maille")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zthe)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zthe)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZPIC", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZPIC", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "zpic Orographie sous-maille")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zpic)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zpic)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZVAL", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZVAL", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "zval Orographie sous-maille")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zval)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zval)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "TANCIEN", NF_DOUBLE, 1, idim3,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "TANCIEN", NF_FLOAT, 1, idim3,nvarid)
#endif
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,t_ancien)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,t_ancien)
#endif
c
      ierr = NF_CLOSE(nid)
c
      RETURN
      END
