!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/phyetat0.F,v 1.2 2004/06/22 11:45:33 lmdzadmin Exp $
!
c
c
      SUBROUTINE readstartphy(fichnom,
     .            rlat,rlon, tsol,tsoil,
     .           albe, solsw, sollw,
     .           fder,radsol,
     .    zmea, zstd, zsig, zgam, zthe, zpic, zval,
     .           tabcntr0)
c======================================================================
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c Objet: Lecture de l'etat initial pour la physique
c======================================================================
      use dimphy
      IMPLICIT none
#include "dimensions.h"
#include "netcdf.inc"
#include "dimsoil.h"
c======================================================================
      integer ngridmx
      parameter (ngridmx=(2+(jjm-1)*iim - 1/jjm))

      CHARACTER*(*) fichnom
      REAL rlat(ngridmx), rlon(ngridmx)
      REAL tsol(ngridmx)
      REAL tsoil(ngridmx,nsoilmx)
      REAL albe(ngridmx)
cIM BEG alblw
      REAL alblw(ngridmx)
cIM END alblw
      REAL radsol(ngridmx)
      REAL sollw(ngridmx)
      real solsw(ngridmx)
      real fder(ngridmx)
      REAL zmea(ngridmx), zstd(ngridmx)
      REAL zsig(ngridmx), zgam(ngridmx), zthe(ngridmx)
      REAL zpic(ngridmx), zval(ngridmx)
      INTEGER length
      PARAMETER (length=100)
      REAL tabcntr0(length)

      REAL xmin, xmax
c
      INTEGER nid, nvarid
      INTEGER ierr, i, nsrf, isoil 
      CHARACTER*2 str2
c
c Ouvrir le fichier contenant l'etat initial:
c
      print*,'fichnom',fichnom
      ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
      IF (ierr.NE.NF_NOERR) THEN
        write(6,*)' Pb d''ouverture du fichier '//fichnom
        write(6,*)' ierr = ', ierr
        CALL ABORT
      ENDIF
c
c Lecture des parametres de controle:
c
      ierr = NF_INQ_VARID (nid, "controle", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <controle> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tabcntr0)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, tabcntr0)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <controle>'
         CALL abort
      ENDIF
c
c
c Lecture des latitudes (coordonnees):
c
      ierr = NF_INQ_VARID (nid, "latitude", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <latitude> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlat)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlat)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <latitude>'
         CALL abort
      ENDIF
c
c Lecture des longitudes (coordonnees):
c
      ierr = NF_INQ_VARID (nid, "longitude", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <longitude> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlon)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlon)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <latitude>'
         CALL abort
      ENDIF
C
c Lecture des temperatures du sol:
c
      ierr = NF_INQ_VARID (nid, "TS", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <TS> est absent'
         PRINT*, "phyetat0: Lecture echouee pour <TS>"
         CALL abort
      ELSE
         PRINT*, 'phyetat0: Le champ <TS> est present'
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tsol)
#else
         ierr = NF_GET_VAR_REAL(nid, nvarid, tsol)
#endif
         IF (ierr.NE.NF_NOERR) THEN
            PRINT*, "phyetat0: Lecture echouee pour <TS>"
            CALL abort
         ENDIF
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, ngridmx
            xmin = MIN(tsol(i),xmin)
            xmax = MAX(tsol(i),xmax)
         ENDDO
         PRINT*,'Temperature du sol <TS>', xmin, xmax
      ENDIF
c
c Lecture des temperatures du sol profond:
c
      DO isoil=1, nsoilmx
      IF (isoil.GT.99) THEN
         PRINT*, "Trop de couches"
         CALL abort
      ENDIF
      WRITE(str2,'(i2.2)') isoil
      ierr = NF_INQ_VARID (nid, 'Tsoil'//str2, nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, "phyetat0: Le champ <Tsoil"//str2//"> est absent"
         PRINT*, "          Il prend donc la valeur de surface"
         DO i=1, ngridmx
             tsoil(i,isoil)=tsol(i)
         ENDDO
      ELSE
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tsoil(1,isoil))
#else
         ierr = NF_GET_VAR_REAL(nid, nvarid, tsoil(1,isoil))
#endif
         IF (ierr.NE.NF_NOERR) THEN
            PRINT*, "Lecture echouee pour <Tsoil"//str2//">"
            CALL abort
         ENDIF
      ENDIF
      ENDDO

c
c Lecture de albedo au sol:
c
      ierr = NF_INQ_VARID (nid, "ALBE", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ALBE> est absent'
         PRINT*, "phyetat0: Lecture echouee pour <ALBE>"
         CALL abort
      ELSE
         PRINT*, 'phyetat0: Le champ <ALBE> est present'
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid, nvarid, albe(1))
#else
         ierr = NF_GET_VAR_REAL(nid, nvarid, albe(1))
#endif
         IF (ierr.NE.NF_NOERR) THEN
            PRINT*, "phyetat0: Lecture echouee pour <ALBE>"
            CALL abort
         ENDIF
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, ngridmx
            xmin = MIN(albe(i),xmin)
            xmax = MAX(albe(i),xmax)
         ENDDO
         PRINT*,'Albedo du sol <ALBE>', xmin, xmax
      ENDIF

c
c Lecture rayonnement solaire au sol:
c
      ierr = NF_INQ_VARID (nid, "solsw", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <solsw> est absent'
         PRINT*, 'mis a zero'
         solsw = 0.
      ELSE
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid, nvarid, solsw)
#else
        ierr = NF_GET_VAR_REAL(nid, nvarid, solsw)
#endif
        IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <solsw>'
          CALL abort
        ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(solsw(i),xmin)
         xmax = MAX(solsw(i),xmax)
      ENDDO
      PRINT*,'Rayonnement solaire au sol solsw:', xmin, xmax
c
c Lecture rayonnement IF au sol:
c
      ierr = NF_INQ_VARID (nid, "sollw", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <sollw> est absent'
         PRINT*, 'mis a zero'
         sollw = 0.
      ELSE
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid, nvarid, sollw)
#else
        ierr = NF_GET_VAR_REAL(nid, nvarid, sollw)
#endif
        IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <sollw>'
          CALL abort
        ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(sollw(i),xmin)
         xmax = MAX(sollw(i),xmax)
      ENDDO
      PRINT*,'Rayonnement IF au sol sollw:', xmin, xmax

c
c Lecture derive des flux:
c
      ierr = NF_INQ_VARID (nid, "fder", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <fder> est absent'
         PRINT*, 'mis a zero'
         fder = 0.
      ELSE
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid, nvarid, fder)
#else
        ierr = NF_GET_VAR_REAL(nid, nvarid, fder)
#endif
        IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <fder>'
          CALL abort
        ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(fder(i),xmin)
         xmax = MAX(fder(i),xmax)
      ENDDO
      PRINT*,'Derive des flux fder:', xmin, xmax

c
c Lecture du rayonnement net au sol:
c
      ierr = NF_INQ_VARID (nid, "RADS", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <RADS> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, radsol)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, radsol)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <RADS>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(radsol(i),xmin)
         xmax = MAX(radsol(i),xmax)
      ENDDO
      PRINT*,'Rayonnement net au sol radsol:', xmin, xmax

c
c Lecture des parametres orographie sous-maille:
c
      ierr = NF_INQ_VARID (nid, "ZMEA", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZMEA> est absent'
         PRINT*, 'mis a zero'
         zmea = 0.
      ELSE
#ifdef NC_DOUBLE
       ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zmea)
#else
       ierr = NF_GET_VAR_REAL(nid, nvarid, zmea)
#endif
       IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZMEA>'
         CALL abort
       ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zmea(i),xmin)
         xmax = MAX(zmea(i),xmax)
      ENDDO
      PRINT*,'zmea:', xmin, xmax
c
      ierr = NF_INQ_VARID (nid, "ZSTD", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZSTD> est absent'
         PRINT*, 'mis a zero'
         zstd = 0.
      ELSE
#ifdef NC_DOUBLE
       ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zstd)
#else
       ierr = NF_GET_VAR_REAL(nid, nvarid, zstd)
#endif
       IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZSTD>'
         CALL abort
       ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zstd(i),xmin)
         xmax = MAX(zstd(i),xmax)
      ENDDO
      PRINT*,'zstd:', xmin, xmax
c
      ierr = NF_INQ_VARID (nid, "ZSIG", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZSIG> est absent'
         PRINT*, 'mis a zero'
         zsig = 0.
      ELSE
#ifdef NC_DOUBLE
       ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zsig)
#else
       ierr = NF_GET_VAR_REAL(nid, nvarid, zsig)
#endif
       IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZSIG>'
         CALL abort
       ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zsig(i),xmin)
         xmax = MAX(zsig(i),xmax)
      ENDDO
      PRINT*,'zsig:', xmin, xmax
c
      ierr = NF_INQ_VARID (nid, "ZGAM", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZGAM> est absent'
         PRINT*, 'mis a zero'
         zgam = 0.
      ELSE
#ifdef NC_DOUBLE
       ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zgam)
#else
       ierr = NF_GET_VAR_REAL(nid, nvarid, zgam)
#endif
       IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZGAM>'
         CALL abort
       ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zgam(i),xmin)
         xmax = MAX(zgam(i),xmax)
      ENDDO
      PRINT*,'zgam:', xmin, xmax
c
      ierr = NF_INQ_VARID (nid, "ZTHE", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZTHE> est absent'
         PRINT*, 'mis a zero'
         zthe = 0.
      ELSE
#ifdef NC_DOUBLE
       ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zthe)
#else
       ierr = NF_GET_VAR_REAL(nid, nvarid, zthe)
#endif
       IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZTHE>'
         CALL abort
       ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zthe(i),xmin)
         xmax = MAX(zthe(i),xmax)
      ENDDO
      PRINT*,'zthe:', xmin, xmax
c
      ierr = NF_INQ_VARID (nid, "ZPIC", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZPIC> est absent'
         PRINT*, 'mis a zero'
         zpic = 0.
      ELSE
#ifdef NC_DOUBLE
       ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zpic)
#else
       ierr = NF_GET_VAR_REAL(nid, nvarid, zpic)
#endif
       IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZPIC>'
         CALL abort
       ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zpic(i),xmin)
         xmax = MAX(zpic(i),xmax)
      ENDDO
      PRINT*,'zpic:', xmin, xmax
c
      ierr = NF_INQ_VARID (nid, "ZVAL", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZVAL> est absent'
         PRINT*, 'mis a zero'
         zval = 0.
      ELSE
#ifdef NC_DOUBLE
       ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zval)
#else
       ierr = NF_GET_VAR_REAL(nid, nvarid, zval)
#endif
       IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZVAL>'
         CALL abort
       ENDIF
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zval(i),xmin)
         xmax = MAX(zval(i),xmax)
      ENDDO
      PRINT*,'zval:', xmin, xmax
c
c Fermer le fichier:
c
      ierr = NF_CLOSE(nid)
c
      RETURN
      END

