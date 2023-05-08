      SUBROUTINE phyetat0 (fichnom,tab0,Lmodif,nsoil,nq,
     .           day_ini,time,
     .           tsurf,tsoil,emis,q2,qsurf,co2ice)
      IMPLICIT none
c======================================================================
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c  Adaptation à Mars : Yann Wanherdrick 
c Objet: Lecture de l etat initial pour la physique
c======================================================================
#include "dimensions.h"
#include "netcdf.inc"
#include "dimphys.h"
#include "comgeomfi.h"
#include "surfdat.h"
#include "planete.h"
#include "dimradmars.h"
#include "yomaer.h"
#include "comcstfi.h"
c======================================================================
      CHARACTER*(*) fichnom
      INTEGER nbsrf !Mars nbsrf a 1 au lieu de 4
      PARAMETER (nbsrf=1) ! nombre de sous-fractions pour une maille

      INTEGER radpas
      REAL co2_ppm
      REAL solaire,time
      REAL tsoil(ngridmx,nsoil,nbsrf)

      REAL xmin, xmax
c
      INTEGER nsoil,nq
      INTEGER ig,iq,lmax
      INTEGER nid, nvarid
      INTEGER ierr, i, nsrf, isoil 
      INTEGER length
      PARAMETER (length=100)
      CHARACTER*7 str7
      CHARACTER*2 str2
      CHARACTER*1 yes
c
      integer   Lmodif,tab0
      REAL      p_rad,p_omeg,p_g,p_mugaz,p_daysec
      REAL tsurf(ngridmx,nbsrf),co2ice(ngridmx),emis(ngridmx)
      REAL q2(ngridmx, llm+1)
      REAL qsurf(ngridmx,nq) 
      INTEGER day_ini,nqold
c
c Ouvrir le fichier contenant l etat initial:
c
 
      ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
      IF (ierr.NE.NF_NOERR) THEN
        write(6,*)' Pb d''ouverture du fichier '//fichnom
        CALL ABORT
      ENDIF


c modifications possibles des variables de tab_cntrl
      PRINT*
      write(*,*) 'TABFI de phyeta0',Lmodif,tab0
      call tabfi (nid,Lmodif,tab0,day_ini,lmax,p_rad,
     .              p_omeg,p_g,p_mugaz,p_daysec,time)
c
c Lecture des latitudes (coordonnees):
c
      ierr = NF_INQ_VARID (nid, "latitude", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <latitude> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, lati)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, lati)
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
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, long)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, long)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <longitude>'
         CALL abort
      ENDIF
c
c Lecture des aires des mailles:
c
      ierr = NF_INQ_VARID (nid, "area", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <area> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, area)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, area)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <area>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(area)
      xmax = MAXVAL(area)
      PRINT*,'Aires des mailles <area>:', xmin, xmax
c
c Lecture du geopotentiel au sol:
c
      ierr = NF_INQ_VARID (nid, "phisfi", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <phisfi> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, phisfi)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, phisfi)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <phisfi>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(phisfi)
      xmax = MAXVAL(phisfi)
      PRINT*,'Geopotentiel au sol <phisfi>:', xmin, xmax
c
c Lecture de l''albedo du sol nu:
c
      ierr = NF_INQ_VARID (nid, "albedodat", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <albedodat> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, albedodat)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, albedodat)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <albedodat>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(albedodat)
      xmax = MAXVAL(albedodat)
      PRINT*,'Albedo du sol nu <albedodat>:', xmin, xmax
c
c Lecture de l''inertie thermique du sol:
c
      ierr = NF_INQ_VARID (nid, "inertiedat", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <inertiedat> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, inertiedat)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, inertiedat)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <inertiedat>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(inertiedat)
      xmax = MAXVAL(inertiedat)
      PRINT*,'Inertie thermique du sol <inertiedat>:', xmin, xmax
c
c ZMEA
c
      ierr = NF_INQ_VARID (nid, "ZMEA", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZMEA> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zmea)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, zmea)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZMEA>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zmea(i),xmin)
         xmax = MAX(zmea(i),xmax)
      ENDDO
      PRINT*,'<zmea>:', xmin, xmax
c
c ZSTD
c
      ierr = NF_INQ_VARID (nid, "ZSTD", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZSTD> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zstd)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, zstd)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZSTD>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zstd(i),xmin)
         xmax = MAX(zstd(i),xmax)
      ENDDO
      PRINT*,'<zstd>:', xmin, xmax
c
c ZSIG
c
      ierr = NF_INQ_VARID (nid, "ZSIG", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZSIG> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zsig)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, zsig)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZSIG>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zsig(i),xmin)
         xmax = MAX(zsig(i),xmax)
      ENDDO
      PRINT*,'<zsig>:', xmin, xmax
c
c ZGAM
c
      ierr = NF_INQ_VARID (nid, "ZGAM", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZGAM> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zgam)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, zgam)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZGAM>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zgam(i),xmin)
         xmax = MAX(zgam(i),xmax)
      ENDDO
      PRINT*,'<zgam>:', xmin, xmax
c
c ZTHE
c
      ierr = NF_INQ_VARID (nid, "ZTHE", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <ZTHE> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, zthe)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, zthe)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZTHE>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zthe(i),xmin)
         xmax = MAX(zthe(i),xmax)
      ENDDO
      PRINT*,'<zthe>:', xmin, xmax
c
c CO2 ice cover
c
      ierr = NF_INQ_VARID (nid, "co2ice", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <co2ice> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, co2ice)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, co2ice)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <co2ice>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(co2ice)
      xmax = MAXVAL(co2ice)
      PRINT*,'CO2 ice cover <co2ice>:', xmin, xmax
c
c Lecture des temperatures du sol:
c
      ierr = NF_INQ_VARID (nid, "tsurf", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <tsurf> est absent'
         PRINT*, '          Mais je vais essayer de lire TS**'
         IF (nbsrf.GT.99) THEN
            PRINT*, "Trop de sous-mailles"
            CALL abort
         ENDIF
         DO nsrf = 1, nbsrf
           WRITE(str2,'(i2.2)') nsrf
           ierr = NF_INQ_VARID (nid, "TS"//str2, nvarid)
           IF (ierr.NE.NF_NOERR) THEN
              PRINT*, "phyetat0: Le champ <TS"//str2//"> est absent"
              CALL abort
           ENDIF
#ifdef NC_DOUBLE
           ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tsurf(1,nsrf))
#else
           ierr = NF_GET_VAR_REAL(nid, nvarid, tsurf(1,nsrf))
#endif
           IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Lecture echouee pour <TS"//str2//">"
             CALL abort
           ENDIF
           xmin = 1.0E+20
           xmax = -1.0E+20
           xmin = MINVAL(tsurf)
           xmax = MAXVAL(tsurf)
           PRINT*,'Temperature du sol TS**:', nsrf, xmin, xmax
         ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <tsurf> est present'
         PRINT*, '          J ignore donc les autres temperatures TS**'
#ifdef NC_DOUBLE
         ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tsurf(1,1))
#else
         ierr = NF_GET_VAR_REAL(nid, nvarid, tsurf(1,1))
#endif
         IF (ierr.NE.NF_NOERR) THEN
            PRINT*, "phyetat0: Lecture echouee pour <TSURF>"
            CALL abort
         ENDIF
         xmin = 1.0E+20
         xmax = -1.0E+20
         xmin = MINVAL(tsurf)
         xmax = MAXVAL(tsurf)
         PRINT*,'Temperature du sol <tsurf>', xmin, xmax
         IF (nbsrf >= 2) THEN
            DO nsrf = 2, nbsrf
               DO i = 1, ngridmx
                  tsurf(i,nsrf) = tsurf(i,1)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
c
c Lecture des temperatures du sol profond:
c
      IF (nsoil.GT.99 .OR. nbsrf.GT.99) THEN
         PRINT*, "Trop de couches ou sous-mailles"
         CALL abort
      ENDIF
      DO nsrf = 1, nbsrf
         DO isoil=1, nsoil
            WRITE(str7,'(i2.2,"srf",i2.2)') isoil, nsrf
            ierr = NF_INQ_VARID (nid, 'tsoil', nvarid)
            IF (ierr.NE.NF_NOERR) THEN
               PRINT*, "phyetat0: Le champ <tsoil> est absent"
               PRINT*, "          Il prend donc la valeur de surface"
               DO i=1, ngridmx
                  tsoil(i,isoil,nsrf)=tsurf(i,nsrf)
               ENDDO
            ELSE
#ifdef NC_DOUBLE
              ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tsoil(1,1,nsrf))
#else
              ierr = NF_GET_VAR_REAL(nid, nvarid, tsoil(1,1,nsrf))
#endif
               IF (ierr.NE.NF_NOERR) THEN
                  PRINT*, "Lecture echouee pour <tsoil>"
                  CALL abort
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(tsoil)
      xmax = MAXVAL(tsoil)
      PRINT*,'Temperatures du sol profond <tsoil>', xmin, xmax
c
c Surface emissivity
c
      ierr = NF_INQ_VARID (nid, "emis", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <emis> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, emis)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, emis)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <emis>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(emis)
      xmax = MAXVAL(emis)
      PRINT*,'Surface emissivity <emis>:', xmin, xmax
c
c pbl wind variance
c
      ierr = NF_INQ_VARID (nid, "q2", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Le champ <q2> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, q2)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, q2)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <q2>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(q2)
      xmax = MAXVAL(q2)
      PRINT*,'pbl wind variance <q2>:', xmin, xmax
c
c tracer on surface
c

      IF(nq.GE.1) THEN
         nqold=nq
         DO iq=1,nq
            str7(1:5)='qsurf'
            WRITE(str7(6:7),'(i2.2)') iq
            ierr = NF_INQ_VARID (nid,str7,nvarid)
            IF (ierr.NE.NF_NOERR) THEN
                write (*,*) ' WARNING : ',str7,' not in the file'
                write (*,*) str7, ' set to 0'
                do ig=1,ngridmx
                    qsurf(ig,iq)=0.
                end do
                nqold=min(iq-1,nqold)
            ELSE
#ifdef NC_DOUBLE
              ierr = NF_GET_VAR_DOUBLE(nid, nvarid,qsurf(1,iq))
#else
              ierr = NF_GET_VAR_REAL(nid, nvarid,qsurf(1,iq))
#endif
              IF (ierr.NE.NF_NOERR) THEN
                PRINT*, 'phyetat0: Lecture echouee pour <',str7,'>'
                CALL abort
              ENDIF
           ENDIF
           xmin = 1.0E+20
           xmax = -1.0E+20
           xmin = MINVAL(qsurf(1:ngridmx,iq))
           xmax = MAXVAL(qsurf(1:ngridmx,iq))
      PRINT*,'tracer on surface <',str7,'>:', xmin, xmax
         ENDDO
         if(nqold.lt.nq) then
c        case when new tracer are added in addition to old ones
             write(*,*)'qsurf 1 to ', nqold,'were already present'
             write(*,*)'qsurf ', nqold+1,' to ', nqmx,'are new'
            yes=' '
            do while ((yes.ne.'y').and.(yes.ne.'n'))
             write(*,*) 'Would you like to reindex qsurf # 1 ->',nqold
             write(*,*) 'to #',nqmx-nqold+1,'->', nqmx,'   (y or n) ?'
             read(*,fmt='(a)') yes
            end do
            if (yes.eq.'y') then
              write(*,*) 'OK, let s reindex qsurf'
                 do ig=1,ngridmx
                    do iq=nqmx,nqmx-nqold+1,-1
                       qsurf(ig,iq)=qsurf(ig,iq-nqmx+nqold)
                    end do
                    do iq=nqmx-nqold,1,-1
                       qsurf(ig,iq)= 0.
                    end do
                 end do
            end if
         end if
      ENDIF
c
c Fermer le fichier:
c
      ierr = NF_CLOSE(nid)
c
      RETURN
      END
