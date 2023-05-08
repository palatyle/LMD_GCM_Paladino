      SUBROUTINE writerestart(fichnom,tab_cntrl,
     .                  phis,vcov,ucov,teta,masse,q,ps)

      USE IOIPSL
      USE infotrac
      USE comconst_mod
      USE comvert_mod, ONLY: ap,bp,presnivs,pa,preff,nivsigs,nivsig,
     .                       aps,bps,scaleheight,pseudoalt,
     .                       disvert_type,pressure_exner

      IMPLICIT NONE

c=======================================================================
c Ecriture du fichier de redemarrage sous format NetCDF (pour newstart)
c=======================================================================
c   Declarations:
c   -------------
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "netcdf.inc"
#include "description.h"

c   Arguments:
c   ----------
      CHARACTER*(*) fichnom
      INTEGER length
      PARAMETER (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      REAL phis(ip1jmp1)
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
      REAL teta(ip1jmp1,llm),q(ip1jmp1,llm,nqtot) 
      REAL ps(ip1jmp1),masse(ip1jmp1,llm)                   

c   Local:
c   ------
      INTEGER iq,l
      INTEGER ierr
      character*20 modname
      character*80 abort_message

c   Variables locales pour NetCDF:
c
      INTEGER dims2(2), dims3(3), dims4(4)
      INTEGER idim_index
      INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
      INTEGER idim_s, idim_sig
      INTEGER idim_tim
      INTEGER nid,nvarid

      character*30 unites


c-----------------------------------------------------------------------
c Creation du fichier:
c-----------------------------------------------------------------------
c
      modname='write_restart_ini'

      ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
      IF (ierr.NE.NF_NOERR) THEN
         WRITE(6,*)" Pb d ouverture du fichier "//fichnom
         WRITE(6,*)' ierr = ', ierr
         CALL ABORT
      ENDIF
c
c Preciser quelques attributs globaux:
c
      ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 27,
     .                       "Fichier demarrage dynamique")
c
c Definir les dimensions du fichiers:
c
      ierr = NF_DEF_DIM (nid, "index", length, idim_index)
      ierr = NF_DEF_DIM (nid, "rlonu", iip1, idim_rlonu)
      ierr = NF_DEF_DIM (nid, "rlatu", jjp1, idim_rlatu)
      ierr = NF_DEF_DIM (nid, "rlonv", iip1, idim_rlonv)
      ierr = NF_DEF_DIM (nid, "rlatv", jjm, idim_rlatv)
      ierr = NF_DEF_DIM (nid, "sigs", llm, idim_s)
      ierr = NF_DEF_DIM (nid, "sig", llmp1, idim_sig)
      ierr = NF_DEF_DIM (nid, "temps", NF_UNLIMITED, idim_tim)
c
      ierr = NF_ENDDEF(nid) ! sortir du mode de definition
c
c Definir et enregistrer certains champs invariants:
c
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
c
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
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"nivsigs",NF_DOUBLE,1,idim_s,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"nivsigs",NF_FLOAT,1,idim_s,nvarid)
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
      ierr = NF_DEF_VAR (nid,"nivsig",NF_DOUBLE,1,idim_sig,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"nivsig",NF_FLOAT,1,idim_sig,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 32,
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
      ierr = NF_DEF_VAR (nid,"ap",NF_DOUBLE,1,idim_sig,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ap",NF_FLOAT,1,idim_sig,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 26,
     .                       "Coefficient A pour hybride")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ap)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"bp",NF_DOUBLE,1,idim_sig,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"bp",NF_FLOAT,1,idim_sig,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 26,
     .                       "Coefficient B pour hybride")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bp)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bp)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"presnivs",NF_DOUBLE,1,idim_s,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"presnivs",NF_FLOAT,1,idim_s,nvarid)
#endif
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,presnivs)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,presnivs)
#endif
c
c Coefficients de passage cov. <-> contra. <--> naturel
c
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
c
c Geopentiel au sol:
c
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
c
c Definir les variables pour pouvoir les enregistrer plus tard:
c
      ierr = NF_REDEF (nid) ! entrer dans le mode de definition
c
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"temps",NF_DOUBLE,1,idim_tim,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"temps",NF_FLOAT,1,idim_tim,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                       "Temps de simulation")
      write(unites,200) 1990,1,1
200   format('days since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "units", 30,
     .                         unites)

c
      dims4(1) = idim_rlonu
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ucov",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ucov",NF_FLOAT,4,dims4,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9,
     .                       "Vitesse U")
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatv
      dims4(3) = idim_s
      dims4(4) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"vcov",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"vcov",NF_FLOAT,4,dims4,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9,
     .                       "Vitesse V")
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"teta",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"teta",NF_FLOAT,4,dims4,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 11,
     .                       "Temperature")
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
      DO iq=1,nqtot
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,tname(iq),NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,tname(iq),NF_FLOAT,4,dims4,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12,ttext(iq))
      ENDDO
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"masse",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"masse",NF_FLOAT,4,dims4,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12,
     .                       "C est quoi ?")
c
      dims3(1) = idim_rlonv
      dims3(2) = idim_rlatu
      dims3(3) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ps",NF_DOUBLE,3,dims3,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ps",NF_FLOAT,3,dims3,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 15,
     .                       "Pression au sol")
c
      ierr = NF_ENDDEF(nid) ! sortir du mode de definition

c-----------------------------------------------------------------------
c Ecriture des variables:
c-----------------------------------------------------------------------
c
      modname='write_restart'

      ierr = NF_INQ_VARID(nid, "temps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         print *, NF_STRERROR(ierr)
         abort_message='Variable temps n est pas definie'
         CALL abort_gcm(modname,abort_message,ierr)
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR1_DOUBLE (nid,nvarid,1,0.)
#else
      ierr = NF_PUT_VAR1_REAL (nid,nvarid,1,0.)
#endif

      ierr = NF_INQ_VARID(nid, "ucov", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable ucov n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ucov)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ucov)
#endif

      ierr = NF_INQ_VARID(nid, "vcov", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable vcov n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,vcov)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,vcov)
#endif

      ierr = NF_INQ_VARID(nid, "teta", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable teta n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,teta)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,teta)
#endif

       do iq=1,nqtot   
        ierr = NF_INQ_VARID(nid, tname(iq), nvarid)
        IF (ierr .NE. NF_NOERR) THEN
           PRINT*, "Variable  tname(iq) n est pas definie"
           CALL abort
        ENDIF
#ifdef NC_DOUBLE
          ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,q(1,1,iq))
#else
          ierr = NF_PUT_VAR_REAL (nid,nvarid,q(1,1,iq))
#endif
      ENDDO
c
      ierr = NF_INQ_VARID(nid, "masse", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable masse n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,masse)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,masse)
#endif
c
      ierr = NF_INQ_VARID(nid, "ps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable ps n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ps)
#endif

      ierr = NF_CLOSE(nid) ! fermer le fichier

c      PRINT*,'iim,jjm,llm,iday_end',iim,jjm,llm,iday_end
c      PRINT*,'rad,omeg,g,cpp,kappa',
c     ,        rad,omeg,g,cpp,kappa

      RETURN
      END

