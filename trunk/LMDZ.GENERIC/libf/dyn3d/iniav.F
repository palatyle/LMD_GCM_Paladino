      SUBROUTINE iniav(fichnom,pday0,pyear0,pfrac0,pperiod,phis,nq)

      USE comvert_mod, ONLY: ap,bp,nivsigs,nivsig,presnivs

      IMPLICIT NONE

c=======================================================================
c
c   pday0    jour de reference pour le stockage
c   pfrac0   decalage du debut du stockage par rapport a pday0 (en jour)
c   pperiod  period de stockage en secondes
c   phis     geopotentiel au sol
c   
c   initialisation de l'ecriture des moyennes temporelles.
c
c=======================================================================

#include "dimensions.h"
#include "paramet.h"
#include "comav.h"
#include "netcdf.inc"
#include "comgeom.h"

      REAL pperiod

      INTEGER*4 pday0,pyear0
      REAL pfrac0
      CHARACTER*(*) fichnom
      REAL phis(ip1jmp1)
c
      CHARACTER str3*3
c
      INTEGER nq
      INTEGER ierr, nid, nvarid, iq
      INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
      INTEGER idim_s, idim_sig
      INTEGER idim_tim
      INTEGER dims2(2), dims3(3), dims4(4)

      period_av=pperiod
      time0_av=pfrac0
c
      ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
      IF (ierr.NE.NF_NOERR) THEN
         WRITE(6,*)" Pb d ouverture du fichier "//fichnom
         WRITE(6,*)' ierr = ', ierr
         CALL ABORT
      ENDIF
c
c Preciser quelques attributs globaux:
c
      ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 32,
     .                       "Fichier moyen histoire dynamique")
c
c Definir les dimensions du fichiers:
c
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
c
      ierr = NF_REDEF (nid)
ccc      ierr = NF_DEF_VAR (nid,"rlonu",NF_DOUBLE,1,idim_rlonu,nvarid)
      ierr = NF_DEF_VAR (nid,"rlonu",NF_FLOAT,1,idim_rlonu,nvarid)
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonu)
#endif
c
      ierr = NF_REDEF (nid)
ccc      ierr = NF_DEF_VAR (nid,"rlatu",NF_DOUBLE,1,idim_rlatu,nvarid)
      ierr = NF_DEF_VAR (nid,"rlatu",NF_FLOAT,1,idim_rlatu,nvarid)
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu)
#endif
c
      ierr = NF_REDEF (nid)
ccc      ierr = NF_DEF_VAR (nid,"rlonv",NF_DOUBLE,1,idim_rlonv,nvarid)
      ierr = NF_DEF_VAR (nid,"rlonv",NF_FLOAT,1,idim_rlonv,nvarid)
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv)
#endif
c
      ierr = NF_REDEF (nid)
ccc      ierr = NF_DEF_VAR (nid,"rlatv",NF_DOUBLE,1,idim_rlatv,nvarid)
      ierr = NF_DEF_VAR (nid,"rlatv",NF_FLOAT,1,idim_rlatv,nvarid)
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatv)
#endif
c
      ierr = NF_REDEF (nid)
ccc      ierr = NF_DEF_VAR (nid,"nivsigs",NF_DOUBLE,1,idim_s,nvarid)
      ierr = NF_DEF_VAR (nid,"nivsigs",NF_FLOAT,1,idim_s,nvarid)
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,nivsigs)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsigs)
#endif
c
      ierr = NF_REDEF (nid)
ccc      ierr = NF_DEF_VAR (nid,"nivsig",NF_DOUBLE,1,idim_sig,nvarid)
      ierr = NF_DEF_VAR (nid,"nivsig",NF_FLOAT,1,idim_sig,nvarid)
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,nivsig)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsig)
#endif
c
      ierr = NF_REDEF (nid)
ccc      ierr = NF_DEF_VAR (nid,"ap",NF_DOUBLE,1,idim_sig,nvarid)
      ierr = NF_DEF_VAR (nid,"ap",NF_FLOAT,1,idim_sig,nvarid)
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ap)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
#endif
c
      ierr = NF_REDEF (nid)
ccc      ierr = NF_DEF_VAR (nid,"bp",NF_DOUBLE,1,idim_sig,nvarid)
      ierr = NF_DEF_VAR (nid,"bp",NF_FLOAT,1,idim_sig,nvarid)
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bp)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bp)
#endif
c
      ierr = NF_REDEF (nid)
ccc      ierr = NF_DEF_VAR (nid,"presnivs",NF_DOUBLE,1,idim_s,nvarid)
      ierr = NF_DEF_VAR (nid,"presnivs",NF_FLOAT,1,idim_s,nvarid)
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,presnivs)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,presnivs)
#endif
c
c Geopentiel au sol:
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatu
ccc      ierr = NF_DEF_VAR (nid,"phisinit",NF_DOUBLE,2,dims2,nvarid)
      ierr = NF_DEF_VAR (nid,"phisinit",NF_FLOAT,2,dims2,nvarid)
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
ccc      ierr = NF_DEF_VAR (nid,"temps",NF_DOUBLE,1,idim_tim,nvarid)
      ierr = NF_DEF_VAR (nid,"temps",NF_FLOAT,1,idim_tim,nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "iniav: Erreur pour definir <temps>"
         CALL abort
      ENDIF
c
      dims4(1) = idim_rlonu
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
ccc      ierr = NF_DEF_VAR (nid,"vitu",NF_DOUBLE,4,dims4,nvarid)
      ierr = NF_DEF_VAR (nid,"vitu",NF_FLOAT,4,dims4,nvarid)
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatv
      dims4(3) = idim_s
      dims4(4) = idim_tim
ccc      ierr = NF_DEF_VAR (nid,"vitv",NF_DOUBLE,4,dims4,nvarid)
      ierr = NF_DEF_VAR (nid,"vitv",NF_FLOAT,4,dims4,nvarid)
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
ccc      ierr = NF_DEF_VAR (nid,"temp",NF_DOUBLE,4,dims4,nvarid)
      ierr = NF_DEF_VAR (nid,"temp",NF_FLOAT,4,dims4,nvarid)
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
      IF(nq.GE.1) THEN
      DO iq=1,nq
      IF (iq.GT.99) THEN
         PRINT*, "Trop de traceurs"
         CALL abort
      ELSE
         str3(1:1)='q'
         WRITE(str3(2:3),'(i2.2)') iq
ccc         ierr = NF_DEF_VAR (nid,str3,NF_DOUBLE,4,dims4,nvarid)
         ierr = NF_DEF_VAR (nid,str3,NF_FLOAT,4,dims4,nvarid)
      ENDIF
      ENDDO
      ENDIF
c
      dims3(1) = idim_rlonv
      dims3(2) = idim_rlatu
      dims3(3) = idim_tim
ccc      ierr = NF_DEF_VAR (nid,"ps",NF_DOUBLE,3,dims3,nvarid)
      ierr = NF_DEF_VAR (nid,"ps",NF_FLOAT,3,dims3,nvarid)
c
      ierr = NF_ENDDEF(nid) ! sortir du mode de definition
      ierr = NF_CLOSE(nid) ! fermer le fichier

      PRINT*,'le fichier hist.moy stoke les variables toutes'
      PRINT*,'les ',period_av,'s'
      PRINT*,'A partir du jour de reference day0'
      RETURN
      END
