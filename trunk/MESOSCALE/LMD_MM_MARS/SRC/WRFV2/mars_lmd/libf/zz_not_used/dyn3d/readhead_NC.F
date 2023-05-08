      SUBROUTINE readhead_NC (fichnom,
     .           day0,
     .           phis,constR)
      IMPLICIT none
c======================================================================
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c  Adaptation à Mars : Yann Wanherdrick 
c Objet: Lecture de l etat initial pour la physique
c======================================================================
#include "netcdf.inc"
c====== includes de l ancien readhead ===
#include "dimensions.h"
#include "paramet.h"
#include "temps.h"
#include "comconst.h"
#include "comvert.h"
#include "comgeom.h"
#include "ener.h"
#include "description.h"

c======================================================================

      CHARACTER*(*) fichnom
      INTEGER nbsrf !Mars nbsrf a 1 au lieu de 4
      PARAMETER (nbsrf=1) ! nombre de sous-fractions pour une maille

      INTEGER radpas

      REAL xmin, xmax
c
      INTEGER  i

c   Variables
c
      INTEGER length,iq
      PARAMETER (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      INTEGER ierr, nid, nvarid
      CHARACTER  str3*3

c
      INTEGER day0
      REAL phis(ip1jmp1),constR
c
c Ouvrir le fichier contenant l etat initial:
c
      ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
      IF (ierr.NE.NF_NOERR) THEN
        write(6,*)' Pb d''ouverture du fichier '//fichnom
        CALL ABORT
      ENDIF
c
c Lecture des parametres de controle:
c
      ierr = NF_INQ_VARID (nid, "controle", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'readhead_NC: Le champ <controle> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'readhead_NC: Lecture echouee pour <controle>'
         CALL abort
      ENDIF


c Info sur la Planete Mars pour la dynamique 
      im         = tab_cntrl(1)
      jm         = tab_cntrl(2)
      lllm       = tab_cntrl(3)
      day_ini    = tab_cntrl(4)
      rad        = tab_cntrl(5)
      omeg       = tab_cntrl(6)
      g          = tab_cntrl(7)
c      mugaz      = tab_cntrl(8)
      cpp        =  744.499
      kappa      = tab_cntrl(9)
      daysec     = tab_cntrl(10)
      dtvr       = tab_cntrl(11)
      etot0      = tab_cntrl(12)
      ptot0      = tab_cntrl(13)
      ztot0      = tab_cntrl(14)
      stot0      = tab_cntrl(15)
      ang0       = tab_cntrl(16)
c pas vrai pour diagfi, seulement pour start      preff      = tab_cntrl(18)
      preff=610.
      WRITE (*,*) 'readhead -     preff ' , preff 
c

      day0=day_ini

      constR=kappa*cpp
      WRITE (*,*) 'constR = ' , constR
      r=constR
      IF(   im.ne.iim           )  THEN
          PRINT 1,im,iim
          STOP
      ELSE  IF( jm.ne.jjm       )  THEN
          PRINT 2,jm,jjm
          STOP
      ELSE  IF( lllm.ne.llm     )  THEN
          PRINT 3,lllm,llm
          STOP
      ENDIF
                                                                       
      ierr = NF_INQ_VARID (nid, "rlonu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "readhead_NC: Le champ <rlonu> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonu)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlonu)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "readhead_NC: Lecture echouee pour <rlonu>"
         CALL abort
      ENDIF
                                                                       
      ierr = NF_INQ_VARID (nid, "rlatv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "readhead_NC: Le champ <rlatv> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatv)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlatv)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "readhead_NC: Lecture echouee pour rlatv"
         CALL abort
      ENDIF

      ierr = NF_GET_VAR_REAL(nid, nvarid, cv)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "readhead_NC: Lecture echouee pour <cv>"
         CALL abort
      ENDIF
c
c Lecture des aires des mailles:
c
      ierr = NF_INQ_VARID (nid, "aire", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'readhead_NC: Le champ <aire> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, aire)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, aire)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'readhead_NC: Lecture echouee pour <aire>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(aire)
      xmax = MAXVAL(aire)
      PRINT*,'Aires des mailles <aire>:', xmin, xmax
c
c Lecture du geopotentiel au sol:
c
      ierr = NF_INQ_VARID (nid, "phisinit", nvarid)
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'readhead_NC: Le champ <phisinit> est absent'
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, phis)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, phis)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'readhead_NC: Lecture echouee pour <phis>'
         CALL abort
      ENDIF
c      PRINT*,'READHEAD_NC  Phis:',phis

      ierr = NF_INQ_VARID (nid, "aps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "readhead_NC: Le champ <aps> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, aps)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, aps)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "readhead_NC: Lecture echouee pour <aps>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "bps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "readhead_NC: Le champ <bps> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, bps)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, bps)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "readhead_NC: Lecture echouee pour <bps>"
         CALL abort
      ENDIF

   1  FORMAT(//10x,'la valeur de im =',i4,2x,'lue sur le fichier de dema
     *rrage est differente de la valeur parametree iim =',i4//)
   2  FORMAT(//10x,'la valeur de jm =',i4,2x,'lue sur le fichier de dema
     *rrage est differente de la valeur parametree jjm =',i4//)
   3  FORMAT(//10x,'la valeur de lmax =',i4,2x,'lue sur le fichier demar
     *rage est differente de la valeur parametree llm =',i4//)
   4  FORMAT(//10x,'la valeur de dtrv =',i4,2x,'lue sur le fichier demar
     *rage est differente de la valeur  dtinteg =',i4//)

      
c Fermer le fichier:
c
      ierr = NF_CLOSE(nid)
c
      RETURN
      END
