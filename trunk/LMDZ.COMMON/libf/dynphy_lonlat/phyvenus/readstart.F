!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/dynetat0.F,v 1.2 2004/06/22 11:45:30 lmdzadmin Exp $
!
      SUBROUTINE readstart(fichnom,nq,vcov,ucov,
     .                    teta,q,masse,ps,phis,tab_cntrl)

      USE infotrac
      USE comconst_mod
      USE comvert_mod, ONLY: ap,bp,presnivs,pa,preff,nivsigs,nivsig,
     .                       aps,bps,scaleheight,pseudoalt,
     .                       disvert_type,pressure_exner
      USE serre_mod, ONLY: clon,clat,grossismx,grossismy
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0
      USE temps_mod, ONLY: day_ref,annee_ref,day_ini,itau_dyn,
     &                     start_time
      USE logic_mod, ONLY: fxyhypb,ysinus,iflag_trac
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van / L.Fairhead
c   -------
c
c   objet:
c   ------
c
c   Lecture de l'etat initial
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "description.h"
#include "netcdf.inc"

c   Arguments:
c   ----------

      CHARACTER*(*) fichnom
      INTEGER nq
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm),teta(ip1jmp1,llm)
      REAL q(ip1jmp1,llm,nq),masse(ip1jmp1,llm)
      REAL ps(ip1jmp1),phis(ip1jmp1)
      REAL time
      INTEGER length
      PARAMETER (length=100)
      REAL tab_cntrl(length) ! tableau des parametres du run

c   Variables 
c
      INTEGER iq,i,j,ij,l
      INTEGER ierr, nid, nvarid

c   local, cas particulier compo.dat
      integer nyread
      real    qy(jjp1,llm,nq)
      character*10 nomy(nq)

c-----------------------------------------------------------------------

c  Ouverture NetCDF du fichier etat initial

      ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
      IF (ierr.NE.NF_NOERR) THEN
        write(6,*)' Pb d''ouverture du fichier start.nc'
        write(6,*)' ierr = ', ierr
        CALL ABORT
      ENDIF

c
      ierr = NF_INQ_VARID (nid, "controle", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <controle> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echoue pour <controle>"
         CALL abort
      ENDIF

      im         = tab_cntrl(1)
      jm         = tab_cntrl(2)
      lllm       = tab_cntrl(3)
      day_ref    = tab_cntrl(4)
      annee_ref  = tab_cntrl(5)
      rad        = tab_cntrl(6)
      omeg       = tab_cntrl(7)
      g          = tab_cntrl(8)
      cpp        = tab_cntrl(9)
      kappa      = tab_cntrl(10)
      daysec     = tab_cntrl(11)
      dtvr       = tab_cntrl(12)
      etot0      = tab_cntrl(13)
      ptot0      = tab_cntrl(14)
      ztot0      = tab_cntrl(15)
      stot0      = tab_cntrl(16)
      ang0       = tab_cntrl(17)
      pa         = tab_cntrl(18)
      preff      = tab_cntrl(19)
c
      clon       = tab_cntrl(20)
      clat       = tab_cntrl(21)
      grossismx  = tab_cntrl(22)
      grossismy  = tab_cntrl(23)
c
      IF ( tab_cntrl(24).EQ.1. )  THEN
        fxyhypb  = . TRUE .
c        dzoomx   = tab_cntrl(25)
c        dzoomy   = tab_cntrl(26)
c        taux     = tab_cntrl(28)
c        tauy     = tab_cntrl(29)
      ELSE
        fxyhypb = . FALSE .
        ysinus  = . FALSE .
        IF( tab_cntrl(27).EQ.1. ) ysinus = . TRUE. 
      ENDIF

      day_ini = tab_cntrl(30)
      itau_dyn = tab_cntrl(31)
      start_time = tab_cntrl(32)
c   .................................................................
c
c
      PRINT*,'rad,omeg,g,cpp,kappa',rad,omeg,g,cpp,kappa

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
         PRINT*, "dynetat0: Le champ <rlonu> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonu)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlonu)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <rlonu>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "rlatu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <rlatu> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatu)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlatu)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <rlatu>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "rlonv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <rlonv> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonv)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlonv)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <rlonv>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "rlatv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <rlatv> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatv)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, rlatv)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour rlatv"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "nivsigs", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <nivsigs> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, nivsigs)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, nivsigs)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <nivsigs>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "nivsig", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <nivsig> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, nivsig)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, nivsig)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <nivsig>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "ap", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <ap> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, ap)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, ap)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <ap>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "bp", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <bp> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, bp)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, bp)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <bp>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "presnivs", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <presnivs> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, presnivs)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, presnivs)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <presnivs>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "cu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <cu> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, cu)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, cu)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <cu>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "cv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <cv> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, cv)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, cv)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <cv>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "aire", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <aire> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, aire)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, aire)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <aire>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "phisinit", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <phisinit> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, phis)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, phis)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <phisinit>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "temps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <temps> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, time)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, time)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee <temps>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "ucov", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <ucov> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, ucov)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, ucov)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <ucov>"
         CALL abort
      ENDIF
 
      ierr = NF_INQ_VARID (nid, "vcov", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <vcov> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, vcov)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, vcov)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <vcov>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "teta", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <teta> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, teta)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, teta)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <teta>"
         CALL abort
      ENDIF

c TNAME: IL EST LU A PARTIR DE traceur.def (mettre l'ancien si
c                changement du nombre de traceurs)

      IF((nq.GE.1).and.(iflag_trac.eq.1)) THEN
      DO iq=1,nq
        ierr =  NF_INQ_VARID (nid, tname(iq), nvarid)
        IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "dynetat0: Le champ <"//tname(iq)//"> est absent"
            PRINT*, "          Il est donc initialise a zero"
            q(:,:,iq)=0.
        ELSE
#ifdef NC_DOUBLE
          ierr = NF_GET_VAR_DOUBLE(nid, nvarid, q(1,1,iq))
#else
          ierr = NF_GET_VAR_REAL(nid, nvarid, q(1,1,iq))
#endif
          IF (ierr .NE. NF_NOERR) THEN
             PRINT*, "dynetat0: Lecture echouee pour "//tname(iq)
             CALL abort
          ENDIF
        ENDIF
      ENDDO
      ENDIF

c--------------------------------------------
c cas particulier: lecture des traceurs 2D dans compo.dat (issu de start 2d)
c 
      if (1.eq.0) then
      OPEN(10,file='compo.dat',status='old',form='formatted',
     . iostat=ierr)
      IF (ierr.ne.0) THEN
       WRITE(6,*)' Pb d''ouverture du fichier de demarrage (compo.dat)'
       WRITE(6,*)' ierr = ', ierr
       CALL exit(1)
      ENDIF
      READ(10,*) nyread
      print*,"nombre de composes chimiques ajoutes:",nyread
      READ(10,*) (((qy(ij,l,iq),ij=1,jjp1),l=1,llm),
     s                   iq=1,nyread)
      do iq=1,nyread
         READ(10,'(1X,A10)') nomy(iq)
         print*,nomy(iq)," = ", tname(iq+10)
         do i=1,iip1
          do j=1,jjp1
           ij = (j-1)*iip1+i
           q(ij,:,iq+10) = qy(j,:,iq)
          enddo
         enddo
      enddo
      CLOSE(10)
      endif
c--------------------------------------------

      ierr = NF_INQ_VARID (nid, "masse", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <masse> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, masse)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, masse)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <masse>"
         CALL abort
      ENDIF

      ierr = NF_INQ_VARID (nid, "ps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Le champ <ps> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, ps)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, ps)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee pour <ps>"
         CALL abort
      ENDIF

      ierr = NF_CLOSE(nid)

       day_ini=day_ini+INT(time)
       time=time-INT(time)

  1   FORMAT(//10x,'la valeur de im =',i4,2x,'lue sur le fichier de dem
     *arrage est differente de la valeur parametree iim =',i4//)
   2  FORMAT(//10x,'la valeur de jm =',i4,2x,'lue sur le fichier de dem
     *arrage est differente de la valeur parametree jjm =',i4//)
   3  FORMAT(//10x,'la valeur de lmax =',i4,2x,'lue sur le fichier dema
     *rrage est differente de la valeur parametree llm =',i4//)
   4  FORMAT(//10x,'la valeur de dtrv =',i4,2x,'lue sur le fichier dema
     *rrage est differente de la valeur  dtinteg =',i4//)

      RETURN
      END
