      SUBROUTINE dynetat0(fichnom,vcov,ucov,
     .                    teta,q,masse,ps,phis,time)
      use infotrac, only: tname, nqtot
      USE comvert_mod, ONLY: pa,preff
      USE comconst_mod, ONLY: im,jm,lllm,daysec,dtvr,rad,omeg,g,cpp,
     .			kappa
      USE logic_mod, ONLY: fxyhypb,ysinus
      USE serre_mod, ONLY: clon,clat,grossismx,grossismy,dzoomx,dzoomy,
     .			taux,tauy
      USE temps_mod, ONLY: day_ini
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0
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
c   Modifs: Oct.2008 read in tracers by name. Ehouarn Millour
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "netcdf.inc"
!#include"advtrac.h"

c   Arguments:
c   ----------

      CHARACTER*(*) fichnom
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm),teta(ip1jmp1,llm)
      REAL q(iip1,jjp1,llm,nqtot),masse(ip1jmp1,llm)
      REAL ps(ip1jmp1),phis(ip1jmp1)

      REAL time

c   Variables 
c
      INTEGER length,iq,i,j,l
      PARAMETER (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      INTEGER ierr, nid, nvarid, nqold
      CHARACTER  str3*3,yes*1


!     added by RW for test
      real pmean,airetot
      integer ij

c-----------------------------------------------------------------------

c  Ouverture NetCDF du fichier etat initial

      ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
      IF (ierr.NE.NF_NOERR) THEN
        write(6,*)' Pb d''ouverture du fichier ',fichnom
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
      day_ini    = tab_cntrl(4)
      rad        = tab_cntrl(5)
      omeg       = tab_cntrl(6)
      g          = tab_cntrl(7)
      cpp        = tab_cntrl(8)
      kappa      = tab_cntrl(9)
      daysec     = tab_cntrl(10)
      dtvr       = tab_cntrl(11)
      etot0      = tab_cntrl(12)
      ptot0      = tab_cntrl(13)
      ztot0      = tab_cntrl(14)
      stot0      = tab_cntrl(15)
      ang0       = tab_cntrl(16)
      pa         = tab_cntrl(17)
      preff      = tab_cntrl(18)
c
      clon       = tab_cntrl(19)
      clat       = tab_cntrl(20)
      grossismx  = tab_cntrl(21)
      grossismy  = tab_cntrl(22)
c
      IF ( tab_cntrl(23).EQ.1. )  THEN
        fxyhypb  = . TRUE .
        dzoomx   = tab_cntrl(24)
        dzoomy   = tab_cntrl(25)
        taux     = tab_cntrl(27)
        tauy     = tab_cntrl(28)
      ELSE
        fxyhypb = . FALSE .
        ysinus  = . FALSE .
        IF( tab_cntrl(26).EQ.1. ) ysinus = . TRUE. 
      ENDIF
c   .................................................................
c
c
      PRINT*,'dynetat0: rad,omeg,g,cpp,kappa',rad,omeg,g,cpp,kappa
 
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

      ierr = NF_INQ_VARID (nid, "Time", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
	     ierr = NF_INQ_VARID (nid, "temps", nvarid)
		 IF (ierr .NE. NF_NOERR) THEN
           PRINT*, "dynetat0: <Time> or <temps> absent"
           CALL abort
         ENDIF
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, time)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, time)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "dynetat0: Lecture echouee <Time>/<temps>"
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


      IF(nqtot.GE.1) THEN
        write(*,*) 'dynetat0: loading tracers'
         IF(nqtot.GT.99) THEN
            PRINT*, "Trop de traceurs"
            CALL abort
         ENDIF
         nqold=nqtot
         DO iq=1,nqtot
!           str3(1:1)='q'
!           WRITE(str3(2:3),'(i2.2)') iq
!           ierr =  NF_INQ_VARID (nid, str3, nvarid)
! NB: tracers are now read in using their name ('tname' from infotrac)
!           write(*,*) "  loading tracer:",trim(tname(iq))
           ierr=NF_INQ_VARID(nid,tname(iq),nvarid)
           IF (ierr .NE. NF_NOERR) THEN
!              PRINT*, "dynetat0: Le champ <"//str3//"> est absent"
              PRINT*, "dynetat0: Le champ <"//trim(tname(iq))//
     &                "> est absent"
              PRINT*, "          Il est donc initialise a zero"
              CALL initial0(ijp1llm,q(1,1,1,iq))
              nqold=min(iq-1,nqold)
           ELSE
#ifdef NC_DOUBLE
           ierr = NF_GET_VAR_DOUBLE(nid, nvarid, q(1,1,1,iq))
#else
           ierr = NF_GET_VAR_REAL(nid, nvarid, q(1,1,1,iq))
#endif
             IF (ierr .NE. NF_NOERR) THEN
!                 PRINT*, "dynetat0: Lecture echouee pour "//str3
               PRINT*,"dynetat0: Lecture echouee pour "//trim(tname(iq))
               CALL abort
             ENDIF
           ENDIF
         ENDDO
         if ((nqold.lt.nqtot).and.(nqold.ge.1)) then   
c        case when new tracer are added in addition to old ones
             write(*,*)'tracers 1 to ', nqold,'were already present'
             write(*,*)'tracers ', nqold+1,' to ', nqtot,'are new'
             write(*,*)' and initialized to zero'
             q(:,:,:,nqold+1:nqtot)=0.0
!             yes=' '
!            do while ((yes.ne.'y').and.(yes.ne.'n'))
!             write(*,*) 'Would you like to reindex tracer # 1 ->',nqold
!             write(*,*) 'to #',nqmx-nqold+1,'->', nqmx,'   (y or n) ?'
!             read(*,fmt='(a)') yes
!            end do
!            if (yes.eq.'y') then
!              write(*,*) 'OK, let s reindex the tracers'
!              do l=1,llm
!                do j=1,jjp1
!                  do i=1,iip1
!                    do iq=nqmx,nqmx-nqold+1,-1
!                       q(i,j,l,iq)=q(i,j,l,iq-nqmx+nqold)   
!                    end do
!                    do iq=nqmx-nqold,1,-1
!                       q(i,j,l,iq)= 0.
!                    end do
!                  end do
!                end do
!              end do
!            end if 
         end if ! of if ((nqold.lt.nqtot).and.(nqold.ge.1))
      ENDIF ! of IF(nqtot.GE.1)

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
