!
! $Header$
!
      SUBROUTINE conemav (dtime,paprs,pplay,t,q,u,v,tra,ntra,
     .             work1,work2,d_t,d_q,d_u,d_v,d_tra,
     .             rain, snow, kbas, ktop,
     .             upwd,dnwd,dnwdbis,Ma,cape,tvp,iflag,
     .             pbase,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr)
 
c
      USE dimphy
      USE infotrac, ONLY : nbtr
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: schema de convection de Emanuel (1991) interface
c======================================================================
c Arguments:
c dtime--input-R-pas d'integration (s)
c s-------input-R-la valeur "s" pour chaque couche
c sigs----input-R-la valeur "sigma" de chaque couche
c sig-----input-R-la valeur de "sigma" pour chaque niveau
c psolpa--input-R-la pression au sol (en Pa)
C pskapa--input-R-exponentiel kappa de psolpa
c h-------input-R-enthalpie potentielle (Cp*T/P**kappa)
c q-------input-R-vapeur d'eau (en kg/kg)
c
c work*: input et output: deux variables de travail,
c                            on peut les mettre a 0 au debut
c ALE-----input-R-energie disponible pour soulevement
c
C d_h-----output-R-increment de l'enthalpie potentielle (h)
c d_q-----output-R-increment de la vapeur d'eau
c rain----output-R-la pluie (mm/s)
c snow----output-R-la neige (mm/s)
c upwd----output-R-saturated updraft mass flux (kg/m**2/s)
c dnwd----output-R-saturated downdraft mass flux (kg/m**2/s)
c dnwd0---output-R-unsaturated downdraft mass flux (kg/m**2/s)
c Cape----output-R-CAPE (J/kg)
c Tvp-----output-R-Temperature virtuelle d'une parcelle soulevee
c                  adiabatiquement a partir du niveau 1 (K)
c deltapb-output-R-distance entre LCL et base de la colonne (<0 ; Pa)
c Ice_flag-input-L-TRUE->prise en compte de la thermodynamique de la glace
c======================================================================
c
#include "dimensions.h"
c
c
       REAL dtime, paprs(klon,klev+1),pplay(klon,klev)
       REAL t(klon,klev),q(klon,klev),u(klon,klev),v(klon,klev)
       REAL tra(klon,klev,nbtr)
       INTEGER ntra
       REAL work1(klon,klev),work2(klon,klev)
c
       REAL d_t(klon,klev),d_q(klon,klev),d_u(klon,klev),d_v(klon,klev)
       REAL d_tra(klon,klev,nbtr)
       REAL rain(klon),snow(klon)
c
       INTEGER kbas(klon),ktop(klon)
       REAL em_ph(klon,klev+1),em_p(klon,klev)
       REAL upwd(klon,klev),dnwd(klon,klev),dnwdbis(klon,klev)
       REAL Ma(klon,klev),cape(klon),tvp(klon,klev)
       INTEGER iflag(klon)
       REAL rflag(klon)
       REAL pbase(klon),bbase(klon)
       REAL dtvpdt1(klon,klev),dtvpdq1(klon,klev)
       REAL dplcldt(klon),dplcldr(klon)
c
       REAL zx_t,zdelta,zx_qs,zcor
c
       INTEGER noff, minorig
       INTEGER i,k,itra
       REAL qs(klon,klev)
       REAL,ALLOCATABLE,SAVE :: cbmf(:)
c$OMP THREADPRIVATE(cbmf)
       INTEGER ifrst
       SAVE ifrst
       DATA ifrst /0/
c$OMP THREADPRIVATE(ifrst)
#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
c
c
      IF (ifrst .EQ. 0) THEN
         ifrst = 1
	 allocate(cbmf(klon))
         DO i = 1, klon
          cbmf(i) = 0.
         ENDDO
      ENDIF

      DO k = 1, klev+1
         DO i=1,klon
         em_ph(i,k) = paprs(i,k) / 100.0
      ENDDO
      ENDDO
c
      DO k = 1, klev
         DO i=1,klon
         em_p(i,k) = pplay(i,k) / 100.0
      ENDDO
      ENDDO

c
      DO k = 1, klev
        DO i = 1, klon
         zx_t = t(i,k)
         zdelta=MAX(0.,SIGN(1.,rtt-zx_t))
         zx_qs= MIN(0.5 , r2es * FOEEW(zx_t,zdelta)/em_p(i,k)/100.0)
         zcor=1./(1.-retv*zx_qs)
         qs(i,k)=zx_qs*zcor
        ENDDO
      ENDDO
c
      noff = 2
      minorig = 2
      CALL convect1(klon,klev,klev+1,noff,minorig,t,q,qs,u,v,
     $              em_p,em_ph,iflag,
     $              d_t,d_q,d_u,d_v,rain,cbmf,dtime,Ma)
c
      DO i = 1,klon
        rain(i) = rain(i)/86400.
        rflag(i)=iflag(i)
      ENDDO
c      call dump2d(iim,jjm-1,rflag(2:klon-1),'FLAG CONVECTION   ')
c     if (klon.eq.1) then
c        print*,'IFLAG ',iflag
c     else
c        write(*,'(96i1)') (iflag(i),i=2,klon-1)
c     endif
      DO k = 1, klev
        DO i = 1, klon
           d_t(i,k) = dtime*d_t(i,k)
           d_q(i,k) = dtime*d_q(i,k)
           d_u(i,k) = dtime*d_u(i,k)
           d_v(i,k) = dtime*d_v(i,k)
        ENDDO
        DO itra = 1,ntra
          DO i = 1, klon
            d_tra(i,k,itra) = 0.
          ENDDO
        ENDDO
      ENDDO
 
c
c
c
      RETURN
      END
 
