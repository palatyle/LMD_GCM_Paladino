      SUBROUTINE concvl (iflag_con,iflag_clos,
     .             dtime,paprs,pplay,
     .             t,q,t_wake,q_wake,s_wake,u,v,tra,ntra,
     .             ALE,ALP,work1,work2,
     .             d_t,d_q,d_u,d_v,d_tra,
     .             rain, snow, kbas, ktop, sigd,
     .             cbmf,upwd,dnwd,dnwdbis,Ma,mip,Vprecip,
     .             cape,cin,tvp,Tconv,iflag,
     .             pbase,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr,
     .             qcondc,wd,pmflxr,pmflxs,
     .             da,phi,mp,dd_t,dd_q,lalim_conv,wght_th)
***************************************************************
*                                                             *
* CONCVL                                                      *
*                                                             *
*                                                             *
* written by   : Sandrine Bony-Lena, 17/05/2003, 11.16.04    *
* modified by :                                               *
***************************************************************
*
c
      USE dimphy
      USE infotrac, ONLY : nbtr
      IMPLICIT none
c======================================================================
c Auteur(s): S. Bony-Lena (LMD/CNRS) date: ???
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
c ALP-----input-R-puissance disponible pour soulevement
c
C d_h-----output-R-increment de l'enthalpie potentielle (h)
c d_q-----output-R-increment de la vapeur d'eau
c rain----output-R-la pluie (mm/s)
c snow----output-R-la neige (mm/s)
c upwd----output-R-saturated updraft mass flux (kg/m**2/s)
c dnwd----output-R-saturated downdraft mass flux (kg/m**2/s)
c dnwd0---output-R-unsaturated downdraft mass flux (kg/m**2/s)
c Ma------output-R-adiabatic ascent mass flux (kg/m2/s)
c mip-----output-R-mass flux shed by adiabatic ascent (kg/m2/s)
c Vprecip-output-R-vertical profile of precipitations (kg/m2/s)
c Tconv---output-R-environment temperature seen by convective scheme (K)
c Cape----output-R-CAPE (J/kg)
c Cin ----output-R-CIN  (J/kg)
c Tvp-----output-R-Temperature virtuelle d'une parcelle soulevee
c                  adiabatiquement a partir du niveau 1 (K)
c deltapb-output-R-distance entre LCL et base de la colonne (<0 ; Pa)
c Ice_flag-input-L-TRUE->prise en compte de la thermodynamique de la glace
c dd_t-----output-R-increment de la temperature du aux descentes precipitantes
c dd_q-----output-R-increment de la vapeur d'eau du aux desc precip
c======================================================================
c
#include "dimensions.h"
c
       INTEGER iflag_con,iflag_clos
c
       REAL dtime, paprs(klon,klev+1),pplay(klon,klev)
       REAL t(klon,klev),q(klon,klev),u(klon,klev),v(klon,klev)
       REAL t_wake(klon,klev),q_wake(klon,klev)
       Real s_wake(klon)
       REAL tra(klon,klev,nbtr)
       INTEGER ntra
       REAL work1(klon,klev),work2(klon,klev),ptop2(klon)
       REAL pmflxr(klon,klev+1),pmflxs(klon,klev+1)
       REAL ALE(klon),ALP(klon)
c
       REAL d_t(klon,klev),d_q(klon,klev),d_u(klon,klev),d_v(klon,klev)
       REAL dd_t(klon,klev),dd_q(klon,klev)
       REAL d_tra(klon,klev,nbtr)
       REAL rain(klon),snow(klon)
c
       INTEGER kbas(klon),ktop(klon)
       REAL em_ph(klon,klev+1),em_p(klon,klev)
       REAL upwd(klon,klev),dnwd(klon,klev),dnwdbis(klon,klev)

!!       REAL Ma(klon,klev), mip(klon,klev),Vprecip(klon,klev)     !jyg
       REAL Ma(klon,klev), mip(klon,klev),Vprecip(klon,klev+1)     !jyg

       real da(klon,klev),phi(klon,klev,klev),mp(klon,klev)
       REAL cape(klon),cin(klon),tvp(klon,klev)
       REAL Tconv(klon,klev)
c
cCR:test: on passe lentr et alim_star des thermiques
       INTEGER lalim_conv(klon)
       REAL wght_th(klon,klev)
       REAL em_sig1feed ! sigma at lower bound of feeding layer
       REAL em_sig2feed ! sigma at upper bound of feeding layer
       REAL em_wght(klev) ! weight density determining the feeding mixture
con enleve le save
c       SAVE em_sig1feed,em_sig2feed,em_wght
c
       INTEGER iflag(klon)
       REAL rflag(klon)
       REAL pbase(klon),bbase(klon)
       REAL dtvpdt1(klon,klev),dtvpdq1(klon,klev)
       REAL dplcldt(klon),dplcldr(klon)
       REAL qcondc(klon,klev)
       REAL wd(klon)
       REAL Plim1(klon),Plim2(klon),asupmax(klon,klev)
       REAL supmax0(klon),asupmaxmin(klon)
c
       REAL sigd(klon)
       REAL zx_t,zdelta,zx_qs,zcor
c
!       INTEGER iflag_mix
!       SAVE iflag_mix
       INTEGER noff, minorig
       INTEGER i,k,itra
       REAL qs(klon,klev),qs_wake(klon,klev)
       REAL cbmf(klon)
cLF       SAVE cbmf
cIM/JYG      REAL, SAVE, ALLOCATABLE :: cbmf(:)
ccc$OMP THREADPRIVATE(cbmf)!       
       REAL cbmflast(klon)
       INTEGER ifrst
       SAVE ifrst
       DATA ifrst /0/
c$OMP THREADPRIVATE(ifrst)

c
C     Variables supplementaires liees au bilan d'energie
c      Real paire(klon)
cLF      Real ql(klon,klev)
c      Save paire
cLF      Save ql
cLF      Real t1(klon,klev),q1(klon,klev)
cLF      Save t1,q1
c      Data paire /1./
       REAL, SAVE, ALLOCATABLE :: ql(:,:), q1(:,:), t1(:,:)
c$OMP THREADPRIVATE(ql, q1, t1)
c
C     Variables liees au bilan d'energie et d'enthalpi
      REAL ztsol(klon)
      REAL      h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot
     $        , h_qs_tot, qw_tot, ql_tot, qs_tot , ec_tot
      SAVE      h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot
     $        , h_qs_tot, qw_tot, ql_tot, qs_tot , ec_tot
c$OMP THREADPRIVATE(h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot)
c$OMP THREADPRIVATE(h_qs_tot, qw_tot, ql_tot, qs_tot , ec_tot)
      REAL      d_h_vcol, d_h_dair, d_qt, d_qw, d_ql, d_qs, d_ec
      REAL      d_h_vcol_phy
      REAL      fs_bound, fq_bound
      SAVE      d_h_vcol_phy
c$OMP THREADPRIVATE(d_h_vcol_phy)
      REAL      zero_v(klon)
      CHARACTER*15 ztit
      INTEGER   ip_ebil  ! PRINT level for energy conserv. diag.
      SAVE      ip_ebil
      DATA      ip_ebil/2/
c$OMP THREADPRIVATE(ip_ebil)
      INTEGER   if_ebil ! level for energy conserv. dignostics
      SAVE      if_ebil
      DATA      if_ebil/2/
c$OMP THREADPRIVATE(if_ebil)
c+jld ec_conser
      REAL d_t_ec(klon,klev)    ! tendance du a la conersion Ec -> E thermique
      REAL ZRCPD
c-jld ec_conser
cLF
      INTEGER nloc
      logical, save :: first=.true.
c$OMP THREADPRIVATE(first)
      INTEGER, SAVE :: itap, igout
c$OMP THREADPRIVATE(itap, igout)
c
#include "YOMCST.h"
#include "YOMCST2.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "iniprint.h"
c
      if (first) then
c Allocate some variables LF 04/2008
c
cIM/JYG allocate(cbmf(klon))
        allocate(ql(klon,klev))
        allocate(t1(klon,klev))
        allocate(q1(klon,klev))
        itap=0
        igout=klon/2+1/klon
      endif
c Incrementer le compteur de la physique
      itap   = itap + 1

c    Copy T into Tconv
      DO k = 1,klev
        DO i = 1,klon
          Tconv(i,k) = T(i,k)
        ENDDO
      ENDDO
c
      IF (if_ebil.ge.1) THEN
        DO i=1,klon
          ztsol(i) = t(i,1)
          zero_v(i)=0.
          Do k = 1,klev
            ql(i,k) = 0.
          ENDDO
        END DO
      END IF
c
cym
      snow(:)=0
      
c      IF (ifrst .EQ. 0) THEN
c         ifrst = 1
       if (first) then
         first=.false.
c
C===========================================================================
C    READ IN PARAMETERS FOR THE CLOSURE AND THE MIXING DISTRIBUTION
C===========================================================================
C
      if (iflag_con.eq.3) then
c     CALL cv3_inicp()
      CALL cv3_inip()
      endif
c
C===========================================================================
C    READ IN PARAMETERS FOR CONVECTIVE INHIBITION BY TROPOS. DRYNESS
C===========================================================================
C
cc$$$         open (56,file='supcrit.data')
cc$$$         read (56,*) Supcrit1, Supcrit2
cc$$$         close (56)
c
         print*, 'supcrit1, supcrit2' ,supcrit1, supcrit2
C
C===========================================================================
C      Initialisation pour les bilans d'eau et d'energie
C===========================================================================
         IF (if_ebil.ge.1) d_h_vcol_phy=0.
c
         DO i = 1, klon
          cbmf(i) = 0.
          sigd(i) = 0.
         ENDDO
      ENDIF   !(ifrst .EQ. 0)

      DO k = 1, klev+1
         DO i=1,klon
         em_ph(i,k) = paprs(i,k) / 100.0
         pmflxr(i,k)=0.
         pmflxs(i,k)=0.
      ENDDO
      ENDDO
c
      DO k = 1, klev
         DO i=1,klon
         em_p(i,k) = pplay(i,k) / 100.0
      ENDDO
      ENDDO
c
!
!  Feeding layer
!
      em_sig1feed = 1.
      em_sig2feed = 0.97
c      em_sig2feed = 0.8
! Relative Weight densities
       do k=1,klev
         em_wght(k)=1.
       end do
cCRtest: couche alim des tehrmiques ponderee par a*
c       DO i = 1, klon
c         do k=1,lalim_conv(i)
c         em_wght(k)=wght_th(i,k)
c         print*,'em_wght=',em_wght(k),wght_th(i,k)
c       end do
c      END DO

      if (iflag_con .eq. 4) then
      DO k = 1, klev
        DO i = 1, klon
         zx_t = t(i,k)
         zdelta=MAX(0.,SIGN(1.,rtt-zx_t))
         zx_qs= MIN(0.5 , r2es * FOEEW(zx_t,zdelta)/em_p(i,k)/100.0)
         zcor=1./(1.-retv*zx_qs)
         qs(i,k)=zx_qs*zcor
        ENDDO
        DO i = 1, klon
         zx_t = t_wake(i,k)
         zdelta=MAX(0.,SIGN(1.,rtt-zx_t))
         zx_qs= MIN(0.5 , r2es * FOEEW(zx_t,zdelta)/em_p(i,k)/100.0)
         zcor=1./(1.-retv*zx_qs)
         qs_wake(i,k)=zx_qs*zcor
        ENDDO
      ENDDO
      else ! iflag_con=3 (modif de puristes qui fait la diffce pour la convergence numerique)
      DO k = 1, klev
        DO i = 1, klon
         zx_t = t(i,k)
         zdelta=MAX(0.,SIGN(1.,rtt-zx_t))
         zx_qs= r2es * FOEEW(zx_t,zdelta)/em_p(i,k)/100.0
         zx_qs= MIN(0.5,zx_qs)
         zcor=1./(1.-retv*zx_qs)
         zx_qs=zx_qs*zcor
         qs(i,k)=zx_qs
        ENDDO
        DO i = 1, klon
         zx_t = t_wake(i,k)
         zdelta=MAX(0.,SIGN(1.,rtt-zx_t))
         zx_qs= r2es * FOEEW(zx_t,zdelta)/em_p(i,k)/100.0
         zx_qs= MIN(0.5,zx_qs)
         zcor=1./(1.-retv*zx_qs)
         zx_qs=zx_qs*zcor
         qs_wake(i,k)=zx_qs
        ENDDO
      ENDDO
      endif ! iflag_con
c
C------------------------------------------------------------------

C Main driver for convection:
C               iflag_con=3 -> nvlle version de KE (JYG)
C		iflag_con = 30  -> equivalent to convect3
C		iflag_con = 4  -> equivalent to convect1/2
c
c
      if (iflag_con.eq.30) then

      CALL cv_driver(klon,klev,klev+1,ntra,iflag_con,
     :              t,q,qs,u,v,tra,
     $              em_p,em_ph,iflag,
     $              d_t,d_q,d_u,d_v,d_tra,rain,
!!     $              pmflxr,cbmf,work1,work2,           !jyg
     $              Vprecip,cbmf,work1,work2,            !jyg
     $              kbas,ktop,
     $              dtime,Ma,upwd,dnwd,dnwdbis,qcondc,wd,cape,
     $              da,phi,mp)
c
      DO i = 1,klon
        cbmf(i) = Ma(i,kbas(i))
      ENDDO
c
      else

cLF   necessary for gathered fields
      nloc=klon
      CALL cva_driver(klon,klev,klev+1,ntra,nloc,
     $              iflag_con,iflag_mix,iflag_clos,dtime,
     :              t,q,qs,t_wake,q_wake,qs_wake,s_wake,u,v,tra,
     $              em_p,em_ph,
     .              ALE,ALP,
     .              em_sig1feed,em_sig2feed,em_wght,
     .              iflag,d_t,d_q,d_u,d_v,d_tra,rain,kbas,ktop,
     $              cbmf,work1,work2,ptop2,sigd,
     $              Ma,mip,Vprecip,upwd,dnwd,dnwdbis,qcondc,wd,
     $              cape,cin,tvp,
     $              dd_t,dd_q,Plim1,Plim2,asupmax,supmax0,
     $              asupmaxmin,lalim_conv)
      endif  
C------------------------------------------------------------------

      DO i = 1,klon
        rain(i) = rain(i)/86400.
        rflag(i)=iflag(i)
      ENDDO

      DO k = 1, klev
        DO i = 1, klon
           d_t(i,k) = dtime*d_t(i,k)
           d_q(i,k) = dtime*d_q(i,k)
           d_u(i,k) = dtime*d_u(i,k)
           d_v(i,k) = dtime*d_v(i,k)
        ENDDO
      ENDDO
c
       if (iflag_con.eq.30) then
       DO itra = 1,ntra
        DO k = 1, klev
         DO i = 1, klon
            d_tra(i,k,itra) =dtime*d_tra(i,k,itra) 
         ENDDO
        ENDDO
       ENDDO 
       endif

      DO k = 1, klev
        DO i = 1, klon
          t1(i,k) = t(i,k)+ d_t(i,k)
          q1(i,k) = q(i,k)+ d_q(i,k)
        ENDDO
      ENDDO
c                                                  !jyg
c--Separation neige/pluie (pour diagnostics)       !jyg
      DO k = 1, klev                               !jyg
      DO i = 1, klon                               !jyg
       IF (t1(i,k).LT.RTT) THEN                    !jyg
         pmflxs(i,k)=Vprecip(i,k)                  !jyg
       ELSE                                        !jyg
         pmflxr(i,k)=Vprecip(i,k)                  !jyg
       ENDIF                                       !jyg
      ENDDO                                        !jyg
      ENDDO                                        !jyg
c
cc      IF (if_ebil.ge.2) THEN
cc        ztit='after convect'
cc        CALL diagetpq(paire,ztit,ip_ebil,2,2,dtime
cc     e      , t1,q1,ql,qs,u,v,paprs,pplay
cc     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
cc         call diagphy(paire,ztit,ip_ebil
cc     e      , zero_v, zero_v, zero_v, zero_v, zero_v
cc     e      , zero_v, rain, zero_v, ztsol
cc     e      , d_h_vcol, d_qt, d_ec
cc     s      , fs_bound, fq_bound )
cc      END IF
C
c
c les traceurs ne sont pas mis dans cette version de convect4:
      if (iflag_con.eq.4) then
       DO itra = 1,ntra
        DO k = 1, klev
         DO i = 1, klon
            d_tra(i,k,itra) = 0.
         ENDDO
        ENDDO
       ENDDO
      endif
c     print*, 'concvl->: dd_t,dd_q ',dd_t(1,1),dd_q(1,1)

        DO k = 1, klev
         DO i = 1, klon
            dtvpdt1(i,k) = 0.
            dtvpdq1(i,k) = 0.
         ENDDO
        ENDDO
        DO i = 1, klon
           dplcldt(i) = 0.
           dplcldr(i) = 0.
        ENDDO
c
       if(prt_level.GE.20) THEN
       DO k=1,klev
!       print*,'physiq apres_add_con i k it d_u d_v d_t d_q qdl0',igout
!    .,k,itap,d_u_con(igout,k) ,d_v_con(igout,k), d_t_con(igout,k),
!    .d_q_con(igout,k),dql0(igout,k)
!      print*,'phys apres_add_con itap Ma cin ALE ALP wak t q undi t q'
!    .,itap,Ma(igout,k),cin(igout),ALE(igout), ALP(igout),
!    . t_wake(igout,k),q_wake(igout,k),t_undi(igout,k),q_undi(igout,k)
!      print*,'phy apres_add_con itap CON rain snow EMA wk1 wk2 Vpp mip'
!    .,itap,rain_con(igout),snow_con(igout),ema_work1(igout,k),
!    .ema_work2(igout,k),Vprecip(igout,k), mip(igout,k)
!      print*,'phy apres_add_con itap upwd dnwd dnwd0 cape tvp Tconv '
!    .,itap,upwd(igout,k),dnwd(igout,k),dnwd0(igout,k),cape(igout),
!    .tvp(igout,k),Tconv(igout,k)
!      print*,'phy apres_add_con itap dtvpdt dtvdq dplcl dplcldr qcondc'
!    .,itap,dtvpdt1(igout,k),dtvpdq1(igout,k),dplcldt(igout),
!    .dplcldr(igout),qcondc(igout,k)
!      print*,'phy apres_add_con itap wd pmflxr Kpmflxr Kp1 Kpmflxs Kp1'
!    .,itap,wd(igout),pmflxr(igout,k),pmflxr(igout,k+1),pmflxs(igout,k)
!    .,pmflxs(igout,k+1)
!      print*,'phy apres_add_con itap da phi mp ftd fqd lalim wgth',
!    .itap,da(igout,k),phi(igout,k,k),mp(igout,k),ftd(igout,k),
!    . fqd(igout,k),lalim_conv(igout),wght_th(igout,k)
      ENDDO
      endif !(prt_level.EQ.20) THEN
c
      RETURN
      END
 
