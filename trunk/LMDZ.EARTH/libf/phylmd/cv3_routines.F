!
! $Id: cv3_routines.F 1467 2010-12-17 11:17:22Z idelkadi $
!
c
c
      SUBROUTINE cv3_param(nd,delt)
      implicit none

c------------------------------------------------------------
c Set parameters for convectL for iflag_con = 3 
c------------------------------------------------------------

C
C   ***  PBCRIT IS THE CRITICAL CLOUD DEPTH (MB) BENEATH WHICH THE ***
C   ***      PRECIPITATION EFFICIENCY IS ASSUMED TO BE ZERO     ***
C   ***  PTCRIT IS THE CLOUD DEPTH (MB) ABOVE WHICH THE PRECIP. ***     
C   ***            EFFICIENCY IS ASSUMED TO BE UNITY            ***
C   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
C   ***  SPFAC IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE ***     
C   ***                        OF CLOUD                         ***
C
C [TAU: CHARACTERISTIC TIMESCALE USED TO COMPUTE ALPHA & BETA]
C   ***    ALPHA AND BETA ARE PARAMETERS THAT CONTROL THE RATE OF ***
C   ***                 APPROACH TO QUASI-EQUILIBRIUM           ***
C   ***    (THEIR STANDARD VALUES ARE 1.0 AND 0.96, RESPECTIVELY) ***
C   ***           (BETA MUST BE LESS THAN OR EQUAL TO 1)        ***
C
C   ***    DTCRIT IS THE CRITICAL BUOYANCY (K) USED TO ADJUST THE ***
C   ***                 APPROACH TO QUASI-EQUILIBRIUM           ***
C   ***                     IT MUST BE LESS THAN 0              ***

#include "cv3param.h"
#include "conema3.h"

      integer nd
      real delt ! timestep (seconds)

      CHARACTER (LEN=20) :: modname='cv3_param'
      CHARACTER (LEN=80) :: abort_message

c noff: integer limit for convection (nd-noff)
c minorig: First level of convection

c -- limit levels for convection:

      noff    = 1
      minorig = 1
      nl=nd-noff
      nlp=nl+1
      nlm=nl-1

c -- "microphysical" parameters:
       sigdz=0.01
c      sigd=0.003
c     sigd   = 0.01
cCR:test sur la fraction des descentes precipitantes
      spfac  = 0.15
      pbcrit = 150.0
      ptcrit = 500.0
cIM lu dans physiq.def via conf_phys.F90     epmax  = 0.993

      omtrain = 45.0 ! used also for snow (no disctinction rain/snow)

c -- misc:

      dtovsh = -0.2 ! dT for overshoot
      dpbase = -40. ! definition cloud base (400m above LCL)
ccc      dttrig = 5.   ! (loose) condition for triggering
      dttrig = 10.   ! (loose) condition for triggering

c -- rate of approach to quasi-equilibrium:

      dtcrit = -2.0
c      dtcrit = -5.0
c      tau    = 3000.
      tau = 1800.
cc      tau=8000.
      beta   = 1.0 - delt/tau
      alpha1 = 1.5e-3
      alpha  = alpha1 * delt/tau
c increase alpha to compensate W decrease:
      alpha  = alpha*1.5

c -- interface cloud parameterization:

      delta=0.01  ! cld

c -- interface with boundary-layer (gust factor): (sb)

      betad=10.0   ! original value (from convect 4.3)

      OPEN(99,file='conv_param.data',status='old',
     $          form='formatted',err=9999)
      READ(99,*,end=9998) dpbase
      READ(99,*,end=9998) pbcrit
      READ(99,*,end=9998) ptcrit
      READ(99,*,end=9998) sigdz
      READ(99,*,end=9998) spfac
9998  Continue
      CLOSE(99)
9999  Continue
      WRITE(*,*)'dpbase=',dpbase
      WRITE(*,*)'pbcrit=',pbcrit
      WRITE(*,*)'ptcrit=',ptcrit
      WRITE(*,*)'sigdz=',sigdz
      WRITE(*,*)'spfac=',spfac

      return
      end

      SUBROUTINE cv3_prelim(len,nd,ndp1,t,q,p,ph
     :                    ,lv,cpn,tv,gz,h,hm,th)
      implicit none

!=====================================================================
! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
! "ori": from convect4.3 (vectorized)
! "convect3": to be exactly consistent with convect3
!=====================================================================

c inputs:
      integer len, nd, ndp1
      real t(len,nd), q(len,nd), p(len,nd), ph(len,ndp1)

c outputs:
      real lv(len,nd), cpn(len,nd), tv(len,nd)
      real gz(len,nd), h(len,nd), hm(len,nd)
      real th(len,nd)

c local variables:
      integer k, i
      real rdcp
      real tvx,tvy ! convect3
      real cpx(len,nd)

#include "cvthermo.h"
#include "cv3param.h"


c ori      do 110 k=1,nlp
! abderr     do 110 k=1,nl ! convect3
       do 110 k=1,nlp
      
        do 100 i=1,len
cdebug          lv(i,k)= lv0-clmcpv*(t(i,k)-t0)
          lv(i,k)= lv0-clmcpv*(t(i,k)-273.15)
          cpn(i,k)=cpd*(1.0-q(i,k))+cpv*q(i,k)
          cpx(i,k)=cpd*(1.0-q(i,k))+cl*q(i,k)
c ori          tv(i,k)=t(i,k)*(1.0+q(i,k)*epsim1)
          tv(i,k)=t(i,k)*(1.0+q(i,k)/eps-q(i,k))
          rdcp=(rrd*(1.-q(i,k))+q(i,k)*rrv)/cpn(i,k)
          th(i,k)=t(i,k)*(1000.0/p(i,k))**rdcp
 100    continue
 110  continue
c
c gz = phi at the full levels (same as p).
c
      do 120 i=1,len
        gz(i,1)=0.0
 120  continue
c ori      do 140 k=2,nlp
      do 140 k=2,nl ! convect3
        do 130 i=1,len
        tvx=t(i,k)*(1.+q(i,k)/eps-q(i,k))       !convect3
        tvy=t(i,k-1)*(1.+q(i,k-1)/eps-q(i,k-1)) !convect3
        gz(i,k)=gz(i,k-1)+0.5*rrd*(tvx+tvy)     !convect3
     &          *(p(i,k-1)-p(i,k))/ph(i,k)      !convect3
c
cc        print *,' gz(',k,')',gz(i,k),' tvx',tvx,' tvy ',tvy
c
c ori         gz(i,k)=gz(i,k-1)+hrd*(tv(i,k-1)+tv(i,k))
c ori    &         *(p(i,k-1)-p(i,k))/ph(i,k)
 130    continue
 140  continue
c
c h  = phi + cpT (dry static energy).
c hm = phi + cp(T-Tbase)+Lq
c
c ori      do 170 k=1,nlp
      do 170 k=1,nl ! convect3
        do 160 i=1,len
          h(i,k)=gz(i,k)+cpn(i,k)*t(i,k)
          hm(i,k)=gz(i,k)+cpx(i,k)*(t(i,k)-t(i,1))+lv(i,k)*q(i,k)
 160    continue
 170  continue

      return
      end

      SUBROUTINE cv3_feed(len,nd,t,q,u,v,p,ph,hm,gz
     :                  ,p1feed,p2feed,wght
     :                  ,wghti,tnk,thnk,qnk,qsnk,unk,vnk
     :                  ,cpnk,hnk,nk,icb,icbmax,iflag,gznk,plcl)
      implicit none

C================================================================
C Purpose: CONVECTIVE FEED
C
C Main differences with cv_feed:
C   - ph added in input
C	- here, nk(i)=minorig
C	- icb defined differently (plcl compared with ph instead of p)
C
C Main differences with convect3:
C 	- we do not compute dplcldt and dplcldr of CLIFT anymore 
C	- values iflag different (but tests identical)
C   - A,B explicitely defined (!...)
C================================================================

#include "cv3param.h"
#include "cvthermo.h"

c inputs:
	  integer len, nd
      real t(len,nd), q(len,nd), p(len,nd)
      real u(len,nd), v(len,nd)
      real hm(len,nd), gz(len,nd)
      real ph(len,nd+1)
      real p1feed(len)
c,  wght(len)
      real wght(nd)
c input-output
      real p2feed(len)
c outputs:
	  integer iflag(len), nk(len), icb(len), icbmax
c      real   wghti(len)
      real wghti(len,nd)
      real   tnk(len), thnk(len), qnk(len), qsnk(len)
      real   unk(len), vnk(len)
      real   cpnk(len), hnk(len), gznk(len)
      real   plcl(len)

c local variables:
      integer i, k, iter, niter
      integer ihmin(len)
      real work(len)
      real pup(len),plo(len),pfeed(len)
      real plclup(len),plcllo(len),plclfeed(len)
      real posit(len)
      logical nocond(len)
!
!-------------------------------------------------------------------
! --- Origin level of ascending parcels for convect3:
!-------------------------------------------------------------------

         do 220 i=1,len
          nk(i)=minorig
          gznk(i)=gz(i,nk(i))
  220    continue
!
!-------------------------------------------------------------------
! --- Adjust feeding layer thickness so that lifting up to the top of
! --- the feeding layer does not induce condensation (i.e. so that
! --- plcl < p2feed).
! --- Method : iterative secant method.
!-------------------------------------------------------------------
!
c 1- First bracketing of the solution : ph(nk+1), p2feed
c
c 1.a- LCL associated to p2feed
      do i = 1,len
        pup(i) = p2feed(i)
      enddo
         call cv3_vertmix(len,nd,iflag,p1feed,pup,p,ph
     i              ,t,q,u,v,wght
     o              ,wghti,nk,tnk,thnk,qnk,qsnk,unk,vnk,plclup)
c 1.b- LCL associated to ph(nk+1)
      do i = 1,len
        plo(i) = ph(i,nk(i)+1)
      enddo
         call cv3_vertmix(len,nd,iflag,p1feed,plo,p,ph
     i              ,t,q,u,v,wght
     o              ,wghti,nk,tnk,thnk,qnk,qsnk,unk,vnk,plcllo)
c 2- Iterations
      niter = 5
      do iter = 1,niter
        do i = 1,len
          plcllo(i) = min(plo(i),plcllo(i))
          plclup(i) = max(pup(i),plclup(i))
          nocond(i) = plclup(i).le.pup(i)
        enddo
        do i = 1,len
          if(nocond(i)) then
             pfeed(i)=pup(i)
          else
             pfeed(i) = (pup(i)*(plo(i)-plcllo(i))+
     :                plo(i)*(plclup(i)-pup(i)))/
     :            (plo(i)-plcllo(i)+plclup(i)-pup(i))
          endif
        enddo
         call cv3_vertmix(len,nd,iflag,p1feed,pfeed,p,ph
     i              ,t,q,u,v,wght
     o              ,wghti,nk,tnk,thnk,qnk,qsnk,unk,vnk,plclfeed)
        do i = 1,len
          posit(i) = (sign(1.,plclfeed(i)-pfeed(i))+1.)*0.5
          if (plclfeed(i) .eq. pfeed(i)) posit(i) = 1.
c- posit = 1 when lcl is below top of feeding layer (plclfeed>pfeed)
c-               => pup=pfeed
c- posit = 0 when lcl is above top of feeding layer (plclfeed<pfeed)
c-               => plo=pfeed
          pup(i) = posit(i)*pfeed(i) + (1.-posit(i))*pup(i)
          plo(i) = (1.-posit(i))*pfeed(i) + posit(i)*plo(i)
          plclup(i) = posit(i)*plclfeed(i) + (1.-posit(i))*plclup(i)
          plcllo(i) = (1.-posit(i))*plclfeed(i) + posit(i)*plcllo(i)
        enddo
      enddo       !  iter
      do i = 1,len
        p2feed(i) = pfeed(i)
        plcl(i) = plclfeed(i)
      enddo
!
      do 175 i=1,len
         cpnk(i)=cpd*(1.0-qnk(i))+cpv*qnk(i)
         hnk(i)=gz(i,1)+cpnk(i)*tnk(i)
 175  continue
!
!-------------------------------------------------------------------
! --- Check whether parcel level temperature and specific humidity
! --- are reasonable
!-------------------------------------------------------------------
       do 250 i=1,len
       if( (     ( tnk(i).lt.250.0    )
     &       .or.( qnk(i).le.0.0      ) )
     &   .and.
     &       ( iflag(i).eq.0) ) iflag(i)=7
 250   continue
c
!-------------------------------------------------------------------
! --- Calculate first level above lcl (=icb)
!-------------------------------------------------------------------

c@      do 270 i=1,len
c@       icb(i)=nlm
c@ 270  continue
c@c
c@      do 290 k=minorig,nl
c@        do 280 i=1,len
c@          if((k.ge.(nk(i)+1)).and.(p(i,k).lt.plcl(i)))
c@     &    icb(i)=min(icb(i),k)
c@ 280    continue
c@ 290  continue
c@c
c@      do 300 i=1,len
c@        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
c@ 300  continue

      do 270 i=1,len
       icb(i)=nlm
 270  continue
c
c la modification consiste a comparer plcl a ph et non a p:
c icb est defini par :  ph(icb)<plcl<ph(icb-1)
c@      do 290 k=minorig,nl
      do 290 k=3,nl-1 ! modif pour que icb soit sup/egal a 2
        do 280 i=1,len
          if( ph(i,k).lt.plcl(i) ) icb(i)=min(icb(i),k)
 280    continue
 290  continue
c

c     print*,'icb dans cv3_feed '
c     write(*,'(64i2)') icb(2:len-1)
c     call dump2d(64,43,'plcl dans cv3_feed ',plcl(2:len-1))

      do 300 i=1,len
c@        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
        if((icb(i).eq.nlm).and.(iflag(i).eq.0))iflag(i)=9
 300  continue

      do 400 i=1,len
        icb(i) = icb(i)-1 ! icb sup ou egal a 2
 400  continue
c
c Compute icbmax.
c
      icbmax=2
      do 310 i=1,len
c!        icbmax=max(icbmax,icb(i))
       if (iflag(i).lt.7) icbmax=max(icbmax,icb(i)) ! sb Jun7th02
 310  continue

      return
      end

      SUBROUTINE cv3_undilute1(len,nd,t,qs,gz,plcl,p,icb,tnk,qnk,gznk
     :                       ,tp,tvp,clw,icbs)
      implicit none

!----------------------------------------------------------------
! Equivalent de TLIFT entre NK et ICB+1 inclus
!
! Differences with convect4:
!		- specify plcl in input
!       - icbs is the first level above LCL (may differ from icb)
!       - in the iterations, used x(icbs) instead x(icb)
!       - many minor differences in the iterations
!		- tvp is computed in only one time
!		- icbs: first level above Plcl (IMIN de TLIFT) in output
!       - if icbs=icb, compute also tp(icb+1),tvp(icb+1) & clw(icb+1)
!----------------------------------------------------------------

#include "cvthermo.h"
#include "cv3param.h"

c inputs:
      integer len, nd
      integer icb(len)
      real t(len,nd), qs(len,nd), gz(len,nd)
      real tnk(len), qnk(len), gznk(len)
      real p(len,nd)
      real plcl(len) ! convect3

c outputs:
      real tp(len,nd), tvp(len,nd), clw(len,nd)

c local variables:
      integer i, k
      integer icb1(len), icbs(len), icbsmax2 ! convect3
      real tg, qg, alv, s, ahg, tc, denom, es, rg
      real ah0(len), cpp(len)
      real ticb(len), gzicb(len)
      real qsicb(len) ! convect3
      real cpinv(len) ! convect3

!-------------------------------------------------------------------
! --- Calculates the lifted parcel virtual temperature at nk,
! --- the actual temperature, and the adiabatic
! --- liquid water content. The procedure is to solve the equation.
!     cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.
!-------------------------------------------------------------------

c
c   ***  Calculate certain parcel quantities, including static energy   ***
c
      do 330 i=1,len
        ah0(i)=(cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i)
     &         +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15))+gznk(i)
        cpp(i)=cpd*(1.-qnk(i))+qnk(i)*cpv
        cpinv(i)=1./cpp(i)
 330  continue
c
c   ***   Calculate lifted parcel quantities below cloud base   ***
c
        do i=1,len                      !convect3
         icb1(i)=MAX(icb(i),2)          !convect3
         icb1(i)=MIN(icb(i),nl)         !convect3
c if icb is below LCL, start loop at ICB+1:
c (icbs est le premier niveau au-dessus du LCL)
         icbs(i)=icb1(i)                !convect3
         if (plcl(i).lt.p(i,icb1(i))) then
             icbs(i)=MIN(icbs(i)+1,nl)  !convect3
         endif
        enddo                           !convect3

        do i=1,len                      !convect3
         ticb(i)=t(i,icbs(i))           !convect3
         gzicb(i)=gz(i,icbs(i))         !convect3
         qsicb(i)=qs(i,icbs(i))         !convect3
        enddo                           !convect3

c
c Re-compute icbsmax (icbsmax2):        !convect3
c                                       !convect3
      icbsmax2=2                        !convect3
      do 310 i=1,len                    !convect3
        icbsmax2=max(icbsmax2,icbs(i))  !convect3
 310  continue                          !convect3

c initialization outputs:

      do k=1,icbsmax2     ! convect3
       do i=1,len         ! convect3
        tp(i,k)  = 0.0    ! convect3
        tvp(i,k) = 0.0    ! convect3
        clw(i,k) = 0.0    ! convect3
       enddo              ! convect3
      enddo               ! convect3

c tp and tvp below cloud base:

        do 350 k=minorig,icbsmax2-1
          do 340 i=1,len
           tp(i,k)=tnk(i)-(gz(i,k)-gznk(i))*cpinv(i)
           tvp(i,k)=tp(i,k)*(1.+qnk(i)/eps-qnk(i)) !whole thing (convect3)
  340     continue
  350   continue
c
c    ***  Find lifted parcel quantities above cloud base    ***
c
        do 360 i=1,len
         tg=ticb(i)
c ori         qg=qs(i,icb(i))
         qg=qsicb(i) ! convect3
cdebug         alv=lv0-clmcpv*(ticb(i)-t0)
         alv=lv0-clmcpv*(ticb(i)-273.15)
c
c First iteration.
c
c ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
          s=cpd*(1.-qnk(i))+cl*qnk(i)         ! convect3
     :      +alv*alv*qg/(rrv*ticb(i)*ticb(i)) ! convect3
          s=1./s
c ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i) ! convect3
          tg=tg+s*(ah0(i)-ahg)
c ori          tg=max(tg,35.0)
cdebug          tc=tg-t0
          tc=tg-273.15
          denom=243.5+tc
          denom=MAX(denom,1.0) ! convect3
c ori          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
c ori          else
c ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
c ori          endif
c ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
          qg=eps*es/(p(i,icbs(i))-es*(1.-eps))
c
c Second iteration.
c

c ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
c ori          s=1./s
c ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i) ! convect3
          tg=tg+s*(ah0(i)-ahg)
c ori          tg=max(tg,35.0)
cdebug          tc=tg-t0
          tc=tg-273.15
          denom=243.5+tc
          denom=MAX(denom,1.0) ! convect3
c ori          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
c ori          else
c ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
c ori          end if
c ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
          qg=eps*es/(p(i,icbs(i))-es*(1.-eps))

         alv=lv0-clmcpv*(ticb(i)-273.15)

c ori c approximation here:
c ori         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
c ori     &   -gz(i,icb(i))-alv*qg)/cpd

c convect3: no approximation:
         tp(i,icbs(i))=(ah0(i)-gz(i,icbs(i))-alv*qg)
     :                /(cpd+(cl-cpd)*qnk(i))

c ori         clw(i,icb(i))=qnk(i)-qg
c ori         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
         clw(i,icbs(i))=qnk(i)-qg
         clw(i,icbs(i))=max(0.0,clw(i,icbs(i)))

         rg=qg/(1.-qnk(i))
c ori         tvp(i,icb(i))=tp(i,icb(i))*(1.+rg*epsi)
c convect3: (qg utilise au lieu du vrai mixing ratio rg)
         tvp(i,icbs(i))=tp(i,icbs(i))*(1.+qg/eps-qnk(i)) !whole thing

  360   continue
c
c ori      do 380 k=minorig,icbsmax2
c ori       do 370 i=1,len
c ori         tvp(i,k)=tvp(i,k)-tp(i,k)*qnk(i)
c ori 370   continue
c ori 380  continue
c

c -- The following is only for convect3:
c
c * icbs is the first level above the LCL:
c    if plcl<p(icb), then icbs=icb+1
c    if plcl>p(icb), then icbs=icb
c
c * the routine above computes tvp from minorig to icbs (included).
c
c * to compute buoybase (in cv3_trigger.F), both tvp(icb) and tvp(icb+1)
c    must be known. This is the case if icbs=icb+1, but not if icbs=icb.
c
c * therefore, in the case icbs=icb, we compute tvp at level icb+1
c   (tvp at other levels will be computed in cv3_undilute2.F)
c

        do i=1,len
         ticb(i)=t(i,icb(i)+1)
         gzicb(i)=gz(i,icb(i)+1)
         qsicb(i)=qs(i,icb(i)+1)
        enddo

        do 460 i=1,len
         tg=ticb(i)
         qg=qsicb(i) ! convect3
cdebug         alv=lv0-clmcpv*(ticb(i)-t0)
         alv=lv0-clmcpv*(ticb(i)-273.15)
c
c First iteration.
c
c ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
          s=cpd*(1.-qnk(i))+cl*qnk(i)         ! convect3
     :      +alv*alv*qg/(rrv*ticb(i)*ticb(i)) ! convect3
          s=1./s
c ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i) ! convect3
          tg=tg+s*(ah0(i)-ahg)
c ori          tg=max(tg,35.0)
cdebug          tc=tg-t0
          tc=tg-273.15
          denom=243.5+tc
          denom=MAX(denom,1.0) ! convect3
c ori          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
c ori          else
c ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
c ori          endif
c ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
          qg=eps*es/(p(i,icb(i)+1)-es*(1.-eps))
c
c Second iteration.
c

c ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
c ori          s=1./s
c ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i) ! convect3
          tg=tg+s*(ah0(i)-ahg)
c ori          tg=max(tg,35.0)
cdebug          tc=tg-t0
          tc=tg-273.15
          denom=243.5+tc
          denom=MAX(denom,1.0) ! convect3
c ori          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
c ori          else
c ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
c ori          end if
c ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
          qg=eps*es/(p(i,icb(i)+1)-es*(1.-eps))

         alv=lv0-clmcpv*(ticb(i)-273.15)

c ori c approximation here:
c ori         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
c ori     &   -gz(i,icb(i))-alv*qg)/cpd

c convect3: no approximation:
         tp(i,icb(i)+1)=(ah0(i)-gz(i,icb(i)+1)-alv*qg)
     :                /(cpd+(cl-cpd)*qnk(i))

c ori         clw(i,icb(i))=qnk(i)-qg
c ori         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
         clw(i,icb(i)+1)=qnk(i)-qg
         clw(i,icb(i)+1)=max(0.0,clw(i,icb(i)+1))

         rg=qg/(1.-qnk(i))
c ori         tvp(i,icb(i))=tp(i,icb(i))*(1.+rg*epsi)
c convect3: (qg utilise au lieu du vrai mixing ratio rg)
         tvp(i,icb(i)+1)=tp(i,icb(i)+1)*(1.+qg/eps-qnk(i)) !whole thing

  460   continue

      return
      end

      SUBROUTINE cv3_trigger(len,nd,icb,plcl,p,th,tv,tvp,thnk,
     o                pbase,buoybase,iflag,sig,w0)
      implicit none

!-------------------------------------------------------------------
! --- TRIGGERING
!
!	- computes the cloud base
!   - triggering (crude in this version)
!	- relaxation of sig and w0 when no convection
!
!	Caution1: if no convection, we set iflag=4 
!              (it used to be 0 in convect3)
!
!	Caution2: at this stage, tvp (and thus buoy) are know up 
!             through icb only!
! -> the buoyancy below cloud base not (yet) set to the cloud base buoyancy
!-------------------------------------------------------------------

#include "cv3param.h"

c input:
      integer len, nd
      integer icb(len)
      real plcl(len), p(len,nd)
      real th(len,nd), tv(len,nd), tvp(len,nd)
      real thnk(len)

c output:
      real pbase(len), buoybase(len)

c input AND output:
      integer iflag(len)
      real sig(len,nd), w0(len,nd)

c local variables:
      integer i,k
      real tvpbase, tvbase, tdif, ath, ath1

c
c ***   set cloud base buoyancy at (plcl+dpbase) level buoyancy
c
      do 100 i=1,len
       pbase(i) = plcl(i) + dpbase
       tvpbase = tvp(i,icb(i))*(pbase(i)-p(i,icb(i)+1))
     :                        /(p(i,icb(i))-p(i,icb(i)+1))
     :         + tvp(i,icb(i)+1)*(p(i,icb(i))-pbase(i))
     :                          /(p(i,icb(i))-p(i,icb(i)+1))
       tvbase = tv(i,icb(i))*(pbase(i)-p(i,icb(i)+1))
     :                      /(p(i,icb(i))-p(i,icb(i)+1))
     :        + tv(i,icb(i)+1)*(p(i,icb(i))-pbase(i))
     :                        /(p(i,icb(i))-p(i,icb(i)+1))
       buoybase(i) = tvpbase - tvbase
100   continue 

c
c   ***   make sure that column is dry adiabatic between the surface  ***
c   ***    and cloud base, and that lifted air is positively buoyant  ***
c   ***                         at cloud base                         ***
c   ***       if not, return to calling program after resetting       ***
c   ***                        sig(i) and w0(i)                       ***
c

c oct3      do 200 i=1,len
c oct3
c oct3       tdif = buoybase(i)
c oct3       ath1 = th(i,1)
c oct3       ath  = th(i,icb(i)-1) - dttrig
c oct3 
c oct3       if (tdif.lt.dtcrit .or. ath.gt.ath1) then
c oct3         do 60 k=1,nl
c oct3            sig(i,k) = beta*sig(i,k) - 2.*alpha*tdif*tdif
c oct3            sig(i,k) = AMAX1(sig(i,k),0.0)
c oct3            w0(i,k)  = beta*w0(i,k)
c oct3   60    continue
c oct3         iflag(i)=4 ! pour version vectorisee
c oct3c convect3         iflag(i)=0
c oct3cccc         return
c oct3       endif
c oct3
c oct3200   continue
 
c -- oct3: on reecrit la boucle 200 (pour la vectorisation)

      do  60 k=1,nl
      do 200 i=1,len

       tdif = buoybase(i)
       ath1 = thnk(i)
       ath  = th(i,icb(i)-1) - dttrig

       if (tdif.lt.dtcrit .or. ath.gt.ath1) then
            sig(i,k) = beta*sig(i,k) - 2.*alpha*tdif*tdif
            sig(i,k) = AMAX1(sig(i,k),0.0)
            w0(i,k)  = beta*w0(i,k)
        iflag(i)=4 ! pour version vectorisee
c convect3         iflag(i)=0
       endif

200   continue
 60   continue

c fin oct3 --

      return
      end

      SUBROUTINE cv3_compress( len,nloc,ncum,nd,ntra
     :    ,iflag1,nk1,icb1,icbs1
     :    ,plcl1,tnk1,qnk1,gznk1,pbase1,buoybase1
     :    ,t1,q1,qs1,u1,v1,gz1,th1
     :    ,tra1
     :    ,h1,lv1,cpn1,p1,ph1,tv1,tp1,tvp1,clw1
     :    ,sig1,w01
     o    ,iflag,nk,icb,icbs
     o    ,plcl,tnk,qnk,gznk,pbase,buoybase
     o    ,t,q,qs,u,v,gz,th
     o    ,tra
     o    ,h,lv,cpn,p,ph,tv,tp,tvp,clw 
     o    ,sig,w0  )
      implicit none

#include "cv3param.h"
      include 'iniprint.h'

c inputs:
      integer len,ncum,nd,ntra,nloc
      integer iflag1(len),nk1(len),icb1(len),icbs1(len)
      real plcl1(len),tnk1(len),qnk1(len),gznk1(len)
      real pbase1(len),buoybase1(len)
      real t1(len,nd),q1(len,nd),qs1(len,nd),u1(len,nd),v1(len,nd)
      real gz1(len,nd),h1(len,nd),lv1(len,nd),cpn1(len,nd)
      real p1(len,nd),ph1(len,nd+1),tv1(len,nd),tp1(len,nd)
      real tvp1(len,nd),clw1(len,nd)
      real th1(len,nd)
      real sig1(len,nd), w01(len,nd)
      real tra1(len,nd,ntra)

c outputs:
c en fait, on a nloc=len pour l'instant (cf cv_driver)
      integer iflag(nloc),nk(nloc),icb(nloc),icbs(nloc)
      real plcl(nloc),tnk(nloc),qnk(nloc),gznk(nloc)
      real pbase(nloc),buoybase(nloc)
      real t(nloc,nd),q(nloc,nd),qs(nloc,nd),u(nloc,nd),v(nloc,nd)
      real gz(nloc,nd),h(nloc,nd),lv(nloc,nd),cpn(nloc,nd)
      real p(nloc,nd),ph(nloc,nd+1),tv(nloc,nd),tp(nloc,nd)
      real tvp(nloc,nd),clw(nloc,nd)
      real th(nloc,nd)
      real sig(nloc,nd), w0(nloc,nd)
      real tra(nloc,nd,ntra)

c local variables:
      integer i,k,nn,j

      CHARACTER (LEN=20) :: modname='cv3_compress'
      CHARACTER (LEN=80) :: abort_message

      do 110 k=1,nl+1
       nn=0
      do 100 i=1,len
      if(iflag1(i).eq.0)then
        nn=nn+1
        sig(nn,k)=sig1(i,k)
        w0(nn,k)=w01(i,k)
        t(nn,k)=t1(i,k)
        q(nn,k)=q1(i,k)
        qs(nn,k)=qs1(i,k)
        u(nn,k)=u1(i,k)
        v(nn,k)=v1(i,k)
        gz(nn,k)=gz1(i,k)
        h(nn,k)=h1(i,k)
        lv(nn,k)=lv1(i,k)
        cpn(nn,k)=cpn1(i,k)
        p(nn,k)=p1(i,k)
        ph(nn,k)=ph1(i,k)
        tv(nn,k)=tv1(i,k)
        tp(nn,k)=tp1(i,k)
        tvp(nn,k)=tvp1(i,k)
        clw(nn,k)=clw1(i,k)
        th(nn,k)=th1(i,k)
      endif
 100    continue
 110  continue

      do 121 j=1,ntra
ccccc      do 111 k=1,nl+1
      do 111 k=1,nd
       nn=0
      do 101 i=1,len
      if(iflag1(i).eq.0)then
       nn=nn+1
       tra(nn,k,j)=tra1(i,k,j)
      endif
 101  continue
 111  continue
 121  continue

      if (nn.ne.ncum) then
         write(lunout,*)'strange! nn not equal to ncum: ',nn,ncum
         abort_message = ''
         CALL abort_gcm (modname,abort_message,1)
      endif

      nn=0
      do 150 i=1,len
      if(iflag1(i).eq.0)then
      nn=nn+1
      pbase(nn)=pbase1(i)
      buoybase(nn)=buoybase1(i)
      plcl(nn)=plcl1(i)
      tnk(nn)=tnk1(i)
      qnk(nn)=qnk1(i)
      gznk(nn)=gznk1(i)
      nk(nn)=nk1(i)
      icb(nn)=icb1(i)
      icbs(nn)=icbs1(i)
      iflag(nn)=iflag1(i)
      endif
 150  continue

      return
      end

      SUBROUTINE cv3_undilute2(nloc,ncum,nd,icb,icbs,nk
     :                       ,tnk,qnk,gznk,hnk,t,q,qs,gz
     :                       ,p,h,tv,lv,pbase,buoybase,plcl
     o                       ,inb,tp,tvp,clw,hp,ep,sigp,buoy)
      implicit none

C---------------------------------------------------------------------
C Purpose:
C     FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
C     &
C     COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
C     FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
C     &
C     FIND THE LEVEL OF NEUTRAL BUOYANCY
C
C Main differences convect3/convect4:
C	- icbs (input) is the first level above LCL (may differ from icb)
C	- many minor differences in the iterations
C	- condensed water not removed from tvp in convect3
C   - vertical profile of buoyancy computed here (use of buoybase)
C   - the determination of inb is different
C   - no inb1, only inb in output
C---------------------------------------------------------------------

#include "cvthermo.h"
#include "cv3param.h"
#include "conema3.h"

c inputs:
      integer ncum, nd, nloc
      integer icb(nloc), icbs(nloc), nk(nloc)
      real t(nloc,nd), q(nloc,nd), qs(nloc,nd), gz(nloc,nd)
      real p(nloc,nd)
      real tnk(nloc), qnk(nloc), gznk(nloc)
      real hnk(nloc)
      real lv(nloc,nd), tv(nloc,nd), h(nloc,nd)
      real pbase(nloc), buoybase(nloc), plcl(nloc)

c outputs:
      integer inb(nloc)
      real tp(nloc,nd), tvp(nloc,nd), clw(nloc,nd)
      real ep(nloc,nd), sigp(nloc,nd), hp(nloc,nd)
      real buoy(nloc,nd)

c local variables:
      integer i, k
      real tg,qg,ahg,alv,s,tc,es,denom,rg,tca,elacrit
      real by, defrac, pden
      real ah0(nloc), cape(nloc), capem(nloc), byp(nloc)
      logical lcape(nloc)
      integer iposit(nloc)

!=====================================================================
! --- SOME INITIALIZATIONS
!=====================================================================

      do 170 k=1,nl
      do 160 i=1,ncum
       ep(i,k)=0.0
       sigp(i,k)=spfac
 160  continue
 170  continue

!=====================================================================
! --- FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
!=====================================================================
c
c ---       The procedure is to solve the equation.
c              cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.
c
c   ***  Calculate certain parcel quantities, including static energy   ***
c
c
      do 240 i=1,ncum
         ah0(i)=(cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i)
cdebug     &         +qnk(i)*(lv0-clmcpv*(tnk(i)-t0))+gznk(i)
     &         +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15))+gznk(i)
 240  continue
c
c
c    ***  Find lifted parcel quantities above cloud base    ***
c
c
	do 300 k=minorig+1,nl
	  do 290 i=1,ncum
c ori	    if(k.ge.(icb(i)+1))then
	    if(k.ge.(icbs(i)+1))then ! convect3
	      tg=t(i,k)
	      qg=qs(i,k)
cdebug	      alv=lv0-clmcpv*(t(i,k)-t0)
	      alv=lv0-clmcpv*(t(i,k)-273.15)
c
c First iteration.
c
c ori	       s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
           s=cpd*(1.-qnk(i))+cl*qnk(i)      ! convect3
     :      +alv*alv*qg/(rrv*t(i,k)*t(i,k)) ! convect3
	       s=1./s
c ori	       ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
           ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gz(i,k) ! convect3
	       tg=tg+s*(ah0(i)-ahg)
c ori	       tg=max(tg,35.0)
cdebug	       tc=tg-t0
	       tc=tg-273.15
	       denom=243.5+tc
           denom=MAX(denom,1.0) ! convect3
c ori	       if(tc.ge.0.0)then
			es=6.112*exp(17.67*tc/denom)
c ori	       else
c ori			es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
c ori	       endif
			qg=eps*es/(p(i,k)-es*(1.-eps))
c
c Second iteration.
c
c ori	       s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
c ori	       s=1./s
c ori	       ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
           ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gz(i,k) ! convect3
	       tg=tg+s*(ah0(i)-ahg)
c ori	       tg=max(tg,35.0)
cdebug	       tc=tg-t0
	       tc=tg-273.15
	       denom=243.5+tc
           denom=MAX(denom,1.0) ! convect3
c ori	       if(tc.ge.0.0)then
			es=6.112*exp(17.67*tc/denom)
c ori	       else
c ori			es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
c ori	       endif
			qg=eps*es/(p(i,k)-es*(1.-eps))
c
cdebug	       alv=lv0-clmcpv*(t(i,k)-t0)
	       alv=lv0-clmcpv*(t(i,k)-273.15)
c      print*,'cpd dans convect2 ',cpd
c      print*,'tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd'
c      print*,tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd

c ori c approximation here:
c ori        tp(i,k)=(ah0(i)-(cl-cpd)*qnk(i)*t(i,k)-gz(i,k)-alv*qg)/cpd

c convect3: no approximation:
           tp(i,k)=(ah0(i)-gz(i,k)-alv*qg)/(cpd+(cl-cpd)*qnk(i))

               clw(i,k)=qnk(i)-qg
               clw(i,k)=max(0.0,clw(i,k))
               rg=qg/(1.-qnk(i))
c ori               tvp(i,k)=tp(i,k)*(1.+rg*epsi)
c convect3: (qg utilise au lieu du vrai mixing ratio rg):
               tvp(i,k)=tp(i,k)*(1.+qg/eps-qnk(i)) ! whole thing
            endif
  290     continue
  300   continue
c
!=====================================================================
! --- SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF
! --- PRECIPITATION FALLING OUTSIDE OF CLOUD
! --- THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)
!=====================================================================
c
c ori      do 320 k=minorig+1,nl
      do 320 k=1,nl ! convect3
        do 310 i=1,ncum
           pden=ptcrit-pbcrit
           ep(i,k)=(plcl(i)-p(i,k)-pbcrit)/pden*epmax
           ep(i,k)=amax1(ep(i,k),0.0)
           ep(i,k)=amin1(ep(i,k),epmax)
           sigp(i,k)=spfac
c ori          if(k.ge.(nk(i)+1))then
c ori            tca=tp(i,k)-t0
c ori            if(tca.ge.0.0)then
c ori              elacrit=elcrit
c ori            else
c ori              elacrit=elcrit*(1.0-tca/tlcrit)
c ori            endif
c ori            elacrit=max(elacrit,0.0)
c ori            ep(i,k)=1.0-elacrit/max(clw(i,k),1.0e-8)
c ori            ep(i,k)=max(ep(i,k),0.0 )
c ori            ep(i,k)=min(ep(i,k),1.0 )
c ori            sigp(i,k)=sigs
c ori          endif
 310    continue
 320  continue
c
!=====================================================================
! --- CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL
! --- VIRTUAL TEMPERATURE
!=====================================================================
c
c dans convect3, tvp est calcule en une seule fois, et sans retirer
c l'eau condensee (~> reversible CAPE)
c
c ori      do 340 k=minorig+1,nl
c ori        do 330 i=1,ncum
c ori        if(k.ge.(icb(i)+1))then
c ori          tvp(i,k)=tvp(i,k)*(1.0-qnk(i)+ep(i,k)*clw(i,k))
c oric         print*,'i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)'
c oric         print*, i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)
c ori        endif
c ori 330    continue
c ori 340  continue

c ori      do 350 i=1,ncum
c ori       tvp(i,nlp)=tvp(i,nl)-(gz(i,nlp)-gz(i,nl))/cpd
c ori 350  continue

      do 350 i=1,ncum       ! convect3
       tp(i,nlp)=tp(i,nl)   ! convect3
 350  continue              ! convect3
c
c=====================================================================
c  --- EFFECTIVE VERTICAL PROFILE OF BUOYANCY (convect3 only):
c=====================================================================

c-- this is for convect3 only:

c first estimate of buoyancy:

      do 500 i=1,ncum
       do 501 k=1,nl
        buoy(i,k)=tvp(i,k)-tv(i,k)
 501   continue
 500  continue

c set buoyancy=buoybase for all levels below base
c for safety, set buoy(icb)=buoybase

      do 505 i=1,ncum
       do 506 k=1,nl
        if((k.ge.icb(i)).and.(k.le.nl).and.(p(i,k).ge.pbase(i)))then
         buoy(i,k)=buoybase(i)
        endif
 506   continue
c       buoy(icb(i),k)=buoybase(i)
      buoy(i,icb(i))=buoybase(i)
 505  continue

c-- end convect3

c=====================================================================
c  --- FIND THE FIRST MODEL LEVEL (INB) ABOVE THE PARCEL'S
c  --- LEVEL OF NEUTRAL BUOYANCY
c=====================================================================
c
c-- this is for convect3 only:

      do 510 i=1,ncum
       inb(i)=nl-1
       iposit(i) = nl
 510  continue

c
c--    iposit(i) = first level, above icb, with positive buoyancy
      do k = 1,nl-1
       do i = 1,ncum
        if (k .ge. icb(i) .and. buoy(i,k) .gt. 0.) then
          iposit(i) = min(iposit(i),k)
        endif
       enddo
      enddo

      do i = 1,ncum
       if (iposit(i) .eq. nl) then
         iposit(i) = icb(i)
       endif
      enddo

      do 530 i=1,ncum
       do 535 k=1,nl-1
        if ((k.ge.iposit(i)).and.(buoy(i,k).lt.dtovsh)) then
         inb(i)=MIN(inb(i),k)
        endif
 535   continue
 530  continue

c-- end convect3

c ori      do 510 i=1,ncum
c ori        cape(i)=0.0
c ori        capem(i)=0.0
c ori        inb(i)=icb(i)+1
c ori        inb1(i)=inb(i)
c ori 510  continue
c
c Originial Code
c
c     do 530 k=minorig+1,nl-1
c       do 520 i=1,ncum
c         if(k.ge.(icb(i)+1))then
c           by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
c           byp=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
c           cape(i)=cape(i)+by
c           if(by.ge.0.0)inb1(i)=k+1
c           if(cape(i).gt.0.0)then
c             inb(i)=k+1
c             capem(i)=cape(i)
c           endif
c         endif
c520    continue
c530  continue
c     do 540 i=1,ncum
c         byp=(tvp(i,nl)-tv(i,nl))*dph(i,nl)/p(i,nl)
c         cape(i)=capem(i)+byp
c         defrac=capem(i)-cape(i)
c         defrac=max(defrac,0.001)
c         frac(i)=-cape(i)/defrac
c         frac(i)=min(frac(i),1.0)
c         frac(i)=max(frac(i),0.0)
c540   continue
c
c K Emanuel fix
c
c     call zilch(byp,ncum)
c     do 530 k=minorig+1,nl-1
c       do 520 i=1,ncum
c         if(k.ge.(icb(i)+1))then
c           by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
c           cape(i)=cape(i)+by
c           if(by.ge.0.0)inb1(i)=k+1
c           if(cape(i).gt.0.0)then
c             inb(i)=k+1
c             capem(i)=cape(i)
c             byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
c           endif
c         endif
c520    continue
c530  continue
c     do 540 i=1,ncum
c         inb(i)=max(inb(i),inb1(i))
c         cape(i)=capem(i)+byp(i)
c         defrac=capem(i)-cape(i)
c         defrac=max(defrac,0.001)
c         frac(i)=-cape(i)/defrac
c         frac(i)=min(frac(i),1.0)
c         frac(i)=max(frac(i),0.0)
c540   continue
c
c J Teixeira fix
c
c ori      call zilch(byp,ncum)
c ori      do 515 i=1,ncum
c ori        lcape(i)=.true.
c ori 515  continue
c ori      do 530 k=minorig+1,nl-1
c ori        do 520 i=1,ncum
c ori          if(cape(i).lt.0.0)lcape(i)=.false.
c ori          if((k.ge.(icb(i)+1)).and.lcape(i))then
c ori            by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
c ori            byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
c ori            cape(i)=cape(i)+by
c ori            if(by.ge.0.0)inb1(i)=k+1
c ori            if(cape(i).gt.0.0)then
c ori              inb(i)=k+1
c ori              capem(i)=cape(i)
c ori            endif
c ori          endif
c ori 520    continue
c ori 530  continue
c ori      do 540 i=1,ncum
c ori          cape(i)=capem(i)+byp(i)
c ori          defrac=capem(i)-cape(i)
c ori          defrac=max(defrac,0.001)
c ori          frac(i)=-cape(i)/defrac
c ori          frac(i)=min(frac(i),1.0)
c ori          frac(i)=max(frac(i),0.0)
c ori 540  continue
c
c=====================================================================
c ---   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL
c=====================================================================
c
      do k = 1,nd
      do i=1,ncum
         hp(i,k)=h(i,k)
      enddo
      enddo

      do 600 k=minorig+1,nl
        do 590 i=1,ncum
        if((k.ge.icb(i)).and.(k.le.inb(i)))then
          hp(i,k)=hnk(i)+(lv(i,k)+(cpd-cpv)*t(i,k))*ep(i,k)*clw(i,k)
        endif
 590    continue
 600  continue

        return
        end

      SUBROUTINE cv3_closure(nloc,ncum,nd,icb,inb
     :                      ,pbase,p,ph,tv,buoy
     o                      ,sig,w0,cape,m,iflag)
      implicit none

!===================================================================
! ---  CLOSURE OF CONVECT3
!
! vectorization: S. Bony
!===================================================================

#include "cvthermo.h"
#include "cv3param.h"

c input:
      integer ncum, nd, nloc
      integer icb(nloc), inb(nloc)
      real pbase(nloc)
      real p(nloc,nd), ph(nloc,nd+1)
      real tv(nloc,nd), buoy(nloc,nd)

c input/output:
      real sig(nloc,nd), w0(nloc,nd)
      integer iflag(nloc)

c output:
      real cape(nloc)
      real m(nloc,nd)

c local variables:
      integer i, j, k, icbmax
      real deltap, fac, w, amu
      real dtmin(nloc,nd), sigold(nloc,nd)
      real cbmflast(nloc)


c -------------------------------------------------------
c -- Initialization
c -------------------------------------------------------

      do k=1,nl
       do i=1,ncum
        m(i,k)=0.0
       enddo
      enddo

c -------------------------------------------------------
c -- Reset sig(i) and w0(i) for i>inb and i<icb
c -------------------------------------------------------

c update sig and w0 above LNB:

      do 100 k=1,nl-1
       do 110 i=1,ncum
        if ((inb(i).lt.(nl-1)).and.(k.ge.(inb(i)+1)))then
         sig(i,k)=beta*sig(i,k)
     :            +2.*alpha*buoy(i,inb(i))*ABS(buoy(i,inb(i)))
         sig(i,k)=AMAX1(sig(i,k),0.0)
         w0(i,k)=beta*w0(i,k)
        endif
 110   continue
 100  continue

c compute icbmax:

      icbmax=2
      do 200 i=1,ncum
        icbmax=MAX(icbmax,icb(i))
 200  continue

c update sig and w0 below cloud base:

      do 300 k=1,icbmax
       do 310 i=1,ncum
        if (k.le.icb(i))then
         sig(i,k)=beta*sig(i,k)-2.*alpha*buoy(i,icb(i))*buoy(i,icb(i))
         sig(i,k)=amax1(sig(i,k),0.0)
         w0(i,k)=beta*w0(i,k)
        endif
310    continue
300    continue

c!      if(inb.lt.(nl-1))then
c!         do 85 i=inb+1,nl-1
c!            sig(i)=beta*sig(i)+2.*alpha*buoy(inb)*
c!     1              abs(buoy(inb))
c!            sig(i)=amax1(sig(i),0.0)
c!            w0(i)=beta*w0(i)
c!   85    continue
c!      end if

c!      do 87 i=1,icb
c!         sig(i)=beta*sig(i)-2.*alpha*buoy(icb)*buoy(icb)
c!         sig(i)=amax1(sig(i),0.0)
c!         w0(i)=beta*w0(i)
c!   87 continue

c -------------------------------------------------------------
c -- Reset fractional areas of updrafts and w0 at initial time
c -- and after 10 time steps of no convection
c -------------------------------------------------------------

      do 400 k=1,nl-1
       do 410 i=1,ncum
        if (sig(i,nd).lt.1.5.or.sig(i,nd).gt.12.0)then
         sig(i,k)=0.0
         w0(i,k)=0.0
        endif
 410   continue
 400  continue

c -------------------------------------------------------------
c -- Calculate convective available potential energy (cape),
c -- vertical velocity (w), fractional area covered by
c -- undilute updraft (sig), and updraft mass flux (m)
c -------------------------------------------------------------

      do 500 i=1,ncum
       cape(i)=0.0
 500  continue

c compute dtmin (minimum buoyancy between ICB and given level k):

      do i=1,ncum
       do k=1,nl
         dtmin(i,k)=100.0
       enddo
      enddo

      do 550 i=1,ncum
       do 560 k=1,nl
         do 570 j=minorig,nl
          if ( (k.ge.(icb(i)+1)).and.(k.le.inb(i)).and.
     :         (j.ge.icb(i)).and.(j.le.(k-1)) )then
           dtmin(i,k)=AMIN1(dtmin(i,k),buoy(i,j))
          endif
 570     continue
 560   continue
 550  continue

c the interval on which cape is computed starts at pbase :

      do 600 k=1,nl
       do 610 i=1,ncum

        if ((k.ge.(icb(i)+1)).and.(k.le.inb(i))) then

         deltap = MIN(pbase(i),ph(i,k-1))-MIN(pbase(i),ph(i,k))
         cape(i)=cape(i)+rrd*buoy(i,k-1)*deltap/p(i,k-1)
         cape(i)=AMAX1(0.0,cape(i))
         sigold(i,k)=sig(i,k)

c         dtmin(i,k)=100.0
c         do 97 j=icb(i),k-1 ! mauvaise vectorisation
c          dtmin(i,k)=AMIN1(dtmin(i,k),buoy(i,j))
c  97     continue

         sig(i,k)=beta*sig(i,k)+alpha*dtmin(i,k)*ABS(dtmin(i,k))
         sig(i,k)=amax1(sig(i,k),0.0)
         sig(i,k)=amin1(sig(i,k),0.01)
         fac=AMIN1(((dtcrit-dtmin(i,k))/dtcrit),1.0)
         w=(1.-beta)*fac*SQRT(cape(i))+beta*w0(i,k)
         amu=0.5*(sig(i,k)+sigold(i,k))*w
         m(i,k)=amu*0.007*p(i,k)*(ph(i,k)-ph(i,k+1))/tv(i,k)
         w0(i,k)=w
        endif

 610   continue
 600  continue

      do 700 i=1,ncum
       w0(i,icb(i))=0.5*w0(i,icb(i)+1)
       m(i,icb(i))=0.5*m(i,icb(i)+1)
     :             *(ph(i,icb(i))-ph(i,icb(i)+1))
     :             /(ph(i,icb(i)+1)-ph(i,icb(i)+2))
       sig(i,icb(i))=sig(i,icb(i)+1)
       sig(i,icb(i)-1)=sig(i,icb(i))
 700  continue
c
cccc 3. Compute final cloud base mass flux and set iflag to 3 if
cccc    cloud base mass flux is exceedingly small and is decreasing (i.e. if
cccc    the final mass flux (cbmflast) is greater than the target mass flux
cccc    (cbmf) ??).
ccc
cc      do i = 1,ncum
cc       cbmflast(i) = 0.
cc      enddo
ccc
cc      do k= 1,nl
cc       do i = 1,ncum
cc        IF (k .ge. icb(i) .and. k .le. inb(i)) THEN
cc         cbmflast(i) = cbmflast(i)+M(i,k)
cc        ENDIF
cc       enddo
cc      enddo
ccc
cc      do i = 1,ncum
cc       IF (cbmflast(i) .lt. 1.e-6) THEN
cc         iflag(i) = 3
cc       ENDIF
cc      enddo
ccc
cc      do k= 1,nl
cc       do i = 1,ncum
cc        IF (iflag(i) .ge. 3) THEN
cc         M(i,k) = 0.
cc         sig(i,k) = 0.
cc         w0(i,k) = 0.
cc        ENDIF
cc       enddo
cc      enddo
ccc
c!      cape=0.0
c!      do 98 i=icb+1,inb
c!         deltap = min(pbase,ph(i-1))-min(pbase,ph(i))
c!         cape=cape+rrd*buoy(i-1)*deltap/p(i-1)
c!         dcape=rrd*buoy(i-1)*deltap/p(i-1)
c!         dlnp=deltap/p(i-1)
c!         cape=amax1(0.0,cape)
c!         sigold=sig(i)

c!         dtmin=100.0
c!         do 97 j=icb,i-1
c!            dtmin=amin1(dtmin,buoy(j))
c!   97    continue

c!         sig(i)=beta*sig(i)+alpha*dtmin*abs(dtmin)
c!         sig(i)=amax1(sig(i),0.0)
c!         sig(i)=amin1(sig(i),0.01)
c!         fac=amin1(((dtcrit-dtmin)/dtcrit),1.0)
c!         w=(1.-beta)*fac*sqrt(cape)+beta*w0(i)
c!         amu=0.5*(sig(i)+sigold)*w
c!         m(i)=amu*0.007*p(i)*(ph(i)-ph(i+1))/tv(i)
c!         w0(i)=w
c!   98 continue
c!      w0(icb)=0.5*w0(icb+1)
c!      m(icb)=0.5*m(icb+1)*(ph(icb)-ph(icb+1))/(ph(icb+1)-ph(icb+2))
c!      sig(icb)=sig(icb+1)
c!      sig(icb-1)=sig(icb)

       return
       end

      SUBROUTINE cv3_mixing(nloc,ncum,nd,na,ntra,icb,nk,inb
     :                    ,ph,t,rr,rs,u,v,tra,h,lv,qnk
     :                    ,unk,vnk,hp,tv,tvp,ep,clw,m,sig
     :   ,ment,qent,uent,vent,nent,sij,elij,ments,qents,traent)
      implicit none

!---------------------------------------------------------------------
! a faire:
!   - vectorisation de la partie normalisation des flux (do 789...)
!---------------------------------------------------------------------

#include "cvthermo.h"
#include "cv3param.h"

c inputs:
      integer ncum, nd, na, ntra, nloc
      integer icb(nloc), inb(nloc), nk(nloc)
      real sig(nloc,nd)
      real qnk(nloc),unk(nloc),vnk(nloc)
      real ph(nloc,nd+1)
      real t(nloc,nd), rr(nloc,nd), rs(nloc,nd)
      real u(nloc,nd), v(nloc,nd)
      real tra(nloc,nd,ntra) ! input of convect3
      real lv(nloc,na), h(nloc,na), hp(nloc,na)
      real tv(nloc,na), tvp(nloc,na), ep(nloc,na), clw(nloc,na)
      real m(nloc,na)        ! input of convect3

c outputs:
      real ment(nloc,na,na), qent(nloc,na,na)
      real uent(nloc,na,na), vent(nloc,na,na)
      real sij(nloc,na,na), elij(nloc,na,na)
      real traent(nloc,nd,nd,ntra)
      real ments(nloc,nd,nd), qents(nloc,nd,nd)
      real sigij(nloc,nd,nd)
      integer nent(nloc,nd)

c local variables:
      integer i, j, k, il, im, jm
      integer num1, num2
      real rti, bf2, anum, denom, dei, altem, cwat, stemp, qp
      real alt, smid, sjmin, sjmax, delp, delm
      real asij(nloc), smax(nloc), scrit(nloc)
      real asum(nloc,nd),bsum(nloc,nd),csum(nloc,nd)
      real wgh
      real zm(nloc,na)
      logical lwork(nloc)

c=====================================================================
c --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
c=====================================================================

c ori        do 360 i=1,ncum*nlp
        do 361 j=1,nl
        do 360 i=1,ncum
          nent(i,j)=0
c in convect3, m is computed in cv3_closure
c ori          m(i,1)=0.0
 360    continue
 361    continue

c ori      do 400 k=1,nlp
c ori       do 390 j=1,nlp
      do 400 j=1,nl
       do 390 k=1,nl
          do 385 i=1,ncum
            qent(i,k,j)=rr(i,j)
            uent(i,k,j)=u(i,j)
            vent(i,k,j)=v(i,j)
            elij(i,k,j)=0.0
cym            ment(i,k,j)=0.0
cym            sij(i,k,j)=0.0
 385      continue
 390    continue
 400  continue

cym
      ment(1:ncum,1:nd,1:nd)=0.0
      sij(1:ncum,1:nd,1:nd)=0.0
      
      do k=1,ntra
       do j=1,nd  ! instead nlp
        do i=1,nd ! instead nlp
         do il=1,ncum
            traent(il,i,j,k)=tra(il,j,k)
         enddo
        enddo
       enddo
      enddo
      zm(:,:)=0.

c=====================================================================
c --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
c --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
c --- FRACTION (sij)
c=====================================================================

      do 750 i=minorig+1, nl

       do 710 j=minorig,nl
        do 700 il=1,ncum
         if( (i.ge.icb(il)).and.(i.le.inb(il)).and.
     :      (j.ge.(icb(il)-1)).and.(j.le.inb(il)))then

          rti=qnk(il)-ep(il,i)*clw(il,i)
          bf2=1.+lv(il,j)*lv(il,j)*rs(il,j)/(rrv*t(il,j)*t(il,j)*cpd)
          anum=h(il,j)-hp(il,i)+(cpv-cpd)*t(il,j)*(rti-rr(il,j))
          denom=h(il,i)-hp(il,i)+(cpd-cpv)*(rr(il,i)-rti)*t(il,j)
          dei=denom
          if(abs(dei).lt.0.01)dei=0.01
          sij(il,i,j)=anum/dei
          sij(il,i,i)=1.0
          altem=sij(il,i,j)*rr(il,i)+(1.-sij(il,i,j))*rti-rs(il,j)
          altem=altem/bf2
          cwat=clw(il,j)*(1.-ep(il,j))
          stemp=sij(il,i,j)
          if((stemp.lt.0.0.or.stemp.gt.1.0.or.altem.gt.cwat)
     :                 .and.j.gt.i)then
           anum=anum-lv(il,j)*(rti-rs(il,j)-cwat*bf2)
           denom=denom+lv(il,j)*(rr(il,i)-rti)
           if(abs(denom).lt.0.01)denom=0.01
           sij(il,i,j)=anum/denom
           altem=sij(il,i,j)*rr(il,i)+(1.-sij(il,i,j))*rti-rs(il,j)
           altem=altem-(bf2-1.)*cwat
          end if
         if(sij(il,i,j).gt.0.0.and.sij(il,i,j).lt.0.95)then
          qent(il,i,j)=sij(il,i,j)*rr(il,i)+(1.-sij(il,i,j))*rti
          uent(il,i,j)=sij(il,i,j)*u(il,i)+(1.-sij(il,i,j))*unk(il)
          vent(il,i,j)=sij(il,i,j)*v(il,i)+(1.-sij(il,i,j))*vnk(il)
c!!!      do k=1,ntra
c!!!      traent(il,i,j,k)=sij(il,i,j)*tra(il,i,k)
c!!!     :      +(1.-sij(il,i,j))*tra(il,nk(il),k)
c!!!      end do
          elij(il,i,j)=altem
          elij(il,i,j)=amax1(0.0,elij(il,i,j))
          ment(il,i,j)=m(il,i)/(1.-sij(il,i,j))
          nent(il,i)=nent(il,i)+1
         end if
         sij(il,i,j)=amax1(0.0,sij(il,i,j))
         sij(il,i,j)=amin1(1.0,sij(il,i,j))
         endif ! new
 700   continue
 710  continue

       do k=1,ntra
        do j=minorig,nl
         do il=1,ncum
          if( (i.ge.icb(il)).and.(i.le.inb(il)).and.
     :       (j.ge.(icb(il)-1)).and.(j.le.inb(il)))then
            traent(il,i,j,k)=sij(il,i,j)*tra(il,i,k)
     :            +(1.-sij(il,i,j))*tra(il,nk(il),k)
          endif
         enddo
        enddo
       enddo

c
c   ***   if no air can entrain at level i assume that updraft detrains  ***
c   ***   at that level and calculate detrained air flux and properties  ***
c

c@      do 170 i=icb(il),inb(il)

      do 740 il=1,ncum
      if ((i.ge.icb(il)).and.(i.le.inb(il)).and.(nent(il,i).eq.0)) then
c@      if(nent(il,i).eq.0)then
      ment(il,i,i)=m(il,i)
      qent(il,i,i)=qnk(il)-ep(il,i)*clw(il,i)
      uent(il,i,i)=unk(il)
      vent(il,i,i)=vnk(il)
      elij(il,i,i)=clw(il,i)
cMAF      sij(il,i,i)=1.0
      sij(il,i,i)=0.0
      end if
 740  continue
 750  continue

      do j=1,ntra
       do i=minorig+1,nl
        do il=1,ncum
         if (i.ge.icb(il) .and. i.le.inb(il) .and. nent(il,i).eq.0) then
          traent(il,i,i,j)=tra(il,nk(il),j)
         endif
        enddo
       enddo
      enddo

      do 100 j=minorig,nl
      do 101 i=minorig,nl
      do 102 il=1,ncum
      if ((j.ge.(icb(il)-1)).and.(j.le.inb(il))
     :    .and.(i.ge.icb(il)).and.(i.le.inb(il)))then
       sigij(il,i,j)=sij(il,i,j)
      endif
 102  continue
 101  continue
 100  continue
c@      enddo

c@170   continue

c=====================================================================
c   ---  NORMALIZE ENTRAINED AIR MASS FLUXES
c   ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
c=====================================================================

      call zilch(asum,nloc*nd)
      call zilch(csum,nloc*nd)
      call zilch(csum,nloc*nd)

      do il=1,ncum
       lwork(il) = .FALSE.
      enddo

      DO 789 i=minorig+1,nl

      num1=0
      do il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) ) num1=num1+1
      enddo
      if (num1.le.0) goto 789


      do 781 il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) ) then
        lwork(il)=(nent(il,i).ne.0)
        qp=qnk(il)-ep(il,i)*clw(il,i)
        anum=h(il,i)-hp(il,i)-lv(il,i)*(qp-rs(il,i))
     :           +(cpv-cpd)*t(il,i)*(qp-rr(il,i))
        denom=h(il,i)-hp(il,i)+lv(il,i)*(rr(il,i)-qp)
     :           +(cpd-cpv)*t(il,i)*(rr(il,i)-qp)
        if(abs(denom).lt.0.01)denom=0.01
        scrit(il)=anum/denom
        alt=qp-rs(il,i)+scrit(il)*(rr(il,i)-qp)
        if(scrit(il).le.0.0.or.alt.le.0.0)scrit(il)=1.0
        smax(il)=0.0
        asij(il)=0.0
       endif
781   continue

      do 175 j=nl,minorig,-1

      num2=0
      do il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) .and.
     :      j.ge.(icb(il)-1) .and. j.le.inb(il)
     :      .and. lwork(il) ) num2=num2+1
      enddo
      if (num2.le.0) goto 175

      do 782 il=1,ncum
      if ( i.ge.icb(il) .and. i.le.inb(il) .and.
     :      j.ge.(icb(il)-1) .and. j.le.inb(il)
     :      .and. lwork(il) ) then

       if(sij(il,i,j).gt.1.0e-16.and.sij(il,i,j).lt.0.95)then
        wgh=1.0
        if(j.gt.i)then
         sjmax=amax1(sij(il,i,j+1),smax(il))
         sjmax=amin1(sjmax,scrit(il))
         smax(il)=amax1(sij(il,i,j),smax(il))
         sjmin=amax1(sij(il,i,j-1),smax(il))
         sjmin=amin1(sjmin,scrit(il))
         if(sij(il,i,j).lt.(smax(il)-1.0e-16))wgh=0.0
         smid=amin1(sij(il,i,j),scrit(il))
        else
         sjmax=amax1(sij(il,i,j+1),scrit(il))
         smid=amax1(sij(il,i,j),scrit(il))
         sjmin=0.0
         if(j.gt.1)sjmin=sij(il,i,j-1)
         sjmin=amax1(sjmin,scrit(il))
        endif
        delp=abs(sjmax-smid)
        delm=abs(sjmin-smid)
        asij(il)=asij(il)+wgh*(delp+delm)
        ment(il,i,j)=ment(il,i,j)*(delp+delm)*wgh
       endif
      endif
782   continue

175   continue

      do il=1,ncum
       if (i.ge.icb(il).and.i.le.inb(il).and.lwork(il)) then
        asij(il)=amax1(1.0e-16,asij(il))
        asij(il)=1.0/asij(il)
        asum(il,i)=0.0
        bsum(il,i)=0.0
        csum(il,i)=0.0
       endif
      enddo

      do 180 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :   .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         ment(il,i,j)=ment(il,i,j)*asij(il)
        endif
       enddo
180   continue

      do 190 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :   .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         asum(il,i)=asum(il,i)+ment(il,i,j)
         ment(il,i,j)=ment(il,i,j)*sig(il,j)
         bsum(il,i)=bsum(il,i)+ment(il,i,j)
        endif
       enddo
190   continue

      do il=1,ncum
       if (i.ge.icb(il).and.i.le.inb(il).and.lwork(il)) then
        bsum(il,i)=amax1(bsum(il,i),1.0e-16)
        bsum(il,i)=1.0/bsum(il,i)
       endif
      enddo

      do 195 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :   .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         ment(il,i,j)=ment(il,i,j)*asum(il,i)*bsum(il,i)
        endif
       enddo
195   continue

      do 197 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :   .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         csum(il,i)=csum(il,i)+ment(il,i,j)
        endif
       enddo
197   continue

      do il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :     .and. csum(il,i).lt.m(il,i) ) then
        nent(il,i)=0
        ment(il,i,i)=m(il,i)
        qent(il,i,i)=qnk(il)-ep(il,i)*clw(il,i)
        uent(il,i,i)=unk(il)
        vent(il,i,i)=vnk(il)
        elij(il,i,i)=clw(il,i)
cMAF        sij(il,i,i)=1.0
        sij(il,i,i)=0.0
       endif
      enddo ! il

      do j=1,ntra
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :     .and. csum(il,i).lt.m(il,i) ) then
         traent(il,i,i,j)=tra(il,nk(il),j)
        endif
       enddo
      enddo
789   continue
c      
c MAF: renormalisation de MENT
      call zilch(zm,nloc*na)
      do jm=1,nd
        do im=1,nd
          do il=1,ncum
          zm(il,im)=zm(il,im)+(1.-sij(il,im,jm))*ment(il,im,jm)
         end do
        end do
      end do
c
      do jm=1,nd
        do im=1,nd
          do il=1,ncum
          if(zm(il,im).ne.0.) then
          ment(il,im,jm)=ment(il,im,jm)*m(il,im)/zm(il,im)
          endif
         end do
       end do
      end do
c
      do jm=1,nd
       do im=1,nd
        do 999 il=1,ncum
         qents(il,im,jm)=qent(il,im,jm)
         ments(il,im,jm)=ment(il,im,jm)
999     continue
       enddo
      enddo

      return
      end

      SUBROUTINE cv3_unsat(nloc,ncum,nd,na,ntra,icb,inb,iflag
     :              ,t,rr,rs,gz,u,v,tra,p,ph
     :              ,th,tv,lv,cpn,ep,sigp,clw
     :              ,m,ment,elij,delt,plcl,coef_clos
     o              ,mp,rp,up,vp,trap,wt,water,evap,b,sigd)
      implicit none


#include "cvthermo.h"
#include "cv3param.h"
#include "cvflag.h"

c inputs:
      integer ncum, nd, na, ntra, nloc
      integer icb(nloc), inb(nloc)
      real delt, plcl(nloc)
      real t(nloc,nd), rr(nloc,nd), rs(nloc,nd),gz(nloc,na)
      real u(nloc,nd), v(nloc,nd)
      real tra(nloc,nd,ntra)
      real p(nloc,nd), ph(nloc,nd+1)
      real ep(nloc,na), sigp(nloc,na), clw(nloc,na)
      real th(nloc,na),tv(nloc,na),lv(nloc,na),cpn(nloc,na)
      real m(nloc,na), ment(nloc,na,na), elij(nloc,na,na)
      real coef_clos(nloc)
c
c input/output
      integer iflag(nloc)
c
c outputs:
      real mp(nloc,na), rp(nloc,na), up(nloc,na), vp(nloc,na)
      real water(nloc,na), evap(nloc,na), wt(nloc,na)
      real trap(nloc,na,ntra)
      real b(nloc,na), sigd(nloc)

c local variables
      integer i,j,k,il,num1,ndp1
      real tinv, delti
      real awat, afac, afac1, afac2, bfac
      real pr1, pr2, sigt, b6, c6, revap, tevap, delth
      real amfac, amp2, xf, tf, fac2, ur, sru, fac, d, af, bf
      real ampmax
      real lvcp(nloc,na)
      real h(nloc,na),hm(nloc,na)
      real wdtrain(nloc)
      logical lwork(nloc)


c------------------------------------------------------

        delti = 1./delt
        tinv=1./3.

        mp(:,:)=0.

        do i=1,nl
         do il=1,ncum
          mp(il,i)=0.0
          rp(il,i)=rr(il,i)
          up(il,i)=u(il,i)
          vp(il,i)=v(il,i)
          wt(il,i)=0.001
          water(il,i)=0.0
          evap(il,i)=0.0
          b(il,i)=0.0
          lvcp(il,i)=lv(il,i)/cpn(il,i)
         enddo
        enddo
        do k=1,ntra
         do i=1,nd
          do il=1,ncum
           trap(il,i,k)=tra(il,i,k)
          enddo
         enddo
        enddo
c
c   ***  check whether ep(inb)=0, if so, skip precipitating    ***
c   ***             downdraft calculation                      ***
c

        do il=1,ncum
          lwork(il)=.TRUE.
          if(ep(il,inb(il)).lt.0.0001)lwork(il)=.FALSE.
        enddo

        call zilch(wdtrain,ncum)

c   ***  Set the fractionnal area sigd of precipitating downdraughts
        do il = 1,ncum
          sigd(il) = sigdz*coef_clos(il)
        enddo

        DO 400 i=nl+1,1,-1

        num1=0
        do il=1,ncum
         if ( i.le.inb(il) .and. lwork(il) ) num1=num1+1
        enddo
        if (num1.le.0) goto 400

c
c   ***  integrate liquid water equation to find condensed water   ***
c   ***                and condensed water flux                    ***
c

c
c    ***                    begin downdraft loop                    ***
c

c
c    ***              calculate detrained precipitation             ***
c
       do il=1,ncum
        if (i.le.inb(il) .and. lwork(il)) then
         if (cvflag_grav) then
          wdtrain(il)=grav*ep(il,i)*m(il,i)*clw(il,i)
         else
          wdtrain(il)=10.0*ep(il,i)*m(il,i)*clw(il,i)
         endif
        endif
       enddo

       if(i.gt.1)then
        do 320 j=1,i-1
         do il=1,ncum
          if (i.le.inb(il) .and. lwork(il)) then
           awat=elij(il,j,i)-(1.-ep(il,i))*clw(il,i)
           awat=amax1(awat,0.0)
           if (cvflag_grav) then
            wdtrain(il)=wdtrain(il)+grav*awat*ment(il,j,i)
           else
            wdtrain(il)=wdtrain(il)+10.0*awat*ment(il,j,i)
           endif
          endif
         enddo
320     continue
       endif

c
c    ***    find rain water and evaporation using provisional   ***
c    ***              estimates of rp(i)and rp(i-1)             ***
c

      do 999 il=1,ncum

       if (i.le.inb(il) .and. lwork(il)) then

      wt(il,i)=45.0

      if(i.lt.inb(il))then
       rp(il,i)=rp(il,i+1)
     :       +(cpd*(t(il,i+1)-t(il,i))+gz(il,i+1)-gz(il,i))/lv(il,i)
       rp(il,i)=0.5*(rp(il,i)+rr(il,i))
      endif
      rp(il,i)=amax1(rp(il,i),0.0)
      rp(il,i)=amin1(rp(il,i),rs(il,i))
      rp(il,inb(il))=rr(il,inb(il))

      if(i.eq.1)then
       afac=p(il,1)*(rs(il,1)-rp(il,1))/(1.0e4+2000.0*p(il,1)*rs(il,1))
      else
       rp(il,i-1)=rp(il,i)
     :          +(cpd*(t(il,i)-t(il,i-1))+gz(il,i)-gz(il,i-1))/lv(il,i)
       rp(il,i-1)=0.5*(rp(il,i-1)+rr(il,i-1))
       rp(il,i-1)=amin1(rp(il,i-1),rs(il,i-1))
       rp(il,i-1)=amax1(rp(il,i-1),0.0)
       afac1=p(il,i)*(rs(il,i)-rp(il,i))/(1.0e4+2000.0*p(il,i)*rs(il,i))
       afac2=p(il,i-1)*(rs(il,i-1)-rp(il,i-1))
     :                /(1.0e4+2000.0*p(il,i-1)*rs(il,i-1))
       afac=0.5*(afac1+afac2)
      endif
      if(i.eq.inb(il))afac=0.0
      afac=amax1(afac,0.0)
      bfac=1./(sigd(il)*wt(il,i))
c
cjyg1
ccc        sigt=1.0
ccc        if(i.ge.icb)sigt=sigp(i)
c prise en compte de la variation progressive de sigt dans
c les couches icb et icb-1:
c 	pour plcl<ph(i+1), pr1=0 & pr2=1
c 	pour plcl>ph(i),   pr1=1 & pr2=0
c 	pour ph(i+1)<plcl<ph(i), pr1 est la proportion a cheval
c    sur le nuage, et pr2 est la proportion sous la base du
c    nuage.
      pr1=(plcl(il)-ph(il,i+1))/(ph(il,i)-ph(il,i+1))
      pr1=max(0.,min(1.,pr1))
      pr2=(ph(il,i)-plcl(il))/(ph(il,i)-ph(il,i+1))
      pr2=max(0.,min(1.,pr2))
      sigt=sigp(il,i)*pr1+pr2
cjyg2
c
cjyg----
c       b6 = bfac*100.*sigd(il)*(ph(il,i)-ph(il,i+1))*sigt*afac
c       c6 = water(il,i+1) + wdtrain(il)*bfac
c        revap=0.5*(-b6+sqrt(b6*b6+4.*c6))
c        evap(il,i)=sigt*afac*revap
c        water(il,i)=revap*revap
cc        print *,' i,b6,c6,revap,evap(il,i),water(il,i),wdtrain(il) ',
cc     $            i,b6,c6,revap,evap(il,i),water(il,i),wdtrain(il)
cc---end jyg---
c
c--------retour à la formulation originale d'Emanuel.
      b6=bfac*50.*sigd(il)*(ph(il,i)-ph(il,i+1))*sigt*afac
      c6=water(il,i+1)+bfac*wdtrain(il)
     :    -50.*sigd(il)*bfac*(ph(il,i)-ph(il,i+1))*evap(il,i+1)
      if(c6.gt.0.0)then
       revap=0.5*(-b6+sqrt(b6*b6+4.*c6))
cjyg    Dans sa formulation originale, Emanuel calcule l'evaporation par:
cc             evap(il,i)=sigt*afac*revap
c     ce qui n'est pas correct. Dans cv_routines, la formulation a été modifiee.
c     Ici,l'evaporation evap est simplement calculee par l'equation de
c     conservation.
       water(il,i)=revap*revap
      else
cjyg----   Correction : si c6 <= 0, water(il,i)=0.
       water(il,i) = 0.
      endif
cJYG/IM : ci-dessous formulation originale de KE
c      evap(il,i)=-evap(il,i+1)
c    :            +(wdtrain(il)+sigd(il)*wt(il,i)*water(il,i+1))
c    :                 /(sigd(il)*(ph(il,i)-ph(il,i+1))*50.)
c
cJYG/IM : ci-dessous modification formulation originale de KE
c        pour eliminer oscillations verticales de pluie se produisant
c        lorsqu'il y a evaporation totale de la pluie
c
c       evap(il,i)= +(wdtrain(il)+sigd(il)*wt(il,i)*water(il,i+1)) !itlmd(jyg)
c     :                 /(sigd(il)*(ph(il,i)-ph(il,i+1))*100.)
c      end if  !itlmd(jyg)
cjyg---   Dans tous les cas, evaporation = [tt ce qui entre dans la couche i]
c                                    moins [tt ce qui sort de la couche i]
       evap(il,i)=
     :       (wdtrain(il)+sigd(il)*wt(il,i)*(water(il,i+1)-water(il,i)))
     :                 /(sigd(il)*(ph(il,i)-ph(il,i+1))*100.)
c
ccc
c    ***  calculate precipitating downdraft mass flux under     ***
c    ***              hydrostatic approximation                 ***
c
      if (i.ne.1) then

      tevap=amax1(0.0,evap(il,i))
      delth=amax1(0.001,(th(il,i)-th(il,i-1)))
      if (cvflag_grav) then
       mp(il,i)=100.*ginv*lvcp(il,i)*sigd(il)*tevap
     :              *(p(il,i-1)-p(il,i))/delth
      else
       mp(il,i)=10.*lvcp(il,i)*sigd(il)*tevap
     :         *(p(il,i-1)-p(il,i))/delth
      endif
c
c    ***           if hydrostatic assumption fails,             ***
c    ***   solve cubic difference equation for downdraft theta  ***
c    ***  and mass flux from two simultaneous differential eqns ***
c
      amfac=sigd(il)*sigd(il)*70.0*ph(il,i)*(p(il,i-1)-p(il,i))
     :          *(th(il,i)-th(il,i-1))/(tv(il,i)*th(il,i))
      amp2=abs(mp(il,i+1)*mp(il,i+1)-mp(il,i)*mp(il,i))
      if(amp2.gt.(0.1*amfac))then
       xf=100.0*sigd(il)*sigd(il)*sigd(il)*(ph(il,i)-ph(il,i+1))
       tf=b(il,i)-5.0*(th(il,i)-th(il,i-1))*t(il,i)
     :               /(lvcp(il,i)*sigd(il)*th(il,i))
       af=xf*tf+mp(il,i+1)*mp(il,i+1)*tinv
       bf=2.*(tinv*mp(il,i+1))**3+tinv*mp(il,i+1)*xf*tf
     :            +50.*(p(il,i-1)-p(il,i))*xf*tevap
       fac2=1.0
       if(bf.lt.0.0)fac2=-1.0
       bf=abs(bf)
       ur=0.25*bf*bf-af*af*af*tinv*tinv*tinv
       if(ur.ge.0.0)then
        sru=sqrt(ur)
        fac=1.0
        if((0.5*bf-sru).lt.0.0)fac=-1.0
        mp(il,i)=mp(il,i+1)*tinv+(0.5*bf+sru)**tinv
     :                  +fac*(abs(0.5*bf-sru))**tinv
       else
        d=atan(2.*sqrt(-ur)/(bf+1.0e-28))
        if(fac2.lt.0.0)d=3.14159-d
        mp(il,i)=mp(il,i+1)*tinv+2.*sqrt(af*tinv)*cos(d*tinv)
       endif
       mp(il,i)=amax1(0.0,mp(il,i))

       if (cvflag_grav) then
Cjyg : il y a vraisemblablement une erreur dans la ligne 2 suivante:
C il faut diviser par (mp(il,i)*sigd(il)*grav) et non par (mp(il,i)+sigd(il)*0.1).
C Et il faut bien revoir les facteurs 100.
        b(il,i-1)=b(il,i)+100.0*(p(il,i-1)-p(il,i))*tevap
     2   /(mp(il,i)+sigd(il)*0.1)
     3 -10.0*(th(il,i)-th(il,i-1))*t(il,i)/(lvcp(il,i)
     : *sigd(il)*th(il,i))
       else
        b(il,i-1)=b(il,i)+100.0*(p(il,i-1)-p(il,i))*tevap
     2   /(mp(il,i)+sigd(il)*0.1)
     3 -10.0*(th(il,i)-th(il,i-1))*t(il,i)/(lvcp(il,i)
     : *sigd(il)*th(il,i))
       endif
       b(il,i-1)=amax1(b(il,i-1),0.0)
      endif
c
c   ***         limit magnitude of mp(i) to meet cfl condition      ***
c
      ampmax=2.0*(ph(il,i)-ph(il,i+1))*delti
      amp2=2.0*(ph(il,i-1)-ph(il,i))*delti
      ampmax=amin1(ampmax,amp2)
      mp(il,i)=amin1(mp(il,i),ampmax)
c
c    ***      force mp to decrease linearly to zero                 ***
c    ***       between cloud base and the surface                   ***
c
c
cc      if(p(il,i).gt.p(il,icb(il)))then
cc       mp(il,i)=mp(il,icb(il))*(p(il,1)-p(il,i))/(p(il,1)-p(il,icb(il)))
cc      endif
      if(ph(il,i) .gt. 0.9*plcl(il)) then
       mp(il,i) = mp(il,i)*(ph(il,1)-ph(il,i))/
     $                     (ph(il,1)-0.9*plcl(il))
      endif

360   continue
      endif ! i.eq.1
c
c    ***       find mixing ratio of precipitating downdraft     ***
c

      if (i.ne.inb(il)) then

      rp(il,i)=rr(il,i)

      if(mp(il,i).gt.mp(il,i+1))then

       if (cvflag_grav) then
        rp(il,i)=rp(il,i+1)*mp(il,i+1)+rr(il,i)*(mp(il,i)-mp(il,i+1))
     :   +100.*ginv*0.5*sigd(il)*(ph(il,i)-ph(il,i+1))
     :                     *(evap(il,i+1)+evap(il,i))
       else
        rp(il,i)=rp(il,i+1)*mp(il,i+1)+rr(il,i)*(mp(il,i)-mp(il,i+1))
     :   +5.*sigd(il)*(ph(il,i)-ph(il,i+1))
     :                      *(evap(il,i+1)+evap(il,i))
       endif
      rp(il,i)=rp(il,i)/mp(il,i)
      up(il,i)=up(il,i+1)*mp(il,i+1)+u(il,i)*(mp(il,i)-mp(il,i+1))
      up(il,i)=up(il,i)/mp(il,i)
      vp(il,i)=vp(il,i+1)*mp(il,i+1)+v(il,i)*(mp(il,i)-mp(il,i+1))
      vp(il,i)=vp(il,i)/mp(il,i)

      do j=1,ntra
      trap(il,i,j)=trap(il,i+1,j)*mp(il,i+1)
     :            +trap(il,i,j)*(mp(il,i)-mp(il,i+1))
      trap(il,i,j)=trap(il,i,j)/mp(il,i)
      end do

      else

       if(mp(il,i+1).gt.1.0e-16)then
        if (cvflag_grav) then
         rp(il,i)=rp(il,i+1)
     :            +100.*ginv*0.5*sigd(il)*(ph(il,i)-ph(il,i+1))
     :            *(evap(il,i+1)+evap(il,i))/mp(il,i+1)
        else
         rp(il,i)=rp(il,i+1)
     :           +5.*sigd(il)*(ph(il,i)-ph(il,i+1))
     :           *(evap(il,i+1)+evap(il,i))/mp(il,i+1)
        endif
       up(il,i)=up(il,i+1)
       vp(il,i)=vp(il,i+1)

       do j=1,ntra
       trap(il,i,j)=trap(il,i+1,j)
       end do

       endif
      endif
      rp(il,i)=amin1(rp(il,i),rs(il,i))
      rp(il,i)=amax1(rp(il,i),0.0)

      endif
      endif
999   continue

400   continue

       return
       end

      SUBROUTINE cv3_yield(nloc,ncum,nd,na,ntra
     :                    ,icb,inb,delt
     :                    ,t,rr,t_wake,rr_wake,s_wake,u,v,tra
     :                    ,gz,p,ph,h,hp,lv,cpn,th,th_wake
     :                    ,ep,clw,m,tp,mp,rp,up,vp,trap
     :                    ,wt,water,evap,b,sigd
     :                    ,ment,qent,hent,iflag_mix,uent,vent
     :                    ,nent,elij,traent,sig
     :                    ,tv,tvp,wghti
     :                    ,iflag,precip,Vprecip,ft,fr,fu,fv,ftra
     :                    ,cbmf,upwd,dnwd,dnwd0,ma,mip
     :                    ,tls,tps,qcondc,wd
     :                    ,ftd,fqd)
      
      implicit none

#include "cvthermo.h"
#include "cv3param.h"
#include "cvflag.h"
#include "conema3.h"

c inputs:
c      print*,'cv3_yield apres include'
      integer iflag_mix
      integer ncum,nd,na,ntra,nloc
      integer icb(nloc), inb(nloc)
      real delt
      real t(nloc,nd), rr(nloc,nd), u(nloc,nd), v(nloc,nd)
      real t_wake(nloc,nd), rr_wake(nloc,nd)
      real s_wake(nloc)
      real tra(nloc,nd,ntra), sig(nloc,nd)
      real gz(nloc,na), ph(nloc,nd+1), h(nloc,na), hp(nloc,na)
      real th(nloc,na), p(nloc,nd), tp(nloc,na)
      real lv(nloc,na), cpn(nloc,na), ep(nloc,na), clw(nloc,na)
      real m(nloc,na), mp(nloc,na), rp(nloc,na), up(nloc,na)
      real vp(nloc,na), wt(nloc,nd), trap(nloc,nd,ntra)
      real water(nloc,na), evap(nloc,na), b(nloc,na), sigd(nloc)
      real ment(nloc,na,na), qent(nloc,na,na), uent(nloc,na,na)
      real hent(nloc,na,na)
cIM bug   real vent(nloc,na,na), nent(nloc,na), elij(nloc,na,na)
      real vent(nloc,na,na), elij(nloc,na,na)
      integer nent(nloc,nd)
      real traent(nloc,na,na,ntra)
      real tv(nloc,nd), tvp(nloc,nd), wghti(nloc,nd)
c      print*,'cv3_yield declarations 1'
c input/output:
      integer iflag(nloc)

c outputs:
      real precip(nloc)
      real ft(nloc,nd), fr(nloc,nd), fu(nloc,nd), fv(nloc,nd)
      real ftd(nloc,nd), fqd(nloc,nd)
      real ftra(nloc,nd,ntra)
      real upwd(nloc,nd), dnwd(nloc,nd), ma(nloc,nd)
      real dnwd0(nloc,nd), mip(nloc,nd)
      real Vprecip(nloc,nd+1)
      real tls(nloc,nd), tps(nloc,nd)
      real qcondc(nloc,nd)                               ! cld
      real wd(nloc)                                      ! gust
      real cbmf(nloc)
c      print*,'cv3_yield declarations 2'
c local variables:
      integer i,k,il,n,j,num1
      real rat, delti
      real ax, bx, cx, dx, ex
      real cpinv, rdcp, dpinv
      real awat(nloc)
      real lvcp(nloc,na), mke(nloc,na)
      real am(nloc), work(nloc), ad(nloc), amp1(nloc)
c!!      real up1(nloc), dn1(nloc)
      real up1(nloc,nd,nd), dn1(nloc,nd,nd)
      real asum(nloc), bsum(nloc), csum(nloc), dsum(nloc)
      real esum(nloc), fsum(nloc), gsum(nloc), hsum(nloc)
      real th_wake(nloc,nd)
      real alpha_qpos(nloc),alpha_qpos1(nloc)
      real qcond(nloc,nd), nqcond(nloc,nd), wa(nloc,nd)  ! cld
      real siga(nloc,nd), sax(nloc,nd), mac(nloc,nd)      ! cld

c      print*,'cv3_yield declarations 3'
c-------------------------------------------------------------

c initialization:

      delti = 1.0/delt
c      print*,'cv3_yield initialisation delt', delt
cprecip,Vprecip,ft,fr,fu,fv,ftra
c     :                    ,cbmf,upwd,dnwd,dnwd0,ma,mip
c     :                    ,tls,tps,qcondc,wd
c     :                    ,ftd,fqd  )
      do il=1,ncum
       precip(il)=0.0
       Vprecip(il,nd+1)=0.0
       wd(il)=0.0     ! gust
      enddo

      do i=1,nd
       do il=1,ncum
         Vprecip(il,i)=0.0
         ft(il,i)=0.0
         fr(il,i)=0.0
         fu(il,i)=0.0
         fv(il,i)=0.0
         upwd(il,i)=0.0
         dnwd(il,i)=0.0
         dnwd0(il,i)=0.0
         mip(il,i)=0.0
         ftd(il,i)=0.0
         fqd(il,i)=0.0
         qcondc(il,i)=0.0                                ! cld
         qcond(il,i)=0.0                                 ! cld
         nqcond(il,i)=0.0                                ! cld
       enddo 
      enddo
c       print*,'cv3_yield initialisation 2'
      do j=1,ntra
       do i=1,nd
        do il=1,ncum
          ftra(il,i,j)=0.0
        enddo
       enddo
      enddo
c       print*,'cv3_yield initialisation 3'
      do i=1,nl
       do il=1,ncum
         lvcp(il,i)=lv(il,i)/cpn(il,i)
       enddo
      enddo


c
c   ***  calculate surface precipitation in mm/day     ***
c
      do il=1,ncum
       if(ep(il,inb(il)).ge.0.0001 .and. iflag(il) .le. 1)then
        if (cvflag_grav) then
           precip(il)=wt(il,1)*sigd(il)*water(il,1)*86400.*1000.
     :               /(rowl*grav)
        else
         precip(il)=wt(il,1)*sigd(il)*water(il,1)*8640.
        endif
       endif
      enddo
c      print*,'cv3_yield apres calcul precip'

C
C   ===  calculate vertical profile of  precipitation in kg/m2/s  ===
C
      do i = 1,nl
      do il=1,ncum
       if(ep(il,inb(il)).ge.0.0001 .and. i.le.inb(il)
     :    .and. iflag(il) .le. 1)then
        if (cvflag_grav) then
           VPrecip(il,i) = wt(il,i)*sigd(il)*water(il,i)/grav
        else
           VPrecip(il,i) = wt(il,i)*sigd(il)*water(il,i)/10.
        endif
       endif
      enddo
      enddo
C
c
c   ***  Calculate downdraft velocity scale    ***
c   ***  NE PAS UTILISER POUR L'INSTANT ***
c
c!      do il=1,ncum
c!        wd(il)=betad*abs(mp(il,icb(il)))*0.01*rrd*t(il,icb(il))
c!     :                                  /(sigd(il)*p(il,icb(il)))
c!      enddo

c
c   ***  calculate tendencies of lowest level potential temperature  ***
c   ***                      and mixing ratio                        ***
c
      do il=1,ncum
       work(il)=1.0/(ph(il,1)-ph(il,2))
       cbmf(il)=0.0
      enddo

      do k=2,nl
       do il=1,ncum
        if (k.ge.icb(il)) then
         cbmf(il)=cbmf(il)+m(il,k)
        endif
       enddo
      enddo

c      print*,'cv3_yield avant ft'
c AM is the part of cbmf taken from the first level
      do il=1,ncum
        am(il)=cbmf(il)*wghti(il,1)
      enddo
c
      do il=1,ncum
        if (iflag(il) .le. 1) then
c convect3      if((0.1*dpinv*am).ge.delti)iflag(il)=4
cjyg  Correction pour conserver l'eau
ccc       ft(il,1)=-0.5*lvcp(il,1)*sigd(il)*(evap(il,1)+evap(il,2)) !precip
       ft(il,1)=-lvcp(il,1)*sigd(il)*evap(il,1)                  !precip

      if (cvflag_grav) then
        ft(il,1)=ft(il,1)-0.009*grav*sigd(il)*mp(il,2)
     :                              *t_wake(il,1)*b(il,1)*work(il)
      else
        ft(il,1)=ft(il,1)-0.09*sigd(il)*mp(il,2)
     :                              *t_wake(il,1)*b(il,1)*work(il)
      endif

      ft(il,1)=ft(il,1)+0.01*sigd(il)*wt(il,1)*(cl-cpd)*water(il,2)
     :     *(t_wake(il,2)-t_wake(il,1))*work(il)/cpn(il,1)

      ftd(il,1) = ft(il,1)                        ! fin precip

      if (cvflag_grav) then                  !sature
      if((0.01*grav*work(il)*am(il)).ge.delti)iflag(il)=1!consist vect
       ft(il,1)=ft(il,1)+0.01*grav*work(il)*am(il)*(t(il,2)-t(il,1)
     :            +(gz(il,2)-gz(il,1))/cpn(il,1))
      else
       if((0.1*work(il)*am(il)).ge.delti)iflag(il)=1 !consistency vect
       ft(il,1)=ft(il,1)+0.1*work(il)*am(il)*(t(il,2)-t(il,1)
     :            +(gz(il,2)-gz(il,1))/cpn(il,1))
      endif
      endif  ! iflag
      enddo


       do j=2,nl
      IF (iflag_mix .gt. 0) then
        do il=1,ncum
c FH WARNING a modifier :
      cpinv=0.
c     cpinv=1.0/cpn(il,1)
         if (j.le.inb(il) .and. iflag(il) .le. 1) then
         if (cvflag_grav) then
          ft(il,1)=ft(il,1)
     :       +0.01*grav*work(il)*ment(il,j,1)*(hent(il,j,1)-h(il,1)
     :       +t(il,1)*(cpv-cpd)*(rr(il,1)-Qent(il,j,1)))*cpinv
         else
          ft(il,1)=ft(il,1)
     :       +0.1*work(il)*ment(il,j,1)*(hent(il,j,1)-h(il,1)
     :       +t(il,1)*(cpv-cpd)*(rr(il,1)-Qent(il,j,1)))*cpinv
         endif  ! cvflag_grav
        endif ! j
       enddo
       ENDIF
        enddo
         ! fin sature


      do il=1,ncum
        if (iflag(il) .le. 1) then
          if (cvflag_grav) then
Cjyg1  Correction pour mieux conserver l'eau (conformite avec CONVECT4.3)
       fr(il,1)=0.01*grav*mp(il,2)*(rp(il,2)-rr_wake(il,1))*work(il)
     :          +sigd(il)*evap(il,1)
ccc     :          +sigd(il)*0.5*(evap(il,1)+evap(il,2))

       fqd(il,1)=fr(il,1)     !precip

       fr(il,1)=fr(il,1)+0.01*grav*am(il)*(rr(il,2)-rr(il,1))*work(il)  !sature

       fu(il,1)=fu(il,1)+0.01*grav*work(il)*(mp(il,2)*(up(il,2)-u(il,1))
     :         +am(il)*(u(il,2)-u(il,1)))
       fv(il,1)=fv(il,1)+0.01*grav*work(il)*(mp(il,2)*(vp(il,2)-v(il,1))
     :         +am(il)*(v(il,2)-v(il,1)))
      else  ! cvflag_grav
       fr(il,1)=0.1*mp(il,2)*(rp(il,2)-rr_wake(il,1))*work(il)
     :          +sigd(il)*evap(il,1)
ccc     :          +sigd(il)*0.5*(evap(il,1)+evap(il,2))
       fqd(il,1)=fr(il,1)  !precip
       fr(il,1)=fr(il,1)+0.1*am(il)*(rr(il,2)-rr(il,1))*work(il)
       fu(il,1)=fu(il,1)+0.1*work(il)*(mp(il,2)*(up(il,2)-u(il,1))
     :         +am(il)*(u(il,2)-u(il,1)))
       fv(il,1)=fv(il,1)+0.1*work(il)*(mp(il,2)*(vp(il,2)-v(il,1))
     :         +am(il)*(v(il,2)-v(il,1)))
         endif ! cvflag_grav
       endif  ! iflag
      enddo ! il


      do j=1,ntra
       do il=1,ncum
        if (iflag(il) .le. 1) then
        if (cvflag_grav) then
         ftra(il,1,j)=ftra(il,1,j)+0.01*grav*work(il)
     :                     *(mp(il,2)*(trap(il,2,j)-tra(il,1,j))
     :             +am(il)*(tra(il,2,j)-tra(il,1,j)))
        else
         ftra(il,1,j)=ftra(il,1,j)+0.1*work(il)
     :                     *(mp(il,2)*(trap(il,2,j)-tra(il,1,j))
     :             +am(il)*(tra(il,2,j)-tra(il,1,j)))
        endif
        endif  ! iflag
       enddo
      enddo

       do j=2,nl
       do il=1,ncum
        if (j.le.inb(il) .and. iflag(il) .le. 1) then
         if (cvflag_grav) then
          fr(il,1)=fr(il,1)
     :       +0.01*grav*work(il)*ment(il,j,1)*(qent(il,j,1)-rr(il,1))
          fu(il,1)=fu(il,1)
     :       +0.01*grav*work(il)*ment(il,j,1)*(uent(il,j,1)-u(il,1))
          fv(il,1)=fv(il,1)
     :       +0.01*grav*work(il)*ment(il,j,1)*(vent(il,j,1)-v(il,1))
         else   ! cvflag_grav
          fr(il,1)=fr(il,1)
     :         +0.1*work(il)*ment(il,j,1)*(qent(il,j,1)-rr(il,1))
          fu(il,1)=fu(il,1)
     :         +0.1*work(il)*ment(il,j,1)*(uent(il,j,1)-u(il,1))
          fv(il,1)=fv(il,1)
     :         +0.1*work(il)*ment(il,j,1)*(vent(il,j,1)-v(il,1))  ! fin sature
         endif  ! cvflag_grav
        endif ! j
       enddo
      enddo

      do k=1,ntra
       do j=2,nl
        do il=1,ncum
         if (j.le.inb(il) .and. iflag(il) .le. 1) then

          if (cvflag_grav) then
           ftra(il,1,k)=ftra(il,1,k)+0.01*grav*work(il)*ment(il,j,1)
     :                *(traent(il,j,1,k)-tra(il,1,k))
          else
           ftra(il,1,k)=ftra(il,1,k)+0.1*work(il)*ment(il,j,1)
     :                *(traent(il,j,1,k)-tra(il,1,k))
          endif

         endif
        enddo
       enddo
      enddo
c      print*,'cv3_yield apres ft'
c
c   ***  calculate tendencies of potential temperature and mixing ratio  ***
c   ***               at levels above the lowest level                   ***
c
c   ***  first find the net saturated updraft and downdraft mass fluxes  ***
c   ***                      through each level                          ***
c

      do 500 i=2,nl+1 ! newvecto: mettre nl au lieu nl+1?

       num1=0
       do il=1,ncum
        if(i.le.inb(il) .and. iflag(il) .le. 1)num1=num1+1
       enddo
       if(num1.le.0)go to 500

       call zilch(amp1,ncum)
       call zilch(ad,ncum)

      do 440 k=1,nl+1
       do 441 il=1,ncum
        if(i.ge.icb(il)) then
          if(k.ge.i+1.and. k.le.(inb(il)+1)) then
            amp1(il)=amp1(il)+m(il,k)
          endif
         else
c AMP1 is the part of cbmf taken from layers I and lower
          if(k.le.i) then
           amp1(il)=amp1(il)+cbmf(il)*wghti(il,k)
          endif
        endif
 441   continue
 440  continue

      do 450 k=1,i
       do 451 j=i+1,nl+1
        do 452 il=1,ncum
         if (i.le.inb(il) .and. j.le.(inb(il)+1)) then
          amp1(il)=amp1(il)+ment(il,k,j)
         endif
452     continue
451    continue
450   continue

      do 470 k=1,i-1
       do 471 j=i,nl+1 ! newvecto: nl au lieu nl+1?
        do 472 il=1,ncum
        if (i.le.inb(il) .and. j.le.inb(il)) then
         ad(il)=ad(il)+ment(il,j,k)
        endif
472     continue
471    continue
470   continue
  
      do 1350 il=1,ncum
      if (i.le.inb(il) .and. iflag(il) .le. 1) then
       dpinv=1.0/(ph(il,i)-ph(il,i+1))
       cpinv=1.0/cpn(il,i)

c convect3      if((0.1*dpinv*amp1).ge.delti)iflag(il)=4
      if (cvflag_grav) then
       if((0.01*grav*dpinv*amp1(il)).ge.delti)iflag(il)=1 ! vecto
      else
       if((0.1*dpinv*amp1(il)).ge.delti)iflag(il)=1 ! vecto
      endif

       ! precip
ccc       ft(il,i)= -0.5*sigd(il)*lvcp(il,i)*(evap(il,i)+evap(il,i+1))
       ft(il,i)= -sigd(il)*lvcp(il,i)*evap(il,i)
        rat=cpn(il,i-1)*cpinv
c
      if (cvflag_grav) then
       ft(il,i)=ft(il,i)-0.009*grav*sigd(il)
     :   *(mp(il,i+1)*t_wake(il,i)*b(il,i)
     :   -mp(il,i)*t_wake(il,i-1)*rat*b(il,i-1))*dpinv
       ft(il,i)=ft(il,i)+0.01*sigd(il)*wt(il,i)*(cl-cpd)*water(il,i+1)
     :           *(t_wake(il,i+1)-t_wake(il,i))*dpinv*cpinv
       ftd(il,i)=ft(il,i)
        ! fin precip
c
           ! sature
       ft(il,i)=ft(il,i)+0.01*grav*dpinv*(amp1(il)*(t(il,i+1)-t(il,i)
     :    +(gz(il,i+1)-gz(il,i))*cpinv)
     :    -ad(il)*(t(il,i)-t(il,i-1)+(gz(il,i)-gz(il,i-1))*cpinv))

c
      IF (iflag_mix .eq. 0) then
       ft(il,i)=ft(il,i)+0.01*grav*dpinv*ment(il,i,i)*(hp(il,i)-h(il,i)
     :    +t(il,i)*(cpv-cpd)*(rr(il,i)-qent(il,i,i)))*cpinv
      ENDIF
c
      else  ! cvflag_grav
       ft(il,i)=ft(il,i)-0.09*sigd(il)
     :   *(mp(il,i+1)*t_wake(il,i)*b(il,i)
     :   -mp(il,i)*t_wake(il,i-1)*rat*b(il,i-1))*dpinv
       ft(il,i)=ft(il,i)+0.01*sigd(il)*wt(il,i)*(cl-cpd)*water(il,i+1)
     :           *(t_wake(il,i+1)-t_wake(il,i))*dpinv*cpinv
       ftd(il,i)=ft(il,i)
        ! fin precip
c
           ! sature
       ft(il,i)=ft(il,i)+0.1*dpinv*(amp1(il)*(t(il,i+1)-t(il,i)
     :    +(gz(il,i+1)-gz(il,i))*cpinv)
     :    -ad(il)*(t(il,i)-t(il,i-1)+(gz(il,i)-gz(il,i-1))*cpinv))

c
      IF (iflag_mix .eq. 0) then
       ft(il,i)=ft(il,i)+0.1*dpinv*ment(il,i,i)*(hp(il,i)-h(il,i)
     :    +t(il,i)*(cpv-cpd)*(rr(il,i)-qent(il,i,i)))*cpinv
      ENDIF
      endif ! cvflag_grav


        if (cvflag_grav) then
c sb: on ne fait pas encore la correction permettant de mieux
c conserver l'eau:
c jyg: correction permettant de mieux conserver l'eau:
ccc         fr(il,i)=0.5*sigd(il)*(evap(il,i)+evap(il,i+1))
         fr(il,i)=sigd(il)*evap(il,i)
     :        +0.01*grav*(mp(il,i+1)*(rp(il,i+1)-rr_wake(il,i))
     :        -mp(il,i)*(rp(il,i)-rr_wake(il,i-1)))*dpinv
         fqd(il,i)=fr(il,i)    ! precip

         fu(il,i)=0.01*grav*(mp(il,i+1)*(up(il,i+1)-u(il,i))
     :             -mp(il,i)*(up(il,i)-u(il,i-1)))*dpinv
         fv(il,i)=0.01*grav*(mp(il,i+1)*(vp(il,i+1)-v(il,i))
     :             -mp(il,i)*(vp(il,i)-v(il,i-1)))*dpinv
        else  ! cvflag_grav
ccc         fr(il,i)=0.5*sigd(il)*(evap(il,i)+evap(il,i+1))
         fr(il,i)=sigd(il)*evap(il,i)
     :        +0.1*(mp(il,i+1)*(rp(il,i+1)-rr_wake(il,i))
     :             -mp(il,i)*(rp(il,i)-rr_wake(il,i-1)))*dpinv
         fqd(il,i)=fr(il,i)    ! precip

         fu(il,i)=0.1*(mp(il,i+1)*(up(il,i+1)-u(il,i))
     :             -mp(il,i)*(up(il,i)-u(il,i-1)))*dpinv
         fv(il,i)=0.1*(mp(il,i+1)*(vp(il,i+1)-v(il,i))
     :             -mp(il,i)*(vp(il,i)-v(il,i-1)))*dpinv
        endif ! cvflag_grav


      if (cvflag_grav) then
       fr(il,i)=fr(il,i)+0.01*grav*dpinv*(amp1(il)*(rr(il,i+1)-rr(il,i))
     :           -ad(il)*(rr(il,i)-rr(il,i-1)))
       fu(il,i)=fu(il,i)+0.01*grav*dpinv*(amp1(il)*(u(il,i+1)-u(il,i))
     :             -ad(il)*(u(il,i)-u(il,i-1)))
       fv(il,i)=fv(il,i)+0.01*grav*dpinv*(amp1(il)*(v(il,i+1)-v(il,i))
     :             -ad(il)*(v(il,i)-v(il,i-1)))
      else  ! cvflag_grav
       fr(il,i)=fr(il,i)+0.1*dpinv*(amp1(il)*(rr(il,i+1)-rr(il,i))
     :           -ad(il)*(rr(il,i)-rr(il,i-1)))
       fu(il,i)=fu(il,i)+0.1*dpinv*(amp1(il)*(u(il,i+1)-u(il,i))
     :             -ad(il)*(u(il,i)-u(il,i-1)))
       fv(il,i)=fv(il,i)+0.1*dpinv*(amp1(il)*(v(il,i+1)-v(il,i))
     :             -ad(il)*(v(il,i)-v(il,i-1)))
      endif ! cvflag_grav

      endif ! i
1350  continue

      do k=1,ntra
       do il=1,ncum
        if (i.le.inb(il) .and. iflag(il) .le. 1) then
         dpinv=1.0/(ph(il,i)-ph(il,i+1))
         cpinv=1.0/cpn(il,i)
         if (cvflag_grav) then
           ftra(il,i,k)=ftra(il,i,k)+0.01*grav*dpinv
     :         *(amp1(il)*(tra(il,i+1,k)-tra(il,i,k))
     :           -ad(il)*(tra(il,i,k)-tra(il,i-1,k)))
         else
           ftra(il,i,k)=ftra(il,i,k)+0.1*dpinv
     :         *(amp1(il)*(tra(il,i+1,k)-tra(il,i,k))
     :           -ad(il)*(tra(il,i,k)-tra(il,i-1,k)))
         endif
        endif
       enddo
      enddo

      do 480 k=1,i-1
c
       do il = 1,ncum
        awat(il)=elij(il,k,i)-(1.-ep(il,i))*clw(il,i)
        awat(il)=amax1(awat(il),0.0)
       enddo
c
      IF (iflag_mix .ne. 0) then
       do il=1,ncum
        if (i.le.inb(il) .and. iflag(il) .le. 1) then
         dpinv=1.0/(ph(il,i)-ph(il,i+1))
         cpinv=1.0/cpn(il,i)
       if (cvflag_grav) then
      ft(il,i)=ft(il,i)
     :       +0.01*grav*dpinv*ment(il,k,i)*(hent(il,k,i)-h(il,i)
     :       +t(il,i)*(cpv-cpd)*(rr(il,i)+awat(il)-Qent(il,k,i)))*cpinv

c
c
       else
      ft(il,i)=ft(il,i)
     :       +0.1*dpinv*ment(il,k,i)*(hent(il,k,i)-h(il,i)
     :       +t(il,i)*(cpv-cpd)*(rr(il,i)+awat(il)-Qent(il,k,i)))*cpinv
       endif  !cvflag_grav
       endif  ! i
       enddo
      ENDIF
c
       do 1370 il=1,ncum
        if (i.le.inb(il) .and. iflag(il) .le. 1) then
         dpinv=1.0/(ph(il,i)-ph(il,i+1))
         cpinv=1.0/cpn(il,i)
       if (cvflag_grav) then
      fr(il,i)=fr(il,i)
     :   +0.01*grav*dpinv*ment(il,k,i)*(qent(il,k,i)-awat(il)-rr(il,i))
      fu(il,i)=fu(il,i)
     :         +0.01*grav*dpinv*ment(il,k,i)*(uent(il,k,i)-u(il,i))
      fv(il,i)=fv(il,i)
     :         +0.01*grav*dpinv*ment(il,k,i)*(vent(il,k,i)-v(il,i))
      else  ! cvflag_grav
      fr(il,i)=fr(il,i)
     :   +0.1*dpinv*ment(il,k,i)*(qent(il,k,i)-awat(il)-rr(il,i))
      fu(il,i)=fu(il,i)
     :   +0.01*grav*dpinv*ment(il,k,i)*(uent(il,k,i)-u(il,i))
      fv(il,i)=fv(il,i)
     :   +0.1*dpinv*ment(il,k,i)*(vent(il,k,i)-v(il,i))
      endif ! cvflag_grav

c (saturated updrafts resulting from mixing)        ! cld
        qcond(il,i)=qcond(il,i)+(elij(il,k,i)-awat(il)) ! cld
        nqcond(il,i)=nqcond(il,i)+1.                ! cld
      endif ! i
1370  continue
480   continue

      do j=1,ntra
       do k=1,i-1
        do il=1,ncum
         if (i.le.inb(il) .and. iflag(il) .le. 1) then
          dpinv=1.0/(ph(il,i)-ph(il,i+1))
          cpinv=1.0/cpn(il,i)
          if (cvflag_grav) then
           ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv*ment(il,k,i)
     :        *(traent(il,k,i,j)-tra(il,i,j))
          else
           ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv*ment(il,k,i)
     :        *(traent(il,k,i,j)-tra(il,i,j))
          endif
         endif
        enddo
       enddo
      enddo

      do 490 k=i,nl+1
c
      IF (iflag_mix .ne. 0) then
       do il=1,ncum
        if (i.le.inb(il) .and. k.le.inb(il)
     $               .and. iflag(il) .le. 1) then
         dpinv=1.0/(ph(il,i)-ph(il,i+1))
         cpinv=1.0/cpn(il,i)
       if (cvflag_grav) then
      ft(il,i)=ft(il,i)
     :       +0.01*grav*dpinv*ment(il,k,i)*(hent(il,k,i)-h(il,i)
     :       +t(il,i)*(cpv-cpd)*(rr(il,i)-Qent(il,k,i)))*cpinv
c
c
       else
      ft(il,i)=ft(il,i)
     :       +0.1*dpinv*ment(il,k,i)*(hent(il,k,i)-h(il,i)
     :       +t(il,i)*(cpv-cpd)*(rr(il,i)-Qent(il,k,i)))*cpinv
       endif  !cvflag_grav
       endif  ! i
       enddo
      ENDIF
c
       do 1380 il=1,ncum
        if (i.le.inb(il) .and. k.le.inb(il)
     $               .and. iflag(il) .le. 1) then
         dpinv=1.0/(ph(il,i)-ph(il,i+1))
         cpinv=1.0/cpn(il,i)

         if (cvflag_grav) then
         fr(il,i)=fr(il,i)
     :         +0.01*grav*dpinv*ment(il,k,i)*(qent(il,k,i)-rr(il,i))
         fu(il,i)=fu(il,i)
     :         +0.01*grav*dpinv*ment(il,k,i)*(uent(il,k,i)-u(il,i))
         fv(il,i)=fv(il,i)
     :         +0.01*grav*dpinv*ment(il,k,i)*(vent(il,k,i)-v(il,i))
         else  ! cvflag_grav
         fr(il,i)=fr(il,i)
     :         +0.1*dpinv*ment(il,k,i)*(qent(il,k,i)-rr(il,i))
         fu(il,i)=fu(il,i)
     :         +0.1*dpinv*ment(il,k,i)*(uent(il,k,i)-u(il,i))
         fv(il,i)=fv(il,i)
     :         +0.1*dpinv*ment(il,k,i)*(vent(il,k,i)-v(il,i))
         endif ! cvflag_grav
        endif ! i and k
1380   continue
490   continue

      do j=1,ntra
       do k=i,nl+1
        do il=1,ncum
         if (i.le.inb(il) .and. k.le.inb(il)
     $                .and. iflag(il) .le. 1) then
          dpinv=1.0/(ph(il,i)-ph(il,i+1))
          cpinv=1.0/cpn(il,i)
          if (cvflag_grav) then
           ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv*ment(il,k,i)
     :         *(traent(il,k,i,j)-tra(il,i,j))
          else
           ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv*ment(il,k,i)
     :             *(traent(il,k,i,j)-tra(il,i,j))
          endif
         endif ! i and k
        enddo
       enddo
      enddo

c sb: interface with the cloud parameterization:          ! cld

      do k=i+1,nl
       do il=1,ncum
        if (k.le.inb(il) .and. i.le.inb(il)
     $               .and. iflag(il) .le. 1) then         ! cld
C (saturated downdrafts resulting from mixing)            ! cld
          qcond(il,i)=qcond(il,i)+elij(il,k,i)            ! cld
          nqcond(il,i)=nqcond(il,i)+1.                    ! cld
        endif                                             ! cld
       enddo                                              ! cld
      enddo                                               ! cld

C (particular case: no detraining level is found)         ! cld
      do il=1,ncum                                        ! cld
       if (i.le.inb(il) .and. nent(il,i).eq.0
     $                 .and. iflag(il) .le. 1) then       ! cld
          qcond(il,i)=qcond(il,i)+(1.-ep(il,i))*clw(il,i) ! cld
          nqcond(il,i)=nqcond(il,i)+1.                    ! cld
       endif                                              ! cld
      enddo                                               ! cld

      do il=1,ncum                                        ! cld
       if (i.le.inb(il) .and. nqcond(il,i).ne.0
     $                   .and. iflag(il) .le. 1) then     ! cld
          qcond(il,i)=qcond(il,i)/nqcond(il,i)            ! cld
       endif                                              ! cld
      enddo

      do j=1,ntra
       do il=1,ncum
        if (i.le.inb(il) .and. iflag(il) .le. 1) then
         dpinv=1.0/(ph(il,i)-ph(il,i+1))
         cpinv=1.0/cpn(il,i)

         if (cvflag_grav) then
          ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv
     :     *(mp(il,i+1)*(trap(il,i+1,j)-tra(il,i,j))
     :     -mp(il,i)*(trap(il,i,j)-trap(il,i-1,j)))
         else
          ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv
     :     *(mp(il,i+1)*(trap(il,i+1,j)-tra(il,i,j))
     :     -mp(il,i)*(trap(il,i,j)-trap(il,i-1,j)))
         endif
        endif ! i
       enddo
      enddo


500   continue


c   ***   move the detrainment at level inb down to level inb-1   ***
c   ***        in such a way as to preserve the vertically        ***
c   ***          integrated enthalpy and water tendencies         ***
c
c Correction bug le 18-03-09
      do 503 il=1,ncum
      IF (iflag(il) .le. 1) THEN
        if (cvflag_grav) then
      ax=0.01*grav*ment(il,inb(il),inb(il))*(hp(il,inb(il))
     : -h(il,inb(il))+t(il,inb(il))*(cpv-cpd)
     : *(rr(il,inb(il))-qent(il,inb(il),inb(il))))
     :  /(cpn(il,inb(il))*(ph(il,inb(il))-ph(il,inb(il)+1)))
      ft(il,inb(il))=ft(il,inb(il))-ax
      ft(il,inb(il)-1)=ft(il,inb(il)-1)+ax*cpn(il,inb(il))
     :    *(ph(il,inb(il))-ph(il,inb(il)+1))/(cpn(il,inb(il)-1)
     :    *(ph(il,inb(il)-1)-ph(il,inb(il))))

      bx=0.01*grav*ment(il,inb(il),inb(il))*(qent(il,inb(il),inb(il))
     :    -rr(il,inb(il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
      fr(il,inb(il))=fr(il,inb(il))-bx
      fr(il,inb(il)-1)=fr(il,inb(il)-1)
     :   +bx*(ph(il,inb(il))-ph(il,inb(il)+1))
     :      /(ph(il,inb(il)-1)-ph(il,inb(il)))

      cx=0.01*grav*ment(il,inb(il),inb(il))*(uent(il,inb(il),inb(il))
     :       -u(il,inb(il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
      fu(il,inb(il))=fu(il,inb(il))-cx
      fu(il,inb(il)-1)=fu(il,inb(il)-1)
     :     +cx*(ph(il,inb(il))-ph(il,inb(il)+1))
     :        /(ph(il,inb(il)-1)-ph(il,inb(il)))

      dx=0.01*grav*ment(il,inb(il),inb(il))*(vent(il,inb(il),inb(il))
     :      -v(il,inb(il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
      fv(il,inb(il))=fv(il,inb(il))-dx
      fv(il,inb(il)-1)=fv(il,inb(il)-1)
     :    +dx*(ph(il,inb(il))-ph(il,inb(il)+1))
     :       /(ph(il,inb(il)-1)-ph(il,inb(il)))
       else
       ax=0.1*ment(il,inb(il),inb(il))*(hp(il,inb(il))
     : -h(il,inb(il))+t(il,inb(il))*(cpv-cpd)
     : *(rr(il,inb(il))-qent(il,inb(il),inb(il))))
     :  /(cpn(il,inb(il))*(ph(il,inb(il))-ph(il,inb(il)+1)))
      ft(il,inb(il))=ft(il,inb(il))-ax
      ft(il,inb(il)-1)=ft(il,inb(il)-1)+ax*cpn(il,inb(il))
     :    *(ph(il,inb(il))-ph(il,inb(il)+1))/(cpn(il,inb(il)-1)
     :    *(ph(il,inb(il)-1)-ph(il,inb(il))))

      bx=0.1*ment(il,inb(il),inb(il))*(qent(il,inb(il),inb(il))
     :    -rr(il,inb(il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
      fr(il,inb(il))=fr(il,inb(il))-bx
      fr(il,inb(il)-1)=fr(il,inb(il)-1)
     :   +bx*(ph(il,inb(il))-ph(il,inb(il)+1))
     :      /(ph(il,inb(il)-1)-ph(il,inb(il)))

      cx=0.1*ment(il,inb(il),inb(il))*(uent(il,inb(il),inb(il))
     :       -u(il,inb(il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
      fu(il,inb(il))=fu(il,inb(il))-cx
      fu(il,inb(il)-1)=fu(il,inb(il)-1)
     :     +cx*(ph(il,inb(il))-ph(il,inb(il)+1))
     :        /(ph(il,inb(il)-1)-ph(il,inb(il)))

      dx=0.1*ment(il,inb(il),inb(il))*(vent(il,inb(il),inb(il))
     :      -v(il,inb(il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
      fv(il,inb(il))=fv(il,inb(il))-dx
      fv(il,inb(il)-1)=fv(il,inb(il)-1)
     :    +dx*(ph(il,inb(il))-ph(il,inb(il)+1))
     :       /(ph(il,inb(il)-1)-ph(il,inb(il)))
       endif
      ENDIF    !iflag
503   continue

      do j=1,ntra
       do il=1,ncum
        IF (iflag(il) .le. 1) THEN
	IF (cvflag_grav) then
        ex=0.01*grav*ment(il,inb(il),inb(il))
     :      *(traent(il,inb(il),inb(il),j)-tra(il,inb(il),j))
     :      /(ph(i l,inb(il))-ph(il,inb(il)+1))
        ftra(il,inb(il),j)=ftra(il,inb(il),j)-ex
        ftra(il,inb(il)-1,j)=ftra(il,inb(il)-1,j)
     :       +ex*(ph(il,inb(il))-ph(il,inb(il)+1))
     :          /(ph(il,inb(il)-1)-ph(il,inb(il)))
	else
        ex=0.1*ment(il,inb(il),inb(il))
     :      *(traent(il,inb(il),inb(il),j)-tra(il,inb(il),j))
     :      /(ph(i l,inb(il))-ph(il,inb(il)+1))
        ftra(il,inb(il),j)=ftra(il,inb(il),j)-ex
        ftra(il,inb(il)-1,j)=ftra(il,inb(il)-1,j)
     :       +ex*(ph(il,inb(il))-ph(il,inb(il)+1))
     :          /(ph(il,inb(il)-1)-ph(il,inb(il)))
        ENDIF   !cvflag grav
        ENDIF    !iflag
       enddo
      enddo

c
c   ***    homogenize tendencies below cloud base    ***
c
c
      do il=1,ncum
       asum(il)=0.0
       bsum(il)=0.0
       csum(il)=0.0
       dsum(il)=0.0
        esum(il)=0.0
        fsum(il)=0.0
        gsum(il)=0.0
        hsum(il)=0.0
      enddo
c
c     do i=1,nl
c      do il=1,ncum
c         th_wake(il,i)=t_wake(il,i)*(1000.0/p(il,i))**rdcp
c      enddo
c     enddo
c
      do i=1,nl
       do il=1,ncum
        if (i.le.(icb(il)-1) .and. iflag(il) .le. 1) then
cjyg  Saturated part : use T profile
      asum(il)=asum(il)+(ft(il,i)-ftd(il,i))*(ph(il,i)-ph(il,i+1))
      bsum(il)=bsum(il)+(fr(il,i)-fqd(il,i))
     :              *(lv(il,i)+(cl-cpd)*(t(il,i)-t(il,1)))
     :                  *(ph(il,i)-ph(il,i+1))
      csum(il)=csum(il)+(lv(il,i)+(cl-cpd)*(t(il,i)-t(il,1)))
     :                      *(ph(il,i)-ph(il,i+1))
      dsum(il)=dsum(il)+t(il,i)*(ph(il,i)-ph(il,i+1))/th(il,i)
cjyg  Unsaturated part : use T_wake profile
      esum(il)=esum(il)+ftd(il,i)*(ph(il,i)-ph(il,i+1))
      fsum(il)=fsum(il)+fqd(il,i)
     :              *(lv(il,i)+(cl-cpd)*(t_wake(il,i)-t_wake(il,1)))
     :                  *(ph(il,i)-ph(il,i+1))
      gsum(il)=gsum(il)+(lv(il,i)+(cl-cpd)*(t_wake(il,i)-t_wake(il,1)))
     :                      *(ph(il,i)-ph(il,i+1))
      hsum(il)=hsum(il)+t_wake(il,i)
     ;                      *(ph(il,i)-ph(il,i+1))/th_wake(il,i)
        endif
       enddo
      enddo

c!!!      do 700 i=1,icb(il)-1
      do i=1,nl
       do il=1,ncum
        if (i.le.(icb(il)-1) .and. iflag(il) .le. 1) then
         ftd(il,i)=esum(il)*t_wake(il,i)/(th_wake(il,i)*hsum(il))
         fqd(il,i)=fsum(il)/gsum(il)
         ft(il,i)=ftd(il,i)+asum(il)*t(il,i)/(th(il,i)*dsum(il))
         fr(il,i)=fqd(il,i)+bsum(il)/csum(il)
        endif
       enddo
      enddo

c
c   ***   Check that moisture stays positive. If not, scale tendencies
c        in order to ensure moisture positivity
      DO il = 1,ncum
      alpha_qpos(il)=1.
       IF (iflag(il) .le. 1) THEN
        if (fr(il,1) .le. 0.) then
            alpha_qpos(il) = max(alpha_qpos(il) ,
     :	   (-delt*fr(il,1))/
     :     (s_wake(il)*rr_wake(il,1)+(1.-s_wake(il))*rr(il,1)))
        end if
       ENDIF
      ENDDO
      DO i = 2,nl
       DO il = 1,ncum
        IF (iflag(il) .le. 1) THEN
          IF (fr(il,i) .le. 0.) THEN
           alpha_qpos1(il)=max(1. , (-delt*fr(il,i))/
     :     (s_wake(il)*rr_wake(il,i)+(1.-s_wake(il))*rr(il,i)))
             IF (alpha_qpos1(il) .ge. alpha_qpos(il))
     :           alpha_qpos(il)=alpha_qpos1(il)
          ENDIF
        ENDIF
       ENDDO
      ENDDO
      DO il = 1,ncum
       IF (iflag(il) .le. 1 .and. alpha_qpos(il) .gt. 1.001) THEN
        alpha_qpos(il) = alpha_qpos(il)*1.1
       ENDIF
      ENDDO
      DO il = 1,ncum
       IF (iflag(il) .le. 1) THEN
        sigd(il) = sigd(il)/alpha_qpos(il)
        precip(il) = precip(il)/alpha_qpos(il)
       ENDIF
      ENDDO
      DO i = 1,nl
       DO il = 1,ncum
        IF (iflag(il) .le. 1) THEN
         fr(il,i) = fr(il,i)/alpha_qpos(il)
         ft(il,i) = ft(il,i)/alpha_qpos(il)
         fqd(il,i) = fqd(il,i)/alpha_qpos(il)
         ftd(il,i) = ftd(il,i)/alpha_qpos(il)
         fu(il,i) = fu(il,i)/alpha_qpos(il)
         fv(il,i) = fv(il,i)/alpha_qpos(il)
         m(il,i) = m(il,i)/alpha_qpos(il)
         mp(il,i) = mp(il,i)/alpha_qpos(il)
         Vprecip(il,i) = Vprecip(il,i)/alpha_qpos(il)
        ENDIF
       ENDDO
      ENDDO
      DO i = 1,nl
      DO j = 1,nl
       DO il = 1,ncum
        IF (iflag(il) .le. 1) THEN
         ment(il,i,j) = ment(il,i,j)/alpha_qpos(il)
        ENDIF
       ENDDO
      ENDDO
      ENDDO
      DO j = 1,ntra
      DO i = 1,nl
       DO il = 1,ncum
        IF (iflag(il) .le. 1) THEN
         ftra(il,i,j) = ftra(il,i,j)/alpha_qpos(il)
        ENDIF
       ENDDO
      ENDDO
      ENDDO

c
c   ***           reset counter and return           ***
c
      do il=1,ncum
       sig(il,nd)=2.0
      enddo


      do i=1,nd
       do il=1,ncum
        upwd(il,i)=0.0
        dnwd(il,i)=0.0
       enddo
      enddo

      do i=1,nl
       do il=1,ncum
        dnwd0(il,i)=-mp(il,i)
       enddo
      enddo
      do i=nl+1,nd
       do il=1,ncum
        dnwd0(il,i)=0.
       enddo
      enddo


      do i=1,nl
       do il=1,ncum
        if (i.ge.icb(il) .and. i.le.inb(il)) then
          upwd(il,i)=0.0
          dnwd(il,i)=0.0
        endif
       enddo
      enddo

      do i=1,nl
       do k=1,nl
        do il=1,ncum
          up1(il,k,i)=0.0
          dn1(il,k,i)=0.0
        enddo
       enddo
      enddo

      do i=1,nl
       do k=i,nl
        do n=1,i-1
         do il=1,ncum
          if (i.ge.icb(il).and.i.le.inb(il).and.k.le.inb(il)) then
             up1(il,k,i)=up1(il,k,i)+ment(il,n,k)
             dn1(il,k,i)=dn1(il,k,i)-ment(il,k,n)
          endif
         enddo
        enddo
       enddo
      enddo

      do i=1,nl
       do k=1,nl
        do il=1,ncum
         if(i.ge.icb(il)) then
          if(k.ge.i.and. k.le.(inb(il))) then
            upwd(il,i)=upwd(il,i)+m(il,k)
          endif
         else
          if(k.lt.i) then
            upwd(il,i)=upwd(il,i)+cbmf(il)*wghti(il,k)
          endif
        endif
cc        print *,'cbmf',il,i,k,cbmf(il),wghti(il,k)
        end do
       end do
      end do

      do i=2,nl
       do k=i,nl
        do il=1,ncum
ctest         if (i.ge.icb(il).and.i.le.inb(il).and.k.le.inb(il)) then
         if (i.le.inb(il).and.k.le.inb(il)) then
            upwd(il,i)=upwd(il,i)+up1(il,k,i)
            dnwd(il,i)=dnwd(il,i)+dn1(il,k,i)
         endif
cc         print *,'upwd',il,i,k,inb(il),upwd(il,i),m(il,k),up1(il,k,i)
        enddo
       enddo
      enddo


c!!!      DO il=1,ncum
c!!!      do i=icb(il),inb(il)
c!!!
c!!!      upwd(il,i)=0.0
c!!!      dnwd(il,i)=0.0
c!!!      do k=i,inb(il)
c!!!      up1=0.0
c!!!      dn1=0.0
c!!!      do n=1,i-1
c!!!      up1=up1+ment(il,n,k)
c!!!      dn1=dn1-ment(il,k,n)
c!!!      enddo
c!!!      upwd(il,i)=upwd(il,i)+m(il,k)+up1
c!!!      dnwd(il,i)=dnwd(il,i)+dn1
c!!!      enddo
c!!!      enddo
c!!!
c!!!      ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        determination de la variation de flux ascendant entre
c        deux niveau non dilue mip
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,nl
       do il=1,ncum
        mip(il,i)=m(il,i)
       enddo
      enddo

      do i=nl+1,nd
       do il=1,ncum
        mip(il,i)=0.
       enddo
      enddo

      do i=1,nd
       do il=1,ncum
        ma(il,i)=0
       enddo
      enddo

      do i=1,nl
       do j=i,nl
        do il=1,ncum
         ma(il,i)=ma(il,i)+m(il,j)
        enddo
       enddo
      enddo

      do i=nl+1,nd
       do il=1,ncum
        ma(il,i)=0.
       enddo
      enddo

      do i=1,nl
       do il=1,ncum
        if (i.le.(icb(il)-1)) then
         ma(il,i)=0
        endif
       enddo
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        icb represente de niveau ou se trouve la
c        base du nuage , et inb le top du nuage
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,nd
       do il=1,ncum
        mke(il,i)=upwd(il,i)+dnwd(il,i)
       enddo
      enddo

      do i=1,nd
       DO 999 il=1,ncum
        rdcp=(rrd*(1.-rr(il,i))-rr(il,i)*rrv)
     :        /(cpd*(1.-rr(il,i))+rr(il,i)*cpv)
        tls(il,i)=t(il,i)*(1000.0/p(il,i))**rdcp
        tps(il,i)=tp(il,i)
999    CONTINUE
      enddo

c
c   *** diagnose the in-cloud mixing ratio   ***            ! cld
c   ***           of condensed water         ***            ! cld
c                                                           ! cld

       do i=1,nd                                            ! cld
        do il=1,ncum                                        ! cld
         mac(il,i)=0.0                                      ! cld
         wa(il,i)=0.0                                       ! cld
         siga(il,i)=0.0                                     ! cld
         sax(il,i)=0.0                                      ! cld
        enddo                                               ! cld
       enddo                                                ! cld

       do i=minorig, nl                                     ! cld
        do k=i+1,nl+1                                       ! cld
         do il=1,ncum                                       ! cld
          if (i.le.inb(il) .and. k.le.(inb(il)+1)
     $                     .and. iflag(il) .le. 1) then     ! cld
            mac(il,i)=mac(il,i)+m(il,k)                     ! cld
          endif                                             ! cld
         enddo                                              ! cld
        enddo                                               ! cld
       enddo                                                ! cld

       do i=1,nl                                            ! cld
        do j=1,i                                            ! cld
         do il=1,ncum                                       ! cld
          if (i.ge.icb(il) .and. i.le.(inb(il)-1)           ! cld
     :      .and. j.ge.icb(il)
     :      .and. iflag(il) .le. 1 ) then                   ! cld
           sax(il,i)=sax(il,i)+rrd*(tvp(il,j)-tv(il,j))     ! cld
     :        *(ph(il,j)-ph(il,j+1))/p(il,j)                ! cld
          endif                                             ! cld
         enddo                                              ! cld
        enddo                                               ! cld
       enddo                                                ! cld

       do i=1,nl                                            ! cld
        do il=1,ncum                                        ! cld
         if (i.ge.icb(il) .and. i.le.(inb(il)-1)            ! cld
     :       .and. sax(il,i).gt.0.0
     :       .and. iflag(il) .le. 1 ) then                  ! cld
           wa(il,i)=sqrt(2.*sax(il,i))                      ! cld
         endif                                              ! cld
        enddo                                               ! cld
       enddo                                                ! cld

       do i=1,nl                                            ! cld
        do il=1,ncum                                        ! cld
         if (wa(il,i).gt.0.0 .and. iflag(il) .le. 1)        ! cld
     :     siga(il,i)=mac(il,i)/wa(il,i)                    ! cld
     :         *rrd*tvp(il,i)/p(il,i)/100./delta            ! cld
          siga(il,i) = min(siga(il,i),1.0)                  ! cld
cIM cf. FH
         if (iflag_clw.eq.0) then
          qcondc(il,i)=siga(il,i)*clw(il,i)*(1.-ep(il,i))   ! cld
     :           + (1.-siga(il,i))*qcond(il,i)              ! cld
         else if (iflag_clw.eq.1) then
          qcondc(il,i)=qcond(il,i)                          ! cld
         endif

        enddo                                               ! cld
       enddo
c        print*,'cv3_yield fin'        
                                              ! cld
        return
        end


      SUBROUTINE cv3_uncompress(nloc,len,ncum,nd,ntra,idcum
     :         ,iflag
     :         ,precip,sig,w0
     :         ,ft,fq,fu,fv,ftra
     :         ,Ma,upwd,dnwd,dnwd0,qcondc,wd,cape
     :         ,iflag1
     :         ,precip1,sig1,w01
     :         ,ft1,fq1,fu1,fv1,ftra1
     :         ,Ma1,upwd1,dnwd1,dnwd01,qcondc1,wd1,cape1
     :                               )
      implicit none

#include "cv3param.h"

c inputs:
      integer len, ncum, nd, ntra, nloc
      integer idcum(nloc)
      integer iflag(nloc)
      real precip(nloc)
      real sig(nloc,nd), w0(nloc,nd)
      real ft(nloc,nd), fq(nloc,nd), fu(nloc,nd), fv(nloc,nd)
      real ftra(nloc,nd,ntra)
      real Ma(nloc,nd)
      real upwd(nloc,nd),dnwd(nloc,nd),dnwd0(nloc,nd)
      real qcondc(nloc,nd)
      real wd(nloc),cape(nloc)

c outputs:
      integer iflag1(len)
      real precip1(len)
      real sig1(len,nd), w01(len,nd)
      real ft1(len,nd), fq1(len,nd), fu1(len,nd), fv1(len,nd)
      real ftra1(len,nd,ntra)
      real Ma1(len,nd)
      real upwd1(len,nd),dnwd1(len,nd),dnwd01(len,nd)
      real qcondc1(nloc,nd)
      real wd1(nloc),cape1(nloc)

c local variables:
      integer i,k,j

        do 2000 i=1,ncum
         precip1(idcum(i))=precip(i)
         iflag1(idcum(i))=iflag(i)
         wd1(idcum(i))=wd(i)
         cape1(idcum(i))=cape(i)
 2000   continue

        do 2020 k=1,nl
          do 2010 i=1,ncum
            sig1(idcum(i),k)=sig(i,k)
            w01(idcum(i),k)=w0(i,k)
            ft1(idcum(i),k)=ft(i,k)
            fq1(idcum(i),k)=fq(i,k)
            fu1(idcum(i),k)=fu(i,k)
            fv1(idcum(i),k)=fv(i,k)
            Ma1(idcum(i),k)=Ma(i,k)
            upwd1(idcum(i),k)=upwd(i,k)
            dnwd1(idcum(i),k)=dnwd(i,k)
            dnwd01(idcum(i),k)=dnwd0(i,k)
            qcondc1(idcum(i),k)=qcondc(i,k)
 2010     continue
 2020   continue

        do 2200 i=1,ncum
          sig1(idcum(i),nd)=sig(i,nd)
2200    continue


        do 2100 j=1,ntra
c oct3         do 2110 k=1,nl
         do 2110 k=1,nd ! oct3
          do 2120 i=1,ncum
            ftra1(idcum(i),k,j)=ftra(i,k,j)
 2120     continue
 2110    continue
 2100   continue
        return
        end

