!
! $Id: cva_driver.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE cva_driver(len,nd,ndp1,ntra,nloc,
     &                   iflag_con,iflag_mix,
     &                   iflag_clos,delt,
     &                   t1,q1,qs1,t1_wake,q1_wake,qs1_wake,s1_wake,
     &                   u1,v1,tra1,
     &                   p1,ph1,
     &                   ALE1,ALP1,
     &                   sig1feed1,sig2feed1,wght1,
     o                   iflag1,ft1,fq1,fu1,fv1,ftra1,
     &                   precip1,kbas1,ktop1,cbmf1,
     &                   sig1,w01,                  !input/output
     &                   ptop21,sigd,
     &                   Ma1,mip1,Vprecip1,upwd1,dnwd1,dnwd01,
     &                   qcondc1,wd1,
     &                   cape1,cin1,tvp1,
     &                   ftd1,fqd1,
     &                   Plim11,Plim21,asupmax1,supmax01,asupmaxmin1
     &                   ,lalim_conv)
***************************************************************
*                                                             *
* CV_DRIVER                                                   *
*                                                             *
*                                                             *
* written by   : Sandrine Bony-Lena , 17/05/2003, 11.19.41    *
* modified by :                                               *
***************************************************************
***************************************************************
C
      USE dimphy
      implicit none
C
C.............................START PROLOGUE............................
C
C PARAMETERS:
C      Name            Type         Usage            Description
C   ----------      ----------     -------  ----------------------------
C
C      len           Integer        Input        first (i) dimension
C      nd            Integer        Input        vertical (k) dimension
C      ndp1          Integer        Input        nd + 1
C      ntra          Integer        Input        number of tracors
C      iflag_con     Integer        Input        version of convect (3/4)
C      iflag_mix     Integer        Input        version of mixing  (0/1/2)
C      iflag_clos    Integer        Input        version of closure (0/1)
C      delt          Real           Input        time step
C      t1            Real           Input        temperature (sat draught envt)
C      q1            Real           Input        specific hum (sat draught envt)
C      qs1           Real           Input        sat specific hum (sat draught envt)
C      t1_wake       Real           Input        temperature (unsat draught envt)
C      q1_wake       Real           Input        specific hum(unsat draught envt)
C      qs1_wake      Real           Input        sat specific hum(unsat draughts envt)
C      s1_wake       Real           Input        fractionnal area covered by wakes
C      u1            Real           Input        u-wind
C      v1            Real           Input        v-wind
C      tra1          Real           Input        tracors
C      p1            Real           Input        full level pressure
C      ph1           Real           Input        half level pressure
C      ALE1          Real           Input        Available lifting Energy
C      ALP1          Real           Input        Available lifting Power
C      sig1feed1     Real           Input        sigma coord at lower bound of feeding layer
C      sig2feed1     Real           Input        sigma coord at upper bound of feeding layer
C      wght1         Real           Input        weight density determining the feeding mixture
C      iflag1        Integer        Output       flag for Emanuel conditions
C      ft1           Real           Output       temp tend
C      fq1           Real           Output       spec hum tend
C      fu1           Real           Output       u-wind tend
C      fv1           Real           Output       v-wind tend
C      ftra1         Real           Output       tracor tend
C      precip1       Real           Output       precipitation
C      kbas1         Integer        Output       cloud base level
C      ktop1         Integer        Output       cloud top level
C      cbmf1         Real           Output       cloud base mass flux
C      sig1          Real           In/Out       section adiabatic updraft
C      w01           Real           In/Out       vertical velocity within adiab updraft
C      ptop21        Real           In/Out       top of entraining zone
C      Ma1           Real           Output       mass flux adiabatic updraft
C      mip1          Real           Output       mass flux shed by the adiabatic updraft
C      Vprecip1      Real           Output       vertical profile of precipitations
C      upwd1         Real           Output       total upward mass flux (adiab+mixed)
C      dnwd1         Real           Output       saturated downward mass flux (mixed)
C      dnwd01        Real           Output       unsaturated downward mass flux
C      qcondc1       Real           Output       in-cld mixing ratio of condensed water
C      wd1           Real           Output       downdraft velocity scale for sfc fluxes
C      cape1         Real           Output       CAPE
C      cin1          Real           Output       CIN
C      tvp1          Real           Output       adiab lifted parcell virt temp
C      ftd1          Real           Output       precip temp tend
C      fqt1          Real           Output       precip spec hum tend
C      Plim11        Real           Output
C      Plim21        Real           Output
C      asupmax1      Real           Output
C      supmax01      Real           Output
C      asupmaxmin1   Real           Output
C S. Bony, Mar 2002:
C 	* Several modules corresponding to different physical processes
C 	* Several versions of convect may be used:
C  		- iflag_con=3: version lmd  (previously named convect3)
C  		- iflag_con=4: version 4.3b (vect. version, previously convect1/2)
C   + tard: 	- iflag_con=5: version lmd with ice (previously named convectg)
C S. Bony, Oct 2002:
C	* Vectorization of convect3 (ie version lmd)
C
C..............................END PROLOGUE.............................
c
c
#include "dimensions.h"
ccccc#include "dimphy.h"
      include 'iniprint.h'

c
c Input
      integer len
      integer nd
      integer ndp1
      integer ntra
      integer iflag_con
      integer iflag_mix
      integer iflag_clos
      real delt
      real t1(len,nd)
      real q1(len,nd)
      real qs1(len,nd)
      real t1_wake(len,nd)
      real q1_wake(len,nd)
      real qs1_wake(len,nd)
      real s1_wake(len)
      real u1(len,nd)
      real v1(len,nd)
      real tra1(len,nd,ntra)
      real p1(len,nd)
      real ph1(len,ndp1)
      real ALE1(len)
      real ALP1(len)
      real sig1feed1 ! pressure at lower bound of feeding layer
      real sig2feed1 ! pressure at upper bound of feeding layer
      real wght1(nd) ! weight density determining the feeding mixture
c
c Output
      integer iflag1(len)
      real ft1(len,nd)
      real fq1(len,nd)
      real fu1(len,nd)
      real fv1(len,nd)
      real ftra1(len,nd,ntra)
      real precip1(len)
      integer kbas1(len)
      integer ktop1(len)
      real cbmf1(len)
!      real cbmflast(len)
      real sig1(len,klev)      !input/output
      real w01(len,klev)       !input/output
      real ptop21(len)
      real Ma1(len,nd)
      real mip1(len,nd)
!      real Vprecip1(len,nd)
      real Vprecip1(len,nd+1)
      real upwd1(len,nd)
      real dnwd1(len,nd)
      real dnwd01(len,nd)
      real qcondc1(len,nd)     ! cld
      real wd1(len)            ! gust
      real cape1(len)
      real cin1(len)
      real tvp1(len,nd)
c
      real ftd1(len,nd)
      real fqd1(len,nd)
      real Plim11(len)
      real Plim21(len)
      real asupmax1(len,nd)
      real supmax01(len)
      real asupmaxmin1(len)
      integer lalim_conv(len)
!-------------------------------------------------------------------
! --- ARGUMENTS
!-------------------------------------------------------------------
! --- On input:
!
!  t:   Array of absolute temperature (K) of dimension ND, with first
!       index corresponding to lowest model level. Note that this array
!       will be altered by the subroutine if dry convective adjustment
!       occurs and if IPBL is not equal to 0.
!
!  q:   Array of specific humidity (gm/gm) of dimension ND, with first
!       index corresponding to lowest model level. Must be defined
!       at same grid levels as T. Note that this array will be altered
!       if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!  qs:  Array of saturation specific humidity of dimension ND, with first
!       index corresponding to lowest model level. Must be defined
!       at same grid levels as T. Note that this array will be altered
!       if dry convective adjustment occurs and if IPBL is not equal to 0.
!
! t_wake: Array of absolute temperature (K), seen by unsaturated draughts,
!       of dimension ND, with first index corresponding to lowest model level.
!
! q_wake: Array of specific humidity (gm/gm), seen by unsaturated draughts,
!       of dimension ND, with first index corresponding to lowest model level.
!       Must be defined at same grid levels as T.
!
!qs_wake: Array of saturation specific humidity, seen by unsaturated draughts,
!       of dimension ND, with first index corresponding to lowest model level.
!       Must be defined at same grid levels as T.
!
!s_wake: Array of fractionnal area occupied by the wakes.
!
!  u:   Array of zonal wind velocity (m/s) of dimension ND, witth first
!       index corresponding with the lowest model level. Defined at
!       same levels as T. Note that this array will be altered if
!       dry convective adjustment occurs and if IPBL is not equal to 0.
!
!  v:   Same as u but for meridional velocity.
!
!  tra: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
!       where NTRA is the number of different tracers. If no
!       convective tracer transport is needed, define a dummy
!       input array of dimension (ND,1). Tracers are defined at
!       same vertical levels as T. Note that this array will be altered
!       if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!  p:   Array of pressure (mb) of dimension ND, with first
!       index corresponding to lowest model level. Must be defined
!       at same grid levels as T.
!
!  ph:  Array of pressure (mb) of dimension ND+1, with first index
!       corresponding to lowest level. These pressures are defined at
!       levels intermediate between those of P, T, Q and QS. The first
!       value of PH should be greater than (i.e. at a lower level than)
!       the first value of the array P.
!
! ALE:  Available lifting Energy
!
! ALP:  Available lifting Power
!
!  nl:  The maximum number of levels to which convection can penetrate, plus 1.
!       NL MUST be less than or equal to ND-1.
!
!  delt: The model time step (sec) between calls to CONVECT
!
!----------------------------------------------------------------------------
! ---   On Output:
!
!  iflag: An output integer whose value denotes the following:
!       VALUE   INTERPRETATION
!       -----   --------------
!         0     Moist convection occurs.
!         1     Moist convection occurs, but a CFL condition
!               on the subsidence warming is violated. This
!               does not cause the scheme to terminate.
!         2     Moist convection, but no precip because ep(inb) lt 0.0001
!         3     No moist convection because new cbmf is 0 and old cbmf is 0.
!         4     No moist convection; atmosphere is not
!               unstable
!         6     No moist convection because ihmin le minorig.
!         7     No moist convection because unreasonable
!               parcel level temperature or specific humidity.
!         8     No moist convection: lifted condensation
!               level is above the 200 mb level.
!         9     No moist convection: cloud base is higher
!               then the level NL-1.
!
!  ft:   Array of temperature tendency (K/s) of dimension ND, defined at same
!        grid levels as T, Q, QS and P.
!
!  fq:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
!        defined at same grid levels as T, Q, QS and P.
!
!  fu:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
!        defined at same grid levels as T.
!
!  fv:   Same as FU, but for forcing of meridional velocity.
!
!  ftra: Array of forcing of tracer content, in tracer mixing ratio per
!        second, defined at same levels as T. Dimensioned (ND,NTRA).
!
!  precip: Scalar convective precipitation rate (mm/day).
!
!  wd:   A convective downdraft velocity scale. For use in surface
!        flux parameterizations. See convect.ps file for details.
!
!  tprime: A convective downdraft temperature perturbation scale (K).
!          For use in surface flux parameterizations. See convect.ps
!          file for details.
!
!  qprime: A convective downdraft specific humidity
!          perturbation scale (gm/gm).
!          For use in surface flux parameterizations. See convect.ps
!          file for details.
!
!  cbmf: The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
!        BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
!        ITS NEXT CALL. That is, the value of CBMF must be "remembered"
!        by the calling program between calls to CONVECT.
!
!  det:   Array of detrainment mass flux of dimension ND.
!
!  ftd:  Array of temperature tendency due to precipitations (K/s) of dimension ND,
!        defined at same grid levels as T, Q, QS and P.
!
!  fqd:  Array of specific humidity tendencies due to precipitations ((gm/gm)/s)
!        of dimension ND, defined at same grid levels as T, Q, QS and P.
!
!-------------------------------------------------------------------
c
c  Local arrays
c

      integer i,k,n,il,j
      integer nword1,nword2,nword3,nword4
      integer icbmax
      integer nk1(klon)
      integer icb1(klon)
      integer icbs1(klon)

      logical ok_inhib  ! True => possible inhibition of convection by dryness
      logical, save :: debut=.true. 
c$OMP THREADPRIVATE(debut)

      real plcl1(klon)
      real tnk1(klon)
      real thnk1(klon)
      real qnk1(klon)
      real gznk1(klon)
      real pnk1(klon)
      real qsnk1(klon)
      real unk1(klon)
      real vnk1(klon)
      real cpnk1(klon)
      real hnk1(klon)
      real pbase1(klon)
      real buoybase1(klon)

      real lv1(klon,klev) ,lv1_wake(klon,klev)
      real cpn1(klon,klev),cpn1_wake(klon,klev)
      real tv1(klon,klev) ,tv1_wake(klon,klev)
      real gz1(klon,klev) ,gz1_wake(klon,klev)
      real hm1(klon,klev) ,hm1_wake(klon,klev)
      real h1(klon,klev)  ,h1_wake(klon,klev)
      real tp1(klon,klev)
      real clw1(klon,klev)
      real th1(klon,klev) ,th1_wake(klon,klev)
c
      real bid(klon,klev)   ! dummy array
c
      integer ncum
c
      integer j1feed(klon)
      integer j2feed(klon)
      real p1feed1(len) ! pressure at lower bound of feeding layer
      real p2feed1(len) ! pressure at upper bound of feeding layer
      real wghti1(len,nd) ! weights of the feeding layers
c
c (local) compressed fields:
c
      integer nloc
c      parameter (nloc=klon) ! pour l'instant

      integer idcum(nloc)
      integer iflag(nloc),nk(nloc),icb(nloc)
      integer nent(nloc,klev)
      integer icbs(nloc)
      integer inb(nloc), inbis(nloc)

      real cbmf(nloc),plcl(nloc)
      real t(nloc,klev),q(nloc,klev),qs(nloc,klev)
      real t_wake(nloc,klev),q_wake(nloc,klev),qs_wake(nloc,klev)
      real s_wake(nloc)
      real u(nloc,klev),v(nloc,klev)
      real gz(nloc,klev),h(nloc,klev)     ,hm(nloc,klev)
      real               h_wake(nloc,klev),hm_wake(nloc,klev)
      real lv(nloc,klev)     ,cpn(nloc,klev)
      real lv_wake(nloc,klev),cpn_wake(nloc,klev)
      real p(nloc,klev),ph(nloc,klev+1),tv(nloc,klev)    ,tp(nloc,klev)
      real                              tv_wake(nloc,klev)
      real clw(nloc,klev)
      real dph(nloc,klev)
      real pbase(nloc), buoybase(nloc), th(nloc,klev)
      real                              th_wake(nloc,klev)
      real tvp(nloc,klev)
      real sig(nloc,klev), w0(nloc,klev), ptop2(nloc)
      real hp(nloc,klev), ep(nloc,klev), sigp(nloc,klev)
      real frac(nloc), buoy(nloc,klev)
      real cape(nloc)
      real cin(nloc)
      real m(nloc,klev)
      real ment(nloc,klev,klev), sij(nloc,klev,klev)
      real qent(nloc,klev,klev)
      real hent(nloc,klev,klev)
      real uent(nloc,klev,klev), vent(nloc,klev,klev)
      real ments(nloc,klev,klev), qents(nloc,klev,klev)
      real elij(nloc,klev,klev)
      real supmax(nloc,klev)
      real ale(nloc),alp(nloc),coef_clos(nloc)
      real sigd(nloc)
!      real mp(nloc,klev), qp(nloc,klev), up(nloc,klev), vp(nloc,klev)
!      real wt(nloc,klev), water(nloc,klev), evap(nloc,klev)
!      real b(nloc,klev), sigd(nloc)
!      save mp,qp,up,vp,wt,water,evap,b
      real, save, allocatable :: mp(:,:),qp(:,:),up(:,:),vp(:,:)
      real, save, allocatable :: wt(:,:),water(:,:),evap(:,:), b(:,:)
c$OMP THREADPRIVATE(mp,qp,up,vp,wt,water,evap,b)
      real  ft(nloc,klev), fq(nloc,klev)
      real ftd(nloc,klev), fqd(nloc,klev)
      real fu(nloc,klev), fv(nloc,klev)
      real upwd(nloc,klev), dnwd(nloc,klev), dnwd0(nloc,klev)
      real Ma(nloc,klev), mip(nloc,klev), tls(nloc,klev)
      real tps(nloc,klev), qprime(nloc), tprime(nloc)
      real precip(nloc)
!      real Vprecip(nloc,klev)
      real Vprecip(nloc,klev+1)
      real tra(nloc,klev,ntra), trap(nloc,klev,ntra)
      real ftra(nloc,klev,ntra), traent(nloc,klev,klev,ntra)
      real qcondc(nloc,klev)  ! cld
      real wd(nloc)           ! gust
      real Plim1(nloc),Plim2(nloc)
      real asupmax(nloc,klev)
      real supmax0(nloc)
      real asupmaxmin(nloc)
c
      real tnk(nloc),qnk(nloc),gznk(nloc)
      real wghti(nloc,nd)
      real hnk(nloc),unk(nloc),vnk(nloc)
      logical, save :: first=.true.
c$OMP THREADPRIVATE(first)
      CHARACTER (LEN=20) :: modname='cva_driver'
      CHARACTER (LEN=80) :: abort_message

c
!      print *, 't1, t1_wake ',(k,t1(1,k),t1_wake(1,k),k=1,klev)
!      print *, 'q1, q1_wake ',(k,q1(1,k),q1_wake(1,k),k=1,klev)

!-------------------------------------------------------------------
! --- SET CONSTANTS AND PARAMETERS
!-------------------------------------------------------------------

       if (first) then
         allocate(mp(nloc,klev), qp(nloc,klev), up(nloc,klev))
         allocate(vp(nloc,klev), wt(nloc,klev), water(nloc,klev))
         allocate(evap(nloc,klev), b(nloc,klev))
         first=.false.
       endif
c -- set simulation flags:
c   (common cvflag)

       CALL cv_flag

c -- set thermodynamical constants:
c 	(common cvthermo)

       CALL cv_thermo(iflag_con)

c -- set convect parameters
c
c 	includes microphysical parameters and parameters that
c  	control the rate of approach to quasi-equilibrium)
c 	(common cvparam)

      if (iflag_con.eq.3) then
       CALL cv3_param(nd,delt)
  
      endif

      if (iflag_con.eq.4) then
       CALL cv_param(nd)
      endif

!---------------------------------------------------------------------
! --- INITIALIZE OUTPUT ARRAYS AND PARAMETERS
!---------------------------------------------------------------------
      nword1=len
      nword2=len*nd
      nword3=len*nd*ntra
      nword4=len*nd*nd
 
!      call izilch(iflag1  ,nword1)
!      call  zilch(iflag1  ,nword1)
      do i=1,len
         iflag1(i)=0
         ktop1(i)=0
         kbas1(i)=0
      enddo
      call  zilch(ft1     ,nword2)
      call  zilch(fq1     ,nword2)
      call  zilch(fu1     ,nword2)
      call  zilch(fv1     ,nword2)
      call  zilch(ftra1   ,nword3)
      call  zilch(precip1 ,nword1)
!      call izilch(kbas1   ,nword1)
!      call  zilch(kbas1   ,nword1)
!      call izilch(ktop1   ,nword1)
!      call  zilch(ktop1   ,nword1)
      call  zilch(cbmf1   ,nword1)
      call  zilch(ptop21  ,nword1)
      call  zilch(Ma1     ,nword2)
      call  zilch(mip1    ,nword2)
!      call  zilch(Vprecip1,nword2)
      Vprecip1=0.
      call  zilch(upwd1   ,nword2)
      call  zilch(dnwd1   ,nword2)
      call  zilch(dnwd01  ,nword2)
      call  zilch(qcondc1 ,nword2)
!test
!      call  zilch(qcondc ,nword2)
      call  zilch(wd1     ,nword1)
      call  zilch(cape1   ,nword1)
      call  zilch(cin1    ,nword1)
      call  zilch(tvp1    ,nword2)
      call  zilch(ftd1    ,nword2)
      call  zilch(fqd1    ,nword2)
      call  zilch(Plim11  ,nword1)
      call  zilch(Plim21  ,nword1)
      call  zilch(asupmax1,nword2)
      call  zilch(supmax01,nword1)
      call  zilch(asupmaxmin1,nword1)
c
      DO il = 1,len
       cin1(il) = -100000.
       cape1(il) = -1.
      ENDDO
c  
      if (iflag_con.eq.3) then
        do il=1,len
         sig1(il,nd)=sig1(il,nd)+1.
         sig1(il,nd)=amin1(sig1(il,nd),12.1)
        enddo
      endif
  
!---------------------------------------------------------------------
! --- INITIALIZE LOCAL ARRAYS AND PARAMETERS
!---------------------------------------------------------------------
!
      do il = 1,nloc
         coef_clos(il)=1.
      enddo

!--------------------------------------------------------------------
! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
!--------------------------------------------------------------------

      if (iflag_con.eq.3) then
  
       if (debut) THEN 
        print*,'Emanuel version 3 nouvelle'
       endif 
!       print*,'t1, q1 ',t1,q1
       CALL cv3_prelim(len,nd,ndp1,t1,q1,p1,ph1      ! nd->na
     o               ,lv1,cpn1,tv1,gz1,h1,hm1,th1)
    
c
       CALL cv3_prelim(len,nd,ndp1,t1_wake,q1_wake,p1,ph1 ! nd->na
     o               ,lv1_wake,cpn1_wake,tv1_wake,gz1_wake
     o               ,h1_wake,bid,th1_wake)
    
      endif
c
      if (iflag_con.eq.4) then
       print*,'Emanuel version 4 '
       CALL cv_prelim(len,nd,ndp1,t1,q1,p1,ph1
     o               ,lv1,cpn1,tv1,gz1,h1,hm1)
      endif

!--------------------------------------------------------------------
! --- CONVECTIVE FEED
!--------------------------------------------------------------------
!
! compute feeding layer potential temperature and mixing ratio :
!
! get bounds of feeding layer
!
c test niveaux couche alimentation KE
       if(sig1feed1.eq.sig2feed1) then
         write(lunout,*)'impossible de choisir sig1feed=sig2feed'
         write(lunout,*)'changer la valeur de sig2feed dans physiq.def'
         abort_message = ''
         CALL abort_gcm (modname,abort_message,1)
       endif
c
       do i=1,len
         p1feed1(i)=sig1feed1*ph1(i,1)
         p2feed1(i)=sig2feed1*ph1(i,1)
ctest maf
c         p1feed1(i)=ph1(i,1)
c         p2feed1(i)=ph1(i,2)
c         p2feed1(i)=ph1(i,3)
ctestCR: on prend la couche alim des thermiques
c          p2feed1(i)=ph1(i,lalim_conv(i)+1)
c          print*,'lentr=',lentr(i),ph1(i,lentr(i)+1),ph1(i,2)
       end do
!
       if (iflag_con.eq.3) then
       endif
      do i=1,len
!      print*,'avant cv3_feed plim',p1feed1(i),p2feed1(i)
      enddo
      if (iflag_con.eq.3) then
 
c     print*, 'IFLAG1 avant cv3_feed'
c     print*,'len,nd',len,nd
c     write(*,'(64i1)') iflag1(2:klon-1)

       CALL cv3_feed(len,nd,t1,q1,u1,v1,p1,ph1,hm1,gz1           ! nd->na
     i         ,p1feed1,p2feed1,wght1
     o         ,wghti1,tnk1,thnk1,qnk1,qsnk1,unk1,vnk1
     o         ,cpnk1,hnk1,nk1,icb1,icbmax,iflag1,gznk1,plcl1)
      endif
    
c     print*, 'IFLAG1 apres cv3_feed'
c     print*,'len,nd',len,nd
c     write(*,'(64i1)') iflag1(2:klon-1)

      if (iflag_con.eq.4) then
       CALL cv_feed(len,nd,t1,q1,qs1,p1,hm1,gz1
     o         ,nk1,icb1,icbmax,iflag1,tnk1,qnk1,gznk1,plcl1)
      endif
c
!      print *, 'cv3_feed-> iflag1, plcl1 ',iflag1(1),plcl1(1)
c
!--------------------------------------------------------------------
! --- UNDILUTE (ADIABATIC) UPDRAFT / 1st part
! (up through ICB for convect4, up through ICB+1 for convect3)
!     Calculates the lifted parcel virtual temperature at nk, the
!     actual temperature, and the adiabatic liquid water content.
!--------------------------------------------------------------------

      if (iflag_con.eq.3) then
    
       CALL cv3_undilute1(len,nd,t1,qs1,gz1,plcl1,p1,icb1,tnk1,qnk1  ! nd->na
     o                    ,gznk1,tp1,tvp1,clw1,icbs1)
      endif
   

      if (iflag_con.eq.4) then
       CALL cv_undilute1(len,nd,t1,q1,qs1,gz1,p1,nk1,icb1,icbmax
     :                        ,tp1,tvp1,clw1)
      endif
c
!-------------------------------------------------------------------
! --- TRIGGERING
!-------------------------------------------------------------------
c
!      print *,' avant triggering, iflag_con ',iflag_con
c
      if (iflag_con.eq.3) then
    
       CALL cv3_trigger(len,nd,icb1,plcl1,p1,th1,tv1,tvp1,thnk1      ! nd->na
     o                 ,pbase1,buoybase1,iflag1,sig1,w01)
    

c     print*, 'IFLAG1 apres cv3_triger'
c     print*,'len,nd',len,nd
c     write(*,'(64i1)') iflag1(2:klon-1)

c     call dump2d(iim,jjm-1,sig1(2)
      endif

      if (iflag_con.eq.4) then
       CALL cv_trigger(len,nd,icb1,cbmf1,tv1,tvp1,iflag1)
      endif
c
c
!=====================================================================
! --- IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY
!=====================================================================

      ncum=0
      do 400 i=1,len
        if(iflag1(i).eq.0)then
           ncum=ncum+1
           idcum(ncum)=i
        endif
 400  continue
c
!       print*,'klon, ncum = ',len,ncum
c
      IF (ncum.gt.0) THEN

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! --- COMPRESS THE FIELDS
!		(-> vectorization over convective gridpoints)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      if (iflag_con.eq.3) then
!       print*,'ncum tv1 ',ncum,tv1 
!       print*,'tvp1 ',tvp1
       CALL cv3a_compress( len,nloc,ncum,nd,ntra
     :    ,iflag1,nk1,icb1,icbs1
     :    ,plcl1,tnk1,qnk1,gznk1,hnk1,unk1,vnk1
     :    ,wghti1,pbase1,buoybase1
     :    ,t1,q1,qs1,t1_wake,q1_wake,qs1_wake,s1_wake
     :    ,u1,v1,gz1,th1,th1_wake
     :    ,tra1
     :    ,h1     ,lv1     ,cpn1   ,p1,ph1,tv1    ,tp1,tvp1,clw1
     :    ,h1_wake,lv1_wake,cpn1_wake     ,tv1_wake
     :    ,sig1,w01,ptop21
     :    ,Ale1,Alp1
     o    ,iflag,nk,icb,icbs
     o    ,plcl,tnk,qnk,gznk,hnk,unk,vnk
     o    ,wghti,pbase,buoybase
     o    ,t,q,qs,t_wake,q_wake,qs_wake,s_wake
     o    ,u,v,gz,th,th_wake
     o    ,tra
     o    ,h     ,lv     ,cpn    ,p,ph,tv    ,tp,tvp,clw
     o    ,h_wake,lv_wake,cpn_wake    ,tv_wake
     o    ,sig,w0,ptop2
     o    ,Ale,Alp  ) 

!       print*,'tv ',tv
!       print*,'tvp ',tvp

      endif

      if (iflag_con.eq.4) then
       CALL cv_compress( len,nloc,ncum,nd
     :    ,iflag1,nk1,icb1
     :    ,cbmf1,plcl1,tnk1,qnk1,gznk1
     :    ,t1,q1,qs1,u1,v1,gz1
     :    ,h1,lv1,cpn1,p1,ph1,tv1,tp1,tvp1,clw1
     o    ,iflag,nk,icb
     o    ,cbmf,plcl,tnk,qnk,gznk
     o    ,t,q,qs,u,v,gz,h,lv,cpn,p,ph,tv,tp,tvp,clw
     o    ,dph )
      endif

!-------------------------------------------------------------------
! --- UNDILUTE (ADIABATIC) UPDRAFT / second part :
! ---   FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
! ---   &
! ---   COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
! ---   FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
! ---   &
! ---   FIND THE LEVEL OF NEUTRAL BUOYANCY
!-------------------------------------------------------------------

      if (iflag_con.eq.3) then
       CALL cv3_undilute2(nloc,ncum,nd,icb,icbs,nk        !na->nd
     :                        ,tnk,qnk,gznk,hnk,t,q,qs,gz
     :                        ,p,h,tv,lv,pbase,buoybase,plcl
     o                        ,inb,tp,tvp,clw,hp,ep,sigp,buoy)
   
      endif

      if (iflag_con.eq.4) then
       CALL cv_undilute2(nloc,ncum,nd,icb,nk
     :                        ,tnk,qnk,gznk,t,q,qs,gz
     :                        ,p,dph,h,tv,lv
     o             ,inb,inbis,tp,tvp,clw,hp,ep,sigp,frac)
      endif
c
!-------------------------------------------------------------------
! --- MIXING(1)   (if iflag_mix .ge. 1)
!-------------------------------------------------------------------
      IF (iflag_con .eq. 3) THEN
       IF (iflag_mix .ge. 1 ) THEN
         CALL zilch(supmax,nloc*klev)    
         CALL cv3p_mixing(nloc,ncum,nd,nd,ntra,icb,nk,inb    ! na->nd
     :                       ,ph,t,q,qs,u,v,tra,h,lv,qnk
     :                       ,unk,vnk,hp,tv,tvp,ep,clw,sig
     :                    ,ment,qent,hent,uent,vent,nent
     :                   ,sij,elij,supmax,ments,qents,traent)
!        print*, 'cv3p_mixing-> supmax ', (supmax(1,k), k=1,nd)
      
       ELSE
        CALL zilch(supmax,nloc*klev)
       ENDIF
      ENDIF
!-------------------------------------------------------------------
! --- CLOSURE
!-------------------------------------------------------------------

c
      if (iflag_con.eq.3) then
       IF (iflag_clos .eq. 0) THEN
        CALL cv3_closure(nloc,ncum,nd,icb,inb              ! na->nd
     :                       ,pbase,p,ph,tv,buoy
     o                       ,sig,w0,cape,m,iflag)
       ENDIF
c
       ok_inhib = iflag_mix .EQ. 2
c
       IF (iflag_clos .eq. 1) THEN
        print *,' pas d appel cv3p_closure'
cc        CALL cv3p_closure(nloc,ncum,nd,icb,inb              ! na->nd
cc    :                       ,pbase,plcl,p,ph,tv,tvp,buoy
cc    :                       ,supmax
cc    o                       ,sig,w0,ptop2,cape,cin,m)
       ENDIF
       IF (iflag_clos .eq. 2) THEN
        CALL cv3p1_closure(nloc,ncum,nd,icb,inb              ! na->nd
     :                       ,pbase,plcl,p,ph,tv,tvp,buoy
     :                       ,supmax,ok_inhib,Ale,Alp
     o                       ,sig,w0,ptop2,cape,cin,m,iflag,coef_clos
     :                       ,Plim1,Plim2,asupmax,supmax0
     :                       ,asupmaxmin,cbmf)
       ENDIF
      endif   ! iflag_con.eq.3
  
      if (iflag_con.eq.4) then
       CALL cv_closure(nloc,ncum,nd,nk,icb
     :                ,tv,tvp,p,ph,dph,plcl,cpn
     o                ,iflag,cbmf)
      endif
c
!      print *,'cv_closure-> cape ',cape(1)
c
!-------------------------------------------------------------------
! --- MIXING(2)
!-------------------------------------------------------------------

      if (iflag_con.eq.3) then
        IF (iflag_mix.eq.0) THEN
         CALL cv3_mixing(nloc,ncum,nd,nd,ntra,icb,nk,inb    ! na->nd
     :                       ,ph,t,q,qs,u,v,tra,h,lv,qnk
     :                       ,unk,vnk,hp,tv,tvp,ep,clw,m,sig
     o   ,ment,qent,uent,vent,nent,sij,elij,ments,qents,traent)
         CALL zilch(hent,nloc*klev*klev)
        ELSE
         CALL cv3_mixscale(nloc,ncum,nd,ment,m)
         if (debut) THEN 
          print *,' cv3_mixscale-> '
         endif !(debut) THEN 
        ENDIF
      endif

      if (iflag_con.eq.4) then
       CALL cv_mixing(nloc,ncum,nd,icb,nk,inb,inbis
     :                     ,ph,t,q,qs,u,v,h,lv,qnk
     :                     ,hp,tv,tvp,ep,clw,cbmf
     o                     ,m,ment,qent,uent,vent,nent,sij,elij)
      endif
c
      if (debut) THEN 
       print *,' cv_mixing ->'
      endif !(debut) THEN 
c      do i = 1,klev
c        print*,'cv_mixing-> i,ment ',i,(ment(1,i,j),j=1,klev)
c      enddo
c
!-------------------------------------------------------------------
! --- UNSATURATED (PRECIPITATING) DOWNDRAFTS
!-------------------------------------------------------------------
      if (iflag_con.eq.3) then
       if (debut) THEN 
        print *,' cva_driver -> cv3_unsat '
       endif !(debut) THEN 
    
       CALL cv3_unsat(nloc,ncum,nd,nd,ntra,icb,inb,iflag    ! na->nd
     :               ,t_wake,q_wake,qs_wake,gz,u,v,tra,p,ph
     :               ,th_wake,tv_wake,lv_wake,cpn_wake
     :               ,ep,sigp,clw
     :               ,m,ment,elij,delt,plcl,coef_clos
     o          ,mp,qp,up,vp,trap,wt,water,evap,b,sigd)
      endif
     
      if (iflag_con.eq.4) then
       CALL cv_unsat(nloc,ncum,nd,inb,t,q,qs,gz,u,v,p,ph
     :                   ,h,lv,ep,sigp,clw,m,ment,elij
     o                   ,iflag,mp,qp,up,vp,wt,water,evap)
      endif
c
      if (debut) THEN 
       print *,'cv_unsat-> '
       debut=.FALSE.
      endif !(debut) THEN
!
c      print *,'cv_unsat-> mp ',mp
c      print *,'cv_unsat-> water ',water
!-------------------------------------------------------------------
! --- YIELD
!     (tendencies, precipitation, variables of interface with other
!      processes, etc)
!-------------------------------------------------------------------

      if (iflag_con.eq.3) then
  
       CALL cv3_yield(nloc,ncum,nd,nd,ntra            ! na->nd
     :                     ,icb,inb,delt
     :                     ,t,q,t_wake,q_wake,s_wake,u,v,tra
     :                     ,gz,p,ph,h,hp,lv,cpn,th,th_wake
     :                     ,ep,clw,m,tp,mp,qp,up,vp,trap
     :                     ,wt,water,evap,b,sigd
     :                    ,ment,qent,hent,iflag_mix,uent,vent
     :                    ,nent,elij,traent,sig
     :                    ,tv,tvp,wghti
     :                    ,iflag,precip,Vprecip,ft,fq,fu,fv,ftra
     :                    ,cbmf,upwd,dnwd,dnwd0,ma,mip
     :                    ,tls,tps,qcondc,wd
     :                    ,ftd,fqd)
!      print *,' cv3_yield -> fqd(1) = ',fqd(1,1)
      endif

      if (iflag_con.eq.4) then
       CALL cv_yield(nloc,ncum,nd,nk,icb,inb,delt
     :              ,t,q,t_wake,q_wake,u,v,tra
     :              ,gz,p,ph,h,hp,lv,cpn,th
     :              ,ep,clw,frac,m,mp,qp,up,vp
     :              ,wt,water,evap
     :              ,ment,qent,uent,vent,nent,elij
     :              ,tv,tvp
     o              ,iflag,wd,qprime,tprime
     o              ,precip,cbmf,ft,fq,fu,fv,Ma,qcondc)
      endif

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! --- UNCOMPRESS THE FIELDS
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


      if (iflag_con.eq.3) then
       CALL cv3a_uncompress(nloc,len,ncum,nd,ntra,idcum
     :          ,iflag,icb,inb
     :          ,precip,cbmf,sig,w0,ptop2
     :          ,ft,fq,fu,fv,ftra
     :          ,Ma,mip,Vprecip,upwd,dnwd,dnwd0
     ;          ,qcondc,wd,cape,cin
     :          ,tvp
     :          ,ftd,fqd
     :          ,Plim1,Plim2,asupmax,supmax0
     :          ,asupmaxmin
     o          ,iflag1,kbas1,ktop1
     o          ,precip1,cbmf1,sig1,w01,ptop21
     o          ,ft1,fq1,fu1,fv1,ftra1
     o          ,Ma1,mip1,Vprecip1,upwd1,dnwd1,dnwd01
     o          ,qcondc1,wd1,cape1,cin1
     o          ,tvp1
     o          ,ftd1,fqd1
     o          ,Plim11,Plim21,asupmax1,supmax01
     o          ,asupmaxmin1)
      endif

      if (iflag_con.eq.4) then
       CALL cv_uncompress(nloc,len,ncum,nd,idcum
     :          ,iflag
     :          ,precip,cbmf
     :          ,ft,fq,fu,fv
     :          ,Ma,qcondc
     o          ,iflag1
     o          ,precip1,cbmf1
     o          ,ft1,fq1,fu1,fv1
     o          ,Ma1,qcondc1 )
      endif

      ENDIF ! ncum>0

9999  continue
      return
      end

