!
! $Header$
!
      subroutine convect1(len,nd,ndp1,noff,minorig,
     &                   t,q,qs,u,v,
     &                   p,ph,iflag,ft,
     &                   fq,fu,fv,precip,cbmf,delt,Ma)
C.............................START PROLOGUE............................
C
C SCCS IDENTIFICATION:  @(#)convect1.f	1.1 04/21/00
C                       19:40:52 /h/cm/library/nogaps4/src/sub/fcst/convect1.f_v
C
C CONFIGURATION IDENTIFICATION:  None
C
C MODULE NAME:  convect1
C
C DESCRIPTION:
C
C convect1     The Emanuel Cumulus Convection Scheme
C
C CONTRACT NUMBER AND TITLE:  None
C
C REFERENCES: Programmers  K. Emanuel (MIT), Timothy F. Hogan, M. Peng (NRL)
C
C CLASSIFICATION:  Unclassified
C
C RESTRICTIONS: None
C
C COMPILER DEPENDENCIES: FORTRAN 77, FORTRAN 90
C
C COMPILE OPTIONS: Fortran 77: -Zu -Wf"-ei -o aggress"
C                  Fortran 90: -O vector3,scalar3,task1,aggress,overindex  -ei -r 2
C
C LIBRARIES OF RESIDENCE: /a/ops/lib/libfcst159.a
C
C USAGE: call convect1(len,nd,noff,minorig,
C    &                   t,q,qs,u,v,
C    &                   p,ph,iflag,ft,
C    &                   fq,fu,fv,precip,cbmf,delt)
C
C PARAMETERS:
C      Name            Type         Usage            Description
C   ----------      ----------     -------  ----------------------------
C
C      len           Integer        Input        first (i) dimension
C      nd            Integer        Input        vertical (k) dimension
C      ndp1          Integer        Input        nd + 1
C      noff          Integer        Input        integer limit for convection (nd-noff)
C      minorig       Integer        Input        First level of convection
C      t             Real           Input        temperature
C      q             Real           Input        specific hum
C      qs            Real           Input        sat specific hum
C      u             Real           Input        u-wind
C      v             Real           Input        v-wind
C      p             Real           Input        full level pressure
C      ph            Real           Input        half level pressure
C      iflag         Integer        Output       iflag on latitude strip
C      ft            Real           Output       temp tend
C      fq            Real           Output       spec hum tend
C      fu            Real           Output       u-wind tend
C      fv            Real           Output       v-wind tend
C      cbmf          Real           In/Out       cumulus mass flux
C      delt          Real           Input        time step
C      iflag         Integer        Output       integer flag for Emanuel conditions
C
C COMMON BLOCKS:
C      Block      Name     Type    Usage              Notes
C     --------  --------   ----    ------   ------------------------
C
C FILES: None
C
C DATA BASES: None
C
C NON-FILE INPUT/OUTPUT: None
C
C ERROR CONDITIONS: None
C
C ADDITIONAL COMMENTS: None
C
C.................MAINTENANCE SECTION................................
C
C MODULES CALLED:
C         Name           Description
C         convect2        Emanuel cumulus convection tendency calculations
C        -------     ----------------------
C LOCAL VARIABLES AND
C          STRUCTURES:
C Name     Type    Description
C -------  ------  -----------
C See Comments Below
C
C i        Integer loop index
C k        Integer loop index
c
C METHOD:
C
C See Emanuel, K. and M. Zivkovic-Rothman, 2000: Development and evaluation of a
C       convective scheme for use in climate models.
C
C FILES: None
C
C INCLUDE FILES: None
C
C MAKEFILE: /a/ops/met/nogaps/src/sub/fcst/fcst159lib.mak
C
C..............................END PROLOGUE.............................
c
c
       USE dimphy
      implicit none
c
cym#include "dimensions.h"
cym#include "dimphy.h"
c
      integer len
      integer nd
      integer ndp1
      integer noff
      real t(len,nd)
      real q(len,nd)
      real qs(len,nd)
      real u(len,nd)
      real v(len,nd)
      real p(len,nd)
      real ph(len,ndp1)
      integer iflag(len)
      real ft(len,nd)
      real fq(len,nd)
      real fu(len,nd)
      real fv(len,nd)
      real precip(len)
      real cbmf(len)
      real Ma(len,nd)
      integer minorig
      real delt,cpd,cpv,cl,rv,rd,lv0,g
      real sigs,sigd,elcrit,tlcrit,omtsnow,dtmax,damp
      real alpha,entp,coeffs,coeffr,omtrain,cu
c
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
!-------------------------------------------------------------------
c
c  Local arrays
c
      integer nl
      integer nlp
      integer nlm
      integer i,k,n
      real delti
      real rowl
      real clmcpv
      real clmcpd
      real cpdmcp
      real cpvmcpd
      real eps
      real epsi
      real epsim1
      real ginv
      real hrd
      real prccon1
      integer icbmax
      real lv(klon,klev)
      real cpn(klon,klev)
      real cpx(klon,klev)
      real tv(klon,klev)
      real gz(klon,klev)
      real hm(klon,klev)
      real h(klon,klev)
      real work(klon)
      integer ihmin(klon)
      integer nk(klon)
      real rh(klon)
      real chi(klon)
      real plcl(klon)
      integer icb(klon)
      real tnk(klon)
      real qnk(klon)
      real gznk(klon)
      real pnk(klon)
      real qsnk(klon)
      real ticb(klon)
      real gzicb(klon)
      real tp(klon,klev)
      real tvp(klon,klev)
      real clw(klon,klev)
c
      real ah0(klon),cpp(klon)
      real tg,qg,s,alv,tc,ahg,denom,es,rg
c
      integer ncum
      integer idcum(klon)
c
      cpd=1005.7
      cpv=1870.0
      cl=4190.0
      rv=461.5
      rd=287.04
      lv0=2.501E6
      g=9.8
C
C   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
C   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
C   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
C   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
C   ***               BETWEEN 0 C AND TLCRIT)                        ***
C   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
C   ***                       FORMULATION                            ***
C   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
C   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
C   ***                        OF CLOUD                              ***
C   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
C   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
C   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
C   ***                          OF RAIN                             ***
C   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
C   ***                          OF SNOW                             ***
C   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
C   ***                         TRANSPORT                            ***
C   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
C   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
C   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
C   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
C   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
C   ***                   (DAMP MUST BE LESS THAN 1)                 ***
c
      sigs=0.12
      sigd=0.05
      elcrit=0.0011
      tlcrit=-55.0
      omtsnow=5.5
      dtmax=0.9
      damp=0.1
      alpha=0.2
      entp=1.5
      coeffs=0.8
      coeffr=1.0
      omtrain=50.0
c
      cu=0.70
      damp=0.1
c
c
c Define nl, nlp, nlm, and delti
c
      nl=nd-noff
      nlp=nl+1
      nlm=nl-1
      delti=1.0/delt
!
!-------------------------------------------------------------------
! --- SET CONSTANTS
!-------------------------------------------------------------------
!
      rowl=1000.0
      clmcpv=cl-cpv
      clmcpd=cl-cpd
      cpdmcp=cpd-cpv
      cpvmcpd=cpv-cpd
      eps=rd/rv
      epsi=1.0/eps
      epsim1=epsi-1.0
      ginv=1.0/g
      hrd=0.5*rd
      prccon1=86400.0*1000.0/(rowl*g)
!
! dtmax is the maximum negative temperature perturbation.
!
!=====================================================================
! --- INITIALIZE OUTPUT ARRAYS AND PARAMETERS
!=====================================================================
!
      do 20 k=1,nd
        do 10 i=1,len
         ft(i,k)=0.0
         fq(i,k)=0.0
         fu(i,k)=0.0
         fv(i,k)=0.0
         tvp(i,k)=0.0
         tp(i,k)=0.0
         clw(i,k)=0.0
         gz(i,k) = 0.
 10     continue
 20   continue
      do 60 i=1,len
        precip(i)=0.0
        iflag(i)=0
 60   continue
c
!=====================================================================
! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
!=====================================================================
      do 110 k=1,nl+1
        do 100 i=1,len
          lv(i,k)= lv0-clmcpv*(t(i,k)-273.15)
          cpn(i,k)=cpd*(1.0-q(i,k))+cpv*q(i,k)
          cpx(i,k)=cpd*(1.0-q(i,k))+cl*q(i,k)
          tv(i,k)=t(i,k)*(1.0+q(i,k)*epsim1)
 100    continue
 110  continue
c
c gz = phi at the full levels (same as p).
c
      do 120 i=1,len
        gz(i,1)=0.0
 120  continue
      do 140 k=2,nlp
        do 130 i=1,len
          gz(i,k)=gz(i,k-1)+hrd*(tv(i,k-1)+tv(i,k))
     &         *(p(i,k-1)-p(i,k))/ph(i,k)
 130    continue
 140  continue
c
c h  = phi + cpT (dry static energy).
c hm = phi + cp(T-Tbase)+Lq
c
      do 170 k=1,nlp
        do 160 i=1,len
          h(i,k)=gz(i,k)+cpn(i,k)*t(i,k)
          hm(i,k)=gz(i,k)+cpx(i,k)*(t(i,k)-t(i,1))+lv(i,k)*q(i,k)
 160    continue
 170  continue
c
!-------------------------------------------------------------------
! --- Find level of minimum moist static energy
! --- If level of minimum moist static energy coincides with
! --- or is lower than minimum allowable parcel origin level,
! --- set iflag to 6.
!-------------------------------------------------------------------
      do 180 i=1,len
       work(i)=1.0e12
       ihmin(i)=nl
 180  continue
      do 200 k=2,nlp
        do 190 i=1,len
         if((hm(i,k).lt.work(i)).and.
     &      (hm(i,k).lt.hm(i,k-1)))then
           work(i)=hm(i,k)
           ihmin(i)=k
         endif
 190    continue
 200  continue
      do 210 i=1,len
        ihmin(i)=min(ihmin(i),nlm)
        if(ihmin(i).le.minorig)then
          iflag(i)=6
        endif
 210  continue
c
!-------------------------------------------------------------------
! --- Find that model level below the level of minimum moist static
! --- energy that has the maximum value of moist static energy
!-------------------------------------------------------------------
 
      do 220 i=1,len
       work(i)=hm(i,minorig)
       nk(i)=minorig
 220  continue
      do 240 k=minorig+1,nl
        do 230 i=1,len
         if((hm(i,k).gt.work(i)).and.(k.le.ihmin(i)))then
           work(i)=hm(i,k)
           nk(i)=k
         endif
 230     continue
 240  continue
!-------------------------------------------------------------------
! --- Check whether parcel level temperature and specific humidity
! --- are reasonable
!-------------------------------------------------------------------
       do 250 i=1,len
       if(((t(i,nk(i)).lt.250.0).or.
     &      (q(i,nk(i)).le.0.0).or.
     &      (p(i,ihmin(i)).lt.400.0)).and.
     &      (iflag(i).eq.0))iflag(i)=7
 250   continue
!-------------------------------------------------------------------
! --- Calculate lifted condensation level of air at parcel origin level
! --- (Within 0.2% of formula of Bolton, MON. WEA. REV.,1980)
!-------------------------------------------------------------------
       do 260 i=1,len
        tnk(i)=t(i,nk(i))
        qnk(i)=q(i,nk(i))
        gznk(i)=gz(i,nk(i))
        pnk(i)=p(i,nk(i))
        qsnk(i)=qs(i,nk(i))
c
        rh(i)=qnk(i)/qsnk(i)
        rh(i)=min(1.0,rh(i))
        chi(i)=tnk(i)/(1669.0-122.0*rh(i)-tnk(i))
        plcl(i)=pnk(i)*(rh(i)**chi(i))
        if(((plcl(i).lt.200.0).or.(plcl(i).ge.2000.0))
     &   .and.(iflag(i).eq.0))iflag(i)=8
 260   continue
!-------------------------------------------------------------------
! --- Calculate first level above lcl (=icb)
!-------------------------------------------------------------------
      do 270 i=1,len
       icb(i)=nlm
 270  continue
c
      do 290 k=minorig,nl
        do 280 i=1,len
          if((k.ge.(nk(i)+1)).and.(p(i,k).lt.plcl(i)))
     &    icb(i)=min(icb(i),k)
 280    continue
 290  continue
c
      do 300 i=1,len
        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
 300  continue
c
c Compute icbmax.
c
      icbmax=2
      do 310 i=1,len
        icbmax=max(icbmax,icb(i))
 310  continue
!
!-------------------------------------------------------------------
! --- Calculates the lifted parcel virtual temperature at nk,
! --- the actual temperature, and the adiabatic
! --- liquid water content. The procedure is to solve the equation.
!     cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.
!-------------------------------------------------------------------
!
      do 320 i=1,len
        tnk(i)=t(i,nk(i))
        qnk(i)=q(i,nk(i))
        gznk(i)=gz(i,nk(i))
        ticb(i)=t(i,icb(i))
        gzicb(i)=gz(i,icb(i))
 320  continue
c
c   ***  Calculate certain parcel quantities, including static energy   ***
c
      do 330 i=1,len
        ah0(i)=(cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i)
     &         +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15))+gznk(i)
        cpp(i)=cpd*(1.-qnk(i))+qnk(i)*cpv
 330  continue
c
c   ***   Calculate lifted parcel quantities below cloud base   ***
c
        do 350 k=minorig,icbmax-1
          do 340 i=1,len
           tp(i,k)=tnk(i)-(gz(i,k)-gznk(i))/cpp(i)
           tvp(i,k)=tp(i,k)*(1.+qnk(i)*epsi)
  340     continue
  350   continue
c
c    ***  Find lifted parcel quantities above cloud base    ***
c
        do 360 i=1,len
         tg=ticb(i)
         qg=qs(i,icb(i))
         alv=lv0-clmcpv*(ticb(i)-273.15)
c
c First iteration.
c
          s=cpd+alv*alv*qg/(rv*ticb(i)*ticb(i))
          s=1./s
          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          tg=tg+s*(ah0(i)-ahg)
          tg=max(tg,35.0)
          tc=tg-273.15
          denom=243.5+tc
          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
          else
           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
          endif
          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
c
c Second iteration.
c
          s=cpd+alv*alv*qg/(rv*ticb(i)*ticb(i))
          s=1./s
          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          tg=tg+s*(ah0(i)-ahg)
          tg=max(tg,35.0)
          tc=tg-273.15
          denom=243.5+tc
          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
          else
           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
          end if
          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
c
         alv=lv0-clmcpv*(ticb(i)-273.15)
         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
     &   -gz(i,icb(i))-alv*qg)/cpd
         clw(i,icb(i))=qnk(i)-qg
         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
         rg=qg/(1.-qnk(i))
         tvp(i,icb(i))=tp(i,icb(i))*(1.+rg*epsi)
  360   continue
c
      do 380 k=minorig,icbmax
       do 370 i=1,len
         tvp(i,k)=tvp(i,k)-tp(i,k)*qnk(i)
 370   continue
 380  continue
c
!-------------------------------------------------------------------
! --- Test for instability.
! --- If there was no convection at last time step and parcel
! --- is stable at icb, then set iflag to 4.
!-------------------------------------------------------------------
 
      do 390 i=1,len
        if((cbmf(i).eq.0.0) .and.(iflag(i).eq.0).and.
     &  (tvp(i,icb(i)).le.(tv(i,icb(i))-dtmax)))iflag(i)=4
 390  continue
 
!=====================================================================
! --- IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY
!=====================================================================
c
      ncum=0
      do 400 i=1,len
        if(iflag(i).eq.0)then
           ncum=ncum+1
           idcum(ncum)=i
        endif
 400  continue
c
c Call convect2, which compresses the points and computes the heating,
c moistening, velocity mixing, and precipiation.
c
c     print*,'cpd avant convect2 ',cpd
      if(ncum.gt.0)then
      call convect2(ncum,idcum,len,nd,ndp1,nl,minorig,
     &              nk,icb,
     &              t,q,qs,u,v,gz,tv,tp,tvp,clw,h,
     &              lv,cpn,p,ph,ft,fq,fu,fv,
     &              tnk,qnk,gznk,plcl,
     &              precip,cbmf,iflag,
     &              delt,cpd,cpv,cl,rv,rd,lv0,g,
     &              sigs,sigd,elcrit,tlcrit,omtsnow,dtmax,damp,
     &              alpha,entp,coeffs,coeffr,omtrain,cu,Ma)
      endif
c
      return
      end
