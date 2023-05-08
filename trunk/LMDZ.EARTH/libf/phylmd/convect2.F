!
! $Id: convect2.F 1403 2010-07-01 09:02:53Z fairhead $
!
      subroutine convect2(ncum,idcum,len,nd,ndp1,nl,minorig,
     &                 nk1,icb1,
     &                 t1,q1,qs1,u1,v1,gz1,tv1,tp1,tvp1,clw1,h1,
     &                 lv1,cpn1,p1,ph1,ft1,fq1,fu1,fv1,
     &                 tnk1,qnk1,gznk1,plcl1,
     &                 precip1,cbmf1,iflag1,
     &                 delt,cpd,cpv,cl,rv,rd,lv0,g,
     &                 sigs,sigd,elcrit,tlcrit,omtsnow,dtmax,damp,
     &                 alpha,entp,coeffs,coeffr,omtrain,cu,Ma)
C.............................START PROLOGUE............................
C
C SCCS IDENTIFICATION:  @(#)convect2.f	1.2 05/18/00
C                       22:06:22 /h/cm/library/nogaps4/src/sub/fcst/convect2.f_v
C
C CONFIGURATION IDENTIFICATION:  None
C
C MODULE NAME:  convect2
C
C DESCRIPTION:
C
C convect1     The Emanuel Cumulus Convection Scheme - compute tendencies
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
C USAGE: call convect2(ncum,idcum,len,nd,nl,minorig,
C    &                 nk1,icb1,
C    &                 t1,q1,qs1,u1,v1,gz1,tv1,tp1,tvp1,clw1,h1,
C    &                 lv1,cpn1,p1,ph1,ft1,fq1,fu1,fv1,
C    &                 tnk1,qnk1,gznk1,plcl1,
C    &                 precip1,cbmf1,iflag1,
C    &                 delt,cpd,cpv,cl,rv,rd,lv0,g,
C    &                 sigs,sigd,elcrit,tlcrit,omtsnow,dtmax,damp,
C    &                 alpha,entp,coeffs,coeffr,omtrain,cu)
C
C PARAMETERS:
C      Name            Type         Usage            Description
C   ----------      ----------     -------  ----------------------------
C
C      ncum          Integer        Input        number of cumulus points
C      idcum         Integer        Input        index of cumulus point
C      len           Integer        Input        first dimension
C      nd            Integer        Input        total vertical dimension
C      ndp1          Integer        Input        nd + 1
C      nl            Integer        Input        vertical dimension for convection
C      minorig       Integer        Input        First level where convection is allow to begin
C      nk1           Integer        Input        First level of convection
C      ncb1          Integer        Input        Level of free convection
C      t1            Real           Input        temperature
C      q1            Real           Input        specific hum
C      qs1           Real           Input        sat specific hum
C      u1            Real           Input        u-wind
C      v1            Real           Input        v-wind
C      gz1           Real           Inout        geop
C      tv1           Real           Input        virtual temp
C      tp1           Real           Input
C      clw1          Real           Inout        cloud liquid water
C      h1            Real           Inout
C      lv1           Real           Inout
C      cpn1          Real           Inout
C      p1            Real           Input        full level pressure
C      ph1           Real           Input        half level pressure
C      ft1           Real           Output       temp tend
C      fq1           Real           Output       spec hum tend
C      fu1           Real           Output       u-wind tend
C      fv1           Real           Output       v-wind tend
C      precip1       Real           Output       prec
C      cbmf1         Real           In/Out       cumulus mass flux
C      iflag1        Integer        Output       iflag on latitude strip
C      delt          Real           Input        time step
C      cpd           Integer        Input        See description below
C      cpv           Integer        Input        See description below
C      cl            Integer        Input        See description below
C      rv            Integer        Input        See description below
C      rd            Integer        Input        See description below
C      lv0           Integer        Input        See description below
C      g             Integer        Input        See description below
C      sigs          Integer        Input        See description below
C      sigd          Integer        Input        See description below
C      elcrit        Integer        Input        See description below
C      tlcrit        Integer        Input        See description below
C      omtsnow       Integer        Input        See description below
C      dtmax         Integer        Input        See description below
C      damp          Integer        Input        See description below
C      alpha         Integer        Input        See description below
C      ent           Integer        Input        See description below
C      coeffs        Integer        Input        See description below
C      coeffr        Integer        Input        See description below
C      omtrain       Integer        Input        See description below
C      cu            Integer        Input        See description below
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
C         zilch        Zero out an array
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
      integer kmax2,imax2,kmin2,imin2
      real ftmax2,ftmin2
      integer kmax,imax,kmin,imin
      real ftmax,ftmin
c
      integer ncum
      integer idcum(len)
      integer len
      integer nd
      integer ndp1
      integer nl
      integer minorig
      integer nk1(len)
      integer icb1(len)
      real t1(len,nd)
      real q1(len,nd)
      real qs1(len,nd)
      real u1(len,nd)
      real v1(len,nd)
      real gz1(len,nd)
      real tv1(len,nd)
      real tp1(len,nd)
      real tvp1(len,nd)
      real clw1(len,nd)
      real h1(len,nd)
      real lv1(len,nd)
      real cpn1(len,nd)
      real p1(len,nd)
      real ph1(len,ndp1)
      real ft1(len,nd)
      real fq1(len,nd)
      real fu1(len,nd)
      real fv1(len,nd)
      real tnk1(len)
      real qnk1(len)
      real gznk1(len)
      real precip1(len)
      real cbmf1(len)
      real plcl1(len)
      integer iflag1(len)
      real delt
      real cpd
      real cpv
      real cl
      real rv
      real rd
      real lv0
      real g
      real sigs    ! SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE
      real sigd    ! SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT
      real elcrit  ! ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm)
      real tlcrit  ! TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-
c                     CONVERSION THRESHOLD IS ASSUMED TO BE ZERO
      real omtsnow ! OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW
      real dtmax   ! DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION
c                    A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC.
      real damp
      real alpha
      real entp    ! ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT FORMULATION
      real coeffs  ! COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION OF SNOW
      real coeffr  ! COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION OF RAIN
      real omtrain ! OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN
      real cu      ! CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM TRANSPORT
c
      real Ma(len,nd)
c
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
c Local arrays.
c
      real work(ncum)
      real t(ncum,klev)
      real q(ncum,klev)
      real qs(ncum,klev)
      real u(ncum,klev)
      real v(ncum,klev)
      real gz(ncum,klev)
      real h(ncum,klev)
      real lv(ncum,klev)
      real cpn(ncum,klev)
      real p(ncum,klev)
      real ph(ncum,klev)
      real ft(ncum,klev)
      real fq(ncum,klev)
      real fu(ncum,klev)
      real fv(ncum,klev)
      real precip(ncum)
      real cbmf(ncum)
      real plcl(ncum)
      real tnk(ncum)
      real qnk(ncum)
      real gznk(ncum)
      real tv(ncum,klev)
      real tp(ncum,klev)
      real tvp(ncum,klev)
      real clw(ncum,klev)
c     real det(ncum,klev)
      real dph(ncum,klev)
c     real wd(ncum)
c     real tprime(ncum)
c     real qprime(ncum)
      real ah0(ncum)
      real ep(ncum,klev)
      real sigp(ncum,klev)
      integer nent(ncum,klev)
      real water(ncum,klev)
      real evap(ncum,klev)
      real mp(ncum,klev)
      real m(ncum,klev)
      real qti
      real wt(ncum,klev)
      real hp(ncum,klev)
      real lvcp(ncum,klev)
      real elij(ncum,klev,klev)
      real ment(ncum,klev,klev)
      real sij(ncum,klev,klev)
      real qent(ncum,klev,klev)
      real uent(ncum,klev,klev)
      real vent(ncum,klev,klev)
      real qp(ncum,klev)
      real up(ncum,klev)
      real vp(ncum,klev)
      real cape(ncum)
      real capem(ncum)
      real frac(ncum)
      real dtpbl(ncum)
      real tvpplcl(ncum)
      real tvaplcl(ncum)
      real dtmin(ncum)
      real w3d(ncum,klev)
      real am(ncum)
      real ents(ncum)
      real uav(ncum)
      real vav(ncum)
c
      integer iflag(ncum)
      integer nk(ncum)
      integer icb(ncum)
      integer inb(ncum)
      integer inb1(ncum)
      integer jtt(ncum)
c
      integer nn,i,k,n,icbmax,nlp,j
      integer ij
      integer nn2,nn3
      real clmcpv
      real clmcpd
      real cpdmcp
      real cpvmcpd
      real eps
      real epsi
      real epsim1
      real tg,qg,s,alv,tc,ahg,denom,es,rg,ginv,rowl
      real delti
      real tca,elacrit
      real by,defrac
c     real byp
      real byp(ncum)
      logical lcape(ncum)
      real dbo
      real bf2,anum,dei,altem,cwat,stemp
      real alt,qp1,smid,sjmax,sjmin
      real delp,delm
      real awat,coeff,afac,revap,dhdp,fac,qstm,rat
      real qsm,sigt,b6,c6
      real dpinv,cpinv
      real fqold,ftold,fuold,fvold
      real wdtrain(ncum),xxx
      real bsum(ncum,klev)
      real asij(ncum)
      real smin(ncum)
      real scrit(ncum)
c     real amp1,ad
      real amp1(ncum),ad(ncum)
      logical lwork(ncum)
      integer num1,num2
c
c     print*,'cpd en entree de convect2 ',cpd
      nlp=nl+1
c
      rowl=1000.0
      ginv=1.0/g
      delti=1.0/delt
c
c Define some thermodynamic variables.
c
      clmcpv=cl-cpv
      clmcpd=cl-cpd
      cpdmcp=cpd-cpv
      cpvmcpd=cpv-cpd
      eps=rd/rv
      epsi=1.0/eps
      epsim1=epsi-1.0
c
c Compress the fields.
c
      do 110 k=1,nl+1
       nn=0
	do 100 i=1,len
	  if(iflag1(i).eq.0)then
	    nn=nn+1
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
	  endif
 100    continue
c       print*,'100 ncum,nn',ncum,nn
 110  continue
      nn=0
      do 150 i=1,len
	if(iflag1(i).eq.0)then
	  nn=nn+1
	  cbmf(nn)=cbmf1(i)
	  plcl(nn)=plcl1(i)
	  tnk(nn)=tnk1(i)
	  qnk(nn)=qnk1(i)
	  gznk(nn)=gznk1(i)
	  nk(nn)=nk1(i)
	  icb(nn)=icb1(i)
	  iflag(nn)=iflag1(i)
	endif
 150  continue
c       print*,'150 ncum,nn',ncum,nn
c
c Initialize the tendencies, det, wd, tprime, qprime.
c
      do 170 k=1,nl
	do 160 i=1,ncum
c         det(i,k)=0.0
	  ft(i,k)=0.0
	  fu(i,k)=0.0
	  fv(i,k)=0.0
	  fq(i,k)=0.0
	  dph(i,k)=ph(i,k)-ph(i,k+1)
	  ep(i,k)=0.0
	  sigp(i,k)=sigs
 160    continue
 170  continue
      do 180 i=1,ncum
c      wd(i)=0.0
c      tprime(i)=0.0
c      qprime(i)=0.0
       precip(i)=0.0
       ft(i,nl+1)=0.0
       fu(i,nl+1)=0.0
       fv(i,nl+1)=0.0
       fq(i,nl+1)=0.0
 180  continue
c
c Compute icbmax.
c
      icbmax=2
      do 230 i=1,ncum
	icbmax=max(icbmax,icb(i))
 230  continue
c
c
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
     &         +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15))+gznk(i)
 240  continue
c
c
c    ***  Find lifted parcel quantities above cloud base    ***
c
c
	do 300 k=minorig+1,nl
	  do 290 i=1,ncum
	    if(k.ge.(icb(i)+1))then
	      tg=t(i,k)
	      qg=qs(i,k)
	      alv=lv0-clmcpv*(t(i,k)-273.15)
c
c First iteration.
c
	       s=cpd+alv*alv*qg/(rv*t(i,k)*t(i,k))
	       s=1./s
	       ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
	       tg=tg+s*(ah0(i)-ahg)
	       tg=max(tg,35.0)
	       tc=tg-273.15
	       denom=243.5+tc
	       if(tc.ge.0.0)then
		es=6.112*exp(17.67*tc/denom)
	       else
		es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
	       endif
		qg=eps*es/(p(i,k)-es*(1.-eps))
c
c Second iteration.
c
	       s=cpd+alv*alv*qg/(rv*t(i,k)*t(i,k))
	       s=1./s
	       ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
	       tg=tg+s*(ah0(i)-ahg)
	       tg=max(tg,35.0)
	       tc=tg-273.15
	       denom=243.5+tc
	       if(tc.ge.0.0)then
		es=6.112*exp(17.67*tc/denom)
	       else
		es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
	       endif
		qg=eps*es/(p(i,k)-es*(1.-eps))
c
	       alv=lv0-clmcpv*(t(i,k)-273.15)
c      print*,'cpd dans convect2 ',cpd
c      print*,'tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd'
c      print*,tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd
	       tp(i,k)=(ah0(i)-(cl-cpd)*qnk(i)*t(i,k)-gz(i,k)-alv*qg)
     &                  /cpd
c              if (.not.cpd.gt.1000.) then
c                  print*,'CPD=',cpd
c                  stop
c              endif
               clw(i,k)=qnk(i)-qg
               clw(i,k)=max(0.0,clw(i,k))
               rg=qg/(1.-qnk(i))
               tvp(i,k)=tp(i,k)*(1.+rg*epsi)
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
      do 320 k=minorig+1,nl
        do 310 i=1,ncum
          if(k.ge.(nk(i)+1))then
            tca=tp(i,k)-273.15
            if(tca.ge.0.0)then
              elacrit=elcrit
            else
              elacrit=elcrit*(1.0-tca/tlcrit)
            endif
            elacrit=max(elacrit,0.0)
            ep(i,k)=1.0-elacrit/max(clw(i,k),1.0e-8)
            ep(i,k)=max(ep(i,k),0.0 )
            ep(i,k)=min(ep(i,k),1.0 )
            sigp(i,k)=sigs
          endif
 310    continue
 320  continue
c
!=====================================================================
! --- CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL
! --- VIRTUAL TEMPERATURE
!=====================================================================
c
      do 340 k=minorig+1,nl
        do 330 i=1,ncum
        if(k.ge.(icb(i)+1))then
          tvp(i,k)=tvp(i,k)*(1.0-qnk(i)+ep(i,k)*clw(i,k))
c         print*,'i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)'
c         print*, i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)
        endif
 330    continue
 340  continue
      do 350 i=1,ncum
       tvp(i,nlp)=tvp(i,nl)-(gz(i,nlp)-gz(i,nl))/cpd
 350  continue
c
c
c=====================================================================
c --- NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
c=====================================================================
c
        do 360 i=1,ncum*nlp
          nent(i,1)=0
          water(i,1)=0.0
          evap(i,1)=0.0
          mp(i,1)=0.0
          m(i,1)=0.0
          wt(i,1)=omtsnow
          hp(i,1)=h(i,1)
c         if(.not.cpn(i,1).gt.900.) then
c         print*,'i,lv(i,1),cpn(i,1)'
c         print*, i,lv(i,1),cpn(i,1)
c         k=(i-1)/ncum+1
c         print*,'i,k',mod(i,ncum),k,'  cpn',cpn(mod(i,ncum),k)
c         stop
c         endif
          lvcp(i,1)=lv(i,1)/cpn(i,1)
 360    continue
c
      do 380 i=1,ncum*nlp*nlp
        elij(i,1,1)=0.0
        ment(i,1,1)=0.0
        sij(i,1,1)=0.0
 380  continue
c
      do 400 k=1,nlp
       do 390 j=1,nlp
          do 385 i=1,ncum
            qent(i,k,j)=q(i,j)
            uent(i,k,j)=u(i,j)
            vent(i,k,j)=v(i,j)
 385      continue
 390    continue
 400  continue
c
      do 420 i=1,ncum
        qp(i,1)=q(i,1)
        up(i,1)=u(i,1)
        vp(i,1)=v(i,1)
 420  continue
      do 440 k=2,nlp
        do 430 i=1,ncum
          qp(i,k)=q(i,k-1)
          up(i,k)=u(i,k-1)
          vp(i,k)=v(i,k-1)
 430    continue
 440  continue
c
c=====================================================================
c  --- FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S
c  --- HIGHEST LEVEL OF NEUTRAL BUOYANCY
c  --- AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)
c=====================================================================
c
      do 510 i=1,ncum
        cape(i)=0.0
        capem(i)=0.0
        inb(i)=icb(i)+1
        inb1(i)=inb(i)
 510  continue
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
      call zilch(byp,ncum)
      do 515 i=1,ncum
        lcape(i)=.true.
 515  continue
      do 530 k=minorig+1,nl-1
        do 520 i=1,ncum
          if(cape(i).lt.0.0)lcape(i)=.false.
          if((k.ge.(icb(i)+1)).and.lcape(i))then
            by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
            byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
            cape(i)=cape(i)+by
            if(by.ge.0.0)inb1(i)=k+1
            if(cape(i).gt.0.0)then
              inb(i)=k+1
              capem(i)=cape(i)
            endif
          endif
 520    continue
 530  continue
      do 540 i=1,ncum
          cape(i)=capem(i)+byp(i)
          defrac=capem(i)-cape(i)
          defrac=max(defrac,0.001)
          frac(i)=-cape(i)/defrac
          frac(i)=min(frac(i),1.0)
          frac(i)=max(frac(i),0.0)
 540  continue
c
c=====================================================================
c ---   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL
c=====================================================================
c
      do 600 k=minorig+1,nl
        do 590 i=1,ncum
        if((k.ge.icb(i)).and.(k.le.inb(i)))then
          hp(i,k)=h(i,nk(i))+(lv(i,k)+(cpd-cpv)*t(i,k))*ep(i,k)*clw(i,k)
        endif
 590    continue
 600  continue
c
c=====================================================================
c ---  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),
c --- AT EACH MODEL LEVEL
c=====================================================================
c
c tvpplcl = parcel temperature lifted adiabatically from level
c           icb-1 to the LCL.
c tvaplcl = virtual temperature at the LCL.
c
      do 610 i=1,ncum
        dtpbl(i)=0.0
        tvpplcl(i)=tvp(i,icb(i)-1)
     &  -rd*tvp(i,icb(i)-1)*(p(i,icb(i)-1)-plcl(i))
     &  /(cpn(i,icb(i)-1)*p(i,icb(i)-1))
        tvaplcl(i)=tv(i,icb(i))
     &  +(tvp(i,icb(i))-tvp(i,icb(i)+1))*(plcl(i)-p(i,icb(i)))
     &  /(p(i,icb(i))-p(i,icb(i)+1))
 610  continue
c
c-------------------------------------------------------------------
c --- Interpolate difference between lifted parcel and
c --- environmental temperatures to lifted condensation level
c-------------------------------------------------------------------
c
c dtpbl = average of tvp-tv in the PBL (k=nk to icb-1).
c
      do 630 k=minorig,icbmax
        do 620 i=1,ncum
        if((k.ge.nk(i)).and.(k.le.(icb(i)-1)))then
          dtpbl(i)=dtpbl(i)+(tvp(i,k)-tv(i,k))*dph(i,k)
        endif
 620    continue
 630  continue
      do 640 i=1,ncum
        dtpbl(i)=dtpbl(i)/(ph(i,nk(i))-ph(i,icb(i)))
        dtmin(i)=tvpplcl(i)-tvaplcl(i)+dtmax+dtpbl(i)
 640  continue
c
c-------------------------------------------------------------------
c --- Adjust cloud base mass flux
c-------------------------------------------------------------------
c
      do 650 i=1,ncum
       work(i)=cbmf(i)
       cbmf(i)=max(0.0,(1.0-damp)*cbmf(i)+0.1*alpha*dtmin(i))
       if((work(i).eq.0.0).and.(cbmf(i).eq.0.0))then
         iflag(i)=3
       endif
 650  continue
c
c-------------------------------------------------------------------
c --- Calculate rates of mixing,  m(i)
c-------------------------------------------------------------------
c
      call zilch(work,ncum)
c
      do 670 j=minorig+1,nl
        do 660 i=1,ncum
          if((j.ge.(icb(i)+1)).and.(j.le.inb(i)))then
             k=min(j,inb1(i))
             dbo=abs(tv(i,k+1)-tvp(i,k+1)-tv(i,k-1)+tvp(i,k-1))
     &       +entp*0.04*(ph(i,k)-ph(i,k+1))
             work(i)=work(i)+dbo
             m(i,j)=cbmf(i)*dbo
          endif
 660    continue
 670  continue
      do 690 k=minorig+1,nl
        do 680 i=1,ncum
          if((k.ge.(icb(i)+1)).and.(k.le.inb(i)))then
            m(i,k)=m(i,k)/work(i)
          endif
 680    continue
 690  continue
c
c
c=====================================================================
c --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
c --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
c --- FRACTION (sij)
c=====================================================================
c
c
       do 750 i=minorig+1,nl
         do 710 j=minorig+1,nl
           do 700 ij=1,ncum
             if((i.ge.(icb(ij)+1)).and.(j.ge.icb(ij))
     &         .and.(i.le.inb(ij)).and.(j.le.inb(ij)))then
               qti=qnk(ij)-ep(ij,i)*clw(ij,i)
               bf2=1.+lv(ij,j)*lv(ij,j)*qs(ij,j)
     &         /(rv*t(ij,j)*t(ij,j)*cpd)
               anum=h(ij,j)-hp(ij,i)+(cpv-cpd)*t(ij,j)*(qti-q(ij,j))
               denom=h(ij,i)-hp(ij,i)+(cpd-cpv)*(q(ij,i)-qti)*t(ij,j)
               dei=denom
               if(abs(dei).lt.0.01)dei=0.01
               sij(ij,i,j)=anum/dei
               sij(ij,i,i)=1.0
               altem=sij(ij,i,j)*q(ij,i)+(1.-sij(ij,i,j))*qti-qs(ij,j)
               altem=altem/bf2
               cwat=clw(ij,j)*(1.-ep(ij,j))
               stemp=sij(ij,i,j)
               if((stemp.lt.0.0.or.stemp.gt.1.0.or.
     1           altem.gt.cwat).and.j.gt.i)then
                 anum=anum-lv(ij,j)*(qti-qs(ij,j)-cwat*bf2)
                 denom=denom+lv(ij,j)*(q(ij,i)-qti)
                 if(abs(denom).lt.0.01)denom=0.01
                 sij(ij,i,j)=anum/denom
                 altem=sij(ij,i,j)*q(ij,i)+(1.-sij(ij,i,j))*qti-qs(ij,j)
                 altem=altem-(bf2-1.)*cwat
               endif
               if(sij(ij,i,j).gt.0.0.and.sij(ij,i,j).lt.0.9)then
                 qent(ij,i,j)=sij(ij,i,j)*q(ij,i)
     &                        +(1.-sij(ij,i,j))*qti
                 uent(ij,i,j)=sij(ij,i,j)*u(ij,i)
     &                        +(1.-sij(ij,i,j))*u(ij,nk(ij))
                 vent(ij,i,j)=sij(ij,i,j)*v(ij,i)
     &                        +(1.-sij(ij,i,j))*v(ij,nk(ij))
                 elij(ij,i,j)=altem
                 elij(ij,i,j)=max(0.0,elij(ij,i,j))
                 ment(ij,i,j)=m(ij,i)/(1.-sij(ij,i,j))
                 nent(ij,i)=nent(ij,i)+1
               endif
             sij(ij,i,j)=max(0.0,sij(ij,i,j))
             sij(ij,i,j)=min(1.0,sij(ij,i,j))
             endif
  700      continue
  710    continue
c
c   ***   If no air can entrain at level i assume that updraft detrains  ***
c   ***   at that level and calculate detrained air flux and properties  ***
c
           do 740 ij=1,ncum
             if((i.ge.(icb(ij)+1)).and.(i.le.inb(ij))
     &       .and.(nent(ij,i).eq.0))then
               ment(ij,i,i)=m(ij,i)
               qent(ij,i,i)=q(ij,nk(ij))-ep(ij,i)*clw(ij,i)
               uent(ij,i,i)=u(ij,nk(ij))
               vent(ij,i,i)=v(ij,nk(ij))
               elij(ij,i,i)=clw(ij,i)
               sij(ij,i,i)=1.0
             endif
 740       continue
 750   continue
c
      do 770 i=1,ncum
        sij(i,inb(i),inb(i))=1.0
 770  continue
c
c=====================================================================
c   ---  NORMALIZE ENTRAINED AIR MASS FLUXES
c   ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
c=====================================================================
c
c
       call zilch(bsum,ncum*nlp)
       do 780 ij=1,ncum
         lwork(ij)=.false.
 780   continue
       do 789 i=minorig+1,nl
c
         num1=0
         do 779 ij=1,ncum
           if((i.ge.icb(ij)+1).and.(i.le.inb(ij)))num1=num1+1
 779     continue
         if(num1.le.0)go to 789
c
           do 781 ij=1,ncum
             if((i.ge.icb(ij)+1).and.(i.le.inb(ij)))then
                lwork(ij)=(nent(ij,i).ne.0)
                qp1=q(ij,nk(ij))-ep(ij,i)*clw(ij,i)
                anum=h(ij,i)-hp(ij,i)-lv(ij,i)*(qp1-qs(ij,i))
                denom=h(ij,i)-hp(ij,i)+lv(ij,i)*(q(ij,i)-qp1)
                if(abs(denom).lt.0.01)denom=0.01
                scrit(ij)=anum/denom
                alt=qp1-qs(ij,i)+scrit(ij)*(q(ij,i)-qp1)
                if(scrit(ij).lt.0.0.or.alt.lt.0.0)scrit(ij)=1.0
                asij(ij)=0.0
                smin(ij)=1.0
             endif
 781       continue
         do 783 j=minorig,nl
c
         num2=0
         do 778 ij=1,ncum
             if((i.ge.icb(ij)+1).and.(i.le.inb(ij))
     &       .and.(j.ge.icb(ij)).and.(j.le.inb(ij))
     &       .and.lwork(ij))num2=num2+1
 778     continue
         if(num2.le.0)go to 783
c
           do 782 ij=1,ncum
             if((i.ge.icb(ij)+1).and.(i.le.inb(ij))
     &       .and.(j.ge.icb(ij)).and.(j.le.inb(ij)).and.lwork(ij))then
                  if(sij(ij,i,j).gt.0.0.and.sij(ij,i,j).lt.0.9)then
                    if(j.gt.i)then
                      smid=min(sij(ij,i,j),scrit(ij))
                      sjmax=smid
                      sjmin=smid
                        if(smid.lt.smin(ij)
     &                  .and.sij(ij,i,j+1).lt.smid)then
                          smin(ij)=smid
                          sjmax=min(sij(ij,i,j+1),sij(ij,i,j),scrit(ij))
                          sjmin=max(sij(ij,i,j-1),sij(ij,i,j))
                          sjmin=min(sjmin,scrit(ij))
                        endif
                    else
                      sjmax=max(sij(ij,i,j+1),scrit(ij))
                      smid=max(sij(ij,i,j),scrit(ij))
                      sjmin=0.0
                      if(j.gt.1)sjmin=sij(ij,i,j-1)
                      sjmin=max(sjmin,scrit(ij))
                    endif
                    delp=abs(sjmax-smid)
                    delm=abs(sjmin-smid)
                    asij(ij)=asij(ij)+(delp+delm)
     &                           *(ph(ij,j)-ph(ij,j+1))
                    ment(ij,i,j)=ment(ij,i,j)*(delp+delm)
     &                           *(ph(ij,j)-ph(ij,j+1))
                  endif
              endif
  782    continue
  783    continue
            do 784 ij=1,ncum
            if((i.ge.icb(ij)+1).and.(i.le.inb(ij)).and.lwork(ij))then
               asij(ij)=max(1.0e-21,asij(ij))
               asij(ij)=1.0/asij(ij)
               bsum(ij,i)=0.0
            endif
 784        continue
            do 786 j=minorig,nl+1
              do 785 ij=1,ncum
                if((i.ge.icb(ij)+1).and.(i.le.inb(ij))
     &          .and.(j.ge.icb(ij)).and.(j.le.inb(ij))
     &          .and.lwork(ij))then
                   ment(ij,i,j)=ment(ij,i,j)*asij(ij)
                   bsum(ij,i)=bsum(ij,i)+ment(ij,i,j)
                endif
 785     continue
 786     continue
             do 787 ij=1,ncum
               if((i.ge.icb(ij)+1).and.(i.le.inb(ij))
     &         .and.(bsum(ij,i).lt.1.0e-18).and.lwork(ij))then
                 nent(ij,i)=0
                 ment(ij,i,i)=m(ij,i)
                 qent(ij,i,i)=q(ij,nk(ij))-ep(ij,i)*clw(ij,i)
                 uent(ij,i,i)=u(ij,nk(ij))
                 vent(ij,i,i)=v(ij,nk(ij))
                 elij(ij,i,i)=clw(ij,i)
                 sij(ij,i,i)=1.0
               endif
  787        continue
  789  continue
c
c=====================================================================
c --- PRECIPITATING DOWNDRAFT CALCULATION
c=====================================================================
c
c   ***  Check whether ep(inb)=0, if so, skip precipitating    ***
c   ***             downdraft calculation                      ***
c
c
c   ***  Integrate liquid water equation to find condensed water   ***
c   ***                and condensed water flux                    ***
c
c
      do 890 i=1,ncum
        jtt(i)=2
        if(ep(i,inb(i)).le.0.0001)iflag(i)=2
        if(iflag(i).eq.0)then
          lwork(i)=.true.
        else
          lwork(i)=.false.
        endif
 890  continue
c
c    ***                    Begin downdraft loop                    ***
c
c
        call zilch(wdtrain,ncum)
        do 899 i=nl+1,1,-1
c
          num1=0
          do 879 ij=1,ncum
            if((i.le.inb(ij)).and.lwork(ij))num1=num1+1
 879      continue
          if(num1.le.0)go to 899
c
c
c    ***        Calculate detrained precipitation             ***
c
          do 891 ij=1,ncum
            if((i.le.inb(ij)).and.(lwork(ij)))then
            wdtrain(ij)=g*ep(ij,i)*m(ij,i)*clw(ij,i)
            endif
 891      continue
c
          if(i.gt.1)then
            do 893 j=1,i-1
              do 892 ij=1,ncum
                if((i.le.inb(ij)).and.(lwork(ij)))then
                  awat=elij(ij,j,i)-(1.-ep(ij,i))*clw(ij,i)
                  awat=max(0.0,awat)
                  wdtrain(ij)=wdtrain(ij)+g*awat*ment(ij,j,i)
                endif
 892          continue
 893      continue
          endif
c
c    ***    Find rain water and evaporation using provisional   ***
c    ***              estimates of qp(i)and qp(i-1)             ***
c
c
c  ***  Value of terminal velocity and coeffecient of evaporation for snow   ***
c
          do 894 ij=1,ncum
            if((i.le.inb(ij)).and.(lwork(ij)))then
            coeff=coeffs
            wt(ij,i)=omtsnow
c
c  ***  Value of terminal velocity and coeffecient of evaporation for rain   ***
c
            if(t(ij,i).gt.273.0)then
              coeff=coeffr
              wt(ij,i)=omtrain
            endif
            qsm=0.5*(q(ij,i)+qp(ij,i+1))
            afac=coeff*ph(ij,i)*(qs(ij,i)-qsm)
     &       /(1.0e4+2.0e3*ph(ij,i)*qs(ij,i))
            afac=max(afac,0.0)
            sigt=sigp(ij,i)
            sigt=max(0.0,sigt)
            sigt=min(1.0,sigt)
            b6=100.*(ph(ij,i)-ph(ij,i+1))*sigt*afac/wt(ij,i)
            c6=(water(ij,i+1)*wt(ij,i+1)+wdtrain(ij)/sigd)/wt(ij,i)
            revap=0.5*(-b6+sqrt(b6*b6+4.*c6))
            evap(ij,i)=sigt*afac*revap
            water(ij,i)=revap*revap
c
c    ***  Calculate precipitating downdraft mass flux under     ***
c    ***              hydrostatic approximation                 ***
c
            if(i.gt.1)then
              dhdp=(h(ij,i)-h(ij,i-1))/(p(ij,i-1)-p(ij,i))
              dhdp=max(dhdp,10.0)
              mp(ij,i)=100.*ginv*lv(ij,i)*sigd*evap(ij,i)/dhdp
              mp(ij,i)=max(mp(ij,i),0.0)
c
c   ***   Add small amount of inertia to downdraft              ***
c
              fac=20.0/(ph(ij,i-1)-ph(ij,i))
              mp(ij,i)=(fac*mp(ij,i+1)+mp(ij,i))/(1.+fac)
c
c    ***      Force mp to decrease linearly to zero                 ***
c    ***      between about 950 mb and the surface                  ***
c
              if(p(ij,i).gt.(0.949*p(ij,1)))then
                 jtt(ij)=max(jtt(ij),i)
                 mp(ij,i)=mp(ij,jtt(ij))*(p(ij,1)-p(ij,i))
     &           /(p(ij,1)-p(ij,jtt(ij)))
              endif
            endif
c
c    ***       Find mixing ratio of precipitating downdraft     ***
c
            if(i.ne.inb(ij))then
              if(i.eq.1)then
                qstm=qs(ij,1)
              else
                qstm=qs(ij,i-1)
              endif
              if(mp(ij,i).gt.mp(ij,i+1))then
                 rat=mp(ij,i+1)/mp(ij,i)
                 qp(ij,i)=qp(ij,i+1)*rat+q(ij,i)*(1.0-rat)+100.*ginv*
     &             sigd*(ph(ij,i)-ph(ij,i+1))*(evap(ij,i)/mp(ij,i))
                 up(ij,i)=up(ij,i+1)*rat+u(ij,i)*(1.-rat)
                 vp(ij,i)=vp(ij,i+1)*rat+v(ij,i)*(1.-rat)
               else
                 if(mp(ij,i+1).gt.0.0)then
                   qp(ij,i)=(gz(ij,i+1)-gz(ij,i)
     &               +qp(ij,i+1)*(lv(ij,i+1)+t(ij,i+1)
     &               *(cl-cpd))+cpd*(t(ij,i+1)-t(ij,i)))
     &               /(lv(ij,i)+t(ij,i)*(cl-cpd))
                   up(ij,i)=up(ij,i+1)
                   vp(ij,i)=vp(ij,i+1)
                 endif
              endif
              qp(ij,i)=min(qp(ij,i),qstm)
              qp(ij,i)=max(qp(ij,i),0.0)
            endif
            endif
 894      continue
 899    continue
c
c   ***  Calculate surface precipitation in mm/day     ***
c
        do 1190 i=1,ncum
          if(iflag(i).le.1)then
cc            precip(i)=precip(i)+wt(i,1)*sigd*water(i,1)*3600.*24000.
cc     &                /(rowl*g)
cc            precip(i)=precip(i)*delt/86400.
            precip(i) = wt(i,1)*sigd*water(i,1)*86400/g
          endif
 1190   continue
c
c
c   ***  Calculate downdraft velocity scale and surface temperature and  ***
c   ***                    water vapor fluctuations                      ***
c
c     wd=beta*abs(mp(icb))*0.01*rd*t(icb)/(sigd*p(icb))
c     qprime=0.5*(qp(1)-q(1))
c     tprime=lv0*qprime/cpd
c
c   ***  Calculate tendencies of lowest level potential temperature  ***
c   ***                      and mixing ratio                        ***
c
        do 1200 i=1,ncum
          work(i)=0.01/(ph(i,1)-ph(i,2))
          am(i)=0.0
 1200   continue
        do 1220 k=2,nl
          do 1210 i=1,ncum
            if((nk(i).eq.1).and.(k.le.inb(i)).and.(nk(i).eq.1))then
              am(i)=am(i)+m(i,k)
            endif
 1210     continue
 1220   continue
        do 1240 i=1,ncum
          if((g*work(i)*am(i)).ge.delti)iflag(i)=1
          ft(i,1)=ft(i,1)+g*work(i)*am(i)*(t(i,2)-t(i,1)
     &    +(gz(i,2)-gz(i,1))/cpn(i,1))
          ft(i,1)=ft(i,1)-lvcp(i,1)*sigd*evap(i,1)
          ft(i,1)=ft(i,1)+sigd*wt(i,2)*(cl-cpd)*water(i,2)*(t(i,2)
     &     -t(i,1))*work(i)/cpn(i,1)
          fq(i,1)=fq(i,1)+g*mp(i,2)*(qp(i,2)-q(i,1))*
     &    work(i)+sigd*evap(i,1)
          fq(i,1)=fq(i,1)+g*am(i)*(q(i,2)-q(i,1))*work(i)
          fu(i,1)=fu(i,1)+g*work(i)*(mp(i,2)*(up(i,2)-u(i,1))
     &    +am(i)*(u(i,2)-u(i,1)))
          fv(i,1)=fv(i,1)+g*work(i)*(mp(i,2)*(vp(i,2)-v(i,1))
     &    +am(i)*(v(i,2)-v(i,1)))
 1240   continue
        do 1260 j=2,nl
           do 1250 i=1,ncum
             if(j.le.inb(i))then
               fq(i,1)=fq(i,1)
     &                 +g*work(i)*ment(i,j,1)*(qent(i,j,1)-q(i,1))
               fu(i,1)=fu(i,1)
     &                 +g*work(i)*ment(i,j,1)*(uent(i,j,1)-u(i,1))
               fv(i,1)=fv(i,1)
     &                 +g*work(i)*ment(i,j,1)*(vent(i,j,1)-v(i,1))
             endif
 1250      continue
 1260   continue
c
c   ***  Calculate tendencies of potential temperature and mixing ratio  ***
c   ***               at levels above the lowest level                   ***
c
c   ***  First find the net saturated updraft and downdraft mass fluxes  ***
c   ***                      through each level                          ***
c
        do 1500 i=2,nl+1
c
          num1=0
          do 1265 ij=1,ncum
            if(i.le.inb(ij))num1=num1+1
 1265     continue
          if(num1.le.0)go to 1500
c
          call zilch(amp1,ncum)
          call zilch(ad,ncum)
c
          do 1280 k=i+1,nl+1
            do 1270 ij=1,ncum
              if((i.ge.nk(ij)).and.(i.le.inb(ij))
     &            .and.(k.le.(inb(ij)+1)))then
                amp1(ij)=amp1(ij)+m(ij,k)
              endif
 1270         continue
 1280     continue
c
          do 1310 k=1,i
            do 1300 j=i+1,nl+1
               do 1290 ij=1,ncum
                 if((j.le.(inb(ij)+1)).and.(i.le.inb(ij)))then
                   amp1(ij)=amp1(ij)+ment(ij,k,j)
                 endif
 1290          continue
 1300       continue
 1310     continue
          do 1340 k=1,i-1
            do 1330 j=i,nl+1
              do 1320 ij=1,ncum
                if((i.le.inb(ij)).and.(j.le.inb(ij)))then
                   ad(ij)=ad(ij)+ment(ij,j,k)
                endif
 1320         continue
 1330       continue
 1340     continue
c
          do 1350 ij=1,ncum
          if(i.le.inb(ij))then
            dpinv=0.01/(ph(ij,i)-ph(ij,i+1))
            cpinv=1.0/cpn(ij,i)
c
            ft(ij,i)=ft(ij,i)
     &       +g*dpinv*(amp1(ij)*(t(ij,i+1)-t(ij,i)
     &       +(gz(ij,i+1)-gz(ij,i))*cpinv)
     &       -ad(ij)*(t(ij,i)-t(ij,i-1)+(gz(ij,i)-gz(ij,i-1))*cpinv))
     &       -sigd*lvcp(ij,i)*evap(ij,i)
            ft(ij,i)=ft(ij,i)+g*dpinv*ment(ij,i,i)*(hp(ij,i)-h(ij,i)+
     &        t(ij,i)*(cpv-cpd)*(q(ij,i)-qent(ij,i,i)))*cpinv
            ft(ij,i)=ft(ij,i)+sigd*wt(ij,i+1)*(cl-cpd)*water(ij,i+1)*
     &        (t(ij,i+1)-t(ij,i))*dpinv*cpinv
            fq(ij,i)=fq(ij,i)+g*dpinv*(amp1(ij)*(q(ij,i+1)-q(ij,i))-
     &        ad(ij)*(q(ij,i)-q(ij,i-1)))
            fu(ij,i)=fu(ij,i)+g*dpinv*(amp1(ij)*(u(ij,i+1)-u(ij,i))-
     &        ad(ij)*(u(ij,i)-u(ij,i-1)))
            fv(ij,i)=fv(ij,i)+g*dpinv*(amp1(ij)*(v(ij,i+1)-v(ij,i))-
     &        ad(ij)*(v(ij,i)-v(ij,i-1)))
         endif
 1350    continue
         do 1370 k=1,i-1
           do 1360 ij=1,ncum
             if(i.le.inb(ij))then
               awat=elij(ij,k,i)-(1.-ep(ij,i))*clw(ij,i)
               awat=max(awat,0.0)
               fq(ij,i)=fq(ij,i)
     &         +g*dpinv*ment(ij,k,i)*(qent(ij,k,i)-awat-q(ij,i))
               fu(ij,i)=fu(ij,i)
     &         +g*dpinv*ment(ij,k,i)*(uent(ij,k,i)-u(ij,i))
               fv(ij,i)=fv(ij,i)
     &         +g*dpinv*ment(ij,k,i)*(vent(ij,k,i)-v(ij,i))
             endif
 1360      continue
 1370    continue
         do 1390 k=i,nl+1
           do 1380 ij=1,ncum
             if((i.le.inb(ij)).and.(k.le.inb(ij)))then
               fq(ij,i)=fq(ij,i)
     &                  +g*dpinv*ment(ij,k,i)*(qent(ij,k,i)-q(ij,i))
               fu(ij,i)=fu(ij,i)
     &                  +g*dpinv*ment(ij,k,i)*(uent(ij,k,i)-u(ij,i))
               fv(ij,i)=fv(ij,i)
     &                  +g*dpinv*ment(ij,k,i)*(vent(ij,k,i)-v(ij,i))
             endif
 1380      continue
 1390    continue
          do 1400 ij=1,ncum
           if(i.le.inb(ij))then
             fq(ij,i)=fq(ij,i)
     &                +sigd*evap(ij,i)+g*(mp(ij,i+1)*
     &                (qp(ij,i+1)-q(ij,i))
     &                -mp(ij,i)*(qp(ij,i)-q(ij,i-1)))*dpinv
             fu(ij,i)=fu(ij,i)
     &                +g*(mp(ij,i+1)*(up(ij,i+1)-u(ij,i))-mp(ij,i)*
     &                (up(ij,i)-u(ij,i-1)))*dpinv
             fv(ij,i)=fv(ij,i)
     &               +g*(mp(ij,i+1)*(vp(ij,i+1)-v(ij,i))-mp(ij,i)*
     &               (vp(ij,i)-v(ij,i-1)))*dpinv
           endif
 1400     continue
 1500   continue
c
c   *** Adjust tendencies at top of convection layer to reflect  ***
c   ***       actual position of the level zero cape             ***
c
        do 503 ij=1,ncum
        fqold=fq(ij,inb(ij))
        fq(ij,inb(ij))=fq(ij,inb(ij))*(1.-frac(ij))
        fq(ij,inb(ij)-1)=fq(ij,inb(ij)-1)
     &   +frac(ij)*fqold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/
     1   (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))*lv(ij,inb(ij))
     &   /lv(ij,inb(ij)-1)
        ftold=ft(ij,inb(ij))
        ft(ij,inb(ij))=ft(ij,inb(ij))*(1.-frac(ij))
        ft(ij,inb(ij)-1)=ft(ij,inb(ij)-1)
     &   +frac(ij)*ftold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/
     1   (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))*cpn(ij,inb(ij))
     &   /cpn(ij,inb(ij)-1)
        fuold=fu(ij,inb(ij))
        fu(ij,inb(ij))=fu(ij,inb(ij))*(1.-frac(ij))
        fu(ij,inb(ij)-1)=fu(ij,inb(ij)-1)
     &   +frac(ij)*fuold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/
     1   (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))
        fvold=fv(ij,inb(ij))
        fv(ij,inb(ij))=fv(ij,inb(ij))*(1.-frac(ij))
        fv(ij,inb(ij)-1)=fv(ij,inb(ij)-1)
     &  +frac(ij)*fvold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/
     1   (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))
 503    continue
c
c   ***   Very slightly adjust tendencies to force exact   ***
c   ***     enthalpy, momentum and tracer conservation     ***
c
        do 682 ij=1,ncum
        ents(ij)=0.0
        uav(ij)=0.0
        vav(ij)=0.0
        do 681 i=1,inb(ij)
         ents(ij)=ents(ij)
     &  +(cpn(ij,i)*ft(ij,i)+lv(ij,i)*fq(ij,i))*(ph(ij,i)-ph(ij,i+1))	
         uav(ij)=uav(ij)+fu(ij,i)*(ph(ij,i)-ph(ij,i+1))
         vav(ij)=vav(ij)+fv(ij,i)*(ph(ij,i)-ph(ij,i+1))
  681	continue
  682   continue
        do 683 ij=1,ncum
        ents(ij)=ents(ij)/(ph(ij,1)-ph(ij,inb(ij)+1))
        uav(ij)=uav(ij)/(ph(ij,1)-ph(ij,inb(ij)+1))
        vav(ij)=vav(ij)/(ph(ij,1)-ph(ij,inb(ij)+1))
 683    continue
        do 642 ij=1,ncum
        do 641 i=1,inb(ij)
         ft(ij,i)=ft(ij,i)-ents(ij)/cpn(ij,i)
         fu(ij,i)=(1.-cu)*(fu(ij,i)-uav(ij))
         fv(ij,i)=(1.-cu)*(fv(ij,i)-vav(ij))
  641	continue
 642    continue
c
        do 1810 k=1,nl+1
          do 1800 i=1,ncum
            if((q(i,k)+delt*fq(i,k)).lt.0.0)iflag(i)=10
 1800     continue
 1810   continue
c
c
        do 1900 i=1,ncum
          if(iflag(i).gt.2)then
          precip(i)=0.0
          cbmf(i)=0.0
          endif
 1900   continue
        do 1920 k=1,nl
         do 1910 i=1,ncum
           if(iflag(i).gt.2)then
             ft(i,k)=0.0
             fq(i,k)=0.0
             fu(i,k)=0.0
             fv(i,k)=0.0
           endif
 1910    continue
 1920   continue
        do 2000 i=1,ncum
         precip1(idcum(i))=precip(i)
         cbmf1(idcum(i))=cbmf(i)
         iflag1(idcum(i))=iflag(i)
 2000   continue
        do 2020 k=1,nl
          do 2010 i=1,ncum
            ft1(idcum(i),k)=ft(i,k)
            fq1(idcum(i),k)=fq(i,k)
            fu1(idcum(i),k)=fu(i,k)
            fv1(idcum(i),k)=fv(i,k)
 2010     continue
 2020   continue
c
      DO k=1,nd
        DO i=1,len
         Ma(i,k) = 0.
        ENDDO
      ENDDO
      DO k=nl,1,-1
        DO i=1,ncum
          Ma(i,k) = Ma(i,k+1)+m(i,k)
        ENDDO
      ENDDO
c
        return
        end

