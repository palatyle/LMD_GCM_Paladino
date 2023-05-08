      SUBROUTINE orosetup
     *         ( nlon   , nlev  , ktest
     *         , kkcrit, kkcrith, kcrit, ksect , kkhlim
     *         , kkenvh, kknu  , kknu2, kkbreak
     *         , paphm1, papm1 , pum1 , pvm1, ptm1, pgeom1, pstab, pstd
     *         , prho  , pri   , ptau, pvph, ppsi, pzdep
     *         , pulow , pvlow  
     *         , ptheta, pgam, pmea, ppic, pval
     *         , pnu  ,  pd1  ,  pd2  ,pdmod  )
C
c**** *gwsetup*
c
c     purpose.
c     --------
c     SET-UP THE ESSENTIAL PARAMETERS OF THE SSO DRAG SCHEME:
C     DEPTH OF LOW WBLOCKED LAYER, LOW-LEVEL FLOW, BACKGROUND
C     STRATIFICATION.....
c
c**   interface.
c     ----------
c          from *orodrag*
c
c        explicit arguments :
c        --------------------
c     ==== inputs ===
c 
c nlon----input-I-Total number of horizontal points that get into physics
c nlev----input-I-Number of vertical levels
c ktest--input-I: Flags to indicate active points
c
c ptsphy--input-R-Time-step (s)
c paphm1--input-R: pressure at model 1/2 layer
c papm1---input-R: pressure at model layer
c pgeom1--input-R: Altitude of layer above ground
c VENUS ATTENTION: CP VARIABLE PSTAB CALCULE EN AMONT DES PARAMETRISATIONS
c pstab-----R-: Brunt-Vaisala freq.^2 at 1/2 layers (input except klev+1)
c ptm1, pum1, pvm1--R-: t, u and v
c pmea----input-R-Mean Orography (m)
C pstd----input-R-SSO standard deviation (m)
c psig----input-R-SSO slope
c pgam----input-R-SSO Anisotropy
c pthe----input-R-SSO Angle
c ppic----input-R-SSO Peacks elevation (m)
c pval----input-R-SSO Valleys elevation (m)

c     ==== outputs ===
c pulow, pvlow -output-R: Low-level wind
c kkcrit----I-: Security value for top of low level flow
c kcrit-----I-: Critical level 
c ksect-----I-: Not used
c kkhlim----I-: Not used
c kkenvh----I-: Top of blocked flow layer
c kknu------I-: Layer that sees mountain peacks
c kknu2-----I-: Layer that sees mountain peacks above mountain mean
c kknub-----I-: Layer that sees mountain mean above valleys
c prho------R-: Density at 1/2 layers
c pri-------R-: Background Richardson Number, Wind shear measured along GW stress
c pvph------R-: Wind in  plan of GW stress, Half levels.
c ppsi------R-: Angle between low level wind and SS0 main axis.
c pd1-------R-| Compared the ratio of the stress
c pd2-------R-| that is along the wind to that Normal to it.
c               pdi define the plane of low level stress
c               compared to the low level wind.
c see p. 108 Lott & Miller (1997).                      
c pdmod-----R-: Norme of pdi

c     === local arrays ===
c
c zvpf------R-: Wind projected in the plan of the low-level stress.

c     ==== outputs ===
c
c        implicit arguments :   none
c        --------------------
c
c     method.
c     -------
c
c
c     externals.
c     ----------
c
c
c     reference.
c     ----------
c
c        see ecmwf research department documentation of the "i.f.s."
c
c     author.
c     -------
c
c     modifications.
c     --------------
c     f.lott  for the new-gwdrag scheme november 1993
c
c-----------------------------------------------------------------------
      use dimphy
      implicit none

#include "YOMCST.h"
#include "YOEGWD.h"

c-----------------------------------------------------------------------
c
c*       0.1   arguments
c              ---------
c
      integer nlon,nlev
      integer kkcrit(nlon),kkcrith(nlon),kcrit(nlon),ksect(nlon),
     *        kkhlim(nlon),ktest(nlon),kkenvh(nlon)
      integer kkbreak(nlon)
c
      real paphm1(nlon,klev+1),papm1(nlon,klev),pum1(nlon,klev),
     *     pvm1(nlon,klev),ptm1(nlon,klev),pgeom1(nlon,klev),
     *     prho(nlon,klev+1),pri(nlon,klev+1),pstab(nlon,klev+1),
     *     ptau(nlon,klev+1),pvph(nlon,klev+1),ppsi(nlon,klev+1),
     *     pzdep(nlon,klev)
       real pulow(nlon),pvlow(nlon),ptheta(nlon),pgam(nlon),pnu(nlon),
     *     pd1(nlon),pd2(nlon),pdmod(nlon)
      real pstd(nlon),pmea(nlon),ppic(nlon),pval(nlon)
c
c-----------------------------------------------------------------------
c
c*       0.2   local arrays
c              ------------
c
c
      integer ilevh ,jl,jk,iii
      real zcons1,zhgeo,zu,zphi
      real zvt1,zvt2,zdwind,zwind,zdelp
      real zstabm,zstabp,zrhom,zrhop
      logical lo 
      logical ll1(klon,klev+1)
      integer kknu(klon),kknu2(klon),kknub(klon),kknul(klon),
     *        kentp(klon),ncount(klon)  
c
      real zhcrit(klon,klev),zvpf(klon,klev),
     *     zdp(klon,klev)
      real znorm(klon),zb(klon),zc(klon),
     *      zulow(klon),zvlow(klon),znup(klon),znum(klon)

c     ------------------------------------------------------------------
c
c*         1.    initialization
c                --------------
c
c       PRINT *,' in orosetup'
 100  continue
c
c     ------------------------------------------------------------------
c
c*         1.1   computational constants
c                -----------------------
c
 110  continue
c
      ilevh =klev/3
c
      zcons1=1./rd
c
c     ------------------------------------------------------------------
c
c*         2.
c                --------------
c
 200  continue
c
c     ------------------------------------------------------------------
c
c*         2.1     define low level wind, project winds in plane of
c*                 low level wind, determine sector in which to take
c*                 the variance and set indicator for critical levels.
c
c
c
      do 2001 jl=kidia,kfdia
      kknu(jl)    =klev
      kknu2(jl)   =klev
      kknub(jl)   =klev
      kknul(jl)   =klev
      kkbreak(jl) =klev + 1
      pgam(jl) =max(pgam(jl),gtsec)
      ll1(jl,klev+1)=.false.
 2001 continue
c
c Ajouter une initialisation (L. Li, le 23fev99):
c
      do jk=klev,ilevh,-1
      do jl=kidia,kfdia
      ll1(jl,jk)= .false.
      ENDDO
      ENDDO

c      VENUS: define break for subcritical stress
c      ----------------------------
      do jk=klev,ilevh,-1
      do jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      !zhgeo=pgeom1(jl,jk)/rg
      !!zhcrit(jl,jk)=ppic(jl)
      !zhcrit(jl,jk)=amin1(ppic(jl)-pmea(jl),pmea(jl)-pval(jl))
      !ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
      !if(ll1(jl,jk) .neqv. ll1(jl,jk+1)) then
      !  kkbreak(jl)=jk
      !endif

      !if (paphm1(jl,jk) .ge. 7.e6) kkbreak(jl)=jk
      !kkbreak(jl) = klev - 2 ! gwd1103
      !kkbreak(jl) = klev - 4 ! gwd1104
      !kkbreak(jl) = klev - 3 ! gwd1105

      endif     
      enddo
      enddo

c
c*      define top of low level flow
c       ----------------------------
      do 2002 jk=klev,ilevh,-1
      do 2003 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      lo=(paphm1(jl,jk)/paphm1(jl,klev+1)).ge.gsigcr
      if(lo) then
        kkcrit(jl)=jk
      endif
      zhcrit(jl,jk)=ppic(jl)-pval(jl)           
      zhgeo=pgeom1(jl,jk)/rg
      ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
C     if(ll1(jl,jk).xor.ll1(jl,jk+1)) then
      if(ll1(jl,jk).neqv.ll1(jl,jk+1)) then
        kknu(jl)=jk
      endif
      if(.not.ll1(jl,ilevh))kknu(jl)=ilevh
      endif
 2003 continue
 2002 continue
      do 2004 jk=klev,ilevh,-1
      do 2005 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
!      zhcrit(jl,jk)=ppic(jl)-pmea(jl)
      zhgeo=pgeom1(jl,jk)/rg
!      ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
      ll1(jl,jk)=(zhgeo.gt.0.5*pstd(jl)) ! TN : do not consider outlier peaks
!      ll1(jl,jk)=(zhgeo.gt.pstd(jl)) ! TN : do not consider outlier peaks
!      ll1(jl,jk)=(zhgeo.gt.2*pstd(jl)) ! TN : do not consider outlier peaks
      if(ll1(jl,jk) .neqv. ll1(jl,jk+1)) then
        kknu2(jl)=jk
      endif
      if(.not.ll1(jl,ilevh))kknu2(jl)=ilevh
      endif
 2005 continue
 2004 continue

!      do 2104 jk=klev,ilevh,-1
!      do 2105 jl=kidia,kfdia
!      if(ktest(jl).eq.1) then
!      zhgeo=pgeom1(jl,jk)/rg
!      ll1(jl,jk)=(zhgeo.gt.pstd(jl)) ! TN : do not consider outlier peaks
!      if(ll1(jl,jk) .neqv. ll1(jl,jk+1)) then
!        kknul(jl)=jk
!      endif
!      if(.not.ll1(jl,ilevh))kknul(jl)=ilevh
!      endif
! 2105 continue
! 2104 continue


      do 2006 jk=klev,ilevh,-1
      do 2007 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      zhcrit(jl,jk)=amin1(ppic(jl)-pmea(jl),pmea(jl)-pval(jl))
      zhgeo=pgeom1(jl,jk)/rg
      ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
c     if(ll1(jl,jk).xor.ll1(jl,jk+1)) then
      if(ll1(jl,jk).neqv.ll1(jl,jk+1)) then
        kknub(jl)=jk
      endif
      if(.not.ll1(jl,ilevh))kknub(jl)=ilevh
      endif
 2007 continue
 2006 continue
c
      do 2010 jl=kidia,kfdia  
      if(ktest(jl).eq.1) then
      kknu(jl)=min(kknu(jl),nktopg)
      kknu2(jl)=min(kknu2(jl),nktopg)
      kknub(jl)=min(kknub(jl),nktopg)
!      kknul(jl)=klev
      endif
 2010 continue      
c
 210  continue
c
cc*     initialize various arrays
c
      do 2107 jl=kidia,kfdia
      prho(jl,klev+1)  =0.0
cym correction en attendant mieux
      prho(jl,1)  =0.0      
      pstab(jl,klev+1) =0.0
      pstab(jl,1)      =0.0
      pri(jl,klev+1)   =9999.0
      ppsi(jl,klev+1)  =0.0
      pri(jl,1)        =0.0
      pvph(jl,1)       =0.0
      pvph(jl,klev+1)  =0.0
cym correction en attendant mieux
cym      pvph(jl,klev)    =0.0
      pulow(jl)        =0.0
      pvlow(jl)        =0.0
      zulow(jl)        =0.0
      zvlow(jl)        =0.0
      kkcrith(jl)      =klev
      kkenvh(jl)       =klev
      kentp(jl)        =klev
      kcrit(jl)        =1
      ncount(jl)       =0
      ll1(jl,klev+1)   =.false.
 2107 continue
c
c*     define flow density and stratification (rho and N2)
c      at semi layers.
c      -------------------------------------------------------
c
      do 223 jk=klev,2,-1
      do 222 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
        zdp(jl,jk)=papm1(jl,jk)-papm1(jl,jk-1)
        prho(jl,jk)=2.*paphm1(jl,jk)*zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
      endif
  222 continue
  223 continue
c      print*,"altitude(m)=",pgeom1(kfdia/2,:)
c      print*,"pression(Pa)=",papm1(kfdia/2,:)
c
c********************************************************************
c
c*     define Low level flow (between ground and peacks-valleys)
c      ---------------------------------------------------------
      do 2115 jk=klev,ilevh,-1
      do 2116 jl=kidia,kfdia
      if(ktest(jl).eq.1)  then
      if(jk.ge.kknu2(jl).and.jk.le.kknul(jl)) then
        pulow(jl)=pulow(jl)+pum1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
        pvlow(jl)=pvlow(jl)+pvm1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
        pstab(jl,klev+1)=pstab(jl,klev+1)
     c                   +pstab(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
        prho(jl,klev+1)=prho(jl,klev+1)
     c                   +prho(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
      end if
      endif
 2116 continue
 2115 continue
      do 2110 jl=kidia,kfdia
      if(ktest(jl).eq.1)  then
      pulow(jl)=pulow(jl)/(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
      pvlow(jl)=pvlow(jl)/(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
      znorm(jl)=max(sqrt(pulow(jl)**2+pvlow(jl)**2),gvsec)
      pvph(jl,klev+1)=znorm(jl)
      pstab(jl,klev+1)=pstab(jl,klev+1)
     c                /(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
      prho(jl,klev+1)=prho(jl,klev+1)
     c                /(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
      endif
 2110 continue

c
c*******  setup orography orientation relative to the low level
C       wind and define parameters of the Anisotropic wave stress.
c
      do 2112 jl=kidia,kfdia
      if(ktest(jl).eq.1)  then
        lo=(pulow(jl).lt.gvsec).and.(pulow(jl).ge.-gvsec)
        if(lo) then
          zu=pulow(jl)+2.*gvsec
        else
          zu=pulow(jl)
        endif
        zphi=atan(pvlow(jl)/zu)
        ppsi(jl,klev+1)=ptheta(jl)*rpi/180.-zphi
        zb(jl)=1.-0.18*pgam(jl)-0.04*pgam(jl)**2
        zc(jl)=0.48*pgam(jl)+0.3*pgam(jl)**2
        pd1(jl)=zb(jl)-(zb(jl)-zc(jl))*(sin(ppsi(jl,klev+1))**2)
        pd2(jl)=(zb(jl)-zc(jl))*sin(ppsi(jl,klev+1))
     *                         *cos(ppsi(jl,klev+1))
        pdmod(jl)=sqrt(pd1(jl)**2+pd2(jl)**2)
      endif
 2112 continue
c
c  ************ projet flow in plane of lowlevel stress *************
C  ************ Find critical levels...                 *************
c
      do 213 jk=1,klev
      do 212 jl=kidia,kfdia
      if(ktest(jl).eq.1)  then
        zvt1       =pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk)
        zvt2       =-pvlow(jl)*pum1(jl,jk)+pulow(jl)*pvm1(jl,jk)
        zvpf(jl,jk)=(zvt1*pd1(jl)+zvt2*pd2(jl))/(znorm(jl)*pdmod(jl))
      endif
      ptau(jl,jk)  =0.0
      pzdep(jl,jk) =0.0
      ppsi(jl,jk)  =0.0
      ll1(jl,jk)   =.false.
  212 continue
  213 continue
      do 215 jk=2,klev
      do 214 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
        zdp(jl,jk)=papm1(jl,jk)-papm1(jl,jk-1)
        pvph(jl,jk)=((paphm1(jl,jk)-papm1(jl,jk-1))*zvpf(jl,jk)+
     *            (papm1(jl,jk)-paphm1(jl,jk))*zvpf(jl,jk-1))
     *            /zdp(jl,jk)
        if(pvph(jl,jk).lt.gvsec) then
          pvph(jl,jk)=gvsec
          kcrit(jl)=jk
        endif
      endif
  214 continue
  215 continue
c
c*         2.3     mean flow richardson number.
c
  230 continue
c
      do 232 jk=2,klev
      do 231 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
        zdwind=max(abs(zvpf(jl,jk)-zvpf(jl,jk-1)),gvsec)
        pri(jl,jk)=pstab(jl,jk)*(zdp(jl,jk)
     *          /(rg*prho(jl,jk)*zdwind))**2
        pri(jl,jk)=max(pri(jl,jk),grcrit)
      endif
  231 continue
  232 continue
  
c
c
c*      define top of 'envelope' layer
c       ----------------------------

      do 233 jl=kidia,kfdia
      pnu (jl)=0.0
      znum(jl)=0.0
 233  continue
      
      do 234 jk=2,klev-1
      do 234 jl=kidia,kfdia
      
      if(ktest(jl).eq.1) then
       
      if (jk.ge.kknu2(jl)) then
          
            znum(jl)=pnu(jl)
            zwind=(pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/
     *            max(sqrt(pulow(jl)**2+pvlow(jl)**2),gvsec)
            zwind=max(sqrt(zwind**2),gvsec)
            zdelp=paphm1(jl,jk+1)-paphm1(jl,jk)
            zstabm=sqrt(max(pstab(jl,jk  ),gssec))
            zstabp=sqrt(max(pstab(jl,jk+1),gssec))
            zrhom=prho(jl,jk  )
            zrhop=prho(jl,jk+1)
            pnu(jl) = pnu(jl) + (zdelp/rg)*
     *            ((zstabp/zrhop+zstabm/zrhom)/2.)/zwind     
            if((znum(jl).le.gfrcrit).and.(pnu(jl).gt.gfrcrit)
     *                          .and.(kkenvh(jl).eq.klev))
     *      kkenvh(jl)=jk
     
      endif    

      endif
      
 234  continue
      
c  calculation of a dynamical mixing height for when the waves
C  BREAK AT LOW LEVEL: The drag will be repartited over
C  a depths that depends on waves vertical wavelength,
C  not just between two adjacent model layers.
c  of gravity waves:

      do 235 jl=kidia,kfdia
      znup(jl)=0.0
      znum(jl)=0.0
 235  continue

      do 236 jk=klev-1,2,-1
      do 236 jl=kidia,kfdia
      
      if(ktest(jl).eq.1) then

            znum(jl)=znup(jl)
            zwind=(pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/
     *            max(sqrt(pulow(jl)**2+pvlow(jl)**2),gvsec)
            zwind=max(sqrt(zwind**2),gvsec)
            zdelp=paphm1(jl,jk+1)-paphm1(jl,jk)
            zstabm=sqrt(max(pstab(jl,jk  ),gssec))
            zstabp=sqrt(max(pstab(jl,jk+1),gssec))
            zrhom=prho(jl,jk  )
            zrhop=prho(jl,jk+1)
            znup(jl) = znup(jl) + (zdelp/rg)*
     *            ((zstabp/zrhop+zstabm/zrhom)/2.)/zwind     
            if((znum(jl).le.rpi/4.).and.(znup(jl).gt.rpi/4.)
     *                          .and.(kkcrith(jl).eq.klev))
     *      kkcrith(jl)=jk
     
      endif
      
 236  continue
 
      do 237 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      kkcrith(jl)=max0(kkcrith(jl),ilevh*2)
      kkcrith(jl)=max0(kkcrith(jl),kknu(jl))
      if(kcrit(jl).ge.kkcrith(jl))kcrit(jl)=1
      endif
 237  continue         
c
c     directional info for flow blocking ************************* 
c
      do 251 jk=1,klev    
      do 252 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      lo=(pum1(jl,jk).lt.gvsec).and.(pum1(jl,jk).ge.-gvsec)
      if(lo) then
        zu=pum1(jl,jk)+2.*gvsec
      else
        zu=pum1(jl,jk)
      endif
       zphi=atan(pvm1(jl,jk)/zu)
       ppsi(jl,jk)=ptheta(jl)*rpi/180.-zphi
      endif
 252  continue
 251  continue

c      forms the vertical 'leakiness' **************************

      do 254  jk=ilevh,klev
      do 253  jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      pzdep(jl,jk)=0
      if(jk.ge.kkenvh(jl).and.kkenvh(jl).ne.klev) then
        pzdep(jl,jk)=(pgeom1(jl,kkenvh(jl)  )-pgeom1(jl,  jk))/
     *               (pgeom1(jl,kkenvh(jl)  )-pgeom1(jl,klev))
      end if
      endif
 253  continue
 254  continue

      return
      end


