      SUBROUTINE orodrag( nlon,nlev 
     i                 , kgwd,  kdx, ktest
     r                 , ptsphy
     r                 , paphm1,papm1,pgeom1,pn2m1,ptm1,pum1,pvm1
     r                 , pmea, pstd, psig, pgam, pthe, ppic, pval
c outputs
     r                 , pulow,pvlow
     r                 , pvom,pvol,pte 
     r                 , blustr,blvstr,pnlow,zeff,ikenvh
c 3D and temporary outputs
     r                 , ztau,iknu2,ikbreak)
      
      use dimphy
      IMPLICIT NONE

c
c
c**** *orodrag* - does the SSO drag  parametrization.
c
c     purpose.
c     --------
c
c     this routine computes the physical tendencies of the
c     prognostic variables u,v  and t due to  vertical transports by
c     subgridscale orographically excited gravity waves, and to
c     low level blocked flow drag.
c
c**   interface.
c     ----------
c          called from *drag_noro*.
c
c          the routine takes its input from the long-term storage:
c          u,v,t and p at t-1.
c
c        explicit arguments :
c        --------------------
c     ==== inputs ===
c nlon----input-I-Total number of horizontal points that get into physics
c nlev----input-I-Number of vertical levels
c
c kgwd- -input-I: Total nb of points where the orography schemes are active
c ktest--input-I: Flags to indicate active points
c kdx----input-I: Locate the physical location of an active point.
c ptsphy--input-R-Time-step (s)
c paphm1--input-R: pressure at model 1/2 layer
c papm1---input-R: pressure at model layer
c pgeom1--input-R: Altitude of layer above ground
c pn2m1---input-R-Brunt-Vaisala freq.^2 at 1/2 layers
c ptm1, pum1, pvm1--R-: t, u and v
c pmea----input-R-Mean Orography (m)
C pstd----input-R-SSO standard deviation (m)
c psig----input-R-SSO slope
c pgam----input-R-SSO Anisotropy
c pthe----input-R-SSO Angle
c ppic----input-R-SSO Peacks elevation (m)
c pval----input-R-SSO Valleys elevation (m)

      integer nlon,nlev,kgwd
      real ptsphy

c     ==== outputs ===
c pulow, pvlow -output-R: Low-level wind
c
c pte -----output-R: T tendency
c pvom-----output-R: U tendency
c pvol-----output-R: V tendency
c
c
c Implicit Arguments:
c ===================
c
c klon-common-I: Number of points seen by the physics
c klev-common-I: Number of vertical layers
c
c     method.
c     -------
c
c     externals.
c     ----------
Coff  integer ismin, ismax
Coff  external ismin, ismax
c
c     reference.
c     ----------
c
c     author.
c     -------
c     m.miller + b.ritter   e.c.m.w.f.     15/06/86.
c
c     f.lott + m. miller    e.c.m.w.f.     22/11/94
c-----------------------------------------------------------------------
c
c
#include "YOMCST.h"
#include "YOEGWD.h"

c-----------------------------------------------------------------------
c
c*       0.1   arguments
c              ---------
c
c
      real  pte(nlon,nlev),
     *      pvol(nlon,nlev),
     *      pvom(nlon,nlev),
     *      pulow(nlon),
     *      pvlow(nlon)
      real  pum1(nlon,nlev),
     *      pvm1(nlon,nlev),
     *      ptm1(nlon,nlev),
     *      pmea(nlon),pstd(nlon),psig(nlon),
     *      pgam(nlon),pthe(nlon),ppic(nlon),pval(nlon),
     *      pgeom1(nlon,nlev),pn2m1(nlon,nlev),
     *      papm1(nlon,nlev),
     *      paphm1(nlon,nlev+1),
     *      pnlow(nlon), ! low level stability
     *      blustr(nlon),blvstr(nlon), ! blocked stress
     *      zeff(nlon) ! effective height

c
      integer  kdx(nlon),ktest(nlon)
c-----------------------------------------------------------------------
c
c*       0.2   local arrays
c              ------------
      integer  isect(klon),
     *         icrit(klon),
     *         ikcrith(klon),
     *         ikenvh(klon),
     *         iknu(klon),
     *         iknu2(klon),
     *         ikbreak(klon),
     *         ikcrit(klon),
     *         ikhlim(klon)
c
      real   ztau(klon,klev+1),
     *       zstab(klon,klev+1),
     *       zvph(klon,klev+1),
     *       zrho(klon,klev+1),
     *       zri(klon,klev+1),
     *       zpsi(klon,klev+1),
     *       zzdep(klon,klev)
      real   zdudt(klon),
     *       zdvdt(klon),
     *       zdtdt(klon),
     *       zdedt(klon),
     *       zvidis(klon),
     *       ztfr(klon),
     *       znu(klon),
     *       zd1(klon),
     *       zd2(klon),
     *       zdmod(klon)


c local quantities:

      integer jl,jk,ji
      real ztmst,zdelp,ztemp,zforc,ztend,rover                
      real zb,zc,zconb,zabsv,zzd1,ratio,zbet,zust,zvst,zdis
   
c
c------------------------------------------------------------------
c
c*         1.    initialization
c                --------------
c
c        print *,' in orodrag'
 100  continue
c
c     ------------------------------------------------------------------
c
c*         1.1   computational constants
c                -----------------------
c
 110  continue
c
c     ztmst=twodt
c     if(nstep.eq.nstart) ztmst=0.5*twodt
      ztmst=ptsphy
c     ------------------------------------------------------------------
c
 120  continue
c
c     ------------------------------------------------------------------
c
c*         1.3   check whether row contains point for printing
c                ---------------------------------------------
c
 130  continue
c
c     ------------------------------------------------------------------
c
c*         2.     precompute basic state variables.
c*                ---------- ----- ----- ----------
c*                define low level wind, project winds in plane of
c*                low level wind, determine sector in which to take
c*                the variance and set indicator for critical levels.
c

  200 continue
c
      do jk=1,klev
       zstab(:,jk) = pn2m1(:,jk)
      enddo
c
      call orosetup
     *     ( nlon, nlev , ktest 
     *     , ikcrit, ikcrith, icrit, isect, ikhlim, ikenvh,iknu,iknu2
     *     , ikbreak
     *     , paphm1, papm1 , pum1   , pvm1 , ptm1 , pgeom1, zstab, pstd
     *     , zrho  , zri   , ztau , zvph , zpsi, zzdep
     *     , pulow, pvlow 
     *     , pthe,pgam,pmea,ppic,pval,znu  ,zd1,  zd2,  zdmod )


      pnlow(:) = sqrt(zstab(:,klev+1))

c
c
c
c***********************************************************
c
c
c*         3.      compute low level stresses using subcritical and
c*                 supercritical forms.computes anisotropy coefficient
c*                 as measure of orographic twodimensionality.
c
  300 continue
c
      call gwstress
     *    ( nlon  , nlev
     *    , ikcrit, isect, ikhlim, ktest, ikcrith, icrit, ikenvh, iknu
     *    , zrho  , zstab, zvph  , pstd,  psig, pmea, ppic, pval
     *    , ztfr   , ztau 
     *    , pgeom1,pgam,zd1,zd2,zdmod,znu,zeff)

c
c
c*         4.      compute stress profile including
c                  trapped waves, wave breaking,
c                  linear decay in stratosphere.
c
  400 continue
c
c

      call gwprofil
     *       (  nlon , nlev
     *       , kgwd   , kdx  , ktest
     *       , ikcrit, ikcrith, icrit , ikenvh, iknu
     *       ,iknu2 , ikbreak, paphm1, zrho , zstab , ztfr , zvph
     *       , zri   , ztau 
 
     *       , zdmod , znu    , psig  , pgam , pstd , ppic , pval)

c
c*         5.      Compute tendencies from waves stress profile.
c                  Compute low level blocked flow drag. 
c*                 --------------------------------------------
c
  500 continue

      
c
c  explicit solution at all levels for the gravity wave
c  implicit solution for the blocked levels

      do 510 jl=kidia,kfdia
      zvidis(jl)=0.0
      zdudt(jl)=0.0
      zdvdt(jl)=0.0
      zdtdt(jl)=0.0
      blustr(jl)=0.0
      blvstr(jl)=0.0
  510 continue
c

      do 524 jk=1,klev
c

C  WAVE STRESS 
C-------------
c
c
      do 523 ji=kidia,kfdia

      if(ktest(ji).eq.1) then

      zdelp=paphm1(ji,jk+1)-paphm1(ji,jk)
      ztemp=-rg*(ztau(ji,jk+1)-ztau(ji,jk))/(zvph(ji,klev+1)*zdelp)

      zdudt(ji)=(pulow(ji)*zd1(ji)-pvlow(ji)*zd2(ji))*ztemp/zdmod(ji)
      zdvdt(ji)=(pvlow(ji)*zd1(ji)+pulow(ji)*zd2(ji))*ztemp/zdmod(ji)
c
c Control Overshoots
c

      if(jk.ge.ntop)then
        rover=0.10
        if(abs(zdudt(ji)).gt.rover*abs(pum1(ji,jk))/ztmst)
     C    zdudt(ji)=rover*abs(pum1(ji,jk))/ztmst*
     C              zdudt(ji)/(abs(zdudt(ji))+1.E-10)
        if(abs(zdvdt(ji)).gt.rover*abs(pvm1(ji,jk))/ztmst)
     C    zdvdt(ji)=rover*abs(pvm1(ji,jk))/ztmst*
     C              zdvdt(ji)/(abs(zdvdt(ji))+1.E-10)
      endif 

      rover=0.25
      zforc=sqrt(zdudt(ji)**2+zdvdt(ji)**2)        
      ztend=sqrt(pum1(ji,jk)**2+pvm1(ji,jk)**2)/ztmst

      if(zforc.ge.rover*ztend)then
        zdudt(ji)=rover*ztend/zforc*zdudt(ji)
        zdvdt(ji)=rover*ztend/zforc*zdvdt(ji)
      endif
c
c BLOCKED FLOW DRAG:
C -----------------
c
      if(jk.gt.ikenvh(ji)) then
         zb=1.0-0.18*pgam(ji)-0.04*pgam(ji)**2
         zc=0.48*pgam(ji)+0.3*pgam(ji)**2
         zconb=2.*ztmst*gkwake*psig(ji)/(4.*pstd(ji))
         zabsv=sqrt(pum1(ji,jk)**2+pvm1(ji,jk)**2)/2.
         zzd1=zb*cos(zpsi(ji,jk))**2+zc*sin(zpsi(ji,jk))**2
         ratio=(cos(zpsi(ji,jk))**2+pgam(ji)*sin(zpsi(ji,jk))**2)/
     *   (pgam(ji)*cos(zpsi(ji,jk))**2+sin(zpsi(ji,jk))**2)
         zbet=max(0.,2.-1./ratio)*zconb*zzdep(ji,jk)*zzd1*zabsv
c
c OPPOSED TO THE WIND
c
         zdudt(ji)=-pum1(ji,jk)/ztmst
         zdvdt(ji)=-pvm1(ji,jk)/ztmst
c
c PERPENDICULAR TO THE SSO MAIN AXIS:
C                            
cmod     zdudt(ji)=-(pum1(ji,jk)*cos(pthe(ji)*rpi/180.)
cmod *              +pvm1(ji,jk)*sin(pthe(ji)*rpi/180.))
cmod *              *cos(pthe(ji)*rpi/180.)/ztmst
cmod     zdvdt(ji)=-(pum1(ji,jk)*cos(pthe(ji)*rpi/180.)
cmod *              +pvm1(ji,jk)*sin(pthe(ji)*rpi/180.))
cmod *              *sin(pthe(ji)*rpi/180.)/ztmst
C
         zdudt(ji)=zdudt(ji)*(zbet/(1.+zbet))
         zdvdt(ji)=zdvdt(ji)*(zbet/(1.+zbet))

c output blocked flow stress
         blustr(ji)  = blustr(ji)
     .              +zdudt(ji)*(paphm1(ji,jk+1)-paphm1(ji,jk))/rg
         blvstr(ji)  = blvstr(ji)
     .              +zdvdt(ji)*(paphm1(ji,jk+1)-paphm1(ji,jk))/rg


      end if
      pvom(ji,jk)=zdudt(ji)
      pvol(ji,jk)=zdvdt(ji)
      zust=pum1(ji,jk)+ztmst*zdudt(ji)
      zvst=pvm1(ji,jk)+ztmst*zdvdt(ji)
      zdis=0.5*(pum1(ji,jk)**2+pvm1(ji,jk)**2-zust**2-zvst**2)
      zdedt(ji)=zdis/ztmst
      zvidis(ji)=zvidis(ji)+zdis*zdelp
c VENUS ATTENTION: CP VARIABLE
      zdtdt(ji)=zdedt(ji)/rcpd
c
c  NO TENDENCIES ON TEMPERATURE .....
c
c  Instead of, pte(ji,jk)=zdtdt(ji), due to mechanical dissipation
c
      pte(ji,jk)=0.0

      endif

  523 continue
  524 continue
c
c
  501 continue

      return
      end

