      SUBROUTINE drag_noro_strato (nlon,nlev,dtime,paprs,pplay,
     e                   pmea,pstd, psig, pgam, pthe,ppic,pval,
     e                   kgwd,kdx,ktest,
     e                   t, u, v,
     s                   pulow, pvlow, pustr, pvstr,
     s                   d_t, d_u, d_v)
c
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): F.Lott (LMD/CNRS) date: 19950201
c Object: Mountain drag interface. Made necessary because:
C 1. in the LMD-GCM Layers are from bottom to top,
C    contrary to most European GCM.
c 2. the altitude above ground of each model layers
c    needs to be known (variable zgeom)
c======================================================================
c Explicit Arguments:
c ==================
c nlon----input-I-Total number of horizontal points that get into physics
c nlev----input-I-Number of vertical levels
c dtime---input-R-Time-step (s)
c paprs---input-R-Pressure in semi layers    (Pa)
c pplay---input-R-Pressure model-layers      (Pa)
c t-------input-R-temperature (K)
c u-------input-R-Horizontal wind (m/s)
c v-------input-R-Meridional wind (m/s)
c pmea----input-R-Mean Orography (m)
C pstd----input-R-SSO standard deviation (m)
c psig----input-R-SSO slope
c pgam----input-R-SSO Anisotropy
c pthe----input-R-SSO Angle
c ppic----input-R-SSO Peacks elevation (m)
c pval----input-R-SSO Valleys elevation (m)
c
c kgwd- -input-I: Total nb of points where the orography schemes are active
c ktest--input-I: Flags to indicate active points
c kdx----input-I: Locate the physical location of an active point.

c pulow, pvlow -output-R: Low-level wind
c pustr, pvstr -output-R: Surface stress due to SSO drag      (Pa)
c
c d_t-----output-R: T increment            
c d_u-----output-R: U increment              
c d_v-----output-R: V increment              
c
c Implicit Arguments:
c ===================
c
c iim--common-I: Number of longitude intervals
c jjm--common-I: Number of latitude intervals
c klon-common-I: Number of points seen by the physics
c                (iim+1)*(jjm+1) for instance
c klev-common-I: Number of vertical layers
c======================================================================
c Local Variables:
c ================
c
c zgeom-----R: Altitude of layer above ground
c pt, pu, pv --R: t u v from top to bottom
c pdtdt, pdudt, pdvdt --R: t u v tendencies (from top to bottom) 
c papmf: pressure at model layer (from top to bottom)
c papmh: pressure at model 1/2 layer (from top to bottom)
c 
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "YOEGWD.h"
c
c  ARGUMENTS
c
      INTEGER nlon,nlev
      REAL dtime
      REAL paprs(nlon,nlev+1)
      REAL pplay(nlon,nlev)
      REAL pmea(nlon),pstd(nlon),psig(nlon),pgam(nlon),pthe(nlon)
      REAL ppic(nlon),pval(nlon)
      REAL pulow(nlon),pvlow(nlon),pustr(nlon),pvstr(nlon)
      REAL t(nlon,nlev), u(nlon,nlev), v(nlon,nlev)
      REAL d_t(nlon,nlev), d_u(nlon,nlev), d_v(nlon,nlev)
c
      INTEGER i, k, kgwd,  kdx(nlon), ktest(nlon)
c
c LOCAL VARIABLES:
c
      REAL zgeom(klon,klev)
      REAL pdtdt(klon,klev), pdudt(klon,klev), pdvdt(klon,klev)
      REAL pt(klon,klev), pu(klon,klev), pv(klon,klev)
      REAL papmf(klon,klev),papmh(klon,klev+1)
      CHARACTER (LEN=20) :: modname='orografi_strato'
      CHARACTER (LEN=80) :: abort_message
c
c INITIALIZE OUTPUT VARIABLES 
c
      DO i = 1,klon
         pulow(i) = 0.0
         pvlow(i) = 0.0
         pustr(i) = 0.0
         pvstr(i) = 0.0
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_u(i,k) = 0.0
         d_v(i,k) = 0.0
         pdudt(i,k)=0.0
         pdvdt(i,k)=0.0
         pdtdt(i,k)=0.0
      ENDDO
      ENDDO
c
c PREPARE INPUT VARIABLES FOR ORODRAG (i.e., ORDERED FROM TOP TO BOTTOM)
C CALCULATE LAYERS HEIGHT ABOVE GROUND)
c
      DO k = 1, klev
      DO i = 1, klon
         pt(i,k) = t(i,klev-k+1) 
         pu(i,k) = u(i,klev-k+1)
         pv(i,k) = v(i,klev-k+1)
         papmf(i,k) = pplay(i,klev-k+1)
      ENDDO
      ENDDO
      DO k = 1, klev+1
      DO i = 1, klon
         papmh(i,k) = paprs(i,klev-k+2)
      ENDDO
      ENDDO
      DO i = 1, klon
         zgeom(i,klev) = RD * pt(i,klev)
     .                  * LOG(papmh(i,klev+1)/papmf(i,klev))
      ENDDO
      DO k = klev-1, 1, -1
      DO i = 1, klon
         zgeom(i,k) = zgeom(i,k+1) + RD * (pt(i,k)+pt(i,k+1))/2.0
     .               * LOG(papmf(i,k+1)/papmf(i,k))
      ENDDO
      ENDDO
c
c CALL SSO DRAG ROUTINES        
c
      CALL orodrag_strato(klon,klev,kgwd,kdx,ktest,
     .            dtime,
     .            papmh, papmf, zgeom,
     .            pt, pu, pv,
     .            pmea, pstd, psig, pgam, pthe, ppic,pval,
     .            pulow,pvlow,
     .            pdudt,pdvdt,pdtdt)
C
C COMPUTE INCREMENTS AND STRESS FROM TENDENCIES

      DO k = 1, klev
      DO i = 1, klon
         d_u(i,klev+1-k) = dtime*pdudt(i,k)
         d_v(i,klev+1-k) = dtime*pdvdt(i,k)
         d_t(i,klev+1-k) = dtime*pdtdt(i,k)
         pustr(i)        = pustr(i)
     .                    +pdudt(i,k)*(papmh(i,k+1)-papmh(i,k))/rg
         pvstr(i)        = pvstr(i)
     .                    +pdvdt(i,k)*(papmh(i,k+1)-papmh(i,k))/rg
      ENDDO
      ENDDO
c
      RETURN
      END

      SUBROUTINE orodrag_strato( nlon,nlev 
     i                 , kgwd,  kdx, ktest
     r                 , ptsphy
     r                 , paphm1,papm1,pgeom1,ptm1,pum1,pvm1
     r                 , pmea, pstd, psig, pgam, pthe, ppic, pval
c outputs
     r                 , pulow,pvlow
     r                 , pvom,pvol,pte )
      
      USE dimphy
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
      integer ismin, ismax
      external ismin, ismax
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
cym#include "dimensions.h"
cym#include "dimphy.h"
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
     *      pgeom1(nlon,nlev),
     *      papm1(nlon,nlev),
     *      paphm1(nlon,nlev+1)
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
c
c
      call orosetup_strato
     *     ( nlon, nlev , ktest 
     *     , ikcrit, ikcrith, icrit, isect, ikhlim, ikenvh,iknu,iknu2
     *     , paphm1, papm1 , pum1   , pvm1 , ptm1 , pgeom1, pstd
     *     , zrho  , zri   , zstab  , ztau , zvph , zpsi, zzdep
     *     , pulow, pvlow 
     *     , pthe,pgam,pmea,ppic,pval,znu  ,zd1,  zd2,  zdmod )


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
      call gwstress_strato
     *    ( nlon  , nlev
     *    , ikcrit, isect, ikhlim, ktest, ikcrith, icrit, ikenvh, iknu
     *    , zrho  , zstab, zvph  , pstd,  psig, pmea, ppic, pval
     *    , ztfr   , ztau 
     *    , pgeom1,pgam,zd1,zd2,zdmod,znu)

c
c
c*         4.      compute stress profile including
c                  trapped waves, wave breaking,
c                  linear decay in stratosphere.
c
  400 continue
c
c

      call gwprofil_strato
     *       (  nlon , nlev
     *       , kgwd   , kdx  , ktest
     *       , ikcrit, ikcrith, icrit  , ikenvh, iknu
     *       ,iknu2 , paphm1, zrho   , zstab , ztfr   , zvph
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

      if(jk.ge.nstra)then
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
      end if
      pvom(ji,jk)=zdudt(ji)
      pvol(ji,jk)=zdvdt(ji)
      zust=pum1(ji,jk)+ztmst*zdudt(ji)
      zvst=pvm1(ji,jk)+ztmst*zdvdt(ji)
      zdis=0.5*(pum1(ji,jk)**2+pvm1(ji,jk)**2-zust**2-zvst**2)
      zdedt(ji)=zdis/ztmst
      zvidis(ji)=zvidis(ji)+zdis*zdelp
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
      SUBROUTINE orosetup_strato
     *         ( nlon   , nlev  , ktest
     *         , kkcrit, kkcrith, kcrit, ksect , kkhlim
     *         , kkenvh, kknu  , kknu2
     *         , paphm1, papm1 , pum1   , pvm1 , ptm1  , pgeom1, pstd
     *         , prho  , pri   , pstab  , ptau , pvph  ,ppsi, pzdep
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
c pstab-----R-: Brunt-Vaisala freq. at 1/2 layers
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
      USE dimphy
      implicit none
c

cym#include "dimensions.h"
cym#include "dimphy.h"
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
      integer ilevh ,jl,jk
      real zcons1,zcons2,zhgeo,zu,zphi
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
c       
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
      zcons2=rg**2/rcpd
c
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
      zhcrit(jl,jk)=ppic(jl)-pmea(jl)
      zhgeo=pgeom1(jl,jk)/rg
      ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
      if(ll1(jl,jk).neqv.ll1(jl,jk+1)) then
        kknu2(jl)=jk
      endif
      if(.not.ll1(jl,ilevh))kknu2(jl)=ilevh
      endif
 2005 continue
 2004 continue
      do 2006 jk=klev,ilevh,-1
      do 2007 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      zhcrit(jl,jk)=amin1(ppic(jl)-pmea(jl),pmea(jl)-pval(jl))
      zhgeo=pgeom1(jl,jk)/rg
      ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
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
      kknul(jl)=klev
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
        pstab(jl,jk)=2.*zcons2/(ptm1(jl,jk)+ptm1(jl,jk-1))*
     *  (1.-rcpd*prho(jl,jk)*(ptm1(jl,jk)-ptm1(jl,jk-1))/zdp(jl,jk))
        pstab(jl,jk)=max(pstab(jl,jk),gssec)
      endif
  222 continue
  223 continue
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
      SUBROUTINE gwstress_strato
     *         (  nlon  , nlev
     *         , kkcrit, ksect, kkhlim, ktest, kkcrith, kcrit, kkenvh
     *         , kknu
     *         , prho  , pstab , pvph  , pstd, psig
     *         , pmea , ppic , pval  , ptfr  , ptau  
     *         , pgeom1 , pgamma , pd1  , pd2   , pdmod , pnu )
c
c**** *gwstress*
c
c     purpose.
c     --------
c  Compute the surface stress due to Gravity Waves, according
c  to the Phillips (1979) theory of 3-D flow above 
c  anisotropic elliptic ridges.

C  The stress is reduced two account for cut-off flow over
C  hill.  The flow only see that part of the ridge located
c  above the blocked layer (see zeff).
c
c**   interface.
c     ----------
c     call *gwstress*  from *gwdrag*
c
c        explicit arguments :
c        --------------------
c     ==== inputs ===
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
c   LOTT and MILLER (1997)  &  LOTT (1999)
c
c     author.
c     -------
c
c     modifications.
c     --------------
c     f. lott put the new gwd on ifs      22/11/93
c
c-----------------------------------------------------------------------
      USE dimphy
      implicit none

cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "YOEGWD.h"

c-----------------------------------------------------------------------
c
c*       0.1   arguments
c              ---------
c
      integer nlon,nlev
      integer kkcrit(nlon),kkcrith(nlon),kcrit(nlon),ksect(nlon),
     *        kkhlim(nlon),ktest(nlon),kkenvh(nlon),kknu(nlon)
c
      real prho(nlon,nlev+1),pstab(nlon,nlev+1),ptau(nlon,nlev+1),
     *     pvph(nlon,nlev+1),ptfr(nlon),
     *     pgeom1(nlon,nlev),pstd(nlon)
c
      real pd1(nlon),pd2(nlon),pnu(nlon),psig(nlon),pgamma(nlon)
      real pmea(nlon),ppic(nlon),pval(nlon)
      real pdmod(nlon)
c
c-----------------------------------------------------------------------
c
c*       0.2   local arrays
c              ------------
c  zeff--real: effective height seen by the flow when there is blocking

      integer jl
      real zeff  
c
c-----------------------------------------------------------------------
c
c*       0.3   functions
c              ---------
c     ------------------------------------------------------------------
c
c*         1.    initialization
c                --------------
c
c      PRINT *,' in gwstress'
 100  continue
c
c*         3.1     gravity wave stress.
c
  300 continue
c
c
      do 301 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      
c  effective mountain height above the blocked flow
  
         zeff=ppic(jl)-pval(jl)
         if(kkenvh(jl).lt.klev)then
         zeff=amin1(GFRCRIT*pvph(jl,klev+1)/sqrt(pstab(jl,klev+1))
     c              ,zeff)
         endif

      
        ptau(jl,klev+1)=gkdrag*prho(jl,klev+1)
     *     *psig(jl)*pdmod(jl)/4./pstd(jl)
     *     *pvph(jl,klev+1)*sqrt(pstab(jl,klev+1))
     *     *zeff**2


c  too small value of stress or  low level flow include critical level
c  or low level flow:  gravity wave stress nul.
                
c       lo=(ptau(jl,klev+1).lt.gtsec).or.(kcrit(jl).ge.kknu(jl))
c    *      .or.(pvph(jl,klev+1).lt.gvcrit)
c       if(lo) ptau(jl,klev+1)=0.0
      
c      print *,jl,ptau(jl,klev+1)

      else
      
          ptau(jl,klev+1)=0.0
          
      endif

  301 continue

c      write(21)(ptau(jl,klev+1),jl=kidia,kfdia)
 
      return
      end

      subroutine gwprofil_strato
     *         ( nlon, nlev
     *         , kgwd ,kdx  , ktest
     *         , kkcrit, kkcrith, kcrit ,  kkenvh, kknu,kknu2
     *         , paphm1, prho   , pstab , ptfr , pvph , pri , ptau
     *         , pdmod   , pnu   , psig ,pgamma, pstd, ppic,pval)

C**** *gwprofil*
C
C     purpose.
C     --------
C
C**   interface.
C     ----------
C          from *gwdrag*
C
C        explicit arguments :
C        --------------------
C     ==== inputs ===
C
C     ==== outputs ===
C
C        implicit arguments :   none
C        --------------------
C
C     method:
C     -------
C     the stress profile for gravity waves is computed as follows:
C     it decreases linearly with heights from the ground 
C     to the low-level indicated by kkcrith,
C     to simulates lee waves or 
C     low-level gravity wave breaking.
C     above it is constant, except when the waves encounter a critical
C     level (kcrit) or when they break.
C     The stress is also uniformly distributed above the level
C     nstra.                                          
C
      USE dimphy
      IMPLICIT NONE

cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "YOEGWD.h"

C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      integer nlon,nlev,kgwd
      integer kkcrit(nlon),kkcrith(nlon),kcrit(nlon)
     *       ,kdx(nlon),ktest(nlon)
     *       ,kkenvh(nlon),kknu(nlon),kknu2(nlon)
C
      real paphm1(nlon,nlev+1), pstab(nlon,nlev+1),
     *     prho  (nlon,nlev+1), pvph (nlon,nlev+1),
     *     pri   (nlon,nlev+1), ptfr (nlon), ptau(nlon,nlev+1)
     
      real pdmod (nlon) , pnu (nlon) , psig(nlon),
     *     pgamma(nlon) , pstd(nlon) , ppic(nlon), pval(nlon)
     
C-----------------------------------------------------------------------
C
C*       0.2   local arrays
C              ------------
C
      integer jl,jk
      real zsqr,zalfa,zriw,zdel,zb,zalpha,zdz2n,zdelp,zdelpt

      real zdz2 (klon,klev) , znorm(klon) , zoro(klon)
      real ztau (klon,klev+1)
C
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
C      print *,' entree gwprofil' 
 100  CONTINUE
C
C
C*    COMPUTATIONAL CONSTANTS.
C     ------------- ----------
C
      do 400 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
      zoro(jl)=psig(jl)*pdmod(jl)/4./pstd(jl)
      ztau(jl,klev+1)=ptau(jl,klev+1)
c     print *,jl,ptau(jl,klev+1)
      ztau(jl,kkcrith(jl))=grahilo*ptau(jl,klev+1)
      endif
  400 continue
  
C
      do 430 jk=klev+1,1,-1
C
C
C*         4.1    constant shear stress until top of the
C                 low-level breaking/trapped layer
  410 CONTINUE
C
      do 411 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
           if(jk.gt.kkcrith(jl)) then
           zdelp=paphm1(jl,jk)-paphm1(jl,klev+1) 
           zdelpt=paphm1(jl,kkcrith(jl))-paphm1(jl,klev+1) 
           ptau(jl,jk)=ztau(jl,klev+1)+zdelp/zdelpt*
     c                 (ztau(jl,kkcrith(jl))-ztau(jl,klev+1))
           else                    
           ptau(jl,jk)=ztau(jl,kkcrith(jl))
           endif
       endif
 411  continue             
C
C*         4.15   constant shear stress until the top of the
C                 low level flow layer.
 415  continue
C        
C
C*         4.2    wave displacement at next level.
C
  420 continue
C
  430 continue

C
C*         4.4    wave richardson number, new wave displacement
C*                and stress:  breaking evaluation and critical 
C                 level
C
                          
      do 440 jk=klev,1,-1

      do 441 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
      znorm(jl)=prho(jl,jk)*sqrt(pstab(jl,jk))*pvph(jl,jk)
      zdz2(jl,jk)=ptau(jl,jk)/amax1(znorm(jl),gssec)/zoro(jl)
      endif
  441 continue

      do 442 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
          if(jk.lt.kkcrith(jl)) then
          if((ptau(jl,jk+1).lt.gtsec).or.(jk.le.kcrit(jl))) then
             ptau(jl,jk)=0.0
          else
               zsqr=sqrt(pri(jl,jk))
               zalfa=sqrt(pstab(jl,jk)*zdz2(jl,jk))/pvph(jl,jk)
               zriw=pri(jl,jk)*(1.-zalfa)/(1+zalfa*zsqr)**2
               if(zriw.lt.grcrit) then
C                 print *,' breaking!!!',ptau(jl,jk)
                  zdel=4./zsqr/grcrit+1./grcrit**2+4./grcrit
                  zb=1./grcrit+2./zsqr
                  zalpha=0.5*(-zb+sqrt(zdel))
                  zdz2n=(pvph(jl,jk)*zalpha)**2/pstab(jl,jk)
                  ptau(jl,jk)=znorm(jl)*zdz2n*zoro(jl)
               endif
                
               ptau(jl,jk)=amin1(ptau(jl,jk),ptau(jl,jk+1))
                  
          endif
          endif
      endif
  442 continue
  440 continue

C  REORGANISATION OF THE STRESS PROFILE AT LOW LEVEL

      do 530 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
         ztau(jl,kkcrith(jl)-1)=ptau(jl,kkcrith(jl)-1)
         ztau(jl,nstra)=ptau(jl,nstra)
      endif
 530  continue      

      do 531 jk=1,klev
      
      do 532 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
                
         if(jk.gt.kkcrith(jl)-1)then

          zdelp=paphm1(jl,jk)-paphm1(jl,klev+1    )
          zdelpt=paphm1(jl,kkcrith(jl)-1)-paphm1(jl,klev+1    )
          ptau(jl,jk)=ztau(jl,klev+1    ) +
     .                (ztau(jl,kkcrith(jl)-1)-ztau(jl,klev+1    ) )*
     .                zdelp/zdelpt
     
        endif
      endif
            
 532  continue    
 
C  REORGANISATION AT THE MODEL TOP....

      do 533 jl=kidia,kfdia
      if(ktest(jl).eq.1)then

         if(jk.lt.nstra)then

          zdelp =paphm1(jl,nstra)
          zdelpt=paphm1(jl,jk)
          ptau(jl,jk)=ztau(jl,nstra)*zdelpt/zdelp 
c         ptau(jl,jk)=ztau(jl,nstra)                

        endif

      endif

 533  continue

 
 531  continue        


 123   format(i4,1x,20(f6.3,1x))


      return
      end
      subroutine lift_noro_strato (nlon,nlev,dtime,paprs,pplay,      
     i                   plat,pmea,pstd, psig, pgam, pthe, ppic,pval,
     i                   kgwd,kdx,ktest,
     i                   t, u, v,
     o                   pulow, pvlow, pustr, pvstr,
     o                   d_t, d_u, d_v)
c
      USE dimphy
      implicit none
c======================================================================
c Auteur(s): F.Lott (LMD/CNRS) date: 19950201
c Object: Mountain lift interface (enhanced vortex stretching).
c         Made necessary because:
C 1. in the LMD-GCM Layers are from bottom to top,
C    contrary to most European GCM.
c 2. the altitude above ground of each model layers
c    needs to be known (variable zgeom)
c======================================================================
c Explicit Arguments:
c ==================
c nlon----input-I-Total number of horizontal points that get into physics
c nlev----input-I-Number of vertical levels
c dtime---input-R-Time-step (s)
c paprs---input-R-Pressure in semi layers    (Pa)
c pplay---input-R-Pressure model-layers      (Pa)
c t-------input-R-temperature (K)
c u-------input-R-Horizontal wind (m/s)
c v-------input-R-Meridional wind (m/s)
c pmea----input-R-Mean Orography (m)
C pstd----input-R-SSO standard deviation (m)
c psig----input-R-SSO slope
c pgam----input-R-SSO Anisotropy
c pthe----input-R-SSO Angle
c ppic----input-R-SSO Peacks elevation (m)
c pval----input-R-SSO Valleys elevation (m)
c
c kgwd- -input-I: Total nb of points where the orography schemes are active
c ktest--input-I: Flags to indicate active points
c kdx----input-I: Locate the physical location of an active point.

c pulow, pvlow -output-R: Low-level wind
c pustr, pvstr -output-R: Surface stress due to SSO drag      (Pa)
c
c d_t-----output-R: T increment
c d_u-----output-R: U increment
c d_v-----output-R: V increment
c
c Implicit Arguments:
c ===================
c
c iim--common-I: Number of longitude intervals
c jjm--common-I: Number of latitude intervals
c klon-common-I: Number of points seen by the physics
c                (iim+1)*(jjm+1) for instance
c klev-common-I: Number of vertical layers
c======================================================================
c Local Variables:
c ================
c
c zgeom-----R: Altitude of layer above ground
c pt, pu, pv --R: t u v from top to bottom
c pdtdt, pdudt, pdvdt --R: t u v tendencies (from top to bottom)
c papmf: pressure at model layer (from top to bottom)
c papmh: pressure at model 1/2 layer (from top to bottom)
c
c======================================================================

cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "YOEGWD.h"
c
c ARGUMENTS
c
      INTEGER nlon,nlev
      REAL dtime
      REAL paprs(klon,klev+1)
      REAL pplay(klon,klev)
      REAL plat(nlon),pmea(nlon)
      REAL pstd(nlon),psig(nlon),pgam(nlon),pthe(nlon)
      REAL ppic(nlon),pval(nlon)
      REAL pulow(nlon),pvlow(nlon),pustr(nlon),pvstr(nlon)
      REAL t(nlon,nlev), u(nlon,nlev), v(nlon,nlev)
      REAL d_t(nlon,nlev), d_u(nlon,nlev), d_v(nlon,nlev)
c
      INTEGER i, k, kgwd,  kdx(nlon), ktest(nlon)
c
c Variables locales:
c
      REAL zgeom(klon,klev)
      REAL pdtdt(klon,klev), pdudt(klon,klev), pdvdt(klon,klev)
      REAL pt(klon,klev), pu(klon,klev), pv(klon,klev)
      REAL papmf(klon,klev),papmh(klon,klev+1)
c
c initialiser les variables de sortie (pour securite)
c

c     print *,'in lift_noro'
      DO i = 1,klon
         pulow(i) = 0.0
         pvlow(i) = 0.0
         pustr(i) = 0.0
         pvstr(i) = 0.0
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_u(i,k) = 0.0
         d_v(i,k) = 0.0
         pdudt(i,k)=0.0
         pdvdt(i,k)=0.0
         pdtdt(i,k)=0.0
      ENDDO
      ENDDO
c
c preparer les variables d'entree (attention: l'ordre des niveaux 
c verticaux augmente du haut vers le bas)
c
      DO k = 1, klev
      DO i = 1, klon
         pt(i,k) = t(i,klev-k+1) 
         pu(i,k) = u(i,klev-k+1)
         pv(i,k) = v(i,klev-k+1)
         papmf(i,k) = pplay(i,klev-k+1)
      ENDDO
      ENDDO
      DO k = 1, klev+1
      DO i = 1, klon
         papmh(i,k) = paprs(i,klev-k+2)
      ENDDO
      ENDDO
      DO i = 1, klon
         zgeom(i,klev) = RD * pt(i,klev)
     .                  * LOG(papmh(i,klev+1)/papmf(i,klev))
      ENDDO
      DO k = klev-1, 1, -1
      DO i = 1, klon
         zgeom(i,k) = zgeom(i,k+1) + RD * (pt(i,k)+pt(i,k+1))/2.0
     .               * LOG(papmf(i,k+1)/papmf(i,k))
      ENDDO
      ENDDO
c
c appeler la routine principale
c

      CALL OROLIFT_strato(klon,klev,kgwd,kdx,ktest,
     .            dtime,
     .            papmh, papmf, zgeom,
     .            pt, pu, pv,
     .            plat,pmea, pstd, psig, pgam, pthe, ppic,pval,
     .            pulow,pvlow,
     .            pdudt,pdvdt,pdtdt)
C
      DO k = 1, klev
      DO i = 1, klon
         d_u(i,klev+1-k) = dtime*pdudt(i,k)
         d_v(i,klev+1-k) = dtime*pdvdt(i,k)
         d_t(i,klev+1-k) = dtime*pdtdt(i,k)
         pustr(i)        = pustr(i)
     .                    +pdudt(i,k)*(papmh(i,k+1)-papmh(i,k))/rg
         pvstr(i)        = pvstr(i)
     .                    +pdvdt(i,k)*(papmh(i,k+1)-papmh(i,k))/rg
      ENDDO
      ENDDO

c     print *,' out lift_noro'
c
      RETURN
      END
      subroutine orolift_strato( nlon,nlev
     I                 , kgwd, kdx, ktest
     R                 , ptsphy
     R                 , paphm1,papm1,pgeom1,ptm1,pum1,pvm1
     R                 , plat
     R                 , pmea, pstd, psig, pgam, pthe,ppic,pval
C OUTPUTS
     R                 , pulow,pvlow
     R                 , pvom,pvol,pte )

C
C**** *OROLIFT: SIMULATE THE GEOSTROPHIC LIFT.
C
C     PURPOSE.
C     --------
C this routine computes the physical tendencies of the
C prognostic variables u,v  when enhanced vortex stretching
C is needed.
C
C**   INTERFACE.
C     ----------
C          CALLED FROM *lift_noro
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
c ptm1, pum1, pvm1--R-: t, u and v
c pmea----input-R-Mean Orography (m)
C pstd----input-R-SSO standard deviation (m)
c psig----input-R-SSO slope
c pgam----input-R-SSO Anisotropy
c pthe----input-R-SSO Angle
c ppic----input-R-SSO Peacks elevation (m)
c pval----input-R-SSO Valleys elevation (m)
c plat----input-R-Latitude (degree)
c
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

C     ----------
C
C     AUTHOR.
C     -------
C     F.LOTT  LMD 22/11/95
C
       USE dimphy
       implicit none
C
C
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "YOEGWD.h"
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
C
      integer nlon,nlev,kgwd
      real ptsphy
      real  pte(nlon,nlev),
     *      pvol(nlon,nlev),
     *      pvom(nlon,nlev),
     *      pulow(nlon),
     *      pvlow(nlon)
      real  pum1(nlon,nlev),
     *      pvm1(nlon,nlev),
     *      ptm1(nlon,nlev),
     *      plat(nlon),pmea(nlon),
     *      pstd(nlon),psig(nlon),pgam(nlon),
     *      pthe(nlon),ppic(nlon),pval(nlon),
     *      pgeom1(nlon,nlev),
     *      papm1(nlon,nlev),
     *      paphm1(nlon,nlev+1)
C
      INTEGER  KDX(NLON),KTEST(NLON)
C-----------------------------------------------------------------------
C
C*       0.2   local arrays

      integer jl,ilevh,jk
      real zhgeo,zdelp,zslow,zsqua,zscav,zbet
C              ------------
      integer  iknub(klon),
     *         iknul(klon)
      logical ll1(klon,klev+1)
C
      real   ztau(klon,klev+1),
     *       ztav(klon,klev+1),
     *       zrho(klon,klev+1)
      real   zdudt(klon),
     *       zdvdt(klon)
      real zhcrit(klon,klev)

      logical lifthigh
      real zcons1,ztmst
      CHARACTER (LEN=20) :: modname='orolift_strato'
      CHARACTER (LEN=80) :: abort_message


C-----------------------------------------------------------------------
C
C*         1.1  initialisations
C               ---------------

      lifthigh=.false.

      if(nlon.ne.klon.or.nlev.ne.klev) then
        abort_message = 'pb dimension'
        CALL abort_gcm (modname,abort_message,1)
      ENDIF
      zcons1=1./rd
      ztmst=ptsphy
C
      do 1001 jl=kidia,kfdia
      zrho(jl,klev+1)  =0.0
      pulow(jl)        =0.0
      pvlow(jl)        =0.0
      iknub(JL)   =klev
      iknul(JL)   =klev
      ilevh=klev/3
      ll1(jl,klev+1)=.false.
      do 1000 jk=1,klev
      pvom(jl,jk)=0.0
      pvol(jl,jk)=0.0
      pte (jl,jk)=0.0
 1000 continue
 1001 continue

C
C*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
C*                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
C*                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.
C
C
C
      do 2006 jk=klev,1,-1
      do 2007 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      zhcrit(jl,jk)=amax1(ppic(jl)-pval(jl),100.)
      zhgeo=pgeom1(jl,jk)/rg
      ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
      if(ll1(jl,jk).neqv.ll1(jl,jk+1)) then
        iknub(jl)=jk
      endif
      endif
 2007 continue
 2006 continue
C

      do 2010 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      iknub(jl)=max(iknub(jl),klev/2)
      iknul(jl)=max(iknul(jl),2*klev/3)
      if(iknub(jl).gt.nktopg) iknub(jl)=nktopg
      if(iknub(jl).eq.nktopg) iknul(jl)=klev
      if(iknub(jl).eq.iknul(jl)) iknub(jl)=iknul(jl)-1
      endif
 2010 continue

      do 223 jk=klev,2,-1
      do 222 jl=kidia,kfdia
        zrho(jl,jk)=2.*paphm1(jl,jk)*zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
  222 continue
  223 continue
c     print *,'  dans orolift: 223'

C********************************************************************
C
c*     define low level flow
C      -------------------
      do 2115 jk=klev,1,-1
      do 2116 jl=kidia,kfdia
      if(ktest(jl).eq.1) THEN
      if(jk.ge.iknub(jl).and.jk.le.iknul(jl)) then
        pulow(JL)=pulow(JL)+pum1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
        pvlow(JL)=pvlow(JL)+pvm1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
        zrho(JL,klev+1)=zrho(JL,klev+1)
     *                 +zrho(JL,JK)*(paphm1(jl,jk+1)-paphm1(jl,jk))
      endif
      endif
 2116 continue
 2115 continue
      do 2110 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      pulow(JL)=pulow(JL)/(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
      pvlow(JL)=pvlow(JL)/(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
      zrho(JL,klev+1)=zrho(Jl,klev+1)
     *               /(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
      endif
 2110 continue


200   continue

C***********************************************************
C
C*         3.      COMPUTE MOUNTAIN LIFT
C
  300 continue
C
      do 301 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
       ztau(jl,klev+1)= - gklift*zrho(jl,klev+1)*2.*romega*
c    *                 (2*pstd(jl)+pmea(jl))*
     *                 2*pstd(jl)*
     *                 sin(rpi/180.*plat(jl))*pvlow(jl)
       ztav(jl,klev+1)=   gklift*zrho(jl,klev+1)*2.*romega*
c    *                 (2*pstd(jl)+pmea(jl))*
     *                 2*pstd(jl)*
     *                 sin(rpi/180.*plat(jl))*pulow(jl)
      else
       ztau(jl,klev+1)=0.0
       ztav(jl,klev+1)=0.0
      endif
301   continue

C
C*         4.      COMPUTE LIFT PROFILE         
C*                 --------------------   
C

  400 continue

      do 401 jk=1,klev
      do 401 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      ztau(jl,jk)=ztau(jl,klev+1)*paphm1(jl,jk)/paphm1(jl,klev+1)
      ztav(jl,jk)=ztav(jl,klev+1)*paphm1(jl,jk)/paphm1(jl,klev+1)
      else
      ztau(jl,jk)=0.0
      ztav(jl,jk)=0.0
      endif
401   continue
C
C
C*         5.      COMPUTE TENDENCIES.
C*                 -------------------
      if(lifthigh)then
C
  500 continue
C
C  EXPLICIT SOLUTION AT ALL LEVELS
C
      do 524 jk=1,klev
      do 523 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      zdelp=paphm1(jl,jk+1)-paphm1(jl,jk)
      zdudt(jl)=-rg*(ztau(jl,jk+1)-ztau(jl,jk))/zdelp
      zdvdt(jl)=-rg*(ztav(jl,jk+1)-ztav(jl,jk))/zdelp
      endif  
  523 continue
  524 continue
C
C  PROJECT PERPENDICULARLY TO U NOT TO DESTROY ENERGY
C
      do 530 jk=1,klev
      do 530 jl=kidia,kfdia
      if(ktest(jl).eq.1) then

        zslow=sqrt(pulow(jl)**2+pvlow(jl)**2)
        zsqua=amax1(sqrt(pum1(jl,jk)**2+pvm1(jl,jk)**2),gvsec)
        zscav=-zdudt(jl)*pvm1(jl,jk)+zdvdt(jl)*pum1(jl,jk)
        if(zsqua.gt.gvsec)then
          pvom(jl,jk)=-zscav*pvm1(jl,jk)/zsqua**2
          pvol(jl,jk)= zscav*pum1(jl,jk)/zsqua**2
        else
          pvom(jl,jk)=0.0
          pvol(jl,jk)=0.0      
        endif  
        zsqua=sqrt(pum1(jl,jk)**2+pum1(jl,jk)**2)               
        if(zsqua.lt.zslow)then
          pvom(jl,jk)=zsqua/zslow*pvom(jl,jk)
          pvol(jl,jk)=zsqua/zslow*pvol(jl,jk)
        endif 

      endif  
530   continue
C
C  6.  LOW LEVEL LIFT, SEMI IMPLICIT:
C  ----------------------------------

      else

        do 601 jl=kidia,kfdia
        if(ktest(jl).eq.1) then
          do jk=klev,iknub(jl),-1
          zbet=gklift*2.*romega*sin(rpi/180.*plat(jl))*ztmst*
     *        (pgeom1(jl,iknub(jl)-1)-pgeom1(jl,  jk))/
     *        (pgeom1(jl,iknub(jl)-1)-pgeom1(jl,klev))
          zdudt(jl)=-pum1(jl,jk)/ztmst/(1+zbet**2)
          zdvdt(jl)=-pvm1(jl,jk)/ztmst/(1+zbet**2)
          pvom(jl,jk)= zbet**2*zdudt(jl) - zbet   *zdvdt(jl)
          pvol(jl,jk)= zbet   *zdudt(jl) + zbet**2*zdvdt(jl)    
          enddo
        endif
 601    continue

      endif

c     print *,' out orolift'

      return
      end
      SUBROUTINE SUGWD_strato(NLON,NLEV,paprs,pplay)
C     
C
C**** *SUGWD* INITIALIZE COMMON YOEGWD CONTROLLING GRAVITY WAVE DRAG
C
C     PURPOSE.
C     --------
C           INITIALIZE YOEGWD, THE COMMON THAT CONTROLS THE
C           GRAVITY WAVE DRAG PARAMETRIZATION.
C    VERY IMPORTANT:
C    ______________
C           THIS ROUTINE SET_UP THE "TUNABLE PARAMETERS" OF THE
C           VARIOUS SSO SCHEMES
C
C**   INTERFACE.
C     ----------
C        CALL *SUGWD* FROM *SUPHEC*
C              -----        ------
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C        PAPRS,PPLAY : Pressure at semi and full model levels
C        NLEV        : number of model levels
c        NLON        : number of points treated in the physics
C
C        IMPLICIT ARGUMENTS :
C        --------------------
C        COMMON YOEGWD
C-GFRCRIT-R:  Critical Non-dimensional mountain Height
C             (HNC in (1),    LOTT 1999)
C-GKWAKE--R:  Bluff-body drag coefficient for low level wake
C             (Cd in (2),     LOTT 1999)
C-GRCRIT--R:  Critical Richardson Number 
C             (Ric, End of first column p791 of LOTT 1999) 
C-GKDRAG--R:  Gravity wave drag coefficient
C             (G in (3),      LOTT 1999)
C-GKLIFT--R:  Mountain Lift coefficient
C             (Cl in (4),     LOTT 1999)
C-GHMAX---R:  Not used
C-GRAHILO-R:  Set-up the trapped waves fraction
C             (Beta , End of first column,  LOTT 1999)
C
C-GSIGCR--R:  Security value for blocked flow depth
C-NKTOPG--I:  Security value for blocked flow level
C-nstra----I:  An estimate to qualify the upper levels of
C             the model where one wants to impose strees
C             profiles
C-GSSECC--R:  Security min value for low-level B-V frequency
C-GTSEC---R:  Security min value for anisotropy and GW stress.
C-GVSEC---R:  Security min value for ulow
C         
C
C     METHOD.
C     -------
C        SEE DOCUMENTATION
C
C     EXTERNALS.
C     ----------
C        NONE
C
C     REFERENCE.
C     ----------
C     Lott, 1999: Alleviation of stationary biases in a GCM through...
C                 Monthly Weather Review, 127, pp 788-801.
C
C     AUTHOR.
C     -------
C        FRANCOIS LOTT        *LMD*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 90-01-01 (MARTIN MILLER, ECMWF)
C        LAST:  99-07-09     (FRANCOIS LOTT,LMD)
C     ------------------------------------------------------------------
      USE dimphy
      USE mod_phys_lmdz_para
      USE mod_grid_phy_lmdz
      IMPLICIT NONE
C
C     -----------------------------------------------------------------
#include "YOEGWD.h"
C      ----------------------------------------------------------------
C
C  ARGUMENTS
      integer nlon,nlev
      REAL paprs(nlon,nlev+1)
      REAL pplay(nlon,nlev)
C
      INTEGER JK
      REAL ZPR,ZTOP,ZSIGT,ZPM1R
      REAL :: pplay_glo(klon_glo,nlev)
      REAL :: paprs_glo(klon_glo,nlev+1)

C
C*       1.    SET THE VALUES OF THE PARAMETERS
C              --------------------------------
C
 100  CONTINUE
C
      PRINT *,' DANS SUGWD NLEV=',NLEV
      GHMAX=10000.
C
      ZPR=100000.
      ZTOP=0.001 
      ZSIGT=0.94
cold  ZPR=80000.
cold  ZSIGT=0.85
C
      CALL gather(pplay,pplay_glo)
      CALL bcast(pplay_glo)
      CALL gather(paprs,paprs_glo)
      CALL bcast(paprs_glo)

      DO 110 JK=1,NLEV
      ZPM1R=pplay_glo(klon_glo/2,jk)/paprs_glo(klon_glo/2,1) 
      IF(ZPM1R.GE.ZSIGT)THEN
         nktopg=JK
      ENDIF
      ZPM1R=pplay_glo(klon_glo/2,jk)/paprs_glo(klon_glo/2,1) 
      IF(ZPM1R.GE.ZTOP)THEN
         nstra=JK
      ENDIF
  110 CONTINUE
c
c  inversion car dans orodrag on compte les niveaux a l'envers
      nktopg=nlev-nktopg+1
      nstra=nlev-nstra
      print *,' DANS SUGWD nktopg=', nktopg
      print *,' DANS SUGWD nstra=', nstra
C
      GSIGCR=0.80
C
      GKDRAG=0.1875
      GRAHILO=0.1   
      GRCRIT=1.00 
      GFRCRIT=1.00
      GKWAKE=0.50
C
      GKLIFT=0.25
      GVCRIT =0.1

      WRITE(UNIT=6,FMT='('' *** SSO essential constants ***'')')
      WRITE(UNIT=6,FMT='('' *** SPECIFIED IN SUGWD ***'')')
      WRITE(UNIT=6,FMT='('' Gravity wave ct '',E13.7,'' '')')GKDRAG
      WRITE(UNIT=6,FMT='('' Trapped/total wave dag '',E13.7,'' '')')
     S      GRAHILO
      WRITE(UNIT=6,FMT='('' Critical Richardson   = '',E13.7,'' '')')
     S                  GRCRIT
      WRITE(UNIT=6,FMT='('' Critical Froude'',e13.7)') GFRCRIT
      WRITE(UNIT=6,FMT='('' Low level Wake bluff cte'',e13.7)') GKWAKE
      WRITE(UNIT=6,FMT='('' Low level lift  cte'',e13.7)') GKLIFT

C
C
C      ----------------------------------------------------------------
C
C*       2.    SET VALUES OF SECURITY PARAMETERS
C              ---------------------------------
C
 200  CONTINUE
C
      GVSEC=0.10
      GSSEC=0.0001
C
      GTSEC=0.00001
C
      RETURN
      END

