!
! $Header$
!
      SUBROUTINE conccm (dtime,paprs,pplay,t,q,conv_q,
     s                   d_t, d_q, rain, snow, kbascm, ktopcm)
c
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: le 14 mars 1996
c Objet: Schema simple (avec flux de masse) pour la convection 
c        (schema standard du modele NCAR CCM2)
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "YOETHF.h"
c
c Entree:
      REAL dtime              ! pas d'integration
      REAL paprs(klon,klev+1) ! pression inter-couche (Pa)
      REAL pplay(klon,klev)   ! pression au milieu de couche (Pa)
      REAL t(klon,klev)       ! temperature (K)
      REAL q(klon,klev)       ! humidite specifique (g/g)
      REAL conv_q(klon,klev)  ! taux de convergence humidite (g/g/s)
c Sortie:
      REAL d_t(klon,klev)     ! incrementation temperature
      REAL d_q(klon,klev)     ! incrementation vapeur
      REAL rain(klon)         ! pluie (mm/s)
      REAL snow(klon)         ! neige (mm/s)
      INTEGER kbascm(klon)    ! niveau du bas de convection
      INTEGER ktopcm(klon)    ! niveau du haut de convection
c
      REAL pt(klon,klev)
      REAL pq(klon,klev)
      REAL pres(klon,klev)
      REAL dp(klon,klev)
      REAL zgeom(klon,klev)
      REAL cmfprs(klon)
      REAL cmfprt(klon)
      INTEGER ntop(klon)
      INTEGER nbas(klon)
      INTEGER i, k
      REAL zlvdcp, zlsdcp, zdelta, zz, za, zb
c
      LOGICAL usekuo ! utiliser convection profonde (schema Kuo)
      PARAMETER (usekuo=.TRUE.)
c
      REAL d_t_bis(klon,klev)
      REAL d_q_bis(klon,klev)
      REAL rain_bis(klon)
      REAL snow_bis(klon)
      INTEGER ibas_bis(klon)
      INTEGER itop_bis(klon)
      REAL d_ql_bis(klon,klev)
      REAL rneb_bis(klon,klev)
c
c initialiser les variables de sortie (pour securite)
      DO i = 1, klon
         rain(i) = 0.0
         snow(i) = 0.0
         kbascm(i) = 0
         ktopcm(i) = 0
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
      ENDDO
      ENDDO
c
c preparer les variables d'entree (attention: l'ordre des niveaux 
c verticaux augmente du haut vers le bas)
      DO k = 1, klev
      DO i = 1, klon
         pt(i,k) = t(i,klev-k+1)
         pq(i,k) = q(i,klev-k+1)
         pres(i,k) = pplay(i,klev-k+1)
         dp(i,k) = paprs(i,klev+1-k)-paprs(i,klev+1-k+1)
      ENDDO
      ENDDO
      DO i = 1, klon
         zgeom(i,klev) = RD * t(i,1) / (0.5*(paprs(i,1)+pplay(i,1)))
     .                      * (paprs(i,1)-pplay(i,1))
      ENDDO
      DO k = 2, klev
      DO i = 1, klon
         zgeom(i,klev+1-k) = zgeom(i,klev+1-k+1)
     .              + RD * 0.5*(t(i,k-1)+t(i,k)) / paprs(i,k)
     .                   * (pplay(i,k-1)-pplay(i,k))
      ENDDO
      ENDDO
c
      CALL cmfmca(dtime, pres, dp, zgeom, pt, pq,
     $                  cmfprt, cmfprs, ntop, nbas)
C
      DO k = 1, klev
      DO i = 1, klon
         d_q(i,klev+1-k) = pq(i,k) - q(i,klev+1-k) 
         d_t(i,klev+1-k) = pt(i,k) - t(i,klev+1-k)
      ENDDO
      ENDDO
c
      DO i = 1, klon
         rain(i) = cmfprt(i) * rhoh2o
         snow(i) = cmfprs(i) * rhoh2o
         kbascm(i) = klev+1 - nbas(i)
         ktopcm(i) = klev+1 - ntop(i)
      ENDDO
c
      IF (usekuo) THEN
      CALL conkuo(dtime, paprs, pplay, t, q, conv_q,
     s            d_t_bis, d_q_bis, d_ql_bis, rneb_bis,
     s            rain_bis, snow_bis, ibas_bis, itop_bis)
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = d_t(i,k) + d_t_bis(i,k)
         d_q(i,k) = d_q(i,k) + d_q_bis(i,k)
      ENDDO
      ENDDO
      DO i = 1, klon
         rain(i) = rain(i) + rain_bis(i)
         snow(i) = snow(i) + snow_bis(i)
         kbascm(i) = MIN(kbascm(i),ibas_bis(i))
         ktopcm(i) = MAX(ktopcm(i),itop_bis(i))
      ENDDO
      DO k = 1, klev ! eau liquide convective est
      DO i = 1, klon ! dispersee dans l'air
         zlvdcp=RLVTT/RCPD/(1.0+RVTMP2*q(i,k))
         zlsdcp=RLSTT/RCPD/(1.0+RVTMP2*q(i,k))
         zdelta = MAX(0.,SIGN(1.,RTT-t(i,k)))
         zz = d_ql_bis(i,k) ! re-evap. de l'eau liquide
         zb = MAX(0.0,zz)
         za = - MAX(0.0,zz) * (zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
         d_t(i,k) = d_t(i,k) + za
         d_q(i,k) = d_q(i,k) + zb
      ENDDO
      ENDDO
      ENDIF
c
      RETURN
      END
      SUBROUTINE cmfmca(deltat, p, dp, gz,
     $                  tb, shb,
     $                  cmfprt, cmfprs, cnt, cnb)
      USE dimphy
      IMPLICIT none
C-----------------------------------------------------------------------
C Moist convective mass flux procedure:
C If stratification is unstable to nonentraining parcel ascent,
C complete an adjustment making use of a simple cloud model
C 
C Code generalized to allow specification of parcel ("updraft")
C properties, as well as convective transport of an arbitrary
C number of passive constituents (see cmrb array).
C----------------------------Code History-------------------------------
C Original version:  J. J. Hack, March 22, 1990
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, G. Taylor, August 1992
c Adaptation au LMD: Z.X. Li, mars 1996 (reference: Hack 1994,
c                    J. Geophys. Res. vol 99, D3, 5551-5568). J'ai
c                    introduit les constantes et les fonctions thermo-
c                    dynamiques du Centre Europeen. J'ai elimine le
c                    re-indicage du code en esperant que cela pourra
c                    simplifier la lecture et la comprehension.
C-----------------------------------------------------------------------
cym#include "dimensions.h"
cym#include "dimphy.h"
      INTEGER pcnst ! nombre de traceurs passifs
      PARAMETER (pcnst=1)
C------------------------------Arguments--------------------------------
C Input arguments
C
      REAL deltat                 ! time step (seconds)
      REAL p(klon,klev)           ! pressure
      REAL dp(klon,klev)          ! delta-p
      REAL gz(klon,klev)          ! geopotential (a partir du sol)
c
      REAL thtap(klon)            ! PBL perturbation theta
      REAL shp(klon)              ! PBL perturbation specific humidity 
      REAL pblh(klon)             ! PBL height (provided by PBL routine)
      REAL cmrp(klon,pcnst)       ! constituent perturbations in PBL
c
c Updated arguments:
c
      REAL tb(klon,klev)         ! temperature (t bar)
      REAL shb(klon,klev)        ! specific humidity (sh bar)
      REAL cmrb(klon,klev,pcnst) ! constituent mixing ratios (cmr bar)
C
C Output arguments
C
      REAL cmfdt(klon,klev)    ! dT/dt due to moist convection
      REAL cmfdq(klon,klev)    ! dq/dt due to moist convection
      REAL cmfmc(klon,klev )   ! moist convection cloud mass flux
      REAL cmfdqr(klon,klev)   ! dq/dt due to convective rainout 
      REAL cmfsl(klon,klev )   ! convective lw static energy flux
      REAL cmflq(klon,klev )   ! convective total water flux
      REAL cmfprt(klon)        ! convective precipitation rate
      REAL cmfprs(klon)        ! convective snowfall rate
      REAL qc(klon,klev)       ! dq/dt due to rainout terms
      INTEGER cnt(klon)        ! top level of convective activity   
      INTEGER cnb(klon)        ! bottom level of convective activity
C------------------------------Parameters-------------------------------
      REAL c0         ! rain water autoconversion coefficient
      PARAMETER (c0=1.0e-4)
      REAL dzmin       ! minimum convective depth for precipitation
      PARAMETER (dzmin=0.0)
      REAL betamn      ! minimum overshoot parameter
      PARAMETER (betamn=0.10)
      REAL cmftau      ! characteristic adjustment time scale
      PARAMETER (cmftau=3600.)
      INTEGER limcnv   ! top interface level limit for convection
      PARAMETER (limcnv=1)
      REAL tpmax       ! maximum acceptable t perturbation (degrees C)
      PARAMETER (tpmax=1.50)
      REAL shpmax      ! maximum acceptable q perturbation (g/g)
      PARAMETER (shpmax=1.50e-3)
      REAL tiny        ! arbitrary small num used in transport estimates
      PARAMETER (tiny=1.0e-36)
      REAL eps         ! convergence criteria (machine dependent)
      PARAMETER (eps=1.0e-13)
      REAL tmelt       ! freezing point of water(req'd for rain vs snow)
      PARAMETER (tmelt=273.15)
      REAL ssfac ! supersaturation bound (detrained air)
      PARAMETER (ssfac=1.001)
C
C---------------------------Local workspace-----------------------------
      REAL gam(klon,klev)     ! L/cp (d(qsat)/dT)
      REAL sb(klon,klev)      ! dry static energy (s bar)
      REAL hb(klon,klev)      ! moist static energy (h bar)
      REAL shbs(klon,klev)    ! sat. specific humidity (sh bar star)
      REAL hbs(klon,klev)     ! sat. moist static energy (h bar star)
      REAL shbh(klon,klev+1)  ! specific humidity on interfaces
      REAL sbh(klon,klev+1)   ! s bar on interfaces
      REAL hbh(klon,klev+1)   ! h bar on interfaces
      REAL cmrh(klon,klev+1)  ! interface constituent mixing ratio 
      REAL prec(klon)         ! instantaneous total precipitation
      REAL dzcld(klon)        ! depth of convective layer (m)
      REAL beta(klon)         ! overshoot parameter (fraction)
      REAL betamx             ! local maximum on overshoot
      REAL eta(klon)          ! convective mass flux (kg/m^2 s)
      REAL etagdt             ! eta*grav*deltat
      REAL cldwtr(klon)       ! cloud water (mass)
      REAL rnwtr(klon)        ! rain water  (mass)
      REAL sc  (klon)         ! dry static energy   ("in-cloud")
      REAL shc (klon)         ! specific humidity   ("in-cloud")
      REAL hc  (klon)         ! moist static energy ("in-cloud")
      REAL cmrc(klon)         ! constituent mix rat ("in-cloud")
      REAL dq1(klon)          ! shb  convective change (lower lvl)
      REAL dq2(klon)          ! shb  convective change (mid level)
      REAL dq3(klon)          ! shb  convective change (upper lvl)
      REAL ds1(klon)          ! sb   convective change (lower lvl)
      REAL ds2(klon)          ! sb   convective change (mid level)
      REAL ds3(klon)          ! sb   convective change (upper lvl)
      REAL dcmr1(klon)        ! cmrb convective change (lower lvl)
      REAL dcmr2(klon)        ! cmrb convective change (mid level)
      REAL dcmr3(klon)        ! cmrb convective change (upper lvl)
      REAL flotab(klon)       ! hc - hbs (mesure d'instabilite)
      LOGICAL ldcum(klon)     ! .true. si la convection existe
      LOGICAL etagt0          ! true if eta > 0.0
      REAL dt                 ! current 2 delta-t (model time step)
      REAL cats     ! modified characteristic adj. time
      REAL rdt      ! 1./dt
      REAL qprime   ! modified specific humidity pert.
      REAL tprime   ! modified thermal perturbation
      REAL pblhgt   ! bounded pbl height (max[pblh,1m])
      REAL fac1     ! intermediate scratch variable
      REAL shprme   ! intermediate specific humidity pert.
      REAL qsattp   ! saturation mixing ratio for 
C                   !  thermally perturbed PBL parcels 
      REAL dz       ! local layer depth
      REAL b1       ! bouyancy measure in detrainment lvl
      REAL b2       ! bouyancy measure in condensation lvl
      REAL g     ! bounded vertical gradient of hb
      REAL tmass ! total mass available for convective exchange
      REAL denom ! intermediate scratch variable
      REAL qtest1! used in negative q test (middle lvl) 
      REAL qtest2! used in negative q test (lower lvl) 
      REAL fslkp ! flux lw static energy (bot interface)
      REAL fslkm ! flux lw static energy (top interface)
      REAL fqlkp ! flux total water (bottom interface)
      REAL fqlkm ! flux total water (top interface)
      REAL botflx! bottom constituent mixing ratio flux
      REAL topflx! top constituent mixing ratio flux
      REAL efac1 ! ratio cmrb to convectively induced change (bot lvl)
      REAL efac2 ! ratio cmrb to convectively induced change (mid lvl)
      REAL efac3 ! ratio cmrb to convectively induced change (top lvl)
c
      INTEGER i,k  ! indices horizontal et vertical
      INTEGER km1  ! k-1 (index offset)
      INTEGER kp1  ! k+1 (index offset)
      INTEGER m    ! constituent index
      INTEGER ktp  ! temporary index used to track top 
      INTEGER is   ! nombre de points a ajuster
C
      REAL tmp1, tmp2, tmp3, tmp4
      REAL zx_t, zx_p, zx_q, zx_qs, zx_gam
      REAL zcor, zdelta, zcvm5
C
      REAL qhalf, sh1, sh2, shbs1, shbs2
#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
      qhalf(sh1,sh2,shbs1,shbs2) = MIN(MAX(sh1,sh2),
     $                            (shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))
C
C-----------------------------------------------------------------------
c pas de traceur pour l'instant
      DO m = 1, pcnst
      DO k = 1, klev
      DO i = 1, klon
         cmrb(i,k,m) = 0.0
      ENDDO
      ENDDO
      ENDDO
c
c Les perturbations de la couche limite sont zero pour l'instant
c
      DO m = 1, pcnst
      DO i = 1, klon
         cmrp(i,m) = 0.0
      ENDDO
      ENDDO
      DO i = 1, klon
         thtap(i) = 0.0
         shp(i) = 0.0
         pblh(i) = 1.0
      ENDDO
C
C Ensure that characteristic adjustment time scale (cmftau) assumed
C in estimate of eta isn't smaller than model time scale (deltat)
C
      dt   = deltat
      cats = MAX(dt,cmftau)
      rdt  = 1.0/dt
C
C Compute sb,hb,shbs,hbs
C
      DO k = 1, klev
      DO i = 1, klon
         zx_t = tb(i,k)
         zx_p = p(i,k)
         zx_q = shb(i,k)
           zdelta=MAX(0.,SIGN(1.,RTT-zx_t))
           zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
           zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*zx_q)
           zx_qs= r2es * FOEEW(zx_t,zdelta)/zx_p
           zx_qs=MIN(0.5,zx_qs)
           zcor=1./(1.-retv*zx_qs)
           zx_qs=zx_qs*zcor
           zx_gam = FOEDE(zx_t,zdelta,zcvm5,zx_qs,zcor)
         shbs(i,k) = zx_qs
         gam(i,k) = zx_gam
      ENDDO
      ENDDO
C
      DO k=limcnv,klev
         DO i=1,klon
            sb (i,k) = RCPD*tb(i,k) + gz(i,k)
            hb (i,k) = sb(i,k) + RLVTT*shb(i,k)
            hbs(i,k) = sb(i,k) + RLVTT*shbs(i,k)
         ENDDO
      ENDDO
C
C Compute sbh, shbh
C
      DO k=limcnv+1,klev
         km1 = k - 1
         DO i=1,klon
            sbh (i,k) =0.5*(sb(i,km1) + sb(i,k))
            shbh(i,k) =qhalf(shb(i,km1),shb(i,k),shbs(i,km1),shbs(i,k))
            hbh (i,k) =sbh(i,k) + RLVTT*shbh(i,k)
         ENDDO
      ENDDO
C
C Specify properties at top of model (not used, but filling anyway)
C
      DO i=1,klon
         sbh (i,limcnv) = sb(i,limcnv)
         shbh(i,limcnv) = shb(i,limcnv)
         hbh (i,limcnv) = hb(i,limcnv)
      ENDDO
C
C Zero vertically independent control, tendency & diagnostic arrays
C
      DO i=1,klon
         prec(i)  = 0.0
         dzcld(i) = 0.0
         cnb(i)   = 0
         cnt(i)   = klev+1
      ENDDO

      DO k = 1, klev
        DO i = 1,klon
         cmfdt(i,k)  = 0.
         cmfdq(i,k)  = 0.
         cmfdqr(i,k) = 0.
         cmfmc(i,k)  = 0.
         cmfsl(i,k)  = 0.
         cmflq(i,k)  = 0.
        ENDDO
      ENDDO
C
C Begin moist convective mass flux adjustment procedure.
C Formalism ensures that negative cloud liquid water can never occur
C
      DO 70 k=klev-1,limcnv+1,-1
         km1 = k - 1
         kp1 = k + 1
         DO 10 i=1,klon
            eta   (i) = 0.0
            beta  (i) = 0.0
            ds1   (i) = 0.0
            ds2   (i) = 0.0
            ds3   (i) = 0.0
            dq1   (i) = 0.0
            dq2   (i) = 0.0
            dq3   (i) = 0.0
C
C Specification of "cloud base" conditions
C
            qprime    = 0.0
            tprime    = 0.0
C
C Assign tprime within the PBL to be proportional to the quantity
C thtap (which will be bounded by tpmax), passed to this routine by 
C the PBL routine.  Don't allow perturbation to produce a dry 
C adiabatically unstable parcel.  Assign qprime within the PBL to be 
C an appropriately modified value of the quantity shp (which will be 
C bounded by shpmax) passed to this routine by the PBL routine.  The 
C quantity qprime should be less than the local saturation value 
C (qsattp=qsat[t+tprime,p]).  In both cases, thtap and shp are
C linearly reduced toward zero as the PBL top is approached.
C
            pblhgt = MAX(pblh(i),1.0)
            IF (gz(i,kp1)/RG.LE.pblhgt .AND. dzcld(i).EQ.0.0) THEN
               fac1   = MAX(0.0,1.0-gz(i,kp1)/RG/pblhgt)
               tprime = MIN(thtap(i),tpmax)*fac1
               qsattp = shbs(i,kp1) + RCPD/RLVTT*gam(i,kp1)*tprime
               shprme = MIN(MIN(shp(i),shpmax)*fac1,
     $                        MAX(qsattp-shb(i,kp1),0.0))
               qprime = MAX(qprime,shprme)
            ELSE
               tprime = 0.0
               qprime = 0.0
            ENDIF
C
C Specify "updraft" (in-cloud) thermodynamic properties
C
            sc (i)    = sb (i,kp1) + RCPD*tprime
            shc(i)    = shb(i,kp1) + qprime
            hc (i)    = sc (i    ) + RLVTT*shc(i)
            flotab(i) = hc(i) - hbs(i,k)
            dz        = dp(i,k)*RD*tb(i,k)/RG/p(i,k)
            IF (flotab(i).gt.0.0) THEN
               dzcld(i) = dzcld(i) + dz
            ELSE
               dzcld(i) = 0.0
            ENDIF
   10    CONTINUE
C
C Check on moist convective instability
C
         is = 0
         DO i = 1, klon
            IF (flotab(i).GT.0.0) THEN
               ldcum(i) = .TRUE.
               is = is + 1
            ELSE
               ldcum(i) = .FALSE.
            ENDIF
         ENDDO
C
         IF (is.EQ.0) THEN
            DO i=1,klon
               dzcld(i) = 0.0
            ENDDO
            GOTO 70
         ENDIF
C
C Current level just below top level => no overshoot
C
         IF (k.le.limcnv+1) THEN
            DO i=1,klon
            IF (ldcum(i)) THEN
               cldwtr(i) = sb(i,k)-sc(i)+flotab(i)/(1.0+gam(i,k))
               cldwtr(i) = MAX(0.0,cldwtr(i))
               beta(i)   = 0.0
            ENDIF
            ENDDO
            GOTO 20
         ENDIF
C
C First guess at overshoot parameter using crude buoyancy closure
C 10% overshoot assumed as a minimum and 1-c0*dz maximum to start
C If pre-existing supersaturation in detrainment layer, beta=0
C cldwtr is temporarily equal to RLVTT*l (l=> liquid water)
C
         DO i=1,klon
         IF (ldcum(i)) THEN
            cldwtr(i) = sb(i,k)-sc(i)+flotab(i)/(1.0+gam(i,k))
            cldwtr(i) = MAX(0.0,cldwtr(i))
            betamx = 1.0 - c0*MAX(0.0,(dzcld(i)-dzmin))
            b1        = (hc(i) - hbs(i,km1))*dp(i,km1)
            b2        = (hc(i) - hbs(i,k  ))*dp(i,k  )
            beta(i)   = MAX(betamn,MIN(betamx, 1.0+b1/b2))
            IF (hbs(i,km1).le.hb(i,km1)) beta(i) = 0.0
         ENDIF
         ENDDO
C
C Bound maximum beta to ensure physically realistic solutions
C
C First check constrains beta so that eta remains positive
C (assuming that eta is already positive for beta equal zero)
c La premiere contrainte de beta est que le flux eta doit etre positif.
C
         DO i=1,klon
         IF (ldcum(i)) THEN
            tmp1 = (1.0+gam(i,k))*(sc(i)-sbh(i,kp1) + cldwtr(i))
     $            - (hbh(i,kp1)-hc(i))*dp(i,k)/dp(i,kp1)
            tmp2 = (1.0+gam(i,k))*(sc(i)-sbh(i,k))
            IF ((beta(i)*tmp2-tmp1).GT.0.0) THEN
               betamx = 0.99*(tmp1/tmp2)
               beta(i) = MAX(0.0,MIN(betamx,beta(i)))
            ENDIF
C
C Second check involves supersaturation of "detrainment layer"
C small amount of supersaturation acceptable (by ssfac factor)
c La 2e contrainte est que la convection ne doit pas sursaturer
c la "detrainment layer", Neanmoins, une petite sursaturation
c est acceptee (facteur ssfac).
C
            IF (hb(i,km1).lt.hbs(i,km1)) THEN
               tmp1 = (1.0+gam(i,k))*(sc(i)-sbh(i,kp1) + cldwtr(i))
     $               - (hbh(i,kp1)-hc(i))*dp(i,k)/dp(i,kp1)
               tmp1 = tmp1/dp(i,k)
               tmp2 = gam(i,km1)*(sbh(i,k)-sc(i) + cldwtr(i)) -
     $                 hbh(i,k) + hc(i) - sc(i) + sbh(i,k)
               tmp3 = (1.0+gam(i,k))*(sc(i)-sbh(i,k))/dp(i,k)
               tmp4 = (dt/cats)*(hc(i)-hbs(i,k))*tmp2
     $               / (dp(i,km1)*(hbs(i,km1)-hb(i,km1))) + tmp3
               IF ((beta(i)*tmp4-tmp1).GT.0.0) THEN
                  betamx = ssfac*(tmp1/tmp4)
                  beta(i)   = MAX(0.0,MIN(betamx,beta(i)))
               ENDIF
            ELSE 
               beta(i) = 0.0
            ENDIF
C
C Third check to avoid introducing 2 delta x thermodynamic
C noise in the vertical ... constrain adjusted h (or theta e)
C so that the adjustment doesn't contribute to "kinks" in h
C
            g = MIN(0.0,hb(i,k)-hb(i,km1))
            tmp3 = (hb(i,k)-hb(i,km1)-g)*(cats/dt) / (hc(i)-hbs(i,k))
            tmp1 = (1.0+gam(i,k))*(sc(i)-sbh(i,kp1) + cldwtr(i))
     $            - (hbh(i,kp1)-hc(i))*dp(i,k)/dp(i,kp1)
            tmp1 = tmp1/dp(i,k)
            tmp1 = tmp3*tmp1 + (hc(i) - hbh(i,kp1))/dp(i,k)
            tmp2 = tmp3*(1.0+gam(i,k))*(sc(i)-sbh(i,k))/dp(i,k)
     $            + (hc(i)-hbh(i,k)-cldwtr(i))
     $             *(1.0/dp(i,k)+1.0/dp(i,kp1))
            IF ((beta(i)*tmp2-tmp1).GT.0.0) THEN
               betamx = 0.0
               IF (tmp2.NE.0.0) betamx = tmp1/tmp2
               beta(i) = MAX(0.0,MIN(betamx,beta(i)))
            ENDIF
         ENDIF
         ENDDO
C
C Calculate mass flux required for stabilization.
C
C Ensure that the convective mass flux, eta, is positive by
C setting negative values of eta to zero..
C Ensure that estimated mass flux cannot move more than the
C minimum of total mass contained in either layer k or layer k+1.
C Also test for other pathological cases that result in non-
C physical states and adjust eta accordingly.
C
   20    CONTINUE
         DO i=1,klon
         IF (ldcum(i)) THEN
            beta(i) = MAX(0.0,beta(i))
            tmp1 = hc(i) - hbs(i,k)
            tmp2 = ((1.0+gam(i,k))*(sc(i)-sbh(i,kp1)+cldwtr(i)) -
     $               beta(i)*(1.0+gam(i,k))*(sc(i)-sbh(i,k)))/dp(i,k) -
     $              (hbh(i,kp1)-hc(i))/dp(i,kp1)
            eta(i) = tmp1/(tmp2*RG*cats)
            tmass = MIN(dp(i,k),dp(i,kp1))/RG
            IF (eta(i).GT.tmass*rdt .OR. eta(i).LE.0.0) eta(i) = 0.0
C
C Check on negative q in top layer (bound beta)
C
            IF(shc(i)-shbh(i,k).LT.0.0 .AND. beta(i)*eta(i).NE.0.0)THEN
               denom = eta(i)*RG*dt*(shc(i) - shbh(i,k))/dp(i,km1)
               beta(i) = MAX(0.0,MIN(-0.999*shb(i,km1)/denom,beta(i)))
            ENDIF
C
C Check on negative q in middle layer (zero eta)
C
            qtest1 = shb(i,k) + eta(i)*RG*dt*((shc(i) - shbh(i,kp1)) -
     $               (1.0 - beta(i))*cldwtr(i)/RLVTT -
     $               beta(i)*(shc(i) - shbh(i,k)))/dp(i,k)
            IF (qtest1.le.0.0) eta(i) = 0.0
C
C Check on negative q in lower layer (bound eta)
C
            fac1 = -(shbh(i,kp1) - shc(i))/dp(i,kp1)
            qtest2 = shb(i,kp1) - eta(i)*RG*dt*fac1
            IF (qtest2 .lt. 0.0) THEN
               eta(i) = 0.99*shb(i,kp1)/(RG*dt*fac1)
            ENDIF
         ENDIF
         ENDDO
C
C
C Calculate cloud water, rain water, and thermodynamic changes
C
         DO 30 i=1,klon
         IF (ldcum(i)) THEN
            etagdt = eta(i)*RG*dt
            cldwtr(i) = etagdt*cldwtr(i)/RLVTT/RG
            rnwtr(i) = (1.0 - beta(i))*cldwtr(i)
            ds1(i) = etagdt*(sbh(i,kp1) - sc(i))/dp(i,kp1)
            dq1(i) = etagdt*(shbh(i,kp1) - shc(i))/dp(i,kp1)
            ds2(i) = (etagdt*(sc(i) - sbh(i,kp1)) +
     $                RLVTT*RG*cldwtr(i) - beta(i)*etagdt*
     $                (sc(i) - sbh(i,k)))/dp(i,k)
            dq2(i) = (etagdt*(shc(i) - shbh(i,kp1)) -
     $                RG*rnwtr(i) - beta(i)*etagdt*
     $                (shc(i) - shbh(i,k)))/dp(i,k)
            ds3(i) = beta(i)*(etagdt*(sc(i) - sbh(i,k)) -
     $               RLVTT*RG*cldwtr(i))/dp(i,km1)
            dq3(i) = beta(i)*etagdt*(shc(i) - shbh(i,k))/dp(i,km1)
C
C Isolate convective fluxes for later diagnostics
C
            fslkp = eta(i)*(sc(i) - sbh(i,kp1))
            fslkm = beta(i)*(eta(i)*(sc(i) - sbh(i,k)) -
     $                       RLVTT*cldwtr(i)*rdt)
            fqlkp = eta(i)*(shc(i) - shbh(i,kp1))
            fqlkm = beta(i)*eta(i)*(shc(i) - shbh(i,k))
C
C
C Update thermodynamic profile (update sb, hb, & hbs later)
C
            tb (i,kp1) = tb(i,kp1)  + ds1(i) / RCPD
            tb (i,k  ) = tb(i,k  )  + ds2(i) / RCPD
            tb (i,km1) = tb(i,km1)  + ds3(i) / RCPD
            shb(i,kp1) = shb(i,kp1) + dq1(i)
            shb(i,k  ) = shb(i,k  ) + dq2(i)
            shb(i,km1) = shb(i,km1) + dq3(i)
            prec(i)    = prec(i)    + rnwtr(i)/rhoh2o
C
C Update diagnostic information for final budget
C Tracking temperature & specific humidity tendencies,
C rainout term, convective mass flux, convective liquid
C water static energy flux, and convective total water flux
C
            cmfdt (i,kp1) = cmfdt (i,kp1) + ds1(i)/RCPD*rdt
            cmfdt (i,k  ) = cmfdt (i,k  ) + ds2(i)/RCPD*rdt
            cmfdt (i,km1) = cmfdt (i,km1) + ds3(i)/RCPD*rdt
            cmfdq (i,kp1) = cmfdq (i,kp1) + dq1(i)*rdt
            cmfdq (i,k  ) = cmfdq (i,k  ) + dq2(i)*rdt
            cmfdq (i,km1) = cmfdq (i,km1) + dq3(i)*rdt
            cmfdqr(i,k  ) = cmfdqr(i,k  ) + (RG*rnwtr(i)/dp(i,k))*rdt
            cmfmc (i,kp1) = cmfmc (i,kp1) + eta(i)
            cmfmc (i,k  ) = cmfmc (i,k  ) + beta(i)*eta(i)
            cmfsl (i,kp1) = cmfsl (i,kp1) + fslkp
            cmfsl (i,k  ) = cmfsl (i,k  ) + fslkm
            cmflq (i,kp1) = cmflq (i,kp1) + RLVTT*fqlkp
            cmflq (i,k  ) = cmflq (i,k  ) + RLVTT*fqlkm
            qc    (i,k  ) =                (RG*rnwtr(i)/dp(i,k))*rdt
         ENDIF
   30    CONTINUE
C
C Next, convectively modify passive constituents
C
         DO 50 m=1,pcnst
         DO 40 i=1,klon
         IF (ldcum(i)) THEN
C
C If any of the reported values of the constituent is negative in
C the three adjacent levels, nothing will be done to the profile
C
            IF ((cmrb(i,kp1,m).LT.0.0) .OR.
     $          (cmrb(i,k,m).LT.0.0) .OR.
     $          (cmrb(i,km1,m).LT.0.0)) GOTO 40
C
C Specify constituent interface values (linear interpolation)
C
            cmrh(i,k  ) = 0.5*(cmrb(i,km1,m) + cmrb(i,k  ,m))
            cmrh(i,kp1) = 0.5*(cmrb(i,k  ,m) + cmrb(i,kp1,m))
C
C Specify perturbation properties of constituents in PBL
C
            pblhgt = MAX(pblh(i),1.0)
            IF (gz(i,kp1)/RG.LE.pblhgt .AND. dzcld(i).EQ.0.) THEN
               fac1 = MAX(0.0,1.0-gz(i,kp1)/RG/pblhgt)
               cmrc(i) = cmrb(i,kp1,m) + cmrp(i,m)*fac1
            ELSE
               cmrc(i) = cmrb(i,kp1,m)
            ENDIF
C
C Determine fluxes, flux divergence => changes due to convection
C Logic must be included to avoid producing negative values. A bit
C messy since there are no a priori assumptions about profiles.
C Tendency is modified (reduced) when pending disaster detected.
C
            etagdt = eta(i)*RG*dt
            botflx   = etagdt*(cmrc(i) - cmrh(i,kp1))
            topflx   = beta(i)*etagdt*(cmrc(i)-cmrh(i,k))
            dcmr1(i) = -botflx/dp(i,kp1)
            efac1    = 1.0
            efac2    = 1.0
            efac3    = 1.0
C
            IF (cmrb(i,kp1,m)+dcmr1(i) .LT. 0.0) THEN
               efac1 = MAX(tiny,ABS(cmrb(i,kp1,m)/dcmr1(i)) - eps)
            ENDIF
C
            IF (efac1.EQ.tiny .OR. efac1.GT.1.0) efac1 = 0.0
            dcmr1(i) = -efac1*botflx/dp(i,kp1)
            dcmr2(i) = (efac1*botflx - topflx)/dp(i,k)
C  
            IF (cmrb(i,k,m)+dcmr2(i) .LT. 0.0) THEN
               efac2 = MAX(tiny,ABS(cmrb(i,k  ,m)/dcmr2(i)) - eps)
            ENDIF
C
            IF (efac2.EQ.tiny .OR. efac2.GT.1.0) efac2 = 0.0
            dcmr2(i) = (efac1*botflx - efac2*topflx)/dp(i,k)
            dcmr3(i) = efac2*topflx/dp(i,km1)
C
            IF (cmrb(i,km1,m)+dcmr3(i) .LT. 0.0) THEN
               efac3 = MAX(tiny,ABS(cmrb(i,km1,m)/dcmr3(i)) - eps)
            ENDIF
C
            IF (efac3.EQ.tiny .OR. efac3.GT.1.0) efac3 = 0.0
            efac3    = MIN(efac2,efac3)
            dcmr2(i) = (efac1*botflx - efac3*topflx)/dp(i,k)
            dcmr3(i) = efac3*topflx/dp(i,km1)
C
            cmrb(i,kp1,m) = cmrb(i,kp1,m) + dcmr1(i)
            cmrb(i,k  ,m) = cmrb(i,k  ,m) + dcmr2(i)
            cmrb(i,km1,m) = cmrb(i,km1,m) + dcmr3(i)
         ENDIF
   40    CONTINUE
   50    CONTINUE              ! end of m=1,pcnst loop
C
         IF (k.EQ.limcnv+1) GOTO 60 ! on ne pourra plus glisser
c
c Dans la procedure de glissage ascendant, les variables thermo-
c dynamiques des couches k et km1 servent au calcul des couches
c superieures. Elles ont donc besoin d'une mise-a-jour.
C
         DO i = 1, klon
         IF (ldcum(i)) THEN
            zx_t = tb(i,k)
            zx_p = p(i,k)
            zx_q = shb(i,k)
              zdelta=MAX(0.,SIGN(1.,RTT-zx_t))
              zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
              zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*zx_q)
              zx_qs= r2es * FOEEW(zx_t,zdelta)/zx_p
              zx_qs=MIN(0.5,zx_qs)
              zcor=1./(1.-retv*zx_qs)
              zx_qs=zx_qs*zcor
              zx_gam = FOEDE(zx_t,zdelta,zcvm5,zx_qs,zcor)
            shbs(i,k) = zx_qs
            gam(i,k) = zx_gam
c
            zx_t = tb(i,km1)
            zx_p = p(i,km1)
            zx_q = shb(i,km1)
              zdelta=MAX(0.,SIGN(1.,RTT-zx_t))
              zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
              zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*zx_q)
              zx_qs= r2es * FOEEW(zx_t,zdelta)/zx_p
              zx_qs=MIN(0.5,zx_qs)
              zcor=1./(1.-retv*zx_qs)
              zx_qs=zx_qs*zcor
              zx_gam = FOEDE(zx_t,zdelta,zcvm5,zx_qs,zcor)
            shbs(i,km1) = zx_qs
            gam(i,km1) = zx_gam
C
            sb (i,k  ) = sb(i,k  ) + ds2(i)
            sb (i,km1) = sb(i,km1) + ds3(i)
            hb (i,k  ) = sb(i,k  ) + RLVTT*shb(i,k)
            hb (i,km1) = sb(i,km1) + RLVTT*shb(i,km1)
            hbs(i,k  ) = sb(i,k  ) + RLVTT*shbs(i,k  )
            hbs(i,km1) = sb(i,km1) + RLVTT*shbs(i,km1)
C
            sbh (i,k) = 0.5*(sb(i,k) + sb(i,km1))
            shbh(i,k) = qhalf(shb(i,km1),shb(i,k)
     $                       ,shbs(i,km1),shbs(i,k))
            hbh (i,k) = sbh(i,k) + RLVTT*shbh(i,k)
            sbh (i,km1) = 0.5*(sb(i,km1) + sb(i,k-2))
            shbh(i,km1) = qhalf(shb(i,k-2),shb(i,km1),
     $                    shbs(i,k-2),shbs(i,km1))
            hbh (i,km1) = sbh(i,km1) + RLVTT*shbh(i,km1)
         ENDIF
         ENDDO
C
C Ensure that dzcld is reset if convective mass flux zero
C specify the current vertical extent of the convective activity
C top of convective layer determined by size of overshoot param.
C
   60    CONTINUE
         DO i=1,klon
            etagt0 = eta(i).gt.0.0
            IF (.not.etagt0) dzcld(i) = 0.0
            IF (etagt0 .and. beta(i).gt.betamn) THEN
               ktp = km1
            ELSE
               ktp = k
            ENDIF
            IF (etagt0) THEN
               cnt(i) = MIN(cnt(i),ktp)
               cnb(i) = MAX(cnb(i),k)
            ENDIF
         ENDDO
   70 CONTINUE        ! end of k loop
C
C determine whether precipitation, prec, is frozen (snow) or not
C
      DO i=1,klon
         IF (tb(i,klev).LT.tmelt .AND. tb(i,klev-1).lt.tmelt) THEN
             cmfprs(i) = prec(i)*rdt
         ELSE
             cmfprt(i) = prec(i)*rdt
         ENDIF
      ENDDO
C
      RETURN  ! we're all done ... return to calling procedure
      END
