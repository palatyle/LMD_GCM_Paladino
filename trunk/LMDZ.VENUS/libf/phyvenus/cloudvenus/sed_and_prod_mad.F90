MODULE SED_AND_PROD_MAD

IMPLICIT NONE

!
! All parameter variable should be modified regarding the sudied system.
!
!
! All vectors arguments of the subroutines defined here should be defined from TOP of
! the atmosphere to the Ground.
!
! There should be one more vertical level than layers. Layer thickness vector should be
! defined on the number of layer. The last item dzlay(nla), where nla is the number of layer,
! is arbitrary and can be set to dzlay(nla-1). Usually dzlay is computed from zlay:
!
!    dzlay(1:nla-1) = zlay(1:nla-1)-zlay(2:nla)
!    dzlay(nla) = dzlay(nla-1)
!
! Level thickness (dzlev) is also defined on the number of layer and can be computed as follow:
!
!    dzlev(1:nla) = zlev(1:nla)-zlev(2:nle)
!
! Where nle is number of vertical level (i.e. nla+1).
!
!
! eta_g and lambda_g functions should be modified regarding the studied atmosphere. The current parameters
! stands for an N2 atmosphere.
!
! eta_g    --> air viscosity at given temperature.
! lambda_g --> air mean free path at given temperature and pressure.
!
! By J.Burgalat
!
! ft = Theoretical Settling velocity at each vertical levels with Fiadero correction

CONTAINS

!============================================================================
! MOMENT RELATED METHODS
!============================================================================
FUNCTION get_r0_1(M0,M3) RESULT(ret)

  use free_param, only: sigma1

  IMPLICIT NONE

  !! Compute the mean radius of the size distribution of the first mode.
  !!
  !! $$ r_{0,1} = \left(\dfrac{M_{3}}{M_{0}}/ \alpha_{1}(3) \right)^{1/3} $$
  REAL, INTENT(in) :: M0
  !! 0th order moment of the size-distribution
  REAL, INTENT(in) :: M3
  !! 3rd order moment of the size-distribution
  REAL             :: ret, alpha_k
  !! Mean radius of the mode 1
  ret = (m3/m0 / alpha_k(3.D0,sigma1))**0.333333D0 !SG: ok+
END FUNCTION get_r0_1


!*****************************************************************************
FUNCTION get_r0_2(M0,M3) RESULT(ret)

  use free_param, only: sigma2

  IMPLICIT NONE

  !! Compute the mean radius of the size distribution of the first mode.
  !!
  !! $$ r_{0,2} = \left(\dfrac{M_{3}}{M_{0}}/ \alpha_{2}(3) \right)^{1/3} $$
  REAL, INTENT(in) :: M0
  !! 0th order moment of the size-distribution
  REAL, INTENT(in) :: M3
  !! 3rd order moment of the size-distribution
  REAL             :: ret, alpha_k
  !! Mean radius of the mode 2
  ret = (m3/m0 / alpha_k(3.D0,sigma2))**0.333333D0
END FUNCTION get_r0_2 !SG: ok ++


!============================================================================
! ATMOSPHERE RELATED METHODS
!============================================================================
FUNCTION effg(z) RESULT(ret)
  !! Compute effective gravitational acceleration.

  use donnees, only: g0, rpla

  IMPLICIT NONE

  REAL, INTENT(in) :: z !! Altitude in meters
  REAL :: ret           !! Effective gravitational acceleration in \(m.s^{-2}\)
  ret = g0 * (rpla/(rpla+z))**2
  RETURN
END FUNCTION effg


!============================================================================
!                 PRODUCTION PROCESS RELATED METHOD
!============================================================================
SUBROUTINE aer_production(dt,nlay,plev,zlev,zlay,dm0_1,dm3_1)

  use free_param, only: sigz, p_aer, r_aer, rho_aer, sig_aer, tx_prod
  use donnees, only: pi

  IMPLICIT NONE

  !! Compute the production of aerosols moments.
  !!
  !! The method computes the tendencies of M0 and M3 of the CCN precursor.
  !! Production values are distributed along a normal law in altitude, centered  at
  !! p_prod pressure level with a fixed sigma of 20km.
  !!
  !! First M3 tendency is computed and M0 is retrieved using the inter-moments relation
  !! of the first mode of drop size-distribution with a fixed mean radius (r_aer).
  REAL, INTENT(in) :: dt
  INTEGER, INTENT(in) :: nlay
  REAL, INTENT(in),  DIMENSION(nlay+1) :: plev    !! Pressure at each vertical level (Pa).
  REAL, INTENT(in),  DIMENSION(nlay+1) :: zlev    !! Altitude of each vertical level (m).
  REAL, INTENT(in),  DIMENSION(nlay)   :: zlay    !! Altitude at the center of each vertical layer (m).
  REAL, INTENT(out), DIMENSION(nlay)   :: dm0_1   !! Tendency of M0 (\(m^{-3}\)).
  REAL, INTENT(out), DIMENSION(nlay)   :: dm3_1   !! Tendency of M3 (\(m^{3}.m^{-3}\)).
  INTEGER               :: i,nla,nle
  REAL, DIMENSION(nlay) :: dzlev
  REAL                  :: zprod,cprod, alpha_k

  REAL, PARAMETER :: fnorm = 1.D0/(dsqrt(2.D0*pi)*sigz)
  REAL, PARAMETER :: znorm = dsqrt(2.D0)*sigz

  nle = SIZE(zlev) 
  nla = nlay+1
  dzlev = zlev(1:nle-1)-zlev(2:nle) !distance between two levels, in meter
  zprod = -1.0D0 ! initialization
  
  ! locate production altitude
  DO i=1, nla-1
     IF (plev(i) .lt. p_aer.AND.plev(i+1) .ge. p_aer) THEN
        zprod = zlay(i) ; EXIT ! the altitude of production is determined here
     ENDIF
  ENDDO
  IF (zprod < 0.D0) THEN ! test
     WRITE(*,'(a)') "In aer_prod, cannot find production altitude"
     call EXIT(11)
  ENDIF

  dm3_1(:)= tx_prod *0.75D0/pi *dt / rho_aer / 2.D0 / dzlev(1:nla) * &
       (erf((zlev(1:nla)-zprod)/znorm) - &
       erf((zlev(2:nla+1)-zprod)/znorm))

  dm0_1(:) = dm3_1(:)/(r_aer**3*alpha_k(3.D0,sig_aer))

END SUBROUTINE aer_production


!============================================================================
! SEDIMENTATION PROCESS RELATED METHODS
!============================================================================
SUBROUTINE let_me_fall_in_peace(dt,nlay,dzlay,mk,ft,dmk)
  !! Compute the tendency of the moment through sedimentation process.
  !!
  !!
  !! The method computes the time evolution of the \(k^{th}\) order moment through sedimentation:
  !!
  !! $$ \dfrac{dM_{k}}{dt} = \dfrac{\Phi_{k}}{dz} $$
  !!
  !! The equation is resolved using a [Crank-Nicolson algorithm](http://en.wikipedia.org/wiki/Crank-Nicolson_method).
  !!
  !! Sedimentation algorithm is quite messy. It appeals to the dark side of the Force and uses evil black magic spells
  !! from ancient times. It is based on \cite{toon1988b,fiadeiro1977,turco1979} and is an update of the algorithm
  !! originally implemented in the LMDZ-Titan 2D GCM.

  IMPLICIT NONE

  REAL, INTENT(in)  :: dt
  INTEGER, INTENT(in) :: nlay
  REAL, INTENT(in), DIMENSION(nlay)  :: dzlay !! Atmospheric thickness between the center of 2 adjacent layers (\(m\))
  REAL, INTENT(in), DIMENSION(nlay)  :: mk    !! \(k^{th}\) order moment to sediment (in \(m^{k}\)).
  REAL, INTENT(in), DIMENSION(nlay+1):: ft    !! Downward sedimentation flux  (effective velocity of the moment).
  REAL, INTENT(out), DIMENSION(nlay) :: dmk   !! Tendency of \(k^{th}\) order moment (in \(m^{k}.m^{-3}\)).
  INTEGER                            :: i,nla,nle
  REAL, DIMENSION(nlay) :: as,bs,cs,mko

  nla = SIZE(mk,DIM=1) ; nle = nla+1

  mko(1:nla) = 0.D0 !initialization
  cs(1:nla) = ft(2:nla+1) - dzlay(1:nla)/dt  
  IF (ANY(cs > 0.D0)) THEN
     ! implicit case
     as(1:nla) = ft(1:nla)
     bs(1:nla) = -(ft(2:nle)+dzlay(1:nla)/dt)
     cs(1:nla) = -dzlay(1:nla)/dt*mk(1:nla)
     ! (Tri)diagonal matrix inversion
     mko(1) = cs(1)/bs(1) ! at the top
     DO i=2,nla 
        mko(i) = (cs(i)-mko(i-1)*as(i))/bs(i) 
     ENDDO
  ELSE
     ! explicit case
     as(1:nla)=-dzlay(1:nla)/dt
     bs(1:nla)=-ft(1:nla)
     ! boundaries
     mko(1)=cs(1)*mk(1)/as(1)
     mko(nla)=(bs(nla)*mk(nla-1)+cs(nla)*mk(nla))/as(nla)
     ! interior
     mko(2:nla-1)=(bs(2:nla-1)*mk(1:nla-2) + &
          cs(2:nla-1)*mk(2:nla-1)   &
          )/as(2:nla-1)
  ENDIF
  dmk = mko - mk

  RETURN
END SUBROUTINE let_me_fall_in_peace


!*****************************************************************************
SUBROUTINE get_weff(dt,nlay,plev,zlev,dzlay,btemp,mk,k,rc,sig,wth,corf)
!     call get_weff(plev,zlev,dzlay,btemp,m0ccn,0.D0,r0,sig,wth,fdcor)
  !! Get the effective settling velocity for aerosols moments.
  !!
  !! This method computes the effective settling velocity of the \(k^{th}\) order moment of aerosol
  !! tracers. The basic settling velocity (\(v^{eff}_{M_{k}}\)) is computed using the following
  !! equation:
  !!
  !! $$
  !! \begin{eqnarray*}
  !! \Phi^{sed}_{M_{k}} &=& \int_{0}^{\infty} n(r) r^{k} \times w(r) dr
  !!                     == M_{k}  \times v^{eff}_{M_{k}} \\
  !!     v^{eff}_{M_{k} &=& \dfrac{2 \rho g r_{c}^{\dfrac{3D_{f}-3}{D_{f}}}}
  !!     {r_{m}^{D_{f}-3}/D_{f}} \times \alpha(k)} \times \left( \alpha \right(
  !!     \frac{D_{f}(k+3)-3}{D_{f}}\left) + \dfrac{A_{kn}\lambda_{g}}{r_{c}^{
  !!     3/D_{f}}} \alpha \right( \frac{D_{f}(k+3)-6}{D_{f}}\left)\left)
  !! \end{eqnarray*}
  !! $$
  !!
  !! \(v^{eff}_{M_{k}\) is then corrected to reduce numerical diffusion of the sedimentation algorithm
  !! as defined in \cite{toon1988b}.
  !!
  !! @warning
  !! Both __df__, __rc__ and __afun__ must be consistent with each other otherwise wrong values will
  !! be computed.

  use free_param, only: rho_aer
  use donnees, only: akn, df, fiadero_max, fiadero_min, no_fiadero_w, rmono
  
  IMPLICIT NONE
  REAL, INTENT(in)  :: dt
  INTEGER, INTENT(in) :: nlay
  REAL, INTENT(in), DIMENSION(nlay+1)          :: plev
  !! Pressure at each level (Pa).
  REAL, INTENT(in), DIMENSION(nlay+1)          :: zlev
  !! Altitude at each level (m).
  REAL, INTENT(in), DIMENSION(nlay)            :: dzlay
  !! Altitude at the center of each layer (Pa).
  REAL, INTENT(in), DIMENSION(nlay)            :: btemp
  !! Temperature at each level (T).
  REAL, INTENT(in), DIMENSION(nlay)            :: mk
  !! Moment of order __k__ (\(m^{k}.m^{-3}\)).
  REAL, INTENT(in), DIMENSION(nlay)            :: rc
  !! Characteristic radius associated to the moment at each vertical __levels__.
  REAL, INTENT(in)                          :: k, sig
  !! Fractal dimension of the aersols.
  REAL, INTENT(out), DIMENSION(nlay+1)           :: wth
  !! Theoretical Settling velocity at each vertical __levels__ (\( wth \times corf = weff\)).
  REAL, INTENT(out), DIMENSION(nlay+1), OPTIONAL :: corf
  !! _Fiadero_ correction factor applied to the theoretical settling velocity at each vertical __levels__.
  INTEGER                             :: i,nla,nle
  REAL                       :: af1,af2,ar1,ar2
  REAL                       :: csto,cslf,ratio,wdt,dzb
  REAL                       :: rb2ra, VISAIR, FPLAIR, alpha_k
  REAL, DIMENSION(nlay+1) :: zcorf
  ! ------------------
  write(*,*)'yupitururu0'
  nla = SIZE(mk,DIM=1) 
  nle = nla+1

  write(*,*)'yupitururu1'
  wth(:) = 0.D0 
  zcorf(:) = 1.D0
  
  ar1 = (3.D0*df -3.D0)/df    
  ar2 = (3.D0*df -6.D0)/df
  af1 = (df*(k+3.D0)-3.D0)/df 
  af2 = (df*(k+3.D0)-6.D0)/df
  rb2ra = rmono**((df-3.D0)/df)
  write(*,*)'yupitururu2'
  DO i=2,nla
     IF (rc(i-1) <= 0.D0) CYCLE
     dzb = (dzlay(i)+dzlay(i-1))/2.D0
     csto = 2.D0*rho_aer*effg(zlev(i))/(9.D0*VISAIR(btemp(i)))
  write(*,*)'yupitururu3'
     ! ***** WARNING ******
     ! The 1st order slip flow correction is hidden here :
     !    (rc(i-1)**ar1 * afun(af1) + cslf/rb2ra * rc(i-1)**ar2 * afun(af2))
     ! wth is simple the classical settling velocity of particle due to gravity (Stokes law).
     ! It may be enhanced using a Taylor series of the slip flow correction around a given (modal) radius
     cslf = akn * FPLAIR(btemp(i),plev(i))
     wth(i) = -csto/(rb2ra*alpha_k(k,sig)) * (rc(i-1)**ar1 * alpha_k(af1,sig) + cslf/rb2ra * rc(i-1)**ar2 * alpha_k(af2,sig))
     ! now correct velocity to reduce numerical diffusion
     IF (.NOT.no_fiadero_w) THEN
        IF (mk(i) <= 0.D0) THEN
           ratio = fiadero_max
        ELSE
           ratio = MAX(MIN(mk(i-1)/mk(i),fiadero_max),fiadero_min)
        ENDIF
        ! apply correction
        IF ((ratio <= 0.9D0 .OR. ratio >= 1.1D0) .AND. wth(i) /= 0.D0) THEN
           wdt = wth(i)*dt
           zcorf(i) = dzb/wdt * (exp(-wdt*log(ratio)/dzb)-1.D0) / (1.D0-ratio)
        ENDIF
     ENDIF
  ENDDO
  ! last value (ground) set to first layer value: arbitrary :)
  wth(i) = wth(i-1)
  zcorf(i) = zcorf(i-1)
  IF (PRESENT(corf)) corf(:) = zcorf(:)
  RETURN
END SUBROUTINE get_weff


!*****************************************************************************
SUBROUTINE aer_sedimentation(dt,nlay,m0,m3,plev,zlev,dzlay,btemp,dm0,dm3,aer_flux)
  !! Interface to sedimentation algorithm.
  !!
  !! The subroutine computes the evolution of each moment of the aerosols tracers
  !! through sedimentation process and returns their tendencies for a timestep.
  !!
  !! It is assumed that the aerosols size-distribution is the same as the mode 1 of cloud drops.
  
  use free_param, only: rho_aer, sig_aer
  use donnees, only: pi, wsed_m0, wsed_m3
  
  IMPLICIT NONE

  REAL, INTENT(in)  :: dt
  INTEGER, INTENT(in) :: nlay
  REAL, INTENT(in), DIMENSION(nlay)   :: m0
  !! 0th order moment of the CCN precursors (\(m^{-3}\)).
  REAL, INTENT(in), DIMENSION(nlay)   :: m3
  !! 3rd order moment of the CCN precursors (\(m^{3}.m^{-3}\))
  REAL, INTENT(in), DIMENSION(nlay+1) :: plev
  !! Pressure at each level (Pa).
  REAL, INTENT(in), DIMENSION(nlay+1) :: zlev
  !! Altitude at each level (m).
  REAL, INTENT(in), DIMENSION(nlay)   :: dzlay
  !! Altitude at the center of each layer (Pa).
  REAL, INTENT(in), DIMENSION(nlay)   :: btemp
  !! Temperature at each level (T).
  REAL, INTENT(out), DIMENSION(nlay)  :: dm0
  !! Tendency of the 0th order moment of the CCN precursors (\(m^{-3}\)).
  REAL, INTENT(out), DIMENSION(nlay)  :: dm3
  !! Tendency of the 3rd order moment of the CCN precursors (\(m^{3}.m^{-3}\)).
  REAL, INTENT(out), DIMENSION(nlay+1):: aer_flux
  !! Aerosol mass (downard) flux
  REAL, DIMENSION(nlay+1):: ft,fdcor,wth
  REAL, DIMENSION(nlay)  :: m0vsed,m3vsed,r0
  REAL                   :: m,n,p
  REAL, PARAMETER        :: fac = 4.D0/3.D0 * pi
  INTEGER                :: nla,nle,i
  
  aer_flux(:) = 0.D0
  nla = size(m0,dim=1) 
  nle = nla+1

  DO i=1,nlay
     r0(i) = get_r0_1(m0(i),m3(i))
  ENDDO

  ! Compute sedimentation process based on control flag
  ! wsed_m0 = true, force all aerosols moments to fall at M0 settling velocity
  ! use M0ccn to compute settling velocity
  IF (wsed_m0) THEN  
     ! M0
     call get_weff(dt,nlay,plev,zlev,dzlay,btemp,m0,0.D0,r0,sig_aer,wth,fdcor)
     ft(:) = wth(:) * fdcor(:)  ; m0vsed(:) = ft(1:nla)
     call let_me_fall_in_peace(dt,nlay,dzlay,m0,-ft,dm0)
     ! M3
     m3vsed(:) = ft(1:nla)
     call let_me_fall_in_peace(dt,nlay,dzlay,m3,-ft,dm3)
     ! get sedimentation flux
     aer_flux(:) = fac * rho_aer * ft(:) * dm3
     write(*,*)'sedimentation- aerosol flux: ', aer_flux

  ! wsed_m3 = true, force all aerosols moments to fall at M3 settling velocity
  ! use M3ccn to compute settling velocity
  ELSEIF (wsed_m3) THEN 
     ! M3
     call get_weff(dt,nlay,plev,zlev,dzlay,btemp,m3,3.D0,r0,sig_aer,wth,fdcor)
     ft(:) = wth(:) * fdcor(:)  ; m3vsed(:) = ft(1:nla)
     call let_me_fall_in_peace(dt,nlay,dzlay,m3,-ft,dm3)
     ! M0
     m0vsed(:) = ft(1:nla)
     call let_me_fall_in_peace(dt,nlay,dzlay,m0,-ft,dm0)
     ! get sedimentation flux
     aer_flux(:) = fac * rho_aer * ft(:) * m3
     write(*,*)'sedimentation- aerosol flux: ', aer_flux

  ! m_wsed_m0 = false and m_wsed_m3 = false, computes the effective settling velocity for each moments
  ELSE 
     ! M0
     call get_weff(dt,nlay,plev,zlev,dzlay,btemp,m0,0.D0,r0,sig_aer,wth,fdcor)
     ft(:) = wth(:) * fdcor(:)  ; m0vsed(:) = ft(1:nla)
     call let_me_fall_in_peace(dt,nlay,dzlay,m0,-ft,dm0)
     ! M3
     call get_weff(dt,nlay,plev,zlev,dzlay,btemp,m3,3.D0,r0,sig_aer,wth,fdcor)
     ft(:) = wth(:) * fdcor(:)  ; m3vsed(:) = ft(1:nla)
     call let_me_fall_in_peace(dt,nlay,dzlay,m3,-ft,dm3)
     ! get sedimentation flux
     aer_flux(:) = fac * rho_aer * ft(:) * m3
     write(*,*)'sedimentation- aerosol flux: ', aer_flux
  ENDIF
  RETURN
END SUBROUTINE aer_sedimentation


!*****************************************************************************
SUBROUTINE drop_sedimentation(dt,nlay,M0mode,M3mode,plev,zlev,dzlay,btemp,mode,dm0ccn,dm0,dm3ccn,dm3liq)
  !! Interface to sedimentation algorithm.
  !!
  !! The subroutine computes the evolution of each moment of the cloud drops
  !! through sedimentation process and returns their tendencies for a timestep.
  !!
  !! m0mode,m3ccn should be vector defined on the vertical structure (top to ground).
  !! m3liq a 2D array with in 1st dimension the vertical distribution (top to ground) and
  !! in the 2nd the list of condensate.
  !!
  !! The global variables wsed_m0 and wsed_m3 control how the settling velocity is 
  !! applied during the sedimentation process:
  !!   - If both are .false., m0mode and the sum of m3 (ccn+liquid) are treated separatly.
  !!   - If mm_wsed_m0 is set to .true., all moments fall at m0mode settling velocity.
  !!   - Otherwise all moments fall at m3sum settling velocity.
  !!
  !! note:
  !!  drop sedimentation flux is not computed here. It can be done as in aer_sedimentation.
  !!  However the overall density of cloud drop should be used (that is the density of 
  !!  CCN + the density of each liquid condensate)
  !!     cld_flux = fac *ft * (m3ccn * rho_aer + SUM(rho_liq(:) + m3liq(:)))
  
  use free_param, only: sigma1, sigma2
  use donnees, only: pi, wsed_m0, wsed_m3
  
  IMPLICIT NONE

  REAL, INTENT(in)  :: dt
  INTEGER, INTENT(in) :: nlay
  REAL, INTENT(in), DIMENSION(nlay,2) :: M0mode
  !! 0th order moment (\(m^{-3}\)).
  REAL, INTENT(in), DIMENSION(nlay,3) :: M3mode
  !! 3rd order moment (\(m^{-3}\)).
  REAL, INTENT(in), DIMENSION(nlay+1) :: plev
  !! Pressure at each level (Pa).
  REAL, INTENT(in), DIMENSION(nlay+1) :: zlev
  !! Altitude at each level (m).
  REAL, INTENT(in), DIMENSION(nlay)   :: dzlay
  !! Altitude at the center of each layer (Pa).
  REAL, INTENT(in), DIMENSION(nlay)   :: btemp
  !! Temperature at each level (T).
  INTEGER, INTENT(in)                 :: mode
  !! Size-distribution mode.
  REAL, INTENT(out), DIMENSION(nlay)  :: dm0
  !! Tendency of the 0th order moment of droplets (\(m^{-3}\)).
  REAL, INTENT(out), DIMENSION(nlay)  :: dm0ccn
  !! Tendency of the 0th order moment of the CCN (\(m^{-3}\)).
  REAL, INTENT(out), DIMENSION(nlay)  :: dm3ccn
  !! Tendency of the 3rd order moment of the CCN (\(m^{3}.m^{-3}\)).
  REAL, INTENT(out), DIMENSION(nlay,2):: dm3liq
  !! Tendency of the 3rd order moments of each condensate (\(m^{3}.m^{-3}\)).
  REAL, DIMENSION(nlay+1):: ft, fdcor, wth
  REAL, DIMENSION(nlay)  :: r0, m3sum
  REAL                   :: m, n, p, sig
  REAL, PARAMETER        :: fac = 4.D0/3.D0 * pi
  INTEGER                :: i, nla, nle, nc

  nla = nlay
  nle = nla + 1 ; nc = nlay

  m3sum(:) = 0.D0
  DO i=1,nla
     m3sum(:) = M3mode(i,1) + M3mode(i,2) + M3mode(i,3) ! liq1 + liq2 + ccn
  ENDDO

  IF (mode == 2) THEN
     sig = sigma2
     DO i=1,nla
        r0(i) = get_r0_2(M0mode(i,1),m3sum(i))
     ENDDO
  ELSE
     sig = sigma1
     DO i=1,nla
        r0(i) = get_r0_1(M0mode(i,1),m3sum(i))
     ENDDO
  ENDIF

  ! Compute sedimentation process based on control flag
  ! wsed_m0 = true, force all aerosols moments to fall at M0 settling velocity
  ! use M0mode to compute settling velocity
  IF (wsed_m0) THEN
     ! use M0mode to compute settling velocity
     call get_weff(dt,nlay,plev,zlev,dzlay,btemp,M0mode(:,1),0.D0,r0,sig,wth,fdcor)
     ft(:) = wth(:) * fdcor(:)

     ! apply the sedimentation process to all other moments 
     call let_me_fall_in_peace(dt,nlay,dzlay,M0mode(:,1),-ft,dm0)
     call let_me_fall_in_peace(dt,nlay,dzlay,M0mode(:,2),-ft,dm0ccn)
     call let_me_fall_in_peace(dt,nlay,dzlay,M3mode(:,3),-ft,dm3ccn)
     call let_me_fall_in_peace(dt,nlay,dzlay,M3mode(:,2),-ft,dm3liq(:,2))
     call let_me_fall_in_peace(dt,nlay,dzlay,M3mode(:,1),-ft,dm3liq(:,1))

  ! Compute sedimentation process based on control flag
  ! wsed_m3 = true, force all aerosols moments to fall at M3 settling velocity
  ! use M3ccn to compute settling velocity
  ELSEIF (wsed_m3) THEN
     ! use M3ccn + M3liq = m3sum to compute settling velocity
     call get_weff(dt,nlay,plev,zlev,dzlay,btemp,m3sum,3.D0,r0,sig,wth,fdcor)
     ft(:) = wth(:) * fdcor(:) 

     ! apply the sedimentation process to all other moments 
     call let_me_fall_in_peace(dt,nlay,dzlay,M0mode(:,1),-ft,dm0)
     call let_me_fall_in_peace(dt,nlay,dzlay,M0mode(:,2),-ft,dm0ccn)
     call let_me_fall_in_peace(dt,nlay,dzlay,M3mode(:,3),-ft,dm3ccn)
     call let_me_fall_in_peace(dt,nlay,dzlay,M3mode(:,2),-ft,dm3liq(:,2))
     call let_me_fall_in_peace(dt,nlay,dzlay,M3mode(:,1),-ft,dm3liq(:,1))

  ELSE
     ! M0mode and M3sum fall independently
     ! M0
     call get_weff(dt,nlay,plev,zlev,dzlay,btemp,m0mode,0.D0,r0,sig,wth,fdcor)
     ft(:) = wth(:) * fdcor(:)
     call let_me_fall_in_peace(dt,nlay,dzlay,M0mode(:,1),-ft,dm0)
     call let_me_fall_in_peace(dt,nlay,dzlay,M0mode(:,2),-ft,dm0ccn)

     ! M3
     call get_weff(dt,nlay,plev,zlev,dzlay,btemp,m3sum,3.D0,r0,sig,wth,fdcor)
     ft(:) = wth(:) * fdcor(:) 
     call let_me_fall_in_peace(dt,nlay,dzlay,M3mode(:,3),-ft,dm3ccn)
     call let_me_fall_in_peace(dt,nlay,dzlay,M3mode(:,2),-ft,dm3liq(:,2))
     call let_me_fall_in_peace(dt,nlay,dzlay,M3mode(:,1),-ft,dm3liq(:,1))
  ENDIF

  RETURN
END SUBROUTINE drop_sedimentation

END MODULE SED_AND_PROD_MAD
