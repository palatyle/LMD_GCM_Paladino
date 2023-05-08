subroutine kcmprof_fn(nlayer,psurf_rcm,qsurf_rcm,Tsurf_rcm,Tstra_rcm,P_rcm,Pl_rcm,z_rcm,T_rcm,q_rcm,m_rcm)

use params_h
use watercommon_h, only : mH2O
use gases_h
use comcstfi_mod, only: mugaz, cpp, g
use callkeys_mod, only: co2cond
implicit none

!     ----------------------------------------------------------------
!     Purpose: create profiles of T, rho_v, rho_n, Pv and Pn following 
!              Kasting 1988
!     Authour: Adapted from a code by E. Marcq by R. Wordsworth (2011)
!     ----------------------------------------------------------------

  integer ilay, nlay
  parameter (nlay=10000) ! number of vertical layers

  ! rcm inputs
  integer nlayer
  real Tsurf_rcm,Tstra_rcm

  ! rcm outputs
  real psurf_rcm,qsurf_rcm
  real P_rcm(1:nlayer)
  real Pl_rcm(1:nlayer+1)
  real z_rcm(1:nlayer)
  real T_rcm(1:nlayer),q_rcm(1:nlayer)
  real m_rcm(1:nlayer+1)

  ! rcm for interpolation (should really use log coords?)
  !double precision p1,p2,pnew,ilay_rcm

  double precision lnp1,lnp2,lnpnew
  real Dp_rcm, dlogp_rcm
  integer ilay_rcm,ilev_rcm,ifinal_rcm

  double precision Dz, Dp
  double precision Ptop, dlogp, Psat_max
  parameter (Ptop=1.0)                             ! Pressure at TOA [Pa]

  double precision T(1:nlay)                       ! temperature [K]
  double precision Ztab(1:nlay)                    ! altitude [m]
  double precision Pv(1:nlay),Pn(1:nlay),P(1:nlay) ! pressure [Pa]
  double precision rho_v(1:nlay), rho_n(1:nlay)    ! density [kg m^-3]
  double precision a_v(1:nlay)                     ! = rho_v/rho_n [kg/kg]
  double precision q_v(1:nlay)                     ! = rho_v/rho_tot [kg/kg]
  double precision mtot(1:nlay)                    ! = (rho_v+rho_n)/(n_v+n_n) [g/mol]

  integer profil_flag(1:nlay) ! 0 = dry, 1 = moist, 2 = isothermal

  ! inputs
  double precision Tsurf   ! surface temperature [K]
  double precision Psurf_v ! surface par. pressure (variable species) [Pa]
  double precision Psurf_n ! surface par. pressure (incondensible species)[Pa]
  double precision Ttop    ! stratospheric temperature [K]

  double precision dTdp        ! [K/Pa]
  double precision dPvdp,dPndp ! [Pa/Pa]
  double precision psat_v      ! local Psat_H2O value
  double precision Tcrit       ! Critical temperature [K]
  double precision rho_vTEMP,rho_nTEMP

  double precision TCO2cond ! for CO2 condensation quasi-hack

  ! variables necessary for steam.f90
  double precision rhol,rhov,nul

  ! for output
  double precision vmr

  logical verbose
  parameter(verbose=.true.)

  logical add_Pvar_to_total
  parameter(add_Pvar_to_total=.true.)

  ! initialise flags
  profil_flag(:) = 0

  !-------------------------------
  ! assign input variables
  m_n     = dble(mugaz/1000.)
  cp_n    = cpp
  ! modify/generalise later??

  Psat_max = 1000000.0 ! maximum vapour pressure [Pa]
                       ! set huge until further notice 

  if(vgas.lt.1)then
     if(psat_max.gt.0.0)then
        print*,'Must have Psat_max=0 if no variable species'
        psat_max=0.0
        !stop
     endif
     print*, 'Assuming pure atmosphere'
     m_v   = 1.0
     tcrit = 1000.0
  elseif(trim(gnom(vgas)).eq.'H2O')then
     m_v   = dble(mH2O/1000.)
     tcrit = 6.47d2
  elseif(trim(gnom(vgas)).eq.'NH3')then
     m_v   = 17.031/1000.
     tcrit = 4.06d2
  elseif(trim(gnom(vgas)).eq.'CH4')then
     m_v   = 16.04/1000.
     tcrit = 1.91d2
     stop
  else
     print*,'Variable gas not recognised!'
     call abort
  endif

  rmn     = rc/m_n
  Ttop    = dble(Tstra_rcm)
  Tsurf   = dble(Tsurf_rcm)

  psat_v  = psat_max
  if(vgas.gt.0)then
     if(trim(gnom(vgas)).eq.'H2O')then
        call Psat_H2O(tsurf,psat_v)
     elseif(trim(gnom(vgas)).eq.'NH3')then
        call Psat_NH3(tsurf,psat_v)
     endif
  endif

  ! Moist adiabat unless greater than or equal to psat_max
  if(psat_v*1d6.lt.psat_max)then
    Psurf_v = Psat_v*1d6
    profil_flag(1) = 1
  else 
    Psurf_v = psat_max
    profil_flag(1) = 0
  endif

  if(add_Pvar_to_total)then
    Psurf_n = dble(psurf_rcm)
    psurf_rcm = real(Psurf_n+Psurf_v)
  else
    Psurf_n = dble(psurf_rcm) - Psurf_v
  endif

  ! include relative humidity option
  !if(satval.lt.1.0)then
  !   Psurf_v = Psurf_v*satval
  !   profil_flag(1) = 0     
  !endif

  if(verbose)then
     print*,'Psat_v  =',psat_v*1d6
     print*,'Tsurf   =',Tsurf,' K'
     print*,'Ttop    =',Ttop,' K'
     print*,'Psurf_v =',Psurf_v,' Pa'
     print*,'Psurf_n =',Psurf_n,' Pa'
     print*,'m_n     =',m_n,' kg/mol'
     print*,'m_v     =',m_v,' kg/mol'
     print*,'rc      =',rc
  endif

  ! define fine pressure grid
  dlogp_rcm = -(log(psurf_rcm)-log(ptop))/nlayer

  P_rcm(1)  = psurf_rcm*exp(dlogp_rcm)
  do ilay_rcm=1,nlayer-1
     P_rcm(ilay_rcm+1) = P_rcm(ilay_rcm)*exp(dlogp_rcm)
  enddo

  Pl_rcm(1) = psurf_rcm
  do ilev_rcm=2,nlayer
     ! log-linear interpolation
     Pl_rcm(ilev_rcm) = exp( log( P_rcm(ilev_rcm)*P_rcm(ilev_rcm-1) )/2 )
  enddo

  !-------------------------------
  ! Layer 1
  T(1)     = Tsurf
  Pv(1)    = Psurf_v
  Pn(1)    = Psurf_n
  rho_n(1) = m_n*Pn(1)/(Rc*Tsurf)
  rho_v(1) = m_v*Pv(1)/(Rc*Tsurf)
  a_v(1)   = rho_v(1)/rho_n(1)

  ! log pressure grid spacing (constant)
  dlogp = -(log(Pn(1)+Pv(1))-log(ptop))/(nlay-1)
  
  call gradients_kcm(profil_flag(1),rho_v(1),rho_n(1),Tsurf,dTdp,dPvdp,dPndp)
  if(verbose)then
     print*, 'dT/dp ground [K/Pa] =',dTdp
  endif

  ! initial delta p, delta z
  Dp = (Pn(1) + Pv(1))*(exp(dlogp) - 1d0)
  Dz = -Dp/(  g*(rho_n(1) + rho_v(1))  )

  !-------------------------------
  ! Layer 2
  T(2)     = tsurf + dTdp*Dp 
  Pv(2)    = Pv(1) + dPvdp*Dp
  Pn(2)    = Pn(1) + dPndp*Dp
  rho_n(2) = m_n*Pn(2)/(Rc*T(2))
  rho_v(2) = m_v*Pv(2)/(Rc*T(2))
  a_v(2)   = rho_v(2)/rho_n(2)

  !-------------------------------
  ! start vertical ascent
  Ztab(1) = 0.
  do ilay=2,nlay-1

     ! calculate altitude levels (for diagnostic only)
     Dz         = -Dp/(  g*(rho_n(ilay) + rho_v(ilay))  )
     Ztab(ilay) = Dz + Ztab(ilay-1)

     ! 1st assume next layer same as last one
     profil_flag(ilay) = profil_flag(ilay-1)
  
     ! update delta p
     Dp = (Pn(ilay)+Pv(ilay))*(exp(dlogp) - 1d0)

     ! intial gradients call to calculate temperature at next level
     call gradients_kcm(profil_flag(ilay),rho_v(ilay),rho_n(ilay),&
                        T(ilay),dTdp,dPvdp,dPndp)

     T(ilay+1) = T(ilay) + dTdp*Dp

     ! test for moist adiabat at next level
     psat_v=psat_max

     if(vgas.gt.0)then
     if(trim(gnom(vgas)).eq.'H2O')then
        call Psat_H2O(T(ilay+1),psat_v)
     elseif(trim(gnom(vgas)).eq.'NH3')then
        call Psat_NH3(T(ilay+1),psat_v)
     endif
     endif

     if (psat_v*1d6 .lt. Pv(ilay)+dPvdp*Dp) then
        profil_flag(ilay)=1
        call gradients_kcm(profil_flag(ilay),rho_v(ilay),rho_n(ilay),&
                         T(ilay),dTdp,dPvdp,dPndp)
     endif

     ! test for stratosphere at next level
     if (T(ilay+1) .le. Ttop) then
        profil_flag(ilay)=2
        T(ilay+1)=Ttop
     endif

     ! calculate pressures at next level
     Pn(ilay+1) = Pn(ilay) + dPndp*Dp
     Pv(ilay+1) = Pv(ilay) + dPvdp*Dp

     if(profil_flag(ilay) .eq. 1)then
       
        psat_v=psat_max 

        if(vgas.gt.0)then
        if(trim(gnom(vgas)).eq.'H2O')then
           call Psat_H2O(T(ilay+1),psat_v)
        elseif(trim(gnom(vgas)).eq.'NH3')then
           call Psat_NH3(T(ilay+1),psat_v)
        endif
        endif

        if(Pv(ilay+1) .lt. psat_v*1e6)then
           Pv(ilay+1)=psat_v*1d6
        endif

     endif

     ! calculate gas densities at next level (assume ideal)
     rho_n(ilay+1) = m_n*Pn(ilay+1)/(rc*T(ilay+1)) 
     select case(profil_flag(ilay)) 
     case(2) ! isothermal
        rho_v(ilay+1) = rho_v(ilay)/rho_n(ilay)*rho_n(ilay+1)
     case(1) ! moist

        ! dont think this is necessary
        !call psat_est(T(ilay+1),psat_v)
        ! modify for ammonia!!!

        rho_v(ilay+1) = m_v*psat_v*1d6/(rc*T(ilay+1))
     case(0) ! dry
        rho_v(ilay+1) = m_v*Pv(ilay+1)/(rc*T(ilay+1)) 
     end select

  enddo

  Ztab(nlay)=Ztab(nlay-1)+Dz

  !-------------------------------
  ! save to kcm1d variables 

  ! surface quantities
  psurf_rcm = Pn(1) + Pv(1)
  qsurf_rcm = rho_v(1)/(rho_v(1) + rho_n(1))

  ! create q_v, mtot for saving
  do ilay=1,nlay
     mtot(ilay) = 1d3*(rho_v(ilay) + rho_n(ilay)) / &
                  (rho_v(ilay)/m_v + rho_n(ilay)/m_n)
     q_v(ilay)  = rho_v(ilay)/(rho_v(ilay) + rho_n(ilay))
     ! CHECK THIS
  enddo


  ! convert to rcm lower-res grid
  z_rcm(:) = 0.0
  T_rcm(:) = 0.0
  q_rcm(:) = 0.0
  m_rcm(:) = 0.0

  m_rcm(1) = real( 1d3*(rho_v(1) + rho_n(1)) / &
                  (rho_v(1)/m_v + rho_n(1)/m_n) )

  ilay_rcm=1
  do ilay=2,nlay

     if(ilay_rcm.le.nlayer)then
     ! interpolate rcm variables

        if(Pn(ilay)+Pv(ilay) .lt. P_rcm(ilay_rcm))then

           if(ilay.eq.1)then
              print*,'Error in create_profils: Psurf here less than Psurf in RCM!'
              call abort
           endif

           lnp1   = log(Pn(ilay-1)+Pv(ilay-1))
           lnp2   = log(Pn(ilay)+Pv(ilay))
           lnpnew = dble(log(P_rcm(ilay_rcm)))

           z_rcm(ilay_rcm) = real(Ztab(ilay-1)*(lnp2-lnpnew)/(lnp2-lnp1) &
                             + Ztab(ilay)*(lnpnew-lnp1)/(lnp2-lnp1))
           T_rcm(ilay_rcm) = real(T(ilay-1)*(lnp2-lnpnew)/(lnp2-lnp1) &
                             + T(ilay)*(lnpnew-lnp1)/(lnp2-lnp1))
           q_rcm(ilay_rcm) = real(q_v(ilay-1)*(lnp2-lnpnew)/(lnp2-lnp1) &
                             + q_v(ilay)*(lnpnew-lnp1)/(lnp2-lnp1))

           m_rcm(ilay_rcm+1) = real(mtot(ilay-1)*(lnp2-lnpnew)/(lnp2-lnp1) &
                             + mtot(ilay)*(lnpnew-lnp1)/(lnp2-lnp1))

           ilay_rcm = ilay_rcm+1
        endif

     endif
  enddo

  ifinal_rcm=ilay_rcm-1
  if(ifinal_rcm.lt.nlayer)then
     if(verbose)then
        print*,'Interpolation in kcmprof stopped at layer',ilay_rcm,'!'
     endif

     do ilay_rcm=ifinal_rcm+1,nlayer

        z_rcm(ilay_rcm) = z_rcm(ilay_rcm-1)
        T_rcm(ilay_rcm) = T_rcm(ilay_rcm-1)
        q_rcm(ilay_rcm) = q_rcm(ilay_rcm-1)
        m_rcm(ilay_rcm+1) = m_rcm(ilay_rcm)

     enddo
  endif

  do ilay=2,nlayer
     if(T_rcm(ilay).lt.Ttop)then
        T_rcm(ilay)=Ttop
     endif
  enddo

!    CO2 condensation 'haircut' of temperature profile if necessary
  if(co2cond)then
     print*,'CO2 condensation haircut - assumes CO2-dominated atmosphere!'
     do ilay=2,nlayer
        if(P_rcm(ilay).lt.518000.)then
           TCO2cond = (-3167.8)/(log(.01*P_rcm(ilay))-23.23) ! Fanale's formula
        else
           TCO2cond = 684.2-92.3*log(P_rcm(ilay))+4.32*log(P_rcm(ilay))**2 
           ! liquid-vapour transition (based on CRC handbook 2003 data)
        endif

        print*,'p=',P_rcm(ilay),', T=',T_rcm(ilay),' Tcond=',TCO2cond
        if(T_rcm(ilay).lt.TCO2cond)then
           T_rcm(ilay)=TCO2cond
        endif
     enddo
  endif

  return
end subroutine kcmprof_fn
