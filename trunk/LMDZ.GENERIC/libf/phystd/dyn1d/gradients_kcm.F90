subroutine gradients_kcm(profil_flag,rho_v,rho_n,T,dTdp,dPvdp,dPndp)

  use params_h
  use gases_h
  implicit none

  ! inputs
  integer profil_flag ! 0 = dry, 1 = moist, 2 = isothermal
  double precision rho_v,rho_n,T

  ! outputs
  double precision dTdp,dPndp,dPvdp
  double precision a_v

  ! internal
  !double precision cp_n,cp_v
  double precision cp_v

  double precision press, rho_plus, rho_minus, dVdT, rho_c
  double precision dlnr,dlna,dpsat,dsv
  double precision s_minus,s_plus
  double precision s_v,s_c,L
  double precision psat_plus,psat_minus,Pn
  double precision nul

  ! functions
  double precision cp_neutral

  !cp_n = cp_neutral(T)

  select case(profil_flag)
  case(2) ! isothermal

     dTdp  = 0.
     a_v   = rho_v/rho_n ! constant here
     dPndp = 1/(1d0+m_n/m_v*a_v)
     dPvdp = 1 - dPndp

  case(1) ! moist

     Pn = rho_n*T*rmn

     if(ngasmx.eq.1)then
        print*,'Cannot have moist adiabat with one gas...'
        stop
     endif

     if(trim(gnom(ngasmx)).eq.'H2O')then

        call psat_H2O(T-2d-1,psat_minus)
        call psat_H2O(T+2d-1,psat_plus)
        call psat_H2O(T,press)

        rho_minus = m_v*psat_minus*1d6/(Rc*(T-2d-1))
        rho_plus  = m_v*psat_plus*1d6/(Rc*(T+2d-1))

        call therm(T-2d-1,rho_minus*1d-3,nul,nul,nul,nul,nul,nul,nul,&
                   nul,nul,press,s_minus,nul)
        call therm(T+2d-1,rho_plus*1d-3,nul,nul,nul,nul,nul,nul,nul,&
                   nul,nul,press,s_plus,nul)
        s_c = 2.06 * log(T/273.15)

        s_plus  = s_plus  * 1d3
        s_minus = s_minus * 1d3
        s_c     = s_c     * 1d3 ! convert to SI

        if(T.lt.280.0)then 
           dpsat = press*1d6 * ( 1730.63*log(10.) / (T-39.714)**2 )
        else
           call tdpsdt(T,dpsat)
           dpsat = dpsat * 1d6 / T
        endif

     elseif(trim(gnom(ngasmx)).eq.'NH3')then

        call psat_NH3(T-2d-1,psat_minus)
        call psat_NH3(T+2d-1,psat_plus)
        call psat_NH3(T,press)

        rho_minus = m_v*psat_minus*1d6/(Rc*(T-2d-1))
        rho_plus  = m_v*psat_plus*1d6/(Rc*(T+2d-1))

        call latheat_NH3(T-2d-1,nul,s_minus)
        call latheat_NH3(T+2d-1,nul,s_plus)
        call latheat_NH3(T,s_c,nul)

        dpsat = press*1d6 * (-2*1.5609d-4*T + 0.1236)

     endif

     dsv  = (s_plus-s_minus)/4d-1  ! dsv*T = ds / d ln[T]
     s_v  = (s_plus+s_minus)/2d0
     dlnr = T/rho_v * (rho_plus-rho_minus)/4d-1  ! d ln[rho_v] / d ln[T]

     if(rho_n/rho_v.lt.1e-5)then
        dlna = -T*dsv/(s_v-s_c)
     else
        a_v  = rho_v/rho_n
        dlna = (rmn*dlnr - cp_n + rmn - a_v*T*dsv)/(a_v*(s_v-s_c)+rmn)  
        ! d ln[alpha_v] / d ln[T]
        ! note cp_n + rmn = cv_n, which is what's required
     endif

     dTdp  = 1d0 / (dpsat + rho_n*rmn*(1d0 + dlnr - dlna)) ! c.f. Marcq S2.2.2
     dPvdp = dTdp * dpsat
     dPndp = 1d0 - dPvdp ! from p = p_v + p_n

  case(0) ! dry

     cp_v=0.0
     if(trim(gnom(ngasmx)).eq.'H2O')then
        cp_v = (32.24+1.923d-3*T+1.055d-5*T**2-3.511d-9*T**3)/m_v
     elseif(trim(gnom(ngasmx)).eq.'NH3')then
        cp_v = 2.058d3
     elseif(trim(gnom(ngasmx)).eq.'CH4')then
        cp_v = 2.226d3
     endif
     
     dTdp  = 1/(rho_n*cp_n+rho_v*cp_v)
     dPndp = 1/(1d0+m_n/m_v*rho_v/rho_n)
     dPvdp = 1/(1d0+m_v/m_n*rho_n/rho_v)
      
  end select


end subroutine gradients_kcm
