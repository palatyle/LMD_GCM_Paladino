subroutine calcenergy_kcm(nlayer,Tsurf,T,Play,Plev,Qsurf,Q,muvar,Eatmtot)


use params_h, only : Rc
use watercommon_h, only : mH2O
use comcstfi_mod, only: g, mugaz
implicit none

!     ----------------------------------------------------------------
!     Purpose: Calculate total energy of the (steam) atmosphere
!     Authour: R. Wordsworth (2011)
!     ----------------------------------------------------------------


  ! inputs
  integer,intent(in) :: nlayer
  real Tsurf,T(1:nlayer)
  real Play(1:nlayer),Plev(1:nlayer+1)
  real Qsurf,Q(1:nlayer)
  real muvar(1,nlayer+1)

  ! internal
  integer il
  real m_n,m_v,cp_n,cp_v,cp_a,Eatmlay
  real s_c,rho_v,L
  double precision p_v, s_v, nul
  real VMR(1:nlayer)

  ! output
  real Eatmtot

  ! functions
  double precision cp_neutral

  m_n = mugaz/1000.
  m_v = mH2O/1000.


  do il=1,nlayer
     VMR(il)=Q(il)*(muvar(1,il+1)/mH2O)
  end do


  Eatmtot = 0.0
  do il=1,nlayer


     cp_n = cp_neutral(dble(T(il)))
     cp_v = (32.24+1.923d-3*T(il) + 1.055d-5*T(il)**2-3.511d-9*T(il)**3)/m_v
     ! / by m_v makes it per kg not per mol

     call psat_H2O(dble(T(il)),p_v)
     p_v   = p_v*1d6
     rho_v = real(p_v)*m_v/(real(Rc)*T(il))

     call therm(dble(T(il)),dble(rho_v)*1d-3,nul,nul,nul,nul,nul,nul,nul,&
                nul,nul,nul,s_v,nul)

     s_c = 2.06 * log(T(il)/273.15)
     s_v = s_v * 1d3
     s_c = s_c * 1d3
     L   = (real(s_v)-s_c)*T(il)

     cp_a    = (1-VMR(il))*cp_n + VMR(il)*cp_v 
     !cp_a    = (1-Q(il))*cp_n + Q(il)*cp_v 

     Eatmlay = ( cp_a*T(il) + L*VMR(il) ) * (Plev(il)-Plev(il+1))/g

     Eatmtot = Eatmtot + Eatmlay 
  enddo

  print*,'WARNING: E-calc currently assumes H2O saturation!'

end subroutine calcenergy_kcm

