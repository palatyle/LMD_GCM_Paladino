










double precision function cp_neutral(T)

  use gases_h

  implicit none

  ! inputs
  double precision T


  ! this function has been disabled in gradients_kcm.F90 because it doesnt
  ! work if you have gaseous mixtures. need to decide whether to generalise
  ! it or simply remove entirely...

  ! Cp_n : cf CO2 dans abe&matsui (1988)
  !cp_neutral = (22.26+5.981d-2*T-3.501d-5*T**2+7.469d-9*T**3)/m_n

  if(trim(gnom(1)).eq.'N2_')then
     cp_neutral = 1040.0
  elseif(trim(gnom(1)).eq.'H2_')then
     cp_neutral = 14310.0
  else
     print*,'Gas not recognised in cp_neutral!'
     call abort
  endif

  
end function cp_neutral









