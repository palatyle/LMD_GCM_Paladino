
SUBROUTINE drop_coagul(TAIR,PAIR,dt,M0_m1,M3_m1,M0_m2,M3_m2)

  implicit none

  real, intent(in) :: TAIR, PAIR
  real, intent(in) :: dt

  real, intent(inout) :: M0_m1, M3_m1, M0_m2, M3_m2

  print*, 'plop coagulation'

END SUBROUTINE drop_coagul
