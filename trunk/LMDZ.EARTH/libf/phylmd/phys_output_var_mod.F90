!
! phys_local_var_mod.F90 1327 2010-03-17 15:33:56Z idelkadi $

      MODULE phys_output_var_mod

      use dimphy
! Variables outputs pour les ecritures des sorties
!======================================================================
!
!
!======================================================================
! Declaration des variables

      REAL, SAVE, ALLOCATABLE :: snow_o(:), zfra_o(:)
!$OMP THREADPRIVATE(snow_o, zfra_o)
      INTEGER, save, ALLOCATABLE ::  itau_con(:)       ! Nombre de pas ou rflag <= 1
!$OMP THREADPRIVATE(itau_con)

CONTAINS

!======================================================================
SUBROUTINE phys_output_var_init
use dimphy

IMPLICIT NONE

      allocate(snow_o(klon), zfra_o(klon))
      allocate(itau_con(klon))

END SUBROUTINE phys_output_var_init

!======================================================================
SUBROUTINE phys_output_var_end
use dimphy
IMPLICIT NONE

      deallocate(snow_o,zfra_o,itau_con)

END SUBROUTINE phys_output_var_end

END MODULE phys_output_var_mod
