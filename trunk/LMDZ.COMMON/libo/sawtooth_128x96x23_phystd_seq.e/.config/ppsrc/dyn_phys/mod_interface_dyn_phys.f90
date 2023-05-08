










! 
! $Id: mod_interface_dyn_phys.F90 2351 2015-08-25 15:14:59Z emillour $
!
MODULE mod_interface_dyn_phys
  INTEGER,SAVE,dimension(:),allocatable :: index_i
  INTEGER,SAVE,dimension(:),allocatable :: index_j
  
  
CONTAINS
  
  SUBROUTINE Init_interface_dyn_phys
  ! dummy routine for seq case
  END SUBROUTINE Init_interface_dyn_phys 
! of #ifdef CPP_PARA
END MODULE mod_interface_dyn_phys
