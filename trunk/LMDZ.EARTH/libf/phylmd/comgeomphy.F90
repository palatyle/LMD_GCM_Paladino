module comgeomphy
   real,save,allocatable :: airephy(:)
   real,save,allocatable :: cuphy(:)
   real,save,allocatable :: cvphy(:)
   real,save,allocatable :: rlatd(:)
   real,save,allocatable :: rlond(:)
!$OMP THREADPRIVATE(airephy,cuphy,cvphy,rlatd,rlond)
contains
  
  subroutine InitComgeomphy
  USE mod_phys_lmdz_para
  implicit none
    
 
    allocate(airephy(klon_omp))
    allocate(cuphy(klon_omp))
    allocate(cvphy(klon_omp))
    allocate(rlatd(klon_omp))
    allocate(rlond(klon_omp))

  end subroutine InitComgeomphy
  
end module comgeomphy
