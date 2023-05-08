MODULE dimphy
  
  INTEGER,SAVE :: klon   ! number of atmospheric columns (for this OpenMP subgrid)
  INTEGER,SAVE :: klev   ! number of atmospheric layers
  INTEGER,SAVE :: klevp1 ! number of atmospheric layers+1
  INTEGER,SAVE :: klevm1 ! number of atmospheric layers-1
!  INTEGER,SAVE :: kflev

!$OMP THREADPRIVATE(klon)

CONTAINS
  
  SUBROUTINE Init_dimphy(klon0,klev0)
  IMPLICIT NONE
  
    INTEGER, INTENT(in) :: klon0
    INTEGER, INTENT(in) :: klev0
    
    klon=klon0
    
!$OMP MASTER 
    klev=klev0
    klevp1=klev+1
    klevm1=klev-1
!    kflev=klev
!$OMP END MASTER    
!$OMP BARRIER
    
  END SUBROUTINE Init_dimphy

  
END MODULE dimphy
