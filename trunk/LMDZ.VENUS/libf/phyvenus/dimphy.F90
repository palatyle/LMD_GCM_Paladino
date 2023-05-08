MODULE dimphy
  
  INTEGER,SAVE :: klon
  INTEGER,SAVE :: kdlon
  INTEGER,SAVE :: kfdia
  INTEGER,SAVE :: kidia
  INTEGER,SAVE :: klev
  INTEGER,SAVE :: klevp1
  INTEGER,SAVE :: klevm1

!$OMP THREADPRIVATE(klon,kfdia,kidia,kdlon)

CONTAINS
  
  SUBROUTINE init_dimphy(klon0,klev0)
  IMPLICIT NONE
  
    INTEGER, INTENT(in) :: klon0
    INTEGER, INTENT(in) :: klev0
    
    klon=klon0
    
    kdlon=klon
    kidia=1
    kfdia=klon
!$OMP MASTER 
    klev=klev0
    klevp1=klev+1
    klevm1=klev-1
!$OMP END MASTER    
    
  END SUBROUTINE init_dimphy

  
END MODULE dimphy
