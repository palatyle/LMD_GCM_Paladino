!
! $Header$
!
MODULE write_field_phy

  CONTAINS 
 
    SUBROUTINE WriteField_phy(name,Field,ll)
    USE dimphy
    USE mod_phys_lmdz_para
    USE mod_grid_phy_lmdz
    USE Write_Field
    
    IMPLICIT NONE
    include 'dimensions.h'
    include 'paramet.h'

    character(len=*)   :: name
    INTEGER :: ll
    real, dimension(klon_omp,ll) :: Field
    real,save,allocatable :: Field_tmp(:,:)
    real, dimension(klon_glo,ll):: New_Field
    real, dimension(iim,jjp1,ll):: Field_2d

    CALL Gather(Field,New_Field)
!$OMP MASTER
    IF (is_mpi_root) THEN	
      CALL Grid1Dto2D_glo(New_Field,Field_2D)
      CALL WriteField(name,Field_2d)
    ENDIF
!$OMP END MASTER

  
   END SUBROUTINE WriteField_phy
 
 END MODULE write_field_phy
