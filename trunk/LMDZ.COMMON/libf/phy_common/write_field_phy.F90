!
! $Id: write_field_phy.F90 2342 2015-08-19 13:21:38Z emillour $
!
MODULE write_field_phy

  ! Dump a field on the global (nbp_lon by nbp_lat) physics grid
  
  CONTAINS 
 
    SUBROUTINE WriteField_phy(name,Field,ll)
    USE mod_phys_lmdz_para, ONLY: klon_omp, is_mpi_root, &
                                  Gather
    USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo, &
                                 Grid1Dto2D_glo
    USE Write_Field, ONLY: WriteField
    
    IMPLICIT NONE

    CHARACTER(len=*),INTENT(IN) :: name
    INTEGER,INTENT(IN) :: ll
    REAL,INTENT(IN) :: Field(klon_omp,ll)

    real, dimension(klon_glo,ll):: New_Field
    real, dimension(nbp_lon,nbp_lat,ll):: Field_2d

    CALL Gather(Field,New_Field)
!$OMP MASTER
    IF (is_mpi_root) THEN	
      CALL Grid1Dto2D_glo(New_Field,Field_2D)
      CALL WriteField(name,Field_2d)
    ENDIF
!$OMP END MASTER

  
   END SUBROUTINE WriteField_phy
 
 END MODULE write_field_phy
