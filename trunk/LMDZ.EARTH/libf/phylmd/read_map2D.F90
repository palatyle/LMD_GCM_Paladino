SUBROUTINE read_map2D(filename, varname, timestep, inverse, varout)
! Open file and read one variable for one timestep.
! Return variable for the given timestep. 
  USE dimphy
  USE netcdf
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para


  IMPLICIT NONE

! Input arguments
  CHARACTER(len=*), INTENT(IN)  :: filename     ! name of file to read
  CHARACTER(len=*), INTENT(IN)  :: varname      ! name of variable in file
  INTEGER, INTENT(IN)           :: timestep     ! actual timestep
  LOGICAL, INTENT(IN)           :: inverse      ! TRUE if latitude needs to be inversed
! Output argument
  REAL, DIMENSION(klon), INTENT(OUT) :: varout  ! The variable read from file for the given timestep

! Local variables
  INTEGER :: j
  INTEGER :: nid, nvarid, ierr
  INTEGER, DIMENSION(3) :: start, count
  CHARACTER(len=20)     :: modname='read_map2D'

  REAL, DIMENSION(nbp_lon,nbp_lat) :: var_glo2D     ! 2D global 
  REAL, DIMENSION(nbp_lon,nbp_lat) :: var_glo2D_tmp ! 2D global
  REAL, DIMENSION(klon_glo)        :: var_glo1D     ! 1D global
  INCLUDE "iniprint.h"

! Read variable from file. Done by master process MPI and master thread OpenMP
  IF (is_mpi_root .AND. is_omp_root) THEN
     ierr = NF90_OPEN(trim(filename), NF90_NOWRITE, nid)
     IF (ierr /= NF90_NOERR) CALL write_err_mess('Problem in opening file')

     ierr = NF90_INQ_VARID(nid, trim(varname), nvarid)
     IF (ierr /= NF90_NOERR) CALL write_err_mess('The variable is absent in file')
     
     start=(/1,1,timestep/)
     count=(/nbp_lon,nbp_lat,1/)
     ierr = NF90_GET_VAR(nid, nvarid, var_glo2D,start,count)
     IF (ierr /= NF90_NOERR) CALL write_err_mess('Problem in reading varaiable')

     ierr = NF90_CLOSE(nid)
     IF (ierr /= NF90_NOERR) CALL write_err_mess('Problem in closing file')

     ! Inverse latitude order
     IF (inverse) THEN
        var_glo2D_tmp(:,:) = var_glo2D(:,:)
        DO j=1, nbp_lat
           var_glo2D(:,j) = var_glo2D_tmp(:,nbp_lat-j+1)
        END DO
     END IF

     ! Transform the global field from 2D to 1D
     CALL grid2Dto1D_glo(var_glo2D,var_glo1D)

     WRITE(lunout,*) 'in read_map2D, filename = ', trim(filename)
     WRITE(lunout,*) 'in read_map2D, varname  = ', trim(varname)
     WRITE(lunout,*) 'in read_map2D, timestep = ', timestep
  ENDIF

! Scatter gloabl 1D variable to all processes
  CALL scatter(var_glo1D, varout)

  CONTAINS
    SUBROUTINE write_err_mess(err_mess)

      CHARACTER(len=*), INTENT(IN) :: err_mess
      INCLUDE "iniprint.h"
      
      WRITE(lunout,*) 'Error in read_map2D, filename = ', trim(filename)
      WRITE(lunout,*) 'Error in read_map2D, varname  = ', trim(varname)
      WRITE(lunout,*) 'Error in read_map2D, timestep = ', timestep

      CALL abort_gcm(modname, err_mess, 1)

    END SUBROUTINE write_err_mess

END SUBROUTINE read_map2D
