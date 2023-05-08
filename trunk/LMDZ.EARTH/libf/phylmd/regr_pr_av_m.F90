! $Id$
module regr_pr_av_m

  ! Author: Lionel GUEZ

  implicit none

contains

  subroutine regr_pr_av(ncid, name, julien, press_in_edg, paprs, v3)

    ! "regr_pr_av" stands for "regrid pressure averaging".
    ! In this procedure:
    ! -- the root process reads 2D latitude-pressure fields from a
    !    NetCDF file, at a given day.
    ! -- the fields are packed to the LMDZ horizontal "physics"
    !    grid and scattered to all threads of all processes;
    ! -- in all the threads of all the processes, the fields are regridded in
    !    pressure to the LMDZ vertical grid.
    ! We assume that, in the input file, the fields have 3 dimensions:
    ! latitude, pressure, julian day.
    ! We assume that the input fields are already on the "rlatu"
    ! latitudes, except that latitudes are in ascending order in the input
    ! file.
    ! We assume that all the inputs fields have the same coordinates.

    ! The target vertical LMDZ grid is the grid of layer boundaries.
    ! Regridding in pressure is done by averaging a step function of pressure.

    ! All the fields are regridded as a single multi-dimensional array
    ! so it saves CPU time to call this procedure once for several NetCDF
    ! variables rather than several times, each time for a single
    ! NetCDF variable.

    use dimphy, only: klon
    use netcdf95, only: nf95_inq_varid, handle_err
    use netcdf, only: nf90_get_var
    use assert_m, only: assert
    use assert_eq_m, only: assert_eq
    use regr1_step_av_m, only: regr1_step_av
    use mod_phys_lmdz_mpi_data, only: is_mpi_root

    use mod_phys_lmdz_transfert_para, only: scatter2d
    ! (pack to the LMDZ horizontal "physics" grid and scatter)

    integer, intent(in):: ncid ! NetCDF ID of the file
    character(len=*), intent(in):: name(:) ! of the NetCDF variables
    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    real, intent(in):: press_in_edg(:)
    ! edges of pressure intervals for input data, in Pa, in strictly
    ! ascending order

    real, intent(in):: paprs(:, :) ! (klon, llm + 1)
    ! (pression pour chaque inter-couche, en Pa)

    real, intent(out):: v3(:, :, :) ! (klon, llm, size(name))
    ! regridded fields on the partial "physics" grid
    ! "v3(i, k, l)" is at longitude "xlon(i)", latitude
    ! "xlat(i)", in pressure interval "[paprs(i, k+1), paprs(i, k)]",
    ! for NetCDF variable "name(l)".

    ! Variables local to the procedure:

    include "dimensions.h"
    integer varid, ncerr ! for NetCDF

    real  v1(iim, jjm + 1, size(press_in_edg) - 1, size(name))
    ! input fields at day "julien", on the global "dynamics" horizontal grid
    ! First dimension is for longitude.
    ! The values are the same for all longitudes.
    ! "v1(:, j, k, l)" is at latitude "rlatu(j)", for
    ! pressure interval "[press_in_edg(k), press_in_edg(k+1)]" and
    ! NetCDF variable "name(l)".

    real v2(klon, size(press_in_edg) - 1, size(name))
    ! fields scattered to the partial "physics" horizontal grid
    ! "v2(i, k, l)" is at longitude "xlon(i)", latitude "xlat(i)",
    ! for pressure interval "[press_in_edg(k), press_in_edg(k+1)]" and
    ! NetCDF variable "name(l)".

    integer i, n_var

    !--------------------------------------------

    call assert(size(v3, 1) == klon, size(v3, 2) == llm, "regr_pr_av v3 klon")
    n_var = assert_eq(size(name), size(v3, 3), "regr_pr_av v3 n_var")
    call assert(shape(paprs) == (/klon, llm+1/), "regr_pr_av paprs")

    !$omp master
    if (is_mpi_root) then
       do i = 1, n_var
          call nf95_inq_varid(ncid, trim(name(i)), varid)
          
          ! Get data at the right day from the input file:
          ncerr = nf90_get_var(ncid, varid, v1(1, :, :, i), &
               start=(/1, 1, julien/))
          call handle_err("regr_pr_av nf90_get_var " // trim(name(i)), ncerr, &
               ncid)
       end do
       
       ! Latitudes are in ascending order in the input file while
       ! "rlatu" is in descending order so we need to invert order:
       v1(1, :, :, :) = v1(1, jjm+1:1:-1, :, :)

       ! Duplicate on all longitudes:
       v1(2:, :, :, :) = spread(v1(1, :, :, :), dim=1, ncopies=iim-1)
    end if
    !$omp end master

    call scatter2d(v1, v2)

    ! Regrid in pressure at each horizontal position:
    do i = 1, klon
       v3(i, llm:1:-1, :) = regr1_step_av(v2(i, :, :), press_in_edg, &
            paprs(i, llm+1:1:-1))
       ! (invert order of indices because "paprs" is in descending order)
    end do

  end subroutine regr_pr_av

end module regr_pr_av_m
