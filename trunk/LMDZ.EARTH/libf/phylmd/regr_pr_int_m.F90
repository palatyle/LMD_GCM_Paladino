! $Id$
module regr_pr_int_m

  ! Author: Lionel GUEZ

  implicit none

contains

  subroutine regr_pr_int(ncid, name, julien, plev, pplay, top_value, v3)

    ! "regr_pr_int" stands for "regrid pressure interpolation".
    ! In this procedure:
    ! -- the root process reads a 2D latitude-pressure field from a
    !    NetCDF file, at a given day.
    ! -- the field is packed to the LMDZ horizontal "physics"
    !    grid and scattered to all threads of all processes;
    ! -- in all the threads of all the processes, the field is regridded in
    !    pressure to the LMDZ vertical grid.
    ! We assume that, in the input file, the field has 3 dimensions:
    ! latitude, pressure, julian day.
    ! We assume that latitudes are in ascending order in the input file.
    ! The target vertical LMDZ grid is the grid of mid-layers.
    ! Regridding is by linear interpolation.

    use dimphy, only: klon
    use netcdf95, only: nf95_inq_varid, handle_err
    use netcdf, only: nf90_get_var
    use assert_m, only: assert
    use regr1_lint_m, only: regr1_lint
    use mod_phys_lmdz_mpi_data, only: is_mpi_root

    use mod_phys_lmdz_transfert_para, only: scatter2d
    ! (pack to the LMDZ horizontal "physics" grid and scatter)

    integer, intent(in):: ncid ! NetCDF ID of the file
    character(len=*), intent(in):: name ! of the NetCDF variable
    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    real, intent(in):: plev(:)
    ! (pressure level of input data, in Pa, in strictly ascending order)

    real, intent(in):: pplay(:, :) ! (klon, llm)
    ! (pression pour le mileu de chaque couche, en Pa)

    real, intent(in):: top_value
    ! (extra value of field at 0 pressure)

    real, intent(out):: v3(:, :) ! (klon, llm)
    ! (regridded field on the partial "physics" grid)
    ! ("v3(i, k)" is at longitude "xlon(i)", latitude
    ! "xlat(i)", middle of layer "k".)

    ! Variables local to the procedure:

    include "dimensions.h"
    integer varid, ncerr ! for NetCDF

    real  v1(iim, jjm + 1, 0:size(plev))
    ! (input field at day "julien", on the global "dynamics" horizontal grid)
    ! (First dimension is for longitude.
    ! The value is the same for all longitudes.
    ! "v1(:, j, k >=1)" is at latitude "rlatu(j)" and pressure "plev(k)".)

    real v2(klon, 0:size(plev))
    ! (field scattered to the partial "physics" horizontal grid)
    ! "v2(i, k >= 1)" is at longitude "xlon(i)", latitude "xlat(i)"
    ! and pressure "plev(k)".)

    integer i

    !--------------------------------------------

    call assert(shape(v3) == (/klon, llm/), "regr_pr_int v3")
    call assert(shape(pplay) == (/klon, llm/), "regr_pr_int pplay")

    !$omp master
    if (is_mpi_root) then
       call nf95_inq_varid(ncid, name, varid)

       ! Get data at the right day from the input file:
       ncerr = nf90_get_var(ncid, varid, v1(1, :, 1:), start=(/1, 1, julien/))
       call handle_err("regr_pr_int nf90_get_var " // name, ncerr, ncid)
       ! Latitudes are in ascending order in the input file while
       ! "rlatu" is in descending order so we need to invert order:
       v1(1, :, 1:) = v1(1, jjm+1:1:-1, 1:)

       ! Complete "v1" with the value at 0 pressure:
       v1(1, :, 0) = top_value

       ! Duplicate on all longitudes:
       v1(2:, :, :) = spread(v1(1, :, :), dim=1, ncopies=iim-1)
    end if
    !$omp end master

    call scatter2d(v1, v2)

    ! Regrid in pressure at each horizontal position:
    do i = 1, klon
       v3(i, llm:1:-1) = regr1_lint(v2(i, :), (/0., plev/), pplay(i, llm:1:-1))
       ! (invert order of indices because "pplay" is in descending order)
    end do

  end subroutine regr_pr_int

end module regr_pr_int_m
