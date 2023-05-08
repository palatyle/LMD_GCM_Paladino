! $Id$
module open_climoz_m

  implicit none

contains

  subroutine open_climoz(ncid, press_in_edg)

    ! This procedure should be called once per "gcm" run, by a single
    ! thread of each MPI process.
    ! The root MPI process opens "climoz_LMDZ.nc", reads the pressure
    ! levels and broadcasts them to the other processes.

    ! We assume that, in "climoz_LMDZ.nc", the pressure levels are in hPa
    ! and in strictly ascending order.

    use netcdf95, only: nf95_open, nf95_close, nf95_gw_var, nf95_inq_varid
    use netcdf, only: nf90_nowrite

    use mod_phys_lmdz_mpi_data, only: is_mpi_root
    use mod_phys_lmdz_mpi_transfert, only: bcast_mpi ! broadcast

    integer, intent(out):: ncid ! of "climoz_LMDZ.nc", OpenMP shared

    real, pointer:: press_in_edg(:)
    ! edges of pressure intervals for ozone climatology, in Pa, in strictly
    ! ascending order, OpenMP shared

    ! Variables local to the procedure:

    real, pointer:: plev(:)
    ! (pressure levels for ozone climatology, converted to Pa, in strictly
    ! ascending order)

    integer varid ! for NetCDF
    integer n_plev ! number of pressure levels in the input data
    integer k

    !---------------------------------------

    print *, "Call sequence information: open_climoz"

    if (is_mpi_root) then
       call nf95_open("climoz_LMDZ.nc", nf90_nowrite, ncid)

       call nf95_inq_varid(ncid, "plev", varid)
       call nf95_gw_var(ncid, varid, plev)
       ! Convert from hPa to Pa because "paprs" and "pplay" are in Pa:
       plev = plev * 100.
       n_plev = size(plev)
    end if

    call bcast_mpi(n_plev)
    if (.not. is_mpi_root) allocate(plev(n_plev))
    call bcast_mpi(plev)
    
    ! Compute edges of pressure intervals:
    allocate(press_in_edg(n_plev + 1))
    if (is_mpi_root) then
       press_in_edg(1) = 0.
       ! We choose edges halfway in logarithm:
       forall (k = 2:n_plev) press_in_edg(k) = sqrt(plev(k - 1) * plev(k))
       press_in_edg(n_plev + 1) = huge(0.)
       ! (infinity, but any value guaranteed to be greater than the
       ! surface pressure would do)
    end if
    call bcast_mpi(press_in_edg)
    deallocate(plev) ! pointer

  end subroutine open_climoz

end module open_climoz_m
