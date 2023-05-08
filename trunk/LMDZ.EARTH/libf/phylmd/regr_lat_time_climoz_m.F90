! $Id$
module regr_lat_time_climoz_m

  ! Author: Lionel GUEZ

  implicit none

  private
  public regr_lat_time_climoz

contains

  subroutine regr_lat_time_climoz(read_climoz)

    ! "regr_lat_time_climoz" stands for "regrid latitude time
    ! climatology ozone".

    ! This procedure reads a climatology of ozone from a NetCDF file,
    ! regrids it in latitude and time, and writes the regridded field
    ! to a new NetCDF file.

    ! The input field depends on time, pressure level and latitude.

    ! If the input field has missing values, they must be signaled by
    ! the "missing_value" attribute.

    ! We assume that the input field is a step function of latitude
    ! and that the input latitude coordinate gives the centers of steps.
    ! Regridding in latitude is made by averaging, with a cosine of
    ! latitude factor.
    ! The target LMDZ latitude grid is the "scalar" grid: "rlatu".
    ! The values of "rlatu" are taken to be the centers of intervals.

    ! We assume that in the input file:

    ! -- Latitude is in degrees.

    ! -- Latitude and pressure are strictly monotonic (as all NetCDF
    ! coordinate variables should be).

    ! -- The time coordinate is in ascending order (even though we do
    ! not use its values).
    ! The input file may contain either values for 12 months or values
    ! for 14 months.
    ! If there are 14 months then we assume that we have (in that order):
    ! December, January, February, ..., November, December, January

    ! -- Missing values are contiguous, at the bottom of
    ! the vertical domain and at the latitudinal boundaries.

    ! If values are all missing at a given latitude and date, then we
    ! replace those missing values by values at the closest latitude,
    ! equatorward, with valid values.
    ! Then, at each latitude and each date, the missing values are replaced
    ! by the lowest valid value above missing values.

    ! Regridding in time is by linear interpolation.
    ! Monthly values are processed to get daily values, on the basis
    ! of a 360-day calendar.
    ! If there are 14 months, we use the first December value to
    ! interpolate values between January 1st and mid-January.
    ! We use the last January value to interpolate values between
    ! mid-December and end of December.
    ! If there are only 12 months in the input file then we assume
    ! periodicity for interpolation at the beginning and at the end of the
    ! year.

    use regr1_step_av_m, only: regr1_step_av
    use regr3_lint_m, only: regr3_lint
    use netcdf95, only: handle_err, nf95_close, nf95_get_att, nf95_gw_var, &
         nf95_inq_dimid, nf95_inq_varid, nf95_inquire_dimension, nf95_open, &
         nf95_put_var
    use netcdf, only: nf90_get_att, nf90_get_var, nf90_noerr, nf90_nowrite
    use assert_m, only: assert

    integer, intent(in):: read_climoz ! read ozone climatology
    ! Allowed values are 1 and 2
    ! 1: read a single ozone climatology that will be used day and night
    ! 2: read two ozone climatologies, the average day and night
    ! climatology and the daylight climatology

    ! Variables local to the procedure:

    include "dimensions.h"
    ! (for "jjm")
    include "paramet.h"
    ! (for the other included files)
    include "comgeom2.h"
    ! (for "rlatv")
    include "comconst.h"
    ! (for "pi")

    integer n_plev ! number of pressure levels in the input data
    integer n_lat ! number of latitudes in the input data
    integer n_month ! number of months in the input data

    real, pointer:: latitude(:)
    ! (of input data, converted to rad, sorted in strictly ascending order)

    real, allocatable:: lat_in_edg(:)
    ! (edges of latitude intervals for input data, in rad, in strictly
    ! ascending order)

    real, pointer:: plev(:)
    ! pressure levels of input data, sorted in strictly ascending
    ! order, converted to hPa

    logical desc_lat ! latitude in descending order in the input file
    logical desc_plev ! pressure levels in descending order in the input file

    real, allocatable:: o3_in(:, :, :, :)
    ! (n_lat, n_plev, n_month, read_climoz)
    ! ozone climatologies from the input file
    ! "o3_in(j, k, :, :)" is at latitude "latitude(j)" and pressure
    ! level "plev(k)".
    ! Third dimension is month index, first value may be December or January.
    ! "o3_in(:, :, :, 1)" is for the day- night average, "o3_in(:, :, :, 2)"
    ! is for daylight.

    real missing_value

    real, allocatable:: o3_regr_lat(:, :, :, :)
    ! (jjm + 1, n_plev, 0:13, read_climoz)
    ! mean of "o3_in" over a latitude interval of LMDZ
    ! First dimension is latitude interval.
    ! The latitude interval for "o3_regr_lat(j,:, :, :)" contains "rlatu(j)".
    ! If "j" is between 2 and "jjm" then the interval is:
    ! [rlatv(j), rlatv(j-1)]
    ! If "j" is 1 or "jjm + 1" then the interval is:
    ! [rlatv(1), pi / 2]
    ! or:
    ! [- pi / 2, rlatv(jjm)]
    ! respectively.
    ! "o3_regr_lat(:, k, :, :)" is for pressure level "plev(k)".
    ! Third dimension is month number, 1 for January.
    ! "o3_regr_lat(:, :, :, 1)" is average day and night,
    ! "o3_regr_lat(:, :, :, 2)" is for daylight.

    real, allocatable:: o3_out(:, :, :, :)
    ! (jjm + 1, n_plev, 360, read_climoz)
    ! regridded ozone climatology
    ! "o3_out(j, k, l, :)" is at latitude "rlatu(j)", pressure
    ! level "plev(k)" and date "January 1st 0h" + "tmidday(l)", in a
    ! 360-day calendar.
    ! "o3_out(:, :, :, 1)" is average day and night,
    ! "o3_out(:, :, :, 2)" is for daylight.

    integer j, k, l,m

    ! For NetCDF:
    integer ncid_in, ncid_out ! IDs for input and output files
    integer varid_plev, varid_time, varid, ncerr, dimid
    character(len=80) press_unit ! pressure unit

    integer varid_in(read_climoz), varid_out(read_climoz)
    ! index 1 is for average ozone day and night, index 2 is for
    ! daylight ozone.

    real, parameter:: tmidmonth(0:13) = (/(-15. + 30. * l, l = 0, 13)/)
    ! (time to middle of month, in days since January 1st 0h, in a
    ! 360-day calendar)
    ! (We add values -15 and 375 so that, for example, day 3 of the year is
    ! interpolated between the December and the January value.)

    real, parameter:: tmidday(360) = (/(l + 0.5, l = 0, 359)/)
    ! (time to middle of day, in days since January 1st 0h, in a
    ! 360-day calendar)

    !---------------------------------

    print *, "Call sequence information: regr_lat_time_climoz"
    call assert(read_climoz == 1 .or. read_climoz == 2, "regr_lat_time_climoz")

    call nf95_open("climoz.nc", nf90_nowrite, ncid_in)

    ! Get coordinates from the input file:

    call nf95_inq_varid(ncid_in, "latitude", varid)
    call nf95_gw_var(ncid_in, varid, latitude)
    ! Convert from degrees to rad, because we will take the sine of latitude:
    latitude = latitude / 180. * pi
    n_lat = size(latitude)
    ! We need to supply the latitudes to "regr1_step_av" in
    ! ascending order, so invert order if necessary:
    desc_lat = latitude(1) > latitude(n_lat)
    if (desc_lat) latitude = latitude(n_lat:1:-1)

    ! Compute edges of latitude intervals:
    allocate(lat_in_edg(n_lat + 1))
    lat_in_edg(1) = - pi / 2
    forall (j = 2:n_lat) lat_in_edg(j) = (latitude(j - 1) + latitude(j)) / 2
    lat_in_edg(n_lat + 1) = pi / 2
    deallocate(latitude) ! pointer

    call nf95_inq_varid(ncid_in, "plev", varid)
    call nf95_gw_var(ncid_in, varid, plev)
    n_plev = size(plev)
    ! We only need the pressure coordinate to copy it to the output file.
    ! The program "gcm" will assume that pressure levels are in
    ! ascending order in the regridded climatology so invert order if
    ! necessary:
    desc_plev = plev(1) > plev(n_plev)
    if (desc_plev) plev = plev(n_plev:1:-1)
    call nf95_get_att(ncid_in, varid, "units", press_unit)
    if (press_unit == "Pa") then
       ! Convert to hPa:
       plev = plev / 100.
    elseif (press_unit /= "hPa") then
       print *, "regr_lat_time_climoz: the only recognized units are Pa " &
            // "and hPa."
       stop 1
    end if

    ! Create the output file and get the variable IDs:
    call prepare_out(ncid_in, n_plev, ncid_out, varid_out, varid_plev, &
         varid_time)

    ! Write remaining coordinate variables:
    call nf95_put_var(ncid_out, varid_plev, plev)
    call nf95_put_var(ncid_out, varid_time, tmidday)

    deallocate(plev) ! pointer

    ! Get the  number of months:
    call nf95_inq_dimid(ncid_in, "time", dimid)
    call nf95_inquire_dimension(ncid_in, dimid, len=n_month)

    allocate(o3_in(n_lat, n_plev, n_month, read_climoz))

    call nf95_inq_varid(ncid_in, "tro3", varid_in(1))
    ncerr = nf90_get_var(ncid_in, varid_in(1), o3_in(:, :, :, 1))
    call handle_err("regr_lat_time_climoz nf90_get_var tro3", ncerr, ncid_in)

    if (read_climoz == 2) then
       call nf95_inq_varid(ncid_in, "tro3_daylight", varid_in(2))
       ncerr = nf90_get_var(ncid_in, varid_in(2), o3_in(:, :, :, 2))
       call handle_err("regr_lat_time_climoz nf90_get_var tro3_daylight", &
            ncerr, ncid_in, varid_in(2))
    end if

    !!!! Aymeric; problem with compilation here.... pb with o3_in
    !AS if (desc_lat) o3_in = o3_in(n_lat:1:-1, :, :, :)
    !AS if (desc_plev) o3_in = o3_in(:, n_plev:1:-1, :, :)

    do m = 1, read_climoz
       ncerr = nf90_get_att(ncid_in, varid_in(m), "missing_value", &
            missing_value)
       if (ncerr == nf90_noerr) then
          do l = 1, n_month
             ! Take care of latitudes where values are all missing:

             ! Next to the south pole:
             j = 1
!AS             do while (o3_in(j, 1, l, m) == missing_value)
!AS                j = j + 1
!AS             end do
!AS             if (j > 1) o3_in(:j-1, :, l, m) = &
!AS                  spread(o3_in(j, :, l, m), dim=1, ncopies=j-1)
             
             ! Next to the north pole:
             j = n_lat
!AS             do while (o3_in(j, 1, l, m) == missing_value)
!AS                j = j - 1
!AS             end do
!AS             if (j < n_lat) o3_in(j+1:, :, l, m) = &
!AS                  spread(o3_in(j, :, l, m), dim=1, ncopies=n_lat-j)

             ! Take care of missing values at high pressure:
             do j = 1, n_lat
                ! Find missing values, starting from top of atmosphere
                ! and going down.
                ! We have already taken care of latitudes full of
                ! missing values so the highest level has a valid value.
                k = 2
!AS                do while  (o3_in(j, k, l, m) /= missing_value .and. k < n_plev)
!AS                   k = k + 1
!AS                end do
                ! Replace missing values with the valid value at the
                ! lowest level above missing values:
!AS                if (o3_in(j, k, l, m) == missing_value) &
!AS                     o3_in(j, k:n_plev, l, m) = o3_in(j, k-1, l, m)
             end do
          end do
       else
          print *, "regr_lat_time_climoz: field ", m, &
               ", no missing value attribute"
       end if
    end do

    call nf95_close(ncid_in)

    allocate(o3_regr_lat(jjm + 1, n_plev, 0:13, read_climoz))
    allocate(o3_out(jjm + 1, n_plev, 360, read_climoz))

    ! Regrid in latitude:
    ! We average with respect to sine of latitude, which is
    ! equivalent to weighting by cosine of latitude:
    if (n_month == 12) then
       print *, &
            "Found 12 months in ozone climatologies, assuming periodicity..."
!AS       o3_regr_lat(jjm+1:1:-1, :, 1:12, :) = regr1_step_av(o3_in, &
!AS            xs=sin(lat_in_edg), xt=sin((/- pi / 2, rlatv(jjm:1:-1), pi / 2/)))
       ! (invert order of indices in "o3_regr_lat" because "rlatu" is
       ! in descending order)

       ! Duplicate January and December values, in preparation of time
       ! interpolation:
       o3_regr_lat(:, :, 0, :) = o3_regr_lat(:, :, 12, :)
       o3_regr_lat(:, :, 13, :) = o3_regr_lat(:, :, 1, :)
    else
       print *, "Using 14 months in ozone climatologies..."
!AS       o3_regr_lat(jjm+1:1:-1, :, :, :) = regr1_step_av(o3_in, &
!AS            xs=sin(lat_in_edg), xt=sin((/- pi / 2, rlatv(jjm:1:-1), pi / 2/)))
       ! (invert order of indices in "o3_regr_lat" because "rlatu" is
       ! in descending order)
    end if

    ! Regrid in time by linear interpolation:
    o3_out = regr3_lint(o3_regr_lat, tmidmonth, tmidday)

    ! Write to file:
    do m = 1, read_climoz
       call nf95_put_var(ncid_out, varid_out(m), o3_out(jjm+1:1:-1, :, :, m))
       ! (The order of "rlatu" is inverted in the output file)
    end do

    call nf95_close(ncid_out)

  end subroutine regr_lat_time_climoz

  !********************************************

  subroutine prepare_out(ncid_in, n_plev, ncid_out, varid_out, varid_plev, &
       varid_time)

    ! This subroutine creates the NetCDF output file, defines
    ! dimensions and variables, and writes one of the coordinate variables.

    use netcdf95, only: nf95_create, nf95_def_dim, nf95_def_var, &
         nf95_put_att, nf95_enddef, nf95_copy_att, nf95_put_var
    use netcdf, only: nf90_clobber, nf90_float, nf90_global

    integer, intent(in):: ncid_in, n_plev
    integer, intent(out):: ncid_out, varid_plev, varid_time

    integer, intent(out):: varid_out(:) ! dim(1 or 2)
    ! "varid_out(1)" is for average ozone day and night,
    ! "varid_out(2)" is for daylight ozone.

    ! Variables local to the procedure:

    include "dimensions.h"
    ! (for "jjm")
    include "paramet.h"
    ! (for the other included files)
    include "comgeom2.h"
    ! (for "rlatu")
    include "comconst.h"
    ! (for "pi")

    integer ncerr
    integer dimid_rlatu, dimid_plev, dimid_time
    integer varid_rlatu

    !---------------------------

    print *, "Call sequence information: prepare_out"

    call nf95_create("climoz_LMDZ.nc", nf90_clobber, ncid_out)

    ! Dimensions:
    call nf95_def_dim(ncid_out, "time", 360, dimid_time)
    call nf95_def_dim(ncid_out, "plev", n_plev, dimid_plev)
    call nf95_def_dim(ncid_out, "rlatu", jjm + 1, dimid_rlatu)

    ! Define coordinate variables:

    call nf95_def_var(ncid_out, "time", nf90_float, dimid_time, varid_time)
    call nf95_put_att(ncid_out, varid_time, "units", "days since 2000-1-1")
    call nf95_put_att(ncid_out, varid_time, "calendar", "360_day")
    call nf95_put_att(ncid_out, varid_time, "standard_name", "time")

    call nf95_def_var(ncid_out, "plev", nf90_float, dimid_plev, varid_plev)
    call nf95_put_att(ncid_out, varid_plev, "units", "millibar")
    call nf95_put_att(ncid_out, varid_plev, "standard_name", "air_pressure")
    call nf95_put_att(ncid_out, varid_plev, "long_name", "air pressure")

    call nf95_def_var(ncid_out, "rlatu", nf90_float, dimid_rlatu, varid_rlatu)
    call nf95_put_att(ncid_out, varid_rlatu, "units", "degrees_north")
    call nf95_put_att(ncid_out, varid_rlatu, "standard_name", "latitude")

    ! Define the primary variables:

    call nf95_def_var(ncid_out, "tro3", nf90_float, &
         (/dimid_rlatu, dimid_plev, dimid_time/), varid_out(1))
    call nf95_put_att(ncid_out, varid_out(1), "long_name", &
         "ozone mole fraction")
    call nf95_put_att(ncid_out, varid_out(1), "standard_name", &
         "mole_fraction_of_ozone_in_air")

    if (size(varid_out) == 2) then
       call nf95_def_var(ncid_out, "tro3_daylight", nf90_float, &
            (/dimid_rlatu, dimid_plev, dimid_time/), varid_out(2))
       call nf95_put_att(ncid_out, varid_out(2), "long_name", &
            "ozone mole fraction in daylight")
    end if

    ! Global attributes:

    ! The following commands, copying attributes, may fail.
    ! That is OK.
    ! It should just mean that the attribute is not defined in the input file.

    call nf95_copy_att(ncid_in, nf90_global, "Conventions", ncid_out, &
         nf90_global, ncerr)
    call handle_err_copy_att("Conventions")

    call nf95_copy_att(ncid_in, nf90_global, "title", ncid_out, nf90_global, &
         ncerr)
    call handle_err_copy_att("title")

    call nf95_copy_att(ncid_in, nf90_global, "institution", ncid_out, &
         nf90_global, ncerr)
    call handle_err_copy_att("institution")

    call nf95_copy_att(ncid_in, nf90_global, "source", ncid_out, nf90_global, &
         ncerr)
    call handle_err_copy_att("source")

    call nf95_put_att(ncid_out, nf90_global, "comment", "Regridded for LMDZ")

    call nf95_enddef(ncid_out)

    ! Write one of the coordinate variables:
    call nf95_put_var(ncid_out, varid_rlatu, rlatu(jjm+1:1:-1) / pi * 180.)
    ! (convert from rad to degrees and sort in ascending order)

  contains

    subroutine handle_err_copy_att(att_name)

      use netcdf, only: nf90_noerr, nf90_strerror

      character(len=*), intent(in):: att_name

      !----------------------------------------

      if (ncerr /= nf90_noerr) then
         print *, "regr_lat_time_climoz_m prepare_out nf95_copy_att " &
              // att_name // " -- " // trim(nf90_strerror(ncerr))
      end if

    end subroutine handle_err_copy_att

  end subroutine prepare_out

end module regr_lat_time_climoz_m
