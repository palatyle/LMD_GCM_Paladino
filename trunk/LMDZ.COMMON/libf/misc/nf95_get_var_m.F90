module nf95_get_var_m

  use netcdf, only: nf90_get_var
  use handle_err_m, only: handle_err

  implicit none

  interface nf95_get_var
     module procedure nf95_get_var_FourByteReal, nf95_get_var_FourByteInt, &
          nf95_get_var_1D_FourByteReal, nf95_get_var_1D_FourByteInt, &
          nf95_get_var_2D_FourByteReal, &
          nf95_get_var_3D_FourByteInt, &
          nf95_get_var_3D_FourByteReal, &
          nf95_get_var_4D_FourByteReal, &
          nf95_get_var_5D_FourByteReal
  end interface

  private
  public nf95_get_var

contains

  subroutine nf95_get_var_FourByteReal(ncid, varid, values, start, ncerr)

    integer, intent(in) :: ncid, varid
    real, intent(out) :: values
    integer, dimension(:), optional, intent(in) :: start
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_FourByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_FourByteReal

  !***********************

  subroutine nf95_get_var_FourByteInt(ncid, varid, values, start, ncerr)

    integer, intent(in) :: ncid, varid
    integer, intent(out) :: values
    integer, dimension(:), optional, intent(in) :: start
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_FourByteInt", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_FourByteInt

  !***********************

  subroutine nf95_get_var_1D_FourByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    integer,                         intent(in) :: ncid, varid
    real, intent(out) :: values(:)
    integer, dimension(:), optional, intent(in) :: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_1D_FourByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_1D_FourByteReal

  !***********************

  subroutine nf95_get_var_1D_FourByteInt(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    integer,                         intent(in) :: ncid, varid
    integer, intent(out) :: values(:)
    integer, dimension(:), optional, intent(in) :: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_1D_FourByteInt", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_1D_FourByteInt

  !***********************

  subroutine nf95_get_var_1D_EightByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    use typesizes, only: eightByteReal

    integer,                         intent(in) :: ncid, varid
    real (kind = EightByteReal),     intent(out) :: values(:)
    integer, dimension(:), optional, intent(in):: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_1D_eightByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_1D_EightByteReal

  !***********************

  subroutine nf95_get_var_2D_FourByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    integer,                         intent(in) :: ncid, varid
    real , intent(out) :: values(:, :)
    integer, dimension(:), optional, intent(in) :: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_2D_FourByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_2D_FourByteReal

  !***********************

  subroutine nf95_get_var_2D_EightByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    use typesizes, only: EightByteReal

    integer,                         intent(in) :: ncid, varid
    real (kind = EightByteReal), intent(out) :: values(:, :)
    integer, dimension(:), optional, intent(in) :: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_2D_EightByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_2D_EightByteReal

  !***********************

  subroutine nf95_get_var_3D_FourByteInt(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    integer, intent(in):: ncid, varid
    integer, intent(out):: values(:, :, :)
    integer, dimension(:), optional, intent(in):: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_3D_FourByteInt", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_3D_FourByteInt

  !***********************

  subroutine nf95_get_var_3D_FourByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    integer,                         intent(in) :: ncid, varid
    real , intent(out) :: values(:, :, :)
    integer, dimension(:), optional, intent(in) :: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_3D_FourByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_3D_FourByteReal

  !***********************

  subroutine nf95_get_var_3D_EightByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    use typesizes, only: eightByteReal

    integer,                         intent(in) :: ncid, varid
    real (kind = EightByteReal),     intent(out) :: values(:, :, :)
    integer, dimension(:), optional, intent(in) :: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_3D_eightByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_3D_EightByteReal

  !***********************

  subroutine nf95_get_var_4D_FourByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    integer,                         intent(in) :: ncid, varid
    real , intent(out) :: values(:, :, :, :)
    integer, dimension(:), optional, intent(in) :: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_4D_FourByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_4D_FourByteReal

  !***********************

  subroutine nf95_get_var_4D_EightByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    use typesizes, only: EightByteReal

    integer, intent(in):: ncid, varid
    real(kind = EightByteReal), intent(out):: values(:, :, :, :)
    integer, dimension(:), optional, intent(in):: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_4D_EightByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_4D_EightByteReal

  !***********************

  subroutine nf95_get_var_5D_FourByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    integer, intent(in):: ncid, varid
    real, intent(out):: values(:, :, :, :, :)
    integer, dimension(:), optional, intent(in):: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_5D_FourByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_5D_FourByteReal

  !***********************

  subroutine nf95_get_var_5D_EightByteReal(ncid, varid, values, start, &
       count_nc, stride, map, ncerr)

    use typesizes, only: EightByteReal

    integer, intent(in):: ncid, varid
    real(kind = EightByteReal), intent(out):: values(:, :, :, :, :)
    integer, dimension(:), optional, intent(in):: start, count_nc, stride, map
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_get_var(ncid, varid, values, start, count_nc, &
         stride, map)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_get_var_5D_EightByteReal", ncerr_not_opt, ncid, &
            varid)
    end if

  end subroutine nf95_get_var_5D_EightByteReal

end module nf95_get_var_m
