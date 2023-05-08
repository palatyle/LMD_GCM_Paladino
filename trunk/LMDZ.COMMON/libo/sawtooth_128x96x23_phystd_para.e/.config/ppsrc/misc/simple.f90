










! $Id$
module simple

  use handle_err_m, only: handle_err
  
  implicit none

  private handle_err

contains

  subroutine nf95_open(path, mode, ncid, chunksize, ncerr)

    use netcdf, only: nf90_open

    character(len=*), intent(in):: path
    integer, intent(in):: mode
    integer, intent(out):: ncid
    integer, intent(inout), optional:: chunksize
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_open(path, mode, ncid, chunksize)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_open " // path, ncerr_not_opt)
    end if

  end subroutine nf95_open

  !************************

  subroutine nf95_inq_dimid(ncid, name, dimid, ncerr)

    use netcdf, only: nf90_inq_dimid

    integer,             intent(in) :: ncid
    character (len = *), intent(in) :: name
    integer,             intent(out) :: dimid
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_inq_dimid(ncid, name, dimid)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_inq_dimid " // name, ncerr_not_opt, ncid)
    end if

  end subroutine nf95_inq_dimid

  !************************

  subroutine nf95_inquire_dimension(ncid, dimid, name, nclen, ncerr)

    use netcdf, only: nf90_inquire_dimension

    integer,                       intent( in) :: ncid, dimid
    character (len = *), optional, intent(out) :: name
    integer,             optional, intent(out) :: nclen
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_inquire_dimension(ncid, dimid, name, nclen)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_inquire_dimension", ncerr_not_opt, ncid)
    end if

  end subroutine nf95_inquire_dimension

  !************************

  subroutine nf95_inq_varid(ncid, name, varid, ncerr)

    use netcdf, only: nf90_inq_varid

    integer,             intent(in) :: ncid
    character(len=*), intent(in):: name
    integer,             intent(out) :: varid
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_inq_varid(ncid, name, varid)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_inq_varid, name = " // name, ncerr_not_opt, ncid)
    end if

  end subroutine nf95_inq_varid

  !************************

  subroutine nf95_inquire_variable(ncid, varid, name, xtype, ndims, dimids, &
       nAtts, ncerr)

    ! In "nf90_inquire_variable", "dimids" is an assumed-size array.
    ! This is not optimal.
    ! We are in the classical case of an array the size of which is
    ! unknown in the calling procedure, before the call.
    ! Here we use a better solution: a pointer argument array.
    ! This procedure associates and defines "dimids" if it is present.

    use netcdf, only: nf90_inquire_variable, nf90_max_var_dims

    integer, intent(in):: ncid, varid
    character(len = *), optional, intent(out):: name
    integer, optional, intent(out) :: xtype, ndims
    integer, dimension(:), optional, pointer :: dimids
    integer, optional, intent(out) :: nAtts
    integer, intent(out), optional :: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt
    integer dimids_local(nf90_max_var_dims)
    integer ndims_not_opt

    !-------------------

    if (present(dimids)) then
       ncerr_not_opt = nf90_inquire_variable(ncid, varid, name, xtype, &
            ndims_not_opt, dimids_local, nAtts)
       allocate(dimids(ndims_not_opt)) ! also works if ndims_not_opt == 0
       dimids = dimids_local(:ndims_not_opt)
       if (present(ndims)) ndims = ndims_not_opt
    else
       ncerr_not_opt = nf90_inquire_variable(ncid, varid, name, xtype, ndims, &
            nAtts=nAtts)
    end if

    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_inquire_variable", ncerr_not_opt, ncid, varid)
    end if

  end subroutine nf95_inquire_variable

  !************************

  subroutine nf95_create(path, cmode, ncid, initialsize, chunksize, ncerr)
    
    use netcdf, only: nf90_create

    character (len = *), intent(in   ) :: path
    integer,             intent(in   ) :: cmode
    integer,             intent(  out) :: ncid
    integer, optional,   intent(in   ) :: initialsize
    integer, optional,   intent(inout) :: chunksize
    integer, intent(out), optional :: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_create(path, cmode, ncid, initialsize, chunksize)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_create " // path, ncerr_not_opt)
    end if

  end subroutine nf95_create

  !************************

  subroutine nf95_def_dim(ncid, name, nclen, dimid, ncerr)

    use netcdf, only: nf90_def_dim

    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent( in) :: nclen
    integer,             intent(out) :: dimid
    integer, intent(out), optional :: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_def_dim(ncid, name, nclen, dimid)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_def_dim " // name, ncerr_not_opt, ncid)
    end if

  end subroutine nf95_def_dim

  !***********************

  subroutine nf95_redef(ncid, ncerr)

    use netcdf, only: nf90_redef

    integer, intent( in) :: ncid
    integer, intent(out), optional :: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_redef(ncid)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_redef", ncerr_not_opt, ncid)
    end if

  end subroutine nf95_redef
  
  !***********************

  subroutine nf95_enddef(ncid, h_minfree, v_align, v_minfree, r_align, ncerr)

    use netcdf, only: nf90_enddef

    integer,           intent( in) :: ncid
    integer, optional, intent( in) :: h_minfree, v_align, v_minfree, r_align
    integer, intent(out), optional :: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_enddef(ncid, h_minfree, v_align, v_minfree, r_align)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_enddef", ncerr_not_opt, ncid)
    end if

  end subroutine nf95_enddef

  !***********************

  subroutine nf95_close(ncid, ncerr)

    use netcdf, only: nf90_close

    integer, intent( in) :: ncid
    integer, intent(out), optional :: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_close(ncid)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_close", ncerr_not_opt)
    end if

  end subroutine nf95_close

  !***********************

  subroutine nf95_copy_att(ncid_in, varid_in, name, ncid_out, varid_out, ncerr)

    use netcdf, only: nf90_copy_att

    integer, intent( in):: ncid_in,  varid_in
    character(len=*), intent( in):: name
    integer, intent( in):: ncid_out, varid_out
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_copy_att(ncid_in, varid_in, name, ncid_out, varid_out)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_copy_att " // name, ncerr_not_opt, ncid_out)
    end if

  end subroutine nf95_copy_att

  !***********************

  subroutine nf95_inquire_attribute(ncid, varid, name, xtype, nclen, attnum, &
       ncerr)

    use netcdf, only: nf90_inquire_attribute

    integer,             intent( in)           :: ncid, varid
    character (len = *), intent( in)           :: name
    integer,             intent(out), optional :: xtype, nclen, attnum
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_inquire_attribute(ncid, varid, name, xtype, nclen, &
         attnum)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_inquire_attribute " // name, ncerr_not_opt, &
            ncid, varid)
    end if

  end subroutine nf95_inquire_attribute

  !***********************

  subroutine nf95_inquire(ncid, nDimensions, nVariables, nAttributes, &
       unlimitedDimId, formatNum, ncerr)

    use netcdf, only: nf90_inquire

    integer,           intent( in) :: ncid
    integer, optional, intent(out) :: nDimensions, nVariables, nAttributes
    integer, optional, intent(out) :: unlimitedDimId, formatNum
    integer, intent(out), optional:: ncerr

    ! Variable local to the procedure:
    integer ncerr_not_opt

    !-------------------

    ncerr_not_opt = nf90_inquire(ncid, nDimensions, nVariables, nAttributes, &
         unlimitedDimId, formatNum)
    if (present(ncerr)) then
       ncerr = ncerr_not_opt
    else
       call handle_err("nf95_inquire", ncerr_not_opt, ncid)
    end if

  end subroutine nf95_inquire

end module simple
