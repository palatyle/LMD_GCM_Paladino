! Copyright Université Reims Champagnne-Ardenne (2010-2015)
! contributor: Jérémie Burgalat
! 
! jeremie.burgalat@univ-reims.fr
! 
! This software is a computer program whose purpose is to compute multi-variate
! linear interpolation.
! 
! This software is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-B
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
! 
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
! 
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-B license and that you accept its terms.

!! file: datasets.F90
!! summary: Dataset module definition file
!! author: J. Burgalat
!! date: 2014,2017

MODULE DATASETS 
  !! dataset definitions module
  !!
  !! This module defines simple derived types that encapsulate data set variables. For a N-dimension 
  !! dataset, N vectors of \(n_{i}\) elements and a single N-dimensional array of 
  !! \(\prod_{i}^{N} n_{i}\) elements are defined.
  !!
  !! If the data set is to be used within [[lintdset(module)]] module then each coordinate values of the 
  !! dataset must be sorted either in ascending or descending order with no duplicates. The module does 
  !! not provide functionnality to sort and to check such requirements.
  !!
  !! The module also provides two interfaces to initialize data sets from input data file which can be 
  !! a NetCDF file (if NetCDF support is enabled at compilation time) or an ASCII file. In the latter 
  !! case, the file must contain a header which is formatted as follows:
  !!
  !! - The first line must contain only one value which is the number of coordinates (__N__)
  !! - The second line must contain all the dimension sizes (that is __N__ values)
  !! - Each other lines should contain __N+1__ columns with, for the __N__ first columns, the 
  !!   coordinates values and finally the point value in the last column.
  !! 
  !! @note
  !! Note that for ASCII file, data must be ordered so first dimensions vary first. Same requirement 
  !! is needed for NetCDF file but in most cases, it is implicitly done (if dimensions are ordered).
  USE NETCDF
  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: read_dset, write_dset, clear_dset, is_in, has_data, debug

  LOGICAL :: debug = .false.  !! A control flag to enable verbose mode

  !> Initialize a data set from either an ASCII or a NetCDF file
  !!
  !! Whatever the kind of file, this interface assumes that you know the dimensionnality of the 
  !! extracted variable. If file content is not compatible with the kind of data set given as 
  !! output argument, or if some I/O error occured, the dataset is cleared.
  !! @note 
  !! Netcdf reader interface is available only if the library has been compiled with NetCDF support.
  INTERFACE read_dset 
    MODULE PROCEDURE ncdf_rd_1d,ncdf_rd_2d,ncdf_rd_3d,ncdf_rd_4d,ncdf_rd_5d
#if HAVE_NC4_FTN
    MODULE PROCEDURE ncdf4_rd_1d,ncdf4_rd_2d,ncdf4_rd_3d,ncdf4_rd_4d,ncdf4_rd_5d
#endif
    MODULE PROCEDURE ascii_rd_1d,ascii_rd_2d,ascii_rd_3d,ascii_rd_4d,ascii_rd_5d
  END INTERFACE

  !> Write a dataset in a netcdf (classic) file.
  INTERFACE write_dset 
    MODULE PROCEDURE nc_wr_1d,nc_wr_2d,nc_wr_3d,nc_wr_4d,nc_wr_5d
  END INTERFACE

  !> Clear the given data set
  INTERFACE clear_dset 
    MODULE PROCEDURE clr_1d_set,clr_2d_set,clr_3d_set,clr_4d_set,clr_5d_set
  END INTERFACE

  !> Check if given point is within the data set
  INTERFACE is_in 
    MODULE PROCEDURE is_in_1d,is_in_2d,is_in_3d,is_in_4d,is_in_5d
  END INTERFACE

  !> Private interface to netcdf informations getters
  INTERFACE get_nc_info
    MODULE PROCEDURE get_nc3_info
#if HAVE_NC4_FTN
    MODULE PROCEDURE get_nc4_info 
#endif
  END INTERFACE

  INTERFACE has_data
    MODULE PROCEDURE has_d_1d,has_d_2d,has_d_3d,has_d_4d,has_d_5d
  END INTERFACE


  TYPE, PUBLIC :: DSET1D
    !! A 1D data set
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: x              !! X coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: data           !! Tabulated function's value at each coordinate
    CHARACTER(len=NF90_MAX_NAME)            :: xname = "X"    !! Name of the X coordinate
    CHARACTER(len=NF90_MAX_NAME)            :: dname = "data" !! Name of the data block
  END TYPE DSET1D

  TYPE, PUBLIC :: DSET2D
    !! A 2D data set
    REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: x              !! X coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: y              !! Y coordinate tabulated values
    REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: data           !! Tabulated function's value at each coordinate
    CHARACTER(len=NF90_MAX_NAME)              :: xname = "X"    !! Name of the X coordinate
    CHARACTER(len=NF90_MAX_NAME)              :: yname = "Y"    !! Name of the Y coordinate
    CHARACTER(len=NF90_MAX_NAME)              :: dname = "data" !! Name of the data block
  END TYPE DSET2D

  TYPE, PUBLIC :: DSET3D
    !! A 3D data set
    REAL(kind=8), DIMENSION(:), ALLOCATABLE     :: x              !! X coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE     :: y              !! Y coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE     :: z              !! Z coordinate tabulated values
    REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: data           !! Tabulated function's value at each coordinate
    CHARACTER(len=NF90_MAX_NAME)                :: xname = "X"    !! Name of the X coordinate
    CHARACTER(len=NF90_MAX_NAME)                :: yname = "Y"    !! Name of the Y coordinate
    CHARACTER(len=NF90_MAX_NAME)                :: zname = "Z"    !! Name of the Z coordinate
    CHARACTER(len=NF90_MAX_NAME)                :: dname = "data" !! Name of the data block
  END TYPE DSET3D

  TYPE, PUBLIC :: DSET4D
    !! A 4D data set
    REAL(kind=8), DIMENSION(:), ALLOCATABLE       :: x              !! X coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE       :: y              !! Y coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE       :: z              !! Z coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE       :: t              !! T coordinate tabulated values
    REAL(kind=8), DIMENSION(:,:,:,:), ALLOCATABLE :: data           !! Tabulated function's value at each coordinate
    CHARACTER(len=NF90_MAX_NAME)                  :: xname = "X"    !! Name of the X coordinate
    CHARACTER(len=NF90_MAX_NAME)                  :: yname = "Y"    !! Name of the Y coordinate
    CHARACTER(len=NF90_MAX_NAME)                  :: zname = "Z"    !! Name of the Z coordinate
    CHARACTER(len=NF90_MAX_NAME)                  :: tname = "T"    !! Name of the T coordinate
    CHARACTER(len=NF90_MAX_NAME)                  :: dname = "data" !! Name of the data block
  END TYPE DSET4D

  TYPE, PUBLIC :: DSET5D
    !! A 5D data set
    REAL(kind=8), DIMENSION(:), ALLOCATABLE         :: x              !! X coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE         :: y              !! Y coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE         :: z              !! Z coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE         :: t              !! T coordinate tabulated values
    REAL(kind=8), DIMENSION(:), ALLOCATABLE         :: w              !! W coordinate tabulated values
    REAL(kind=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: data           !! Tabulated function's value at each coordinate
    CHARACTER(len=NF90_MAX_NAME)                    :: xname = "X"    !! Name of the X coordinate
    CHARACTER(len=NF90_MAX_NAME)                    :: yname = "Y"    !! Name of the Y coordinate
    CHARACTER(len=NF90_MAX_NAME)                    :: zname = "Z"    !! Name of the Z coordinate
    CHARACTER(len=NF90_MAX_NAME)                    :: tname = "T"    !! Name of the T coordinate
    CHARACTER(len=NF90_MAX_NAME)                    :: wname = "W"    !! Name of the W coordinate
    CHARACTER(len=NF90_MAX_NAME)                    :: dname = "data" !! Name of the data block
  END TYPE DSET5D

  CONTAINS


  FUNCTION is_in_1d(set,x) RESULT(ret)
    !! Check if point is in the 1D data set
    TYPE(DSET1D), INTENT(in) :: set !! Dataset object to search in
    REAL(kind=8), INTENT(in) :: x   !! coordinate of the point to check
    LOGICAL :: ret                   !! .true. if the point is in the data set, .false. otherwise
    REAL(kind=8) :: l,u
    ret=.true.
    l  = set%x(1) ; u= set%x(size(set%x))
    IF ((x>=l.EQV.x<=u).OR.(x<=l.EQV.x>=u)) ret=.false.
    RETURN
  END FUNCTION is_in_1d

  FUNCTION is_in_2d(set,x,y) RESULT(ret)
    !! Check if point is in the 2D data set
    TYPE(DSET2D), INTENT(in) :: set !! Dataset object to search in
    REAL(kind=8), INTENT(in) :: x   !! X coordinate of the point to check
    REAL(kind=8), INTENT(in) :: y   !! Y coordinate of the point to check
    LOGICAL :: ret                   !! .true. if the point is in the data set, .false. otherwise
    REAL(kind=8) :: l,u
    ret=.false.
    l  = set%x(1) ; u= set%x(size(set%x))
    IF ((x>l.AND.x>u).OR.(x<l.AND.x<u)) RETURN
    IF ((x>l.AND.x>u).OR.(x<l.AND.x<u)) RETURN 
    l  = set%y(1) ; u= set%y(size(set%y))
    IF ((y>l.AND.y>u).OR.(y<l.AND.y<u)) RETURN
    ret=.true.
    RETURN
  END FUNCTION is_in_2d

  FUNCTION is_in_3d(set,x,y,z) RESULT(ret)
    !! Check if point is in the 3D data set
    TYPE(DSET3D), INTENT(in) :: set !! Dataset object to search in
    REAL(kind=8), INTENT(in) :: x   !! X coordinate of the point to check
    REAL(kind=8), INTENT(in) :: y   !! Y coordinate of the point to check
    REAL(kind=8), INTENT(in) :: z   !! Z coordinate of the point to check
    LOGICAL :: ret                   !! .true. if the point is in the data set, .false. otherwise
    REAL(kind=8) :: l,u
    ret=.false.
    l  = set%x(1) ; u= set%x(size(set%x))
    IF ((x>l.AND.x>u).OR.(x<l.AND.x<u)) RETURN 
    l  = set%y(1) ; u= set%y(size(set%y))
    IF ((y>l.AND.y>u).OR.(y<l.AND.y<u)) RETURN
    l  = set%z(1) ; u= set%z(size(set%z))
    IF ((z>l.AND.z>u).OR.(z<l.AND.z<u)) RETURN
    ret=.true.
    RETURN
  END FUNCTION is_in_3d

  FUNCTION is_in_4d(set,x,y,z,t) RESULT(ret)
    !! Check if point is in the 4D data set
    TYPE(DSET4D), INTENT(in) :: set !! Dataset object to search in
    REAL(kind=8), INTENT(in) :: x   !! X coordinate of the point to check
    REAL(kind=8), INTENT(in) :: y   !! Y coordinate of the point to check
    REAL(kind=8), INTENT(in) :: z   !! Z coordinate of the point to check
    REAL(kind=8), INTENT(in) :: t   !! T coordinate of the point to check
    LOGICAL :: ret                   !! .true. if the point is in the data set, .false. otherwise
    REAL(kind=8) :: l,u
    ret=.false.
    l  = set%x(1) ; u= set%x(size(set%x))
    IF ((x>l.AND.x>u).OR.(x<l.AND.x<u)) RETURN 
    l  = set%y(1) ; u= set%y(size(set%y))
    IF ((y>l.AND.y>u).OR.(y<l.AND.y<u)) RETURN
    l  = set%z(1) ; u= set%z(size(set%z))
    IF ((z>l.AND.z>u).OR.(z<l.AND.z<u)) RETURN
    l  = set%t(1) ; u= set%t(size(set%t))
    IF ((t>l.AND.t>u).OR.(t<l.AND.t<u)) RETURN
    ret=.true.
    RETURN
  END FUNCTION is_in_4d

  FUNCTION is_in_5d(set,x,y,z,t,w) RESULT(ret)
    !! Check if point is in the 4D data set
    TYPE(DSET5D), INTENT(in) :: set !! Dataset object to search in
    REAL(kind=8), INTENT(in) :: x   !! X coordinate of the point to check
    REAL(kind=8), INTENT(in) :: y   !! Y coordinate of the point to check
    REAL(kind=8), INTENT(in) :: z   !! Z coordinate of the point to check
    REAL(kind=8), INTENT(in) :: t   !! T coordinate of the point to check
    REAL(kind=8), INTENT(in) :: w   !! W coordinate of the point to check
    LOGICAL :: ret                   !! .true. if the point is in the data set, .false. otherwise
    REAL(kind=8) :: l,u
    ret=.false.
    l  = set%x(1) ; u= set%x(size(set%x))
    IF ((x>l.AND.x>u).OR.(x<l.AND.x<u)) RETURN 
    l  = set%y(1) ; u= set%y(size(set%y))
    IF ((y>l.AND.y>u).OR.(y<l.AND.y<u)) RETURN
    l  = set%z(1) ; u= set%z(size(set%z))
    IF ((z>l.AND.z>u).OR.(z<l.AND.z<u)) RETURN
    l  = set%t(1) ; u= set%t(size(set%t))
    IF ((t>l.AND.t>u).OR.(t<l.AND.t<u)) RETURN
    l  = set%w(1) ; u= set%w(size(set%w))
    IF ((w>l.AND.w>u).OR.(w<l.AND.w<u)) RETURN
    ret=.true.
    RETURN
  END FUNCTION is_in_5d

  FUNCTION ascii_header(path,ds,dp) RESULT(ret)
    !! Read ASCII file header
    !! 
    !! The method assumes the header is on two lines :
    !!
    !! - the first line must contain a single value which is the number of
    !!   dimensions.
    !! - the second must contain N values with the size of each dimensions (N
    !!   referring to the first line number).
    !!
    !! The file remains open in the logical unit 666 unless some error occured. 
    CHARACTER(len=*), INTENT(in)       :: path !! Path of the ASCII file to read
    INTEGER, INTENT(out), DIMENSION(:) :: ds   !! Size of each dimensions
    INTEGER, INTENT(out), DIMENSION(:) :: dp   !! Product of each lower dimension size (1 for 1st)
    LOGICAL :: ret                             !! .true. if no error occured, .false. otherwise
    INTEGER           :: nd,i,e
    CHARACTER(len=15) :: i2s
    INQUIRE(666,OPENED=ret)
    IF (ret) THEN
      WRITE(*,*) 'ERROR: LUN 666 already used...'
      ret = .false. ; RETURN
    ENDIF
    ret = .false.
    OPEN(666,file=TRIM(path))
    READ(666,*) nd
    IF (nd /= SIZE(ds)) THEN
      WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
      WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be "//TRIM(i2s)//"D"
      RETURN
    ENDIF
    READ(666,*,iostat=e) ds
    IF (e /= 0) THEN
      WRITE(*,'(a)') 'ERROR: Cannot get dimensions size'
      CLOSE(666) ; RETURN
    ENDIF
    dp(1) = 1
    DO i=2,nd
      dp(i)=PRODUCT(ds(:i-1))
    ENDDO
    ret = .true.
  END FUNCTION ascii_header

  FUNCTION nc_get_dim_by_name(fid,dn,values,ds,did) RESULT(ret)
    !! Get informations about a dimension.
    !!
    !! The method gets the values and size of a given dimension.
    !!
    !! Errors occur if:
    !!
    !! - the given name is not a dimension name.
    !! - the size of the dimension is less than 2 (the method assumes dimension has an afferent variable with values).
    !! - the method can not retrieve the values of the afferent variable.
    INTEGER, INTENT(in)                                  :: fid
      !! Id of the NetCDF file
    CHARACTER(len=NF90_MAX_NAME), INTENT(in)             :: dn
      !! Name of the dimension.
    REAL(kind=8), INTENT(out), DIMENSION(:), ALLOCATABLE :: values
      !! Values of the dimension
    INTEGER, INTENT(out)                                 :: ds 
      !! Size of __values__
    INTEGER, INTENT(out)                                 :: did
      !! Id of the NetCDF dimension 
    LOGICAL :: ret
      !! .true. if no error(s) occured, .false. otherwise
    INTEGER                      :: vid,err
    CHARACTER(len=15)            :: i2s
    ret = .false.
    ! --- Get dimension informations
    err = NF90_INQ_DIMID(fid,dn,did)
    IF (err /= NF90_NOERR) RETURN
    err = NF90_INQUIRE_DIMENSION(fid,did,len=ds)
    IF (err /= NF90_NOERR) RETURN
    IF (ds < 2) RETURN
    err = NF90_INQ_VARID(fid,TRIM(dn),vid)
    IF (err /= NF90_NOERR) RETURN
    ALLOCATE(values(ds))
    err = NF90_GET_VAR(fid,vid,values)
    IF (err /= NF90_NOERR) RETURN
    ret = .true.
  END FUNCTION nc_get_dim_by_name

  FUNCTION nc_get_dim_by_id(fid,did,values,ds,dn) RESULT(ret)
    !! Get informations about a dimension.
    !!
    !! The method gets the values and size of a given dimension.
    !!
    !! Errors occur if:
    !!
    !! - the given id is not a dimension identifier.
    !! - the size of the dimension is less than 2 (the method assumes dimension has an afferent variable with values).
    !! - the method can not retrieve the values of the afferent variable.
    INTEGER, INTENT(in)                                  :: fid
      !! Id of the NetCDF file
    INTEGER, INTENT(in)                                  :: did
      !! Id of the NetCDF dimension 
    REAL(kind=8), INTENT(out), DIMENSION(:), ALLOCATABLE :: values
      !! Values of the dimension
    INTEGER, INTENT(out)                                 :: ds 
      !! Size of __values__
    CHARACTER(len=NF90_MAX_NAME), INTENT(out)            :: dn
      !! Name of the dimension.
    LOGICAL :: ret
      !! .true. if no error(s) occured, .false. otherwise
    INTEGER                      :: vid,err
    CHARACTER(len=15)            :: i2s
    ret = .false.
    ! --- Get dimension informations
    IF (NF90_INQUIRE_DIMENSION(fid,did,dn,ds) /= NF90_NOERR) RETURN
    IF (ds < 2) RETURN
    IF (NF90_INQ_VARID(fid,TRIM(dn),vid) /= NF90_NOERR) RETURN
    ALLOCATE(values(ds))
    IF (NF90_GET_VAR(fid,vid,values) /= NF90_NOERR) RETURN
    ret = .true.
  END FUNCTION nc_get_dim_by_id

  FUNCTION get_nc3_info(path,variable,ncid,vid,dimids,verbose) RESULT(ret)
    !! Get variable informations from NetCDF file
    !! 
    !! The method attempts to read the given NetCDF file header and retrieves
    !! variable's informations from it:
    !!
    !! - parent file ID
    !! - variable ID
    !! - variable's dimensions ID
    !! - variable's dimensions sizes
    !!
    !! The method always opens and closes the file. If the file cannot be opened,
    !! __ncid__ is set to -1.
    CHARACTER(len=*), INTENT(in)                    :: path
      !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in)                    :: variable
      !! Name of the NetCDF variable to extract 
    INTEGER, INTENT(out)                            :: ncid
      !! NetCDF file ID 
    INTEGER, INTENT(out)                            :: vid
      !! NetCDF Variable ID
    INTEGER, INTENT(out), DIMENSION(:), ALLOCATABLE :: dimids
      !! Id of the variable's dimensions. __dimids__ is not allocated on error.
    LOGICAL, INTENT(in), OPTIONAL                   :: verbose
      !! True to print out message on error (default to False).
    LOGICAL :: ret
      !! .true. if no errors occured, .false. otherwise
    INTEGER :: fid,ty,nc,err
    LOGICAL :: zlog
    zlog = .false. ; IF (PRESENT(verbose)) zlog = verbose
    ret = .false.
    ncid = -1
    ! Opens file
    IF (NF90_OPEN(TRIM(path),NF90_NOWRITE,fid) /= NF90_NOERR) THEN
      IF (zlog) WRITE(*,'(a)') 'ERROR: Cannot open '//trim(path)
      RETURN
    ENDIF
    ncid = fid
    ! Searches for variable
    IF (NF90_INQ_VARID(ncid,TRIM(variable),vid) /= NF90_NOERR) THEN
      IF (zlog) WRITE(*,'(a)') 'Cannot find '//TRIM(variable)//' in '//trim(path)
      nc = NF90_CLOSE(fid)
      RETURN
    ENDIF
    ! Get variable infos
    ! 1st call to get type and number of dimensions)
    IF (NF90_INQUIRE_VARIABLE(ncid,vid,xtype=ty,ndims=nc) /= NF90_NOERR) THEN
      IF (zlog) WRITE(*,'(a)') 'Cannot access to '//TRIM(variable)//' informations'
      nc = NF90_CLOSE(fid)
      RETURN
    ELSE
      ! Checks type
      IF (ty == NF90_CHAR) THEN
        IF (zlog) WRITE(*,'(a)') 'Inconsistent variable type (should be numeric)'
        nc = NF90_CLOSE(fid)
        RETURN
      ENDIF
      ALLOCATE(dimids(nc))
    ENDIF
    ! Gets variable's dimensions informations
    ! first get dimensions id
    IF (NF90_INQUIRE_VARIABLE(ncid,vid,dimids=dimids) /= NF90_NOERR) THEN
      IF (zlog) WRITE(*,'(a)') 'Cannot access to '//TRIM(variable)//' informations'
      nc = NF90_CLOSE(fid)
      DEALLOCATE(dimids)
      RETURN
    ENDIF
    nc = NF90_CLOSE(fid)
    ret = .true. 
  END FUNCTION get_nc3_info

#if HAVE_NC4_FTN 
  FUNCTION get_nc4_info(path,variable,group,ncid,vid,dimids,verbose) RESULT(ret)
    !! Get variable informations from NetCDF4 file
    !!
    !! The method attempts to read the given NetCDF file header and retrieves
    !! variable's informations from it :
    !!
    !! - parent group/file ID
    !! - variable ID
    !! - variable's dimensions ID
    !! - variable's dimensions sizes
    !!
    !! The method always opens and closes the file. If the file cannot be opened,
    !! __ncid__ is set to -1.
    CHARACTER(len=*), INTENT(in)                    :: path
      !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in)                    :: variable
      !! Name of the NetCDF variable to extract 
    CHARACTER(len=*), INTENT(in)                    :: group
      !! Fullname of the variable's parent group (should be set to empty string for root group)
    INTEGER, INTENT(out)                            :: ncid
      !! NetCDF group ID 
    INTEGER, INTENT(out)                            :: vid
      !! NetCDF Variable ID
    INTEGER, INTENT(out), DIMENSION(:), ALLOCATABLE :: dimids
      !! Id of the variable's dimensions. On error, __dimids__ is not allocated. 
    LOGICAL, INTENT(in), OPTIONAL                   :: verbose
      !! True to print out message on error (default to False).
    LOGICAL :: ret
      !! .true. if no errors occured, .false. otherwise
    INTEGER :: fid,ty,nc,err
    LOGICAL :: zlog
    zlog = .false. ; IF (PRESENT(verbose)) zlog = verbose
    ret = .false.
    ncid = -1
    ! Opens file
    IF (NF90_OPEN(TRIM(path),NF90_NOWRITE,fid) /= NF90_NOERR) THEN
      IF (zlog) WRITE(*,'(a)') 'ERROR: Cannot open '//trim(path)
      RETURN
    ENDIF
    ! Searches for group based on its fullname
    IF (LEN_TRIM(group) == 0) THEN
      nc = NF90_NOERR ; ncid = fid
    ELSE
      nc = NF90_INQ_GRP_FULL_NCID(fid,TRIM(group),ncid)
    ENDIF
    SELECT CASE(nc)
      ! NF90_ENOGRP is missing from Netcdf-Fortran4.2 : its value=-125
      CASE(-125)
        IF (zlog) WRITE(*,'(a)') TRIM(group)//' does not exist in '//TRIM(path)
        nc = NF90_CLOSE(fid) ; RETURN
      CASE(NF90_ENOTNC4,NF90_ESTRICTNC3)
        IF (zlog) WRITE(*,'(a)') TRIM(path)//' is not a NetCDF-4 file with "group" feature'
        nc = NF90_CLOSE(fid) ; RETURN
      CASE(NF90_EHDFERR)
        IF (zlog) WRITE(*,'(a)') "Too bad, an HDF5 error has been reported..."
        nc = NF90_CLOSE(fid) ; RETURN
    END SELECT
    ! Searches for variable
    IF (NF90_INQ_VARID(ncid,TRIM(variable),vid) /= NF90_NOERR) THEN
      IF (zlog) WRITE(*,'(a)') 'Cannot find '//TRIM(variable)//' in '//trim(path)
      nc = NF90_CLOSE(fid)
      RETURN
    ENDIF
    ! Get variable infos
    ! 1st call to get type and number of dimensions)
    IF (NF90_INQUIRE_VARIABLE(ncid,vid,xtype=ty,ndims=nc) /= NF90_NOERR) THEN
      IF (zlog) WRITE(*,'(a)') 'Cannot access to '//TRIM(variable)//' informations'
      nc = NF90_CLOSE(fid)
      RETURN
    ELSE
      ! Checks type
      IF (ty == NF90_CHAR) THEN
        IF (zlog) WRITE(*,'(a)') 'Inconsistent variable type (should be numeric)'
        nc = NF90_CLOSE(fid)
        RETURN
      ENDIF
      ALLOCATE(dimids(nc))
    ENDIF
    ! Gets variable's dimensions informations
    ! first get dimensions id
    IF (NF90_INQUIRE_VARIABLE(ncid,vid,dimids=dimids) /= NF90_NOERR) THEN
      IF (zlog) WRITE(*,'(a)') 'Cannot access to '//TRIM(variable)//' informations'
      nc = NF90_CLOSE(fid)
      DEALLOCATE(dimids)
      RETURN
    ENDIF
    nc = NF90_CLOSE(fid)
    ret = .true. 
  END FUNCTION get_nc4_info
#endif

  !-------------------------
  ! NetCDF data file readers
  !-------------------------

  FUNCTION ncdf_rd_1d(path,variable,set) RESULT(ret)
    !! Read a 1D data set from a NetCDF file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET1D), INTENT(out)    :: set      !! Output dataset object
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=NF90_MAX_NAME)       :: dn
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,variable,fi,vi,di)) RETURN
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 1) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get coordinates values 
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf_rd_1d

  FUNCTION ncdf_rd_2d(path,variable,set) RESULT(ret)
    !! Read a 2D data set from a NetCDF file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET2D), INTENT(out)    :: set      !! Output dataset object
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,variable,fi,vi,di)) RETURN 
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 2) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get dimensions informations
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Y coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(2),set%y,ds(2),set%yname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1),ds(2)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf_rd_2d

  FUNCTION ncdf_rd_3d(path,variable,set) RESULT(ret)
    !! Read a 3D data set from a NetCDF file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET3D), INTENT(out)    :: set      !! Output dataset object
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,variable,fi,vi,di)) RETURN 
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 3) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get dimensions informations
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Y coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(2),set%y,ds(2),set%yname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Z coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(3),set%z,ds(3),set%zname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1),ds(2),ds(3)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf_rd_3d

  FUNCTION ncdf_rd_4d(path,variable,set) RESULT(ret)
    !! Read a 4D data set from a NetCDF file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET4D), INTENT(out)    :: set      !! Output dataset object
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,variable,fi,vi,di)) RETURN 
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 4) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get dimensions informations
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Y coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(2),set%y,ds(2),set%yname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Z coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(3),set%z,ds(3),set%zname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ T coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(4),set%t,ds(4),set%tname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1),ds(2),ds(3),ds(4)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf_rd_4d

  FUNCTION ncdf_rd_5d(path,variable,set) RESULT(ret)
    !! Read a 5D data set from a NetCDF file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET5D), INTENT(out)    :: set      !! Output dataset object
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,variable,fi,vi,di)) RETURN
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 5) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get dimensions informations
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Y coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(2),set%y,ds(2),set%yname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Z coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(3),set%z,ds(3),set%zname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ T coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(4),set%t,ds(4),set%tname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ W coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(5),set%w,ds(5),set%wname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1),ds(2),ds(3),ds(4),ds(5)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf_rd_5d

  ! NetCDF4 with groups !
#if HAVE_NC4_FTN 

  FUNCTION ncdf4_rd_1d(path,group,variable,set) RESULT(ret)
    !! Read a 1D data set from a NetCDF4 file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: group    !! Parent group of the variable (empty for root)
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET1D), INTENT(out)    :: set      !! Output dataset
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=NF90_MAX_NAME)       :: dn
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,group,variable,fi,vi,di)) RETURN 
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 1) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get coordinates values 
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf4_rd_1d

  FUNCTION ncdf4_rd_2d(path,group,variable,set) RESULT(ret)
    !! Read a 2D data set from a NetCDF4 file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: group    !! Parent group of the variable (empty for root)
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET2D), INTENT(out)    :: set      !! Output dataset
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,group,variable,fi,vi,di)) rETURN
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 2) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get dimensions informations
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Y coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(2),set%y,ds(2),set%yname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1),ds(2)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf4_rd_2d

  FUNCTION ncdf4_rd_3d(path,group,variable,set) RESULT(ret)
    !! Read a 3D data set from a NetCDF4 file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: group    !! Parent group of the variable (empty for root)
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET3D), INTENT(out)    :: set      !! Output dataset
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,group,variable,fi,vi,di)) RETURN 
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 3) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get dimensions informations
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Y coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(2),set%y,ds(2),set%yname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Z coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(3),set%z,ds(3),set%zname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1),ds(2),ds(3)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf4_rd_3d

  FUNCTION ncdf4_rd_4d(path,group,variable,set) RESULT(ret)
    !! Read a 4D data set from a NetCDF4 file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: group    !! Parent group of the variable (empty for root)
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET4D), INTENT(out)    :: set      !! Output dataset
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,group,variable,fi,vi,di)) RETURN 
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 4) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get dimensions informations
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Y coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(2),set%y,ds(2),set%yname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Z coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(3),set%z,ds(3),set%zname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ T coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(4),set%t,ds(4),set%tname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1),ds(2),ds(3),ds(4)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf4_rd_4d

  FUNCTION ncdf4_rd_5d(path,group,variable,set) RESULT(ret)
    !! Read a 5D data set from a NetCDF4 file
    CHARACTER(len=*), INTENT(in) :: path     !! Path of the NetCDF file
    CHARACTER(len=*), INTENT(in) :: group    !! Parent group of the variable (empty for root)
    CHARACTER(len=*), INTENT(in) :: variable !! Name of the NetCDF variable to extract 
    TYPE(DSET5D), INTENT(out)    :: set      !! Output dataset
    LOGICAL :: ret                           !! .true. if no errors occured, .false. otherwise
    INTEGER                            :: fi,vi,nd,iret
    INTEGER, DIMENSION(:), ALLOCATABLE :: di,ds
    CHARACTER(len=15)                  :: i2s
    ret = .false.
    ! --- Reset the data set
    CALL clear_dset(set)
    ! --- Check NetCDF file info
    IF (.NOT.get_nc_info(path,group,variable,fi,vi,di)) RETURN
    iret = NF90_OPEN(path,NF90_NOWRITE,fi)
    ! --- Check dimension size
    nd = SIZE(di)
    IF (nd /= 5) THEN
      IF (debug) THEN
        WRITE(i2s,*) nd ; i2s=ADJUSTL(i2s)
        WRITE(*,'(a)') "ERROR: Incompatible size, DSET should be"//TRIM(i2s)//"D"
      ENDIF
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Get dimensions informations
    ALLOCATE(ds(nd))
    ! ------ X coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(1),set%x,ds(1),set%xname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Y coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(2),set%y,ds(2),set%yname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ Z coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(3),set%z,ds(3),set%zname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ T coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(4),set%t,ds(4),set%tname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    ! ------ W coordinate
    IF (.NOT.nc_get_dim_by_id(fi,di(5),set%w,ds(5),set%wname)) THEN
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
      CALL clear_dset(set) ; RETURN
    ENDIF
    ! --- Read data
    ALLOCATE(set%data(ds(1),ds(2),ds(3),ds(4),ds(5)))
    IF (NF90_GET_VAR(fi,vi,set%data) /= NF90_NOERR) THEN
      IF (debug) WRITE(*,'(a)') "ERROR:"//TRIM(variable)//": Cannot get values"
      nd = NF90_CLOSE(fi) ; CALL clear_dset(set) ; RETURN
    ENDIF
    nd = NF90_CLOSE(fi)
    set%dname = variable
    ret = .true.
    RETURN
  END FUNCTION ncdf4_rd_5d
#endif

  !------------------------
  ! ASCII data file readers
  !------------------------

  FUNCTION ascii_rd_1d(path,set) RESULT(ret)
    !! Read a 1D data set from a ASCII file
    CHARACTER(len=*), INTENT(in) :: path !! Path of the ASCII file to read
    TYPE(dset1d), INTENT(out)    :: set  !! output 1D dataset
    LOGICAL :: ret                       !! .true. if no error occured, .false. otherwise
    INTEGER                    :: i,e
    REAL(kind=8), DIMENSION(2) :: vl
    INTEGER, DIMENSION(1)      :: cc,ds,dp
    ret = .false.
    CALL clear_dset(set)
    IF (.NOT.ascii_header(path,ds,dp)) RETURN
    ALLOCATE(set%x(ds(1)), &
             set%data(ds(1)))
    cc(:) = 1
    DO i=1,ds(1)
      READ(666,*) vl(:)
      set%x(cc(1)) = vl(1)
      set%data(cc(1)) = vl(2)
      cc(1) = cc(1) + 1
    ENDDO
    READ(666,*,iostat=e) vl(1) 
    IF (e == 0) THEN
      IF (debug) WRITE(*,*) 'ERROR: Extra value found...'
      ret = .false. ; call clear_dset(set) ; CLOSE(666) ; RETURN
    ENDIF
    ret = .true.
    CLOSE(666)
    RETURN
  END FUNCTION ascii_rd_1d

  FUNCTION ascii_rd_2d(path,set) RESULT(ret)
    !! Read a 2D data set from a ASCII file
    CHARACTER(len=*), INTENT(in) :: path !! Path of the ASCII file to read
    TYPE(dset2d), INTENT(out)    :: set  !! output 2D dataset
    LOGICAL :: ret                       !! .true. if no error occured, .false. otherwise
    INTEGER                    :: i,e
    REAL(kind=8), DIMENSION(3) :: vl 
    INTEGER, DIMENSION(2)      :: cc,ds,dp
    ret = .false.
    CALL clear_dset(set)
    IF (.NOT.ascii_header(path,ds,dp)) RETURN
    ALLOCATE(set%x(ds(1)),set%y(ds(2)), &
             set%data(ds(1),ds(2)))
    cc(:) = 1
    DO i=1,PRODUCT(ds) 
      READ(666,*,iostat=e) vl     ! Reads the line
      IF (e /= 0) THEN
        IF (debug) WRITE(*,'(a)') 'ERROR: Malformed file, line: ',i+2
        ret = .false. ; call clear_dset(set) ; CLOSE(666) ; RETURN
      ENDIF
      ! --- X coordinate
      set%x(cc(1)) = vl(1)
      ! --- Y coordinate
      set%y(cc(2)) = vl(2)
      ! --- Data
      set%data(cc(1),cc(2)) = vl(3)
      ! - Update counters
      IF (mod(i,dp(1))==0) cc(1) = cc(1)+1 ; IF (cc(1) > ds(1)) cc(1)=1
      IF (mod(i,dp(2))==0) cc(2) = cc(2)+1 ; IF (cc(2) > ds(2)) cc(2)=1
    ENDDO
    READ(666,*,iostat=e) vl(1)
    IF (e == 0) THEN
      IF (debug) WRITE(*,'(a)') 'ERROR: Extra value found'
      ret = .false. ; call clear_dset(set) ; CLOSE(666) ; RETURN
    ENDIF
    ret = .true.
    CLOSE(666)
    RETURN
  END FUNCTION ascii_rd_2d

  FUNCTION ascii_rd_3d(path,set) RESULT(ret)
    !! Read a 3D data set from a ASCII file
    CHARACTER(len=*), INTENT(in) :: path !! Path of the ASCII file to read
    TYPE(dset3d), INTENT(out)    :: set  !! output 3D dataset
    LOGICAL :: ret                       !! .true. if no error occured, .false. otherwise
    INTEGER                    :: i,e
    REAL(kind=8), DIMENSION(4) :: vl 
    INTEGER, DIMENSION(3)      :: cc,ds,dp
    ret = .false.
    CALL clear_dset(set)
    IF (.NOT.ascii_header(path,ds,dp)) RETURN
    ALLOCATE(set%x(ds(1)),set%y(ds(2)),set%z(ds(3)), &
             set%data(ds(1),ds(2),ds(3)))
    cc(:) = 1
    DO i=1,PRODUCT(ds) 
      READ(666,*,iostat=e) vl     ! Reads the line
      IF (e /= 0) THEN
        IF (debug) WRITE(*,'(a)') 'ERROR: Malformed file, line: ',i+2
        ret = .false. ; call clear_dset(set) ; CLOSE(666) ; RETURN
      ENDIF
      ! --- X coordinate
      set%x(cc(1)) = vl(1)
      ! --- Y coordinate
      set%y(cc(2)) = vl(2)
      ! --- Z coordinate
      set%z(cc(3)) = vl(3)
      ! --- Data
      set%data(cc(1),cc(2),cc(3)) = vl(4)
      ! - Update counters
      IF (mod(i,dp(1))==0) cc(1) = cc(1)+1 ; IF (cc(1) > ds(1)) cc(1)=1
      IF (mod(i,dp(2))==0) cc(2) = cc(2)+1 ; IF (cc(2) > ds(2)) cc(2)=1
      IF (mod(i,dp(3))==0) cc(3) = cc(3)+1 ; IF (cc(3) > ds(3)) cc(3)=1
    ENDDO
    READ(666,*,iostat=e) vl(1)
    IF (e == 0) THEN
      IF (debug) WRITE(*,'(a)') 'ERROR: Extra value found'
      ret = .false. ; call clear_dset(set) ; CLOSE(666) ; RETURN
    ENDIF
    ret = .true.
    CLOSE(666)
    RETURN
  END FUNCTION ascii_rd_3d

  FUNCTION ascii_rd_4d(path,set) RESULT(ret)
    !! Read a 4D data set from a ASCII file
    CHARACTER(len=*), INTENT(in) :: path !! Path of the ASCII file to read
    TYPE(dset4d), INTENT(out)    :: set  !! output 4D dataset
    LOGICAL :: ret                       !! .true. if no error occured, .false. otherwise
    INTEGER                    :: i,e
    REAL(kind=8), DIMENSION(5) :: vl 
    INTEGER, DIMENSION(4)      :: cc,ds,dp
    ret = .false.
    CALL clear_dset(set)
    IF (.NOT.ascii_header(path,ds,dp)) RETURN
    ALLOCATE(set%x(ds(1)),set%y(ds(2)),set%z(ds(3)),set%t(ds(4)), &
             set%data(ds(1),ds(2),ds(3),ds(4)))
    cc(:) = 1
    DO i=1,PRODUCT(ds) 
      READ(666,*,iostat=e) vl     ! Reads the line
      IF (e /= 0) THEN
        IF (debug) WRITE(*,'(a)') 'ERROR: Malformed file, line: ',i+2
        ret = .false. ; call clear_dset(set) ; CLOSE(666) ; RETURN
      ENDIF
      ! --- X coordinate
      set%x(cc(1)) = vl(1)
      ! --- Y coordinate
      set%y(cc(2)) = vl(2)
      ! --- Z coordinate
      set%z(cc(3)) = vl(3)
      ! --- T coordinate
      set%t(cc(4)) = vl(4)
      ! --- Data
      set%data(cc(1),cc(2),cc(3),cc(4)) = vl(5)
      ! - Update counters
      IF (mod(i,dp(1))==0) cc(1) = cc(1)+1 ; IF (cc(1) > ds(1)) cc(1)=1
      IF (mod(i,dp(2))==0) cc(2) = cc(2)+1 ; IF (cc(2) > ds(2)) cc(2)=1
      IF (mod(i,dp(3))==0) cc(3) = cc(3)+1 ; IF (cc(3) > ds(3)) cc(3)=1
      IF (mod(i,dp(4))==0) cc(4) = cc(4)+1 ; IF (cc(4) > ds(4)) cc(4)=1
    ENDDO
    READ(666,*,iostat=e) vl(1)
    IF (e == 0) THEN
      IF (debug) WRITE(*,'(a)') 'ERROR: Extra value found'
      ret = .false. ; call clear_dset(set) ; CLOSE(666) ; RETURN
    ENDIF
    CLOSE(666)
    ret = .true.
    RETURN
  END FUNCTION ascii_rd_4d

  FUNCTION ascii_rd_5d(path,set) RESULT(ret)
    !! Read a 5D data set from a ASCII file
    CHARACTER(len=*), INTENT(in) :: path !! Path of the ASCII file to read
    TYPE(dset5d), INTENT(out)    :: set  !! output 5D dataset
    LOGICAL :: ret                       !! .true. if no error occured, .false. otherwise
    INTEGER                    :: i,e
    REAL(kind=8), DIMENSION(6) :: vl 
    INTEGER, DIMENSION(5)      :: cc,ds,dp
    ret = .false.
    CALL clear_dset(set)
    IF (.NOT.ascii_header(path,ds,dp)) RETURN
    ALLOCATE(set%x(ds(1)),set%y(ds(2)),set%z(ds(3)),set%t(ds(4)),set%w(ds(5)), &
             set%data(ds(1),ds(2),ds(3),ds(4),ds(5)))
    cc(:) = 1
    DO i=1,PRODUCT(ds) 
      READ(666,*,iostat=e) vl     ! Reads the line
      IF (e /= 0) THEN
        IF (debug) WRITE(*,'(a)') 'ERROR: Malformed file, line: ',i+2
        ret = .false. ; call clear_dset(set) ; CLOSE(666) ; RETURN
      ENDIF
      ! --- X coordinate
      set%x(cc(1)) = vl(1)
      ! --- Y coordinate
      set%y(cc(2)) = vl(2)
      ! --- Z coordinate
      set%z(cc(3)) = vl(3)
      ! --- T coordinate
      set%t(cc(4)) = vl(4)
      ! --- W coordinate
      set%w(cc(5)) = vl(5)
      ! --- Data
      set%data(cc(1),cc(2),cc(3),cc(4),cc(5)) = vl(6)
      ! - Update counters
      IF (mod(i,dp(1)) == 0) cc(1)=cc(1)+1 ; IF (cc(1) > ds(1)) cc(1)=1
      IF (mod(i,dp(2)) == 0) cc(2)=cc(2)+1 ; IF (cc(2) > ds(2)) cc(2)=1
      IF (mod(i,dp(3)) == 0) cc(3)=cc(3)+1 ; IF (cc(3) > ds(3)) cc(3)=1
      IF (mod(i,dp(4)) == 0) cc(4)=cc(4)+1 ; IF (cc(4) > ds(4)) cc(4)=1
      IF (mod(i,dp(5)) == 0) cc(5)=cc(5)+1 ; IF (cc(5) > ds(5)) cc(5)=1
    ENDDO
    READ(666,*,iostat=e) vl(1)
    IF (e == 0) THEN
      IF (debug) WRITE(*,'(a)') 'ERROR: Extra value found'
      ret = .false. ; call clear_dset(set) ; CLOSE(666) ; RETURN
    ENDIF
    CLOSE(666)
    ret = .true.
    RETURN
  END FUNCTION ascii_rd_5d

  SUBROUTINE clr_1d_set(set)
    !! Clear the given 1D data set
    TYPE(DSET1D), INTENT(inout) :: set !! dataset object to clear
    IF (ALLOCATED(set%data)) DEALLOCATE(set%data)
    IF (ALLOCATED(set%x))    DEALLOCATE(set%x)
    RETURN
  END SUBROUTINE clr_1d_set

  SUBROUTINE clr_2d_set(set)
    !! Clear the given 2D data set
    TYPE(DSET2D), INTENT(inout) :: set !! dataset object to clear
    IF (ALLOCATED(set%data)) DEALLOCATE(set%data)
    IF (ALLOCATED(set%x))    DEALLOCATE(set%x)
    IF (ALLOCATED(set%y))    DEALLOCATE(set%y)
    RETURN
  END SUBROUTINE clr_2d_set

  SUBROUTINE clr_3d_set(set)
    !! Clear the given 3D data set
    TYPE(DSET3D), INTENT(inout) :: set !! dataset object to clear
    IF (ALLOCATED(set%data)) DEALLOCATE(set%data)
    IF (ALLOCATED(set%x))    DEALLOCATE(set%x)
    IF (ALLOCATED(set%y))    DEALLOCATE(set%y)
    IF (ALLOCATED(set%z))    DEALLOCATE(set%z)
    RETURN
  END SUBROUTINE clr_3d_set

  SUBROUTINE clr_4d_set(set)
    !! Clear the given 4D data set
    TYPE(DSET4D), INTENT(inout) :: set !! dataset object to clear
    IF (ALLOCATED(set%data)) DEALLOCATE(set%data)
    IF (ALLOCATED(set%x))    DEALLOCATE(set%x)
    IF (ALLOCATED(set%y))    DEALLOCATE(set%y)
    IF (ALLOCATED(set%z))    DEALLOCATE(set%z)
    IF (ALLOCATED(set%t))    DEALLOCATE(set%t)
    RETURN
  END SUBROUTINE clr_4d_set

  SUBROUTINE clr_5d_set(set)
    !! Clear the given 5D data set
    TYPE(DSET5D), INTENT(inout) :: set !! dataset object to clear 
    IF (ALLOCATED(set%data)) DEALLOCATE(set%data)
    IF (ALLOCATED(set%x))    DEALLOCATE(set%x)
    IF (ALLOCATED(set%y))    DEALLOCATE(set%y)
    IF (ALLOCATED(set%z))    DEALLOCATE(set%z)
    IF (ALLOCATED(set%t))    DEALLOCATE(set%t)
    IF (ALLOCATED(set%w))    DEALLOCATE(set%w)
    RETURN
  END SUBROUTINE clr_5d_set

  FUNCTION has_d_1d(dset) RESULT(yes)
    !! Check whether or not the dataset has data.
    TYPE(DSET1D), INTENT(in)     :: dset !! Dataset to check
    LOGICAL                      :: yes  !! return status
    yes =  (ALLOCATED(dset%data).AND. &
            ALLOCATED(dset%x))
  END FUNCTION has_d_1d

  FUNCTION has_d_2d(dset) RESULT(yes)
    !! Check whether or not the dataset has data.
    TYPE(DSET2D), INTENT(in)     :: dset !! Dataset to check
    LOGICAL                      :: yes  !! return status
    yes =  (ALLOCATED(dset%data).AND. &
            ALLOCATED(dset%x)   .AND. &   
            ALLOCATED(dset%y))
  END FUNCTION has_d_2d

  FUNCTION has_d_3d(dset) RESULT(yes)
    !! Check whether or not the dataset has data.
    TYPE(DSET3D), INTENT(in)     :: dset !! Dataset to check
    LOGICAL                      :: yes  !! return status
    yes =  (ALLOCATED(dset%data).AND. &
            ALLOCATED(dset%x)   .AND. &   
            ALLOCATED(dset%y)   .AND. &
            ALLOCATED(dset%z))
  END FUNCTION has_d_3d

  FUNCTION has_d_4d(dset) RESULT(yes)
    !! Check whether or not the dataset has data.
    TYPE(DSET4D), INTENT(in)     :: dset !! Dataset to check
    LOGICAL                      :: yes  !! return status
    yes =  (ALLOCATED(dset%data).AND. &
            ALLOCATED(dset%x)   .AND. &   
            ALLOCATED(dset%y)   .AND. &
            ALLOCATED(dset%z)   .AND. &
            ALLOCATED(dset%t))
  END FUNCTION has_d_4d

  FUNCTION has_d_5d(dset) RESULt(yes)
    !! Check whether or not the dataset has data.
    TYPE(DSET5D), INTENT(in)     :: dset !! Dataset to check
    LOGICAL                      :: yes  !! return status
    yes =  (ALLOCATED(dset%data).AND. &
            ALLOCATED(dset%x)   .AND. &   
            ALLOCATED(dset%y)   .AND. &
            ALLOCATED(dset%z)   .AND. &
            ALLOCATED(dset%t)   .AND. &
            ALLOCATED(dset%w)) 
  END FUNCTION has_d_5d

  FUNCTION nc_wr_1d(path,dset) RESULT(ret)
    CHARACTER(len=*), INTENT(in) :: path !! Path of the netcdf file.
    TYPE(DSET1D), INTENT(in)     :: dset !! Dataset to write
    LOGICAL                      :: ret  !! Return status
    INTEGER                                 :: fi,vi,nd,ds,iret
    INTEGER, DIMENSION(:), ALLOCATABLE      :: di
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tv
    CHARACTER(len=15)                       :: i2s
    LOGICAL                                 :: hxd,ok
    INTEGER                                 :: xid,vid 
    ret = has_data(dset)
    IF (.NOT.ret) RETURN
    ret = .false.
    INQUIRE(FILE=TRIM(path),EXIST=ok)
    IF (ok) THEN
      IF (NF90_OPEN(TRIM(path), NF90_WRITE, fi) /= NF90_NOERR) RETURN
    ELSE
      IF (NF90_CREATE(TRIM(path), NF90_NOCLOBBER, fi) /= NF90_NOERR) RETURN
    ENDIF
    iret = NF90_CLOSE(fi)
    ! if variable already exist: get out of here !
    IF (get_nc_info(path,dset%dname,fi,vi,di)) RETURN

    iret = NF90_OPEN(path,NF90_WRITE,fi) 

    hxd = nc_get_dim_by_name(fi,dset%xname,tv,ds,xid) 
    IF (hxd) THEN
      IF(.NOT.compare(dset%x,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF

    IF (.NOT.hxd) THEN
      ret = nc_set_dim(fi,dset%xname,dset%x,xid)
      IF (.NOT.ret) THEN ; iret= NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    ret = .false.
    iret = nf90_redef(fi)
    iret = nf90_def_var(fi, TRIM(dset%dname), nf90_double,(/xid/),vid) 
    IF (iret /= 0) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    iret = nf90_enddef(fi)
    iret = nf90_put_var(fi,vid,dset%data)
    ret = (iret == 0)
    iret=NF90_CLOSE(fi)
    ret = .true.
    RETURN
  END FUNCTION nc_wr_1d

  FUNCTION nc_wr_2d(path,dset) RESULT(ret)
    CHARACTER(len=*), INTENT(in) :: path !! Path of the netcdf file.
    TYPE(DSET2D), INTENT(in)     :: dset !! Dataset to write
    LOGICAL                      :: ret  !! Return status
    INTEGER                                 :: fi,vi,nd,ds,iret,nferr
    INTEGER, DIMENSION(:), ALLOCATABLE      :: di
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tv
    CHARACTER(len=15)                       :: i2s
    LOGICAL                                 :: hxd,hyd,ok
    INTEGER                                 :: xid,yid,vid 
    ret = has_data(dset)
    IF (.NOT.ret) RETURN
    ret = .false.
    INQUIRE(FILE=TRIM(path),EXIST=ok)
    IF (ok) THEN
      IF (NF90_OPEN(TRIM(path), NF90_WRITE, fi) /= NF90_NOERR) RETURN
    ELSE
      IF (NF90_CREATE(TRIM(path), NF90_NOCLOBBER, fi) /= NF90_NOERR) RETURN
    ENDIF
    iret = NF90_CLOSE(fi)

    ! if variable already exist: get out of here !
    IF (get_nc_info(path,dset%dname,fi,vi,di)) RETURN

    iret = NF90_OPEN(path,NF90_WRITE,fi) 

    hxd = nc_get_dim_by_name(fi,dset%xname,tv,ds,xid) 
    IF (hxd) THEN
      IF(.NOT.compare(dset%x,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    hyd = nc_get_dim_by_name(fi,dset%yname,tv,ds,yid) 
    IF (hyd) THEN
      IF(.NOT.compare(dset%y,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF

    IF (.NOT.hxd) THEN
      ret = nc_set_dim(fi,dset%xname,dset%x,xid)
      IF (.NOT.ret) THEN ; iret= NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.hyd) THEN
      ret = nc_set_dim(fi,dset%yname,dset%y,yid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    ret = .false.
    iret = nf90_redef(fi)
    iret = nf90_def_var(fi, TRIM(dset%dname), nf90_double,(/xid,yid/),vid) 
    IF (iret /= 0) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    iret = nf90_enddef(fi)
    iret = nf90_put_var(fi,vid,dset%data)
    ret = (iret == 0)
    iret=NF90_CLOSE(fi)
    ret = .true.
    RETURN
  END FUNCTION nc_wr_2d

  FUNCTION nc_wr_3d(path,dset) RESULT(ret)
    CHARACTER(len=*), INTENT(in) :: path !! Path of the netcdf file.
    TYPE(DSET3D), INTENT(in)     :: dset !! Dataset to write
    LOGICAL                      :: ret  !! Return status
    INTEGER                                 :: fi,vi,nd,ds,iret
    INTEGER, DIMENSION(:), ALLOCATABLE      :: di
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tv
    CHARACTER(len=15)                       :: i2s
    LOGICAL                                 :: hxd,hyd,hzd,ok
    INTEGER                                 :: xid,yid,zid,vid 
    ret = has_data(dset)
    IF (.NOT.ret) RETURN
    ret = .false.
    INQUIRE(FILE=TRIM(path),EXIST=ok)
    IF (ok) THEN
      IF (NF90_OPEN(TRIM(path), NF90_WRITE, fi) /= NF90_NOERR) RETURN
    ELSE
      IF (NF90_CREATE(TRIM(path), NF90_NOCLOBBER, fi) /= NF90_NOERR) RETURN
    ENDIF
    iret = NF90_CLOSE(fi)

    ! if variable already exist: get out of here !
    IF (get_nc_info(path,dset%dname,fi,vi,di)) RETURN

    iret = NF90_OPEN(path,NF90_WRITE,fi) 

    hxd = nc_get_dim_by_name(fi,dset%xname,tv,ds,xid) 
    IF (hxd) THEN
      IF(.NOT.compare(dset%x,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    hyd = nc_get_dim_by_name(fi,dset%yname,tv,ds,yid) 
    IF (hyd) THEN
      IF(.NOT.compare(dset%y,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    hzd = nc_get_dim_by_name(fi,dset%zname,tv,ds,zid) 
    IF (hzd) THEN 
      IF(.NOT.compare(dset%z,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF

    IF (.NOT.hxd) THEN
      ret = nc_set_dim(fi,dset%xname,dset%x,xid)
      IF (.NOT.ret) THEN ; iret= NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.hyd) THEN
      ret = nc_set_dim(fi,dset%yname,dset%y,yid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.hzd) THEN
      ret = nc_set_dim(fi,dset%zname,dset%z,zid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    ret = .false.
    iret = nf90_redef(fi)
    iret = nf90_def_var(fi, TRIM(dset%dname), nf90_double,(/xid,yid,zid/),vid) 
    IF (iret /= 0) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    iret = nf90_enddef(fi)
    iret = nf90_put_var(fi,vid,dset%data)
    ret = (iret == 0)
    iret=NF90_CLOSE(fi)
    ret = .true.
    RETURN
  END FUNCTION nc_wr_3d

  FUNCTION nc_wr_4d(path,dset) RESULT(ret)
    CHARACTER(len=*), INTENT(in) :: path !! Path of the netcdf file.
    TYPE(DSET4D), INTENT(in)     :: dset !! Dataset to write
    LOGICAL                      :: ret  !! Return status
    INTEGER                                 :: fi,vi,nd,ds,iret
    INTEGER, DIMENSION(:), ALLOCATABLE      :: di
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tv
    CHARACTER(len=15)                       :: i2s
    LOGICAL                                 :: hxd,hyd,hzd,htd,ok
    INTEGER                                 :: xid,yid,zid,tid,vid 
    ret = has_data(dset)
    IF (.NOT.ret) RETURN
    ret = .false.
    INQUIRE(FILE=TRIM(path),EXIST=ok)
    IF (ok) THEN
      IF (NF90_OPEN(TRIM(path), NF90_WRITE, fi) /= NF90_NOERR) RETURN
    ELSE
      IF (NF90_CREATE(TRIM(path), NF90_NOCLOBBER, fi) /= NF90_NOERR) RETURN
    ENDIF
    iret = NF90_CLOSE(fi)

    ! if variable already exist: get out of here !
    IF (get_nc_info(path,dset%dname,fi,vi,di)) RETURN

    iret = NF90_OPEN(path,NF90_WRITE,fi) 

    hxd = nc_get_dim_by_name(fi,dset%xname,tv,ds,xid) 
    IF (hxd) THEN
      IF(.NOT.compare(dset%x,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    hyd = nc_get_dim_by_name(fi,dset%yname,tv,ds,yid) 
    IF (hyd) THEN
      IF(.NOT.compare(dset%y,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    hzd = nc_get_dim_by_name(fi,dset%zname,tv,ds,zid) 
    IF (hzd) THEN 
      IF(.NOT.compare(dset%z,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    htd = nc_get_dim_by_name(fi,dset%tname,tv,ds,tid) 
    IF (htd) THEN
      IF(.NOT.compare(dset%t,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF

    IF (.NOT.hxd) THEN
      ret = nc_set_dim(fi,dset%xname,dset%x,xid)
      IF (.NOT.ret) THEN ; iret= NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.hyd) THEN
      ret = nc_set_dim(fi,dset%yname,dset%y,yid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.hzd) THEN
      ret = nc_set_dim(fi,dset%zname,dset%z,zid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.htd) THEN
      ret = nc_set_dim(fi,dset%tname,dset%t,tid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    ret = .false.
    iret = nf90_redef(fi)
    iret = nf90_def_var(fi, TRIM(dset%dname), nf90_double,(/xid,yid,zid,tid/),vid) 
    IF (iret /= 0) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    iret = nf90_enddef(fi)
    iret = nf90_put_var(fi,vid,dset%data)
    ret = (iret == 0)
    iret=NF90_CLOSE(fi)
    ret = .true.
    RETURN
  END FUNCTION nc_wr_4d

  FUNCTION nc_wr_5d(path,dset) RESULT(ret)
    CHARACTER(len=*), INTENT(in) :: path !! Path of the netcdf file.
    TYPE(DSET5D), INTENT(in)     :: dset !! Dataset to write
    LOGICAL                      :: ret  !! Return status
    INTEGER                                 :: fi,vi,nd,ds,iret
    INTEGER, DIMENSION(:), ALLOCATABLE      :: di
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tv
    CHARACTER(len=15)                       :: i2s
    LOGICAL                                 :: hxd,hyd,hzd,htd,hwd,ok 
    INTEGER                                 :: xid,yid,zid,tid,wid,vid 
    ret = has_data(dset)
    IF (.NOT.ret) RETURN
    ret = .false.
    INQUIRE(FILE=TRIM(path),EXIST=ok)
    IF (ok) THEN
      IF (NF90_OPEN(TRIM(path), NF90_WRITE, fi) /= NF90_NOERR) RETURN
    ELSE
      IF (NF90_CREATE(TRIM(path), NF90_NOCLOBBER, fi) /= NF90_NOERR) RETURN
    ENDIF
    iret = NF90_CLOSE(fi)

    ! if variable already exist: get out of here !
    IF (get_nc_info(path,dset%dname,fi,vi,di)) RETURN

    iret = NF90_OPEN(path,NF90_WRITE,fi) 

    hxd = nc_get_dim_by_name(fi,dset%xname,tv,ds,xid) 
    IF (hxd) THEN
      IF(.NOT.compare(dset%x,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    hyd = nc_get_dim_by_name(fi,dset%yname,tv,ds,yid) 
    IF (hyd) THEN
      IF(.NOT.compare(dset%y,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    hzd = nc_get_dim_by_name(fi,dset%zname,tv,ds,zid) 
    IF (hzd) THEN 
      IF(.NOT.compare(dset%z,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    htd = nc_get_dim_by_name(fi,dset%tname,tv,ds,tid) 
    IF (htd) THEN
      IF(.NOT.compare(dset%t,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF
    hwd = nc_get_dim_by_name(fi,dset%wname,tv,ds,wid) 
    IF (hwd) THEN 
      IF (.NOT.compare(dset%w,tv)) THEN ; ret = .false. ; RETURN ; ENDIF
    ENDIF

    IF (.NOT.hxd) THEN
      ret = nc_set_dim(fi,dset%xname,dset%x,xid)
      IF (.NOT.ret) THEN ; iret= NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.hyd) THEN
      ret = nc_set_dim(fi,dset%yname,dset%y,yid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.hzd) THEN
      ret = nc_set_dim(fi,dset%zname,dset%z,zid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.htd) THEN
      ret = nc_set_dim(fi,dset%tname,dset%t,tid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    IF (.NOT.hwd) THEN
      ret = nc_set_dim(fi,dset%wname,dset%w,wid)
      IF (.NOT.ret) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    ENDIF
    ret = .false.
    iret = nf90_redef(fi)
    iret = nf90_def_var(fi, TRIM(dset%dname), nf90_double,(/xid,yid,zid,tid,wid/),vid) 
    IF (iret /= 0) THEN ; iret = NF90_CLOSE(fi) ; RETURN ; ENDIF
    iret = nf90_enddef(fi)
    iret = nf90_put_var(fi,vid,dset%data)
    ret = (iret == 0)
    iret=NF90_CLOSE(fi)
    ret = .true.
    RETURN
  END FUNCTION nc_wr_5d

  FUNCTION nc_set_dim(fid,name,values,dimid) RESULT(ret)
    !! Set a new dimension (and its associated values) in a NetCDF file.
    INTEGER, INTENT(in)                    :: fid    !! Netcdf file id.
    CHARACTER(len=*), INTENT(in)           :: name   !! Name of the dimension to set.
    REAL(kind=8), DIMENSION(:), INTENT(in) :: values !! Associated values.
    INTEGER, INTENT(out)                   :: dimid  !! Output dimension id.
    LOGICAL                                :: ret    !! Return status
    INTEGER :: iret,vid
    ret = .false.
    iret = nf90_redef(fid)
    iret = nf90_def_dim(fid, TRIM(name), size(values), dimid)
    WRITE(*,'(a)') "nc_set_dim: "//TRIM(name)//" --> "//TRIM(nf90_strerror(iret))
    IF (iret /=0) RETURN
    iret = nf90_def_var(fid, TRIM(name), nf90_double,(/dimid/),vid) 
    IF (iret /=0) RETURN
    iret = nf90_enddef(fid)
    iret = nf90_put_var(fid, vid,values)
    ret = .true.
  END FUNCTION nc_set_dim

  FUNCTION compare(vec1,vec2,tol) RESULT(ret)
    !! Compare two vector of double.
    !!
    !! The test is checked against the difference (element wise) of the two vector compared to
    !! to a given tolerance.
    REAL(kind=8), DIMENSION(:), INTENT(in) :: vec1 !! First vector to compare
    REAL(kind=8), DIMENSION(:), INTENT(in) :: vec2 !! Second vector to compare
    REAL(kind=8), INTENT(in), OPTIONAL     :: tol  !! Tolerance to apply for the comparison
    LOGICAL :: ret                                 !! .true. if both vectors are equal, .false. otherwise.
    REAL(kind=8) :: ztol
    INTEGER      :: i
    ztol = 1d-10 ; IF (PRESENT(tol)) ztol = abs(tol)
    ret = .false.
    IF (size(vec1) /= size(vec2)) RETURN

    DO i=1,size(vec1) ; IF (abs(vec1(i)-vec2(i)) > ztol) RETURN ; ENDDO
    ret = .true.
  END FUNCTION compare



END MODULE DATASETS
