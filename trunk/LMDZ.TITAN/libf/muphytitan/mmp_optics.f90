! Copyright 2017 Université de Reims Champagne-Ardenne 
! Contributor: J. Burgalat (GSMA, URCA)
! email of the author : jeremie.burgalat@univ-reims.fr
! 
! This software is a computer program whose purpose is to compute
! microphysics processes using a two-moments scheme.
! 
! This library is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
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
! knowledge of the CeCILL license and that you accept its terms.

!! file: mmp_optics.f90
!! summary: Interface for YAMMS aerosols optical properties calculations.
!! author: J. Burgalat
!! date: 2017
MODULE MMP_OPTICS
  !! Optical properties of spherical/fractal aerosols using moments
  !!
  !!
  !! The module contains an initialization function, [mmp_optics(module):mmp_init_aer_optics(function)],
  !! that must be called before any calls of the other methods. On failure, it returns .false. and
  !! consequently, all calls to the other methods will fail !
  !!
  !! If openMP is enabled the call to [mmp_optics(module):mmp_init_aer_optics(function)] should be
  !! done by a single thread.
  !!
  !! Then the module provides 4 four public methods to compute optical properties in infrared and
  !! visible channels as a function of moments of the size-distribution:
  !! 
  !! - EXT, the total extinction opacity.
  !! - SSA, the single scattering albedo.
  !! - ASF, the asymetry factor.
  !!
  !! Fractals and spherical aerosols are calculated sperately, but each EXT, SSA and ASF should be added
  !! to get the global optical properties.
  USE MMP_GLOBALS
  USE DATASETS

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mmp_optic_file ! from mmp_globals :)
  PUBLIC :: mmp_initialize_optics
  PUBLIC :: mmp_sph_optics_vis,mmp_sph_optics_ir
  PUBLIC :: mmp_fra_optics_vis,mmp_fra_optics_ir

  ! OPTICAL PROPERTIES !

  !> Extinction opacty table (spherical,IR).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: ext_s_i
  !> Single scattering albedo table (spherical,IR).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: sca_s_i
  !> Asymetry factor table (spherical,IR).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: asf_s_i
  !> Extinction opacty table (fractal,IR).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: ext_f_i
  !> Single scattering albedo table (fractal,IR).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: sca_f_i
  !> Asymetry factor table (fractal,IR).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: asf_f_i

  !> Extinction opacty table (spherical,VIS).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: ext_s_v
  !> Single scattering albedo table (spherical,VIS).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: sca_s_v
  !> Asymetry factor table (spherical,VIS).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: asf_s_v
  !> Extinction opacty table (fractal,VIS).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: ext_f_v
  !> Single scattering albedo table (fractal,VIS).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: sca_f_v
  !> Asymetry factor table (fractal,VIS).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: asf_f_v


  INTEGER, SAVE      :: mmp_nrc    = -1  !! Size of radius grid.

  !> Characteristic radius grid
  REAL(kind=8), DIMENSION(:), ALLOCATABLE, SAVE :: mmp_rc

  CONTAINS

  SUBROUTINE mmp_initialize_optics(path)
    !! Initialize optics data for aerosols optical properties computation.
    !!
    !! @note
    !! If the subroutine fails to initialize parameters, the run is aborted.
    !!
    !! @warning
    !! The method assumes YAMMS model has been already intialized correctly !
    !!
    !! @warning
    !! If OpenMP is activated, this subroutine must be called in an $OMP SINGLE statement as it 
    !! initializes global variables that are not thread private.
    !!
    !! '''
    !! !$OMP SINGLE
    !! call mmp_initialize(...)
    !! !$OMP END SINGLE
    CHARACTER(len=*), INTENT(in) :: path !! Path of NetCDF look-up tables file.
    LOGICAL :: ret
    WRITE(*,'(a)') "*** mmp_init_aer_optics speaking ***"
    WRITE(*,'(a)') "I'm about to initialize look-up tables of aerosols optical properties."
    WRITE(*,'(a)') "If something's wrong... I will abort the program !"
    IF (.NOT.mm_ini) THEN
      call abort_program(error("[mmp_init_aer_optics] Too bad mmp_initialize has not been called yet !",2))
    ENDIF
    !  look-up tables
    ret = read_lookup_tables(path)
    IF (.NOT.ret) &
      call abort_program(error("[mmp_init_aer_optics] Failed to retrieve data.",2))
  END SUBROUTINE mmp_initialize_optics

  FUNCTION mmp_sph_optics_vis(M0,M3,iwn,ext,sca,ssa,asf) RESULT(ret)
    !! Compute optical properties of the spherical mode in the visible range.
    REAL(kind=mm_wp), INTENT(in)  :: M0  !! 0th order moment of the spherical mode (m-2).
    REAL(kind=mm_wp), INTENT(in)  :: M3  !! 3rd order moment of the spherical mode (m3.m-2).
    INTEGER, INTENT(in)           :: iwn !! Index of the wavenumber to compute
    REAL(kind=mm_wp), INTENT(out) :: ext !! Extinction opacity.
    REAL(kind=mm_wp), INTENT(out) :: sca !! Scattering.
    REAL(kind=mm_wp), INTENT(out) :: ssa !! Single scattering albedo. 
    REAL(kind=mm_wp), INTENT(out) :: asf !! Asymetry factor.
    LOGICAL                       :: ret !! true on success, false otherwise.
    INTEGER          :: i,ridx
    REAL(kind=mm_wp) :: rc1,rc2,rx,rc
    ret = .false.
    IF (mmp_nrc == -1) RETURN
    ret = .true.
    rc = mm_get_rcs(M0,M3)
    ridx = mmp_nrc
    DO i=1, mmp_nrc
      IF (rc < mmp_rc(i)) THEN
        ridx = i-1
        EXIT
      ENDIF
    ENDDO
    IF (ridx == 0) THEN
      ! out of range lower bound
      ext = ext_s_v(1,iwn)
      sca = sca_s_v(1,iwn)
      asf = asf_s_v(1,iwn)
      ssa = sca/ext
    ELSE IF (ridx == mmp_nrc) THEN
      ! out of range upper bound
      ext = ext_s_v(mmp_nrc,iwn)
      sca = sca_s_v(mmp_nrc,iwn)
      asf = asf_s_v(mmp_nrc,iwn)
      ssa = sca/ext
    ELSE
      ! in range: interpolate
      rc1 = mmp_rc(ridx) ; rc2 = mmp_rc(ridx+1)
      rx = (rc-rc1)/(rc2-rc1)
      ext = exp(log(ext_s_v(ridx,iwn))*(1d0-rx) + log(ext_s_v(ridx+1,iwn))*rx)
      sca = exp(log(sca_s_v(ridx,iwn))*(1d0-rx) + log(sca_s_v(ridx+1,iwn))*rx)
      asf = asf_s_v(ridx,iwn)*(1d0-rx) + asf_s_v(ridx+1,iwn)*rx
      ssa = sca/ext
    ENDIF
    ! scale by M0
    ext = ext * M0
    RETURN
  END FUNCTION mmp_sph_optics_vis

  FUNCTION mmp_sph_optics_ir(M0,M3,iwn,ext,sca,ssa,asf) RESULT(ret)
    !! Compute optical properties of the spherical mode in the infra-red range.
    REAL(kind=mm_wp), INTENT(in)  :: M0  !! 0th order moment of the spherical mode (m-2).
    REAL(kind=mm_wp), INTENT(in)  :: M3  !! 3rd order moment of the spherical mode (m3.m-2).
    INTEGER, INTENT(in)           :: iwn !! Index of the wavenumber to compute
    REAL(kind=mm_wp), INTENT(out) :: ext !! Extinction opacity.
    REAL(kind=mm_wp), INTENT(out) :: sca !! Scattering.
    REAL(kind=mm_wp), INTENT(out) :: ssa !! Single scattering albedo. 
    REAL(kind=mm_wp), INTENT(out) :: asf !! Asymetry factor.
    LOGICAL                       :: ret !! true on success, false otherwise.
    INTEGER          :: i,ridx
    REAL(kind=mm_wp) :: rc1,rc2,rx,rc
    ret = .false.
    IF (mmp_nrc == -1) RETURN
    ret = .true.
    rc = mm_get_rcs(M0,M3)
    ridx = mmp_nrc
    DO i=1, mmp_nrc
      IF (rc < mmp_rc(i)) THEN
        ridx = i-1
        EXIT
      ENDIF
    ENDDO
    IF (ridx == 0) THEN
      ! out of range lower bound
      ext = ext_s_i(1,iwn)
      sca = sca_s_i(1,iwn)
      asf = asf_s_i(1,iwn)
      ssa = sca/ext
    ELSE IF (ridx == mmp_nrc) THEN
      ! out of range upper bound
      ext = ext_s_i(mmp_nrc,iwn)
      sca = sca_s_i(mmp_nrc,iwn)
      asf = asf_s_i(mmp_nrc,iwn)
      ssa = sca/ext
    ELSE
      ! in range: interpolate
      rc1 = mmp_rc(ridx) ; rc2 = mmp_rc(ridx+1)
      rx = (rc-rc1)/(rc2-rc1)
      ext = exp(log(ext_s_i(ridx,iwn))*(1d0-rx) + log(ext_s_i(ridx+1,iwn))*rx)
      sca = exp(log(sca_s_i(ridx,iwn))*(1d0-rx) + log(sca_s_i(ridx+1,iwn))*rx)
      asf = asf_s_i(ridx,iwn)*(1d0-rx) + asf_s_i(ridx+1,iwn)*rx
      ssa = sca/ext
    ENDIF
    ! scale by M0
    ext = ext * M0
    RETURN
  END FUNCTION mmp_sph_optics_ir

  FUNCTION mmp_fra_optics_vis(M0,M3,iwn,ext,sca,ssa,asf) RESULT(ret)
    !! Compute optical properties of the spherical mode in the visible range.
    REAL(kind=mm_wp), INTENT(in)  :: M0  !! 0th order moment of the fractal mode (m-2).
    REAL(kind=mm_wp), INTENT(in)  :: M3  !! 3rd order moment of the fractal mode (m3.m-2).
    INTEGER, INTENT(in)           :: iwn !! Index of the wavenumber to compute.
    REAL(kind=mm_wp), INTENT(out) :: ext !! Extinction opacity.
    REAL(kind=mm_wp), INTENT(out) :: sca !! Scattering.
    REAL(kind=mm_wp), INTENT(out) :: ssa !! Single scattering albedo. 
    REAL(kind=mm_wp), INTENT(out) :: asf !! Asymetry factor.
    LOGICAL                       :: ret !! true on success, false otherwise.
    INTEGER          :: i,ridx
    REAL(kind=mm_wp) :: rc1,rc2,rx,rc
    ret = .false.
    IF (mmp_nrc == -1) RETURN
    ret = .true.
    rc = mm_get_rcs(M0,M3)
    ridx = mmp_nrc
    DO i=1, mmp_nrc
      IF (rc < mmp_rc(i)) THEN
        ridx = i-1
        EXIT
      ENDIF
    ENDDO
    IF (ridx == 0) THEN
      ! out of range lower bound
      ext = ext_f_v(1,iwn)
      sca = sca_f_v(1,iwn)
      asf = asf_f_v(1,iwn)
      ssa = sca/ext
    ELSE IF (ridx == mmp_nrc) THEN
      ! out of range upper bound
      ext = ext_f_v(mmp_nrc,iwn)
      sca = sca_f_v(mmp_nrc,iwn)
      asf = asf_f_v(mmp_nrc,iwn)
      ssa = sca/ext
    ELSE
      ! in range: interpolate
      rc1 = mmp_rc(ridx) ; rc2 = mmp_rc(ridx+1)
      rx = (rc-rc1)/(rc2-rc1)
      ext = exp(log(ext_f_v(ridx,iwn))*(1d0-rx) + log(ext_f_v(ridx+1,iwn))*rx)
      sca = exp(log(sca_f_v(ridx,iwn))*(1d0-rx) + log(sca_f_v(ridx+1,iwn))*rx)
      asf = asf_f_v(ridx,iwn)*(1d0-rx) + asf_f_v(ridx+1,iwn)*rx
      ssa = sca/ext
    ENDIF
    ! scale by M0
    ext = ext * M0
    RETURN
  END FUNCTION mmp_fra_optics_vis

  FUNCTION mmp_fra_optics_ir(M0,M3,iwn,ext,sca,ssa,asf) RESULT(ret)
    !! Compute optical properties of the spherical mode in the infra-red range.
    REAL(kind=mm_wp), INTENT(in)  :: M0  !! 0th order moment of the spherical mode (m-2).
    REAL(kind=mm_wp), INTENT(in)  :: M3  !! 3rd order moment of the spherical mode (m3.m-2).
    INTEGER, INTENT(in)           :: iwn !! Index of the wavenumber to compute
    REAL(kind=mm_wp), INTENT(out) :: ext !! Extinction opacity.
    REAL(kind=mm_wp), INTENT(out) :: sca !! Scattering.
    REAL(kind=mm_wp), INTENT(out) :: ssa !! Single scattering albedo. 
    REAL(kind=mm_wp), INTENT(out) :: asf !! Asymetry factor.
    LOGICAL                       :: ret !! true on success, false otherwise.
    INTEGER          :: i,ridx
    REAL(kind=mm_wp) :: rc1,rc2,rx,rc
    ret = .false.
    IF (mmp_nrc == -1) RETURN
    ret = .true.
    rc = mm_get_rcs(M0,M3)
    ridx = mmp_nrc
    DO i=1, mmp_nrc
      IF (rc < mmp_rc(i)) THEN
        ridx = i-1
        EXIT
      ENDIF
    ENDDO
    IF (ridx == 0) THEN
      ! out of range lower bound
      ext = ext_f_i(1,iwn)
      sca = sca_f_i(1,iwn)
      asf = asf_f_i(1,iwn)
      ssa = sca/ext
    ELSE IF (ridx == mmp_nrc) THEN
      ! out of range upper bound
      ext = ext_f_i(mmp_nrc,iwn)
      sca = sca_f_i(mmp_nrc,iwn)
      asf = asf_f_i(mmp_nrc,iwn)
      ssa = sca/ext
    ELSE
      ! in range: interpolate
      rc1 = mmp_rc(ridx) ; rc2 = mmp_rc(ridx+1)
      rx = (rc-rc1)/(rc2-rc1)
      ext = exp(log(ext_f_i(ridx,iwn))*(1d0-rx) + log(ext_f_i(ridx+1,iwn))*rx)
      sca = exp(log(sca_f_i(ridx,iwn))*(1d0-rx) + log(sca_f_i(ridx+1,iwn))*rx)
      asf = asf_f_i(ridx,iwn)*(1d0-rx) + asf_f_i(ridx+1,iwn)*rx
      ssa = sca/ext
    ENDIF
    ! scale by M0
    ext = ext * M0
    RETURN
  END FUNCTION mmp_fra_optics_ir

  FUNCTION read_lookup_tables(path) RESULT(ret)
    !! Read look-up tables.
    CHARACTER(len=*), INTENT(in) :: path !! Path of the look-up tables netcdf file.
    LOGICAL                      :: ret  !! .true. on success, .false. otherwise.
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: sigmas_vi,sigmas_ir

    TYPE(DSET2D) :: dset
    ! data(nrc,ni|nv)
    ! INFRARED
    ret = read_dset(path,"ext_s_i",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'ext_s_i' table" ; RETURN ; ENDIF 
    ext_s_i = dset%data ; sigmas_ir = dset%y ; mmp_rc = dset%x
    ret = read_dset(path,"ext_f_i",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'ext_f_i' table" ; RETURN ; ENDIF 
    ext_f_i = dset%data 
    ret = read_dset(path,"sca_s_i",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'sca_s_i' table" ; RETURN ; ENDIF 
    sca_s_i = dset%data 
    ret = read_dset(path,"sca_f_i",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'sca_f_i' table" ; RETURN ; ENDIF 
    sca_f_i = dset%data 
    ret = read_dset(path,"asf_s_i",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'asf_s_i' table" ; RETURN ; ENDIF 
    asf_s_i = dset%data 
    ret = read_dset(path,"asf_f_i",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'asf_f_i' table" ; RETURN ; ENDIF 
    asf_f_i = dset%data 
    ! VISIBLE
    ret = read_dset(path,"ext_s_v",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'ext_s_v' table" ; RETURN ; ENDIF 
    ext_s_v = dset%data ; sigmas_vi = dset%y
    ret = read_dset(path,"ext_f_v",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'ext_f_v' table" ; RETURN ; ENDIF 
    ext_f_v = dset%data 
    ret = read_dset(path,"sca_s_v",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'sca_s_v' table" ; RETURN ; ENDIF 
    sca_s_v = dset%data 
    ret = read_dset(path,"sca_f_v",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'sca_f_v' table" ; RETURN ; ENDIF 
    sca_f_v = dset%data 
    ret = read_dset(path,"asf_s_v",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'asf_s_v' table" ; RETURN ; ENDIF 
    asf_s_v = dset%data 
    ret = read_dset(path,"asf_f_v",dset)
    IF (.NOT.ret) THEN ; WRITE(*,'(a)') "[read_tables] cannot read 'asf_f_v' table" ; RETURN ; ENDIF 
    asf_f_v = dset%data 
    mmp_nrc = SIZE(mmp_rc)
    ret = .true.
    RETURN
  END FUNCTION read_lookup_tables

END MODULE MMP_OPTICS

