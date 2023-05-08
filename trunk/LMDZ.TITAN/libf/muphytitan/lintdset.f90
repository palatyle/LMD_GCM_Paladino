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

!! file: lintdset.f90
!! summary: Linear interpolation generic interfaces
!! author: burgalat
!! date: 2010-2014 


MODULE LINTDSET
  !! Generic linear interpolation using simple [[datasets(module)]]
  USE LINT_PREC
  USE LINTCORE
  USE DATASETS
  IMPLICIT NONE

  PUBLIC  :: lint_dset, hdcd_lint_dset,              &
             dset1d, dset2d, dset3d, dset4d, dset5d, &
             read_dset, clear_dset, is_in
  PRIVATE 

  !> Interface to multi-variate linear interpolation using data set
  !!
  !! For __n__, the number of variables, user must provide:
  !!
  !! - __n__ scalar(s)/vector(s) with the coordinate(s) of the point(s) to 
  !!   interpolate.
  !! - A single __n__-D data set with the tabulated function.
  !! - A function, satisfying the [[lintcore(module):locate(interface)]] interface, that
  !!   should locate the neighbourhood of the point to interpolate in the given data
  !!   set. 
  !! - Either a scalar or vector which stores the interpolated values.
  !! 
  !! This interface also provides vectorized version of the linear
  !! interpolation. In this case, input coordinates should be vector of the same
  !! size as well as the output argument.
  INTERFACE lint_dset
    MODULE PROCEDURE l1d_dset_sc, l2d_dset_sc, l3d_dset_sc, l4d_dset_sc, &
                     l5d_dset_sc 
    MODULE PROCEDURE l1d_dset_ve, l2d_dset_ve, l3d_dset_ve, l4d_dset_ve, &
                     l5d_dset_ve 
  END INTERFACE 

  !> Interface to multi-variate hard-coded linear interpolation using data set
  !!
  !! Same remarks as in [[lintdset(module):lint_dset(interface)]] apply here.
  INTERFACE hdcd_lint_dset
    MODULE PROCEDURE hdcdl1d_dset_sc, hdcdl2d_dset_sc, hdcdl3d_dset_sc, &
                     hdcdl4d_dset_sc, hdcdl5d_dset_sc
    MODULE PROCEDURE hdcdl1d_dset_ve, hdcdl2d_dset_ve, hdcdl3d_dset_ve, &
                     hdcdl4d_dset_ve, hdcdl5d_dset_ve
  END INTERFACE

  CONTAINS

  FUNCTION l1d_dset_sc(x,set,locator,res) RESULT(ret)
    !! Wrapper interface for 1D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    TYPE(dset1d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(2,2) :: g1d 
    INTEGER                       :: i,ix
    ret = .false.
    res = HUGE(res)
    ix = locator(x,set%x) ; IF (ix == 0) RETURN
    DO i=0,1
      g1d(i+1,1) = set%x(ix+i)
      g1d(i+1,2) = set%data(ix+i)
    ENDDO
    res = lintc_((/x/),g1d) 
    ret = .true.
  END FUNCTION l1d_dset_sc

  FUNCTION l2d_dset_sc(x,y,set,locator,res) RESULT(ret)
    !! Wrapper interface for 2D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: y       !! Y coordinate of the point to interpolate
    TYPE(dset2d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(4,3) :: g2d 
    INTEGER                       :: i,ix0,iy0,a,b
    ret = .false.
    res = HUGE(res)
    ix0 = locator(x,set%x) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,set%y) ; IF (iy0 == 0) RETURN
    DO i=1,4
      a=ix0+MOD((i-1),2)   ; g2d(i,1) = set%x(a)
      b=iy0+MOD((i-1)/2,2) ; g2d(i,2) = set%y(b) 
      g2d(i,3) = set%data(a,b)
    ENDDO
    res = lintc_((/x,y/),g2d) 
    ret = .true.
  END FUNCTION l2d_dset_sc

  FUNCTION l3d_dset_sc(x,y,z,set,locator,res) RESULT(ret)
    !! Wrapper interface for 3D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: z       !! Z coordinate of the point to interpolate
    TYPE(dset3d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(8,4) :: g3d 
    INTEGER                       :: i,ix0,iy0,iz0,a,b,c
    ret = .false.
    res = HUGE(res)
    ix0 = locator(x,set%x) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,set%y) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,set%z) ; IF (iz0 == 0) RETURN
    DO i=1,8
      a=ix0+MOD((i-1),2)   ; g3d(i,1) = set%x(a)
      b=iy0+MOD((i-1)/2,2) ; g3d(i,2) = set%y(b) 
      c=iz0+MOD((i-1)/4,2) ; g3d(i,3) = set%z(c)
      g3d(i,4) = set%data(a,b,c)
    ENDDO
    res = lintc_((/x,y,z/),g3d) 
    ret = .true.
  END FUNCTION l3d_dset_sc

  FUNCTION l4d_dset_sc(x,y,z,t,set,locator,res) RESULT(ret)
    !! Wrapper interface for 4D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: t       !! T coordinate of the point to interpolate
    TYPE(dset4d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(16,5) :: g4d 
    INTEGER                        :: i,ix0,iy0,iz0,it0,a,b,c,d
    ret = .false.
    res = HUGE(res)
    ix0 = locator(x,set%x) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,set%y) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,set%z) ; IF (iz0 == 0) RETURN
    it0 = locator(t,set%t) ; IF (it0 == 0) RETURN
    DO i=1,16
      a=ix0+MOD((i-1),2)   ; g4d(i,1) = set%x(a)
      b=iy0+MOD((i-1)/2,2) ; g4d(i,2) = set%y(b) 
      c=iz0+MOD((i-1)/4,2) ; g4d(i,3) = set%z(c)
      d=it0+MOD((i-1)/8,2) ; g4d(i,4) = set%t(d)
      g4d(i,5) = set%data(a,b,c,d)
    ENDDO
    res = lintc_((/x,y,z,t/),g4d) 
    ret = .true.
  END FUNCTION l4d_dset_sc

  FUNCTION l5d_dset_sc(x,y,z,t,w,set,locator,res) RESULT(ret)
    !! Wrapper interface for 5D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: t       !! T coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: w       !! W coordinate of the point to interpolate
    TYPE(dset5d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(32,6) :: g5d 
    INTEGER                        :: i,ix0,iy0,iz0,it0,iw0,a,b,c,d,e
    ret = .false.
    res = HUGE(res)
    ix0 = locator(x,set%x) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,set%y) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,set%z) ; IF (iz0 == 0) RETURN
    it0 = locator(t,set%t) ; IF (it0 == 0) RETURN
    iw0 = locator(w,set%w) ; IF (iw0 == 0) RETURN
    DO i=1,32
      a=ix0+MOD((i-1),2)    ; g5d(i,1) = set%x(a)
      b=iy0+MOD((i-1)/2,2)  ; g5d(i,2) = set%y(b) 
      c=iz0+MOD((i-1)/4,2)  ; g5d(i,3) = set%z(c)
      d=it0+MOD((i-1)/8,2)  ; g5d(i,4) = set%t(d)
      e=iw0+MOD((i-1)/16,2) ; g5d(i,5) = set%w(e)
      g5d(i,6) = set%data(a,b,c,d,e)
    ENDDO
    res = lintc_((/x,y,z,t,w/),g5d) 
    ret = .true.
  END FUNCTION l5d_dset_sc

  FUNCTION l1d_dset_ve(x,set,locator,res) RESULT(ret)
    !! Wrapper interface for 1D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    TYPE(dset1d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(2,2) :: g1d 
    INTEGER                       :: nv,i,j,ix0
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x) 
      IF (ix0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      DO i=0,1
        g1d(i+1,1) = set%x(ix0+i)
        g1d(i+1,2) = set%data(ix0+i)
      ENDDO
      res(j) = lintc_((/x(j)/),g1d) 
    ENDDO
    ret = .true.
  END FUNCTION l1d_dset_ve

  FUNCTION l2d_dset_ve(x,y,set,locator,res) RESULT(ret)
    !! Wrapper interface for 2D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    TYPE(dset2d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(4,3) :: g2d 
    INTEGER                       :: nv,i,j,ix0,iy0,a,b
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x) 
      iy0 = locator(y(j),set%y)
      IF (ix0 == 0 .OR. iy0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      DO i=1,4
        a=ix0+MOD((i-1),2)   ; g2d(i,1) = set%x(a)
        b=iy0+MOD((i-1)/2,2) ; g2d(i,2) = set%y(b) 
        g2d(i,3) = set%data(a,b)
      ENDDO
      res(j) = lintc_((/x(j),y(j)/),g2d) 
    ENDDO
    ret = .true.
  END FUNCTION l2d_dset_ve

  FUNCTION l3d_dset_ve(x,y,z,set,locator,res) RESULT(ret)
    !! Wrapper interface for 3D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    TYPE(dset3d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(8,4) :: g3d 
    INTEGER                       :: nv,i,j,ix0,iy0,iz0,a,b,c
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x) 
      iy0 = locator(y(j),set%y)
      iz0 = locator(z(j),set%z)
      IF (ix0==0 .OR. iy0==0 .OR. iz0==0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      DO i=1,8
        a=ix0+MOD((i-1),2)   ; g3d(i,1) = set%x(a)
        b=iy0+MOD((i-1)/2,2) ; g3d(i,2) = set%y(b) 
        c=iz0+MOD((i-1)/4,2) ; g3d(i,3) = set%z(c)
        g3d(i,4) = set%data(a,b,c)
      ENDDO
      res(j) = lintc_((/x(j),y(j),z(j)/),g3d) 
    ENDDO
    ret = .true.
  END FUNCTION l3d_dset_ve

  FUNCTION l4d_dset_ve(x,y,z,t,set,locator,res) RESULT(ret)
    !! Wrapper interface for 4D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: t       !! T coordinate of the points to interpolate
    TYPE(dset4d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(16,5) :: g4d 
    INTEGER                        :: nv,i,j,ix0,iy0,iz0,it0,a,b,c,d
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x)
      iy0 = locator(y(j),set%y)
      iz0 = locator(z(j),set%z)
      it0 = locator(t(j),set%t)
      IF (ix0==0 .OR. iy0==0 .OR. iz0==0 .OR. it0==0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      DO i=1,16
        a=ix0+MOD((i-1),2)   ; g4d(i,1) = set%x(a)
        b=iy0+MOD((i-1)/2,2) ; g4d(i,2) = set%y(b) 
        c=iz0+MOD((i-1)/4,2) ; g4d(i,3) = set%z(c)
        d=it0+MOD((i-1)/8,2) ; g4d(i,4) = set%t(d)
        g4d(i,5) = set%data(a,b,c,d)
      ENDDO
      res(j) = lintc_((/x(j),y(j),z(j),t(j)/),g4d) 
    ENDDO
    ret = .true.
  END FUNCTION l4d_dset_ve

  FUNCTION l5d_dset_ve(x,y,z,t,w,set,locator,res) RESULT(ret)
    !! Wrapper interface for 5D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: t       !! T coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: w       !! W coordinate of the points to interpolate
    TYPE(dset5d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(32,6) :: g5d 
    INTEGER                        :: nv,i,j,ix0,iy0,iz0,it0,iw0,a,b,c,d,e
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x)
      iy0 = locator(y(j),set%y)
      iz0 = locator(z(j),set%z)
      it0 = locator(t(j),set%t)
      iw0 = locator(w(j),set%w)
      IF (ix0==0 .OR. iy0==0 .OR. iz0==0 .OR. it0==0 .OR. iw0==0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      DO i=1,32
        a=ix0+MOD((i-1),2)    ; g5d(i,1) = set%x(a)
        b=iy0+MOD((i-1)/2,2)  ; g5d(i,2) = set%y(b) 
        c=iz0+MOD((i-1)/4,2)  ; g5d(i,3) = set%z(c)
        d=it0+MOD((i-1)/8,2)  ; g5d(i,4) = set%t(d)
        e=iw0+MOD((i-1)/16,2) ; g5d(i,5) = set%w(e)
        g5d(i,6) = set%data(a,b,c,d,e)
      ENDDO
      res(j) = lintc_((/x(j),y(j),z(j),t(j),w(j)/),g5d) 
    ENDDO
    ret = .true.
  END FUNCTION l5d_dset_ve

  !--------------------!
  ! HARD CODED VERSION !
  !--------------------!

  FUNCTION hdcdl1d_dset_sc(x,set, locator,res) RESULT(ret)
    !! Hard-coded for 1D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    TYPE(dset1d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd
    INTEGER       :: ix0,ix1
    ret = .false.
    res = HUGE(res)
    ix0 = locator(x,set%x) ; IF (ix0 == 0) RETURN
    ix1 = ix0+1
    xd = (x-set%x(ix0))/(set%x(ix1)-set%x(ix0))
    res = set%data(ix0)*(1d0-xd)+set%data(ix1)*xd
    ret = .true.
  END FUNCTION hdcdl1d_dset_sc

  FUNCTION hdcdl2d_dset_sc(x,y,set,locator,res) RESULT(ret)
    !! Hard-coded for 2D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: y       !! Y coordinate of the point to interpolate
    TYPE(dset2d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd
    REAL(kind=wp) :: f0,f1
    INTEGER       :: ix0,iy0,ix1,iy1
    ret = .false.
    res = HUGE(res)
    ix0 = locator(x,set%x) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,set%y) ; IF (iy0 == 0) RETURN
    ix1 = ix0 + 1 ; iy1 = iy0 + 1
    xd  = (x-set%x(ix0))/(set%x(ix1)-set%x(ix0))
    yd  = (y-set%y(iy0))/(set%y(iy1)-set%y(iy0))
    f0  = set%data(ix0,iy0)*(1d0-xd)+set%data(ix1,iy0)*xd
    f1  = set%data(ix0,iy1)*(1d0-xd)+set%data(ix1,iy1)*xd
    res = f0*(1d0-yd)+f1*yd
    ret = .true.
  END FUNCTION hdcdl2d_dset_sc 

  FUNCTION hdcdl3d_dset_sc(x,y,z,set,locator,res) RESULT(ret)
    !! Hard-coded for 3D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: z       !! Z coordinate of the point to interpolate
    TYPE(dset3d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd
    REAL(kind=wp) :: f00,f10,f01,f11,f0,f1
    INTEGER       :: ix0,iy0,iz0,ix1,iy1,iz1
    ret = .false.
    res = HUGE(res)
    ix0 = locator(x,set%x) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,set%y) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,set%z) ; IF (iz0 == 0) RETURN
    ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1
    xd = (x-set%x(ix0))/(set%x(ix1)-set%x(ix0))
    yd = (y-set%y(iy0))/(set%y(iy1)-set%y(iy0))
    zd = (z-set%z(iz0))/(set%z(iz1)-set%z(iz0))

    f00 = set%data(ix0,iy0,iz0)*(1d0-xd)+set%data(ix1,iy0,iz0)*xd
    f10 = set%data(ix0,iy1,iz0)*(1d0-xd)+set%data(ix1,iy1,iz0)*xd
    f01 = set%data(ix0,iy0,iz1)*(1d0-xd)+set%data(ix1,iy0,iz1)*xd
    f11 = set%data(ix0,iy1,iz1)*(1d0-xd)+set%data(ix1,iy1,iz1)*xd

    f0 = f00 *(1d0-yd)+f10*yd
    f1 = f00 *(1d0-yd)+f10*yd

    res = f0*(1d0-zd)+f1*zd
    ret = .true.
  END FUNCTION hdcdl3d_dset_sc

  FUNCTION hdcdl4d_dset_sc(x,y,z,t,set,locator,res) RESULT(ret)
    !! Hard-coded for 4D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: t       !! T coordinate of the point to interpolate
    TYPE(dset4d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd,td
    REAL(kind=wp) :: f000,f100,f010,f110,f001,f101,f011,f111, &
                     f00,f10,f01,f11,f0,f1
    INTEGER       :: ix0,iy0,iz0,it0,ix1,iy1,iz1,it1
    ret = .false.
    res = HUGE(res)
    ix0 = locator(x,set%x) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,set%y) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,set%z) ; IF (iz0 == 0) RETURN
    it0 = locator(t,set%t) ; IF (it0 == 0) RETURN
    
    ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1 ; it1=it0+1
    xd = (x-set%x(ix0))/(set%x(ix1)-set%x(ix0))
    yd = (y-set%y(iy0))/(set%y(iy1)-set%y(iy0))
    zd = (z-set%z(iz0))/(set%z(iz1)-set%z(iz0))
    td = (t-set%t(it0))/(set%t(it1)-set%t(it0))

    f000 = set%data(ix0,iy0,iz0,it0)*(1d0-xd)+set%data(ix1,iy0,iz0,it0)*xd
    f100 = set%data(ix0,iy1,iz0,it0)*(1d0-xd)+set%data(ix1,iy1,iz0,it0)*xd
    f010 = set%data(ix0,iy0,iz1,it0)*(1d0-xd)+set%data(ix1,iy0,iz1,it0)*xd
    f110 = set%data(ix0,iy1,iz1,it0)*(1d0-xd)+set%data(ix1,iy1,iz1,it0)*xd
    f001 = set%data(ix0,iy0,iz0,it1)*(1d0-xd)+set%data(ix1,iy0,iz0,it1)*xd
    f101 = set%data(ix0,iy1,iz0,it1)*(1d0-xd)+set%data(ix1,iy1,iz0,it1)*xd
    f011 = set%data(ix0,iy0,iz1,it1)*(1d0-xd)+set%data(ix1,iy0,iz1,it1)*xd
    f111 = set%data(ix0,iy1,iz1,it1)*(1d0-xd)+set%data(ix1,iy1,iz1,it1)*xd

    f00 = f000*(1d0-yd)+f100*yd
    f10 = f010*(1d0-yd)+f110*yd 
    f01 = f001*(1d0-yd)+f101*yd 
    f11 = f011*(1d0-yd)+f111*yd 

    f0 = f00 *(1d0-zd)+f10*zd
    f1 = f01 *(1d0-zd)+f11*zd

    res = f0*(1d0-td)+f1*td
    ret = .true.
  END FUNCTION hdcdl4d_dset_sc

  FUNCTION hdcdl5d_dset_sc(x,y,z,t,w,set,locator,res) RESULT(ret)
    !! Hard-coded for 5D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)  :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: t       !! T coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)  :: w       !! W coordinate of the point to interpolate
    TYPE(dset5d), INTENT(in)   :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out) :: res     !! Interpolated value
    PROCEDURE(locate)          :: locator !! Locator function
    LOGICAL :: ret                        !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd,td,wd
    REAL(kind=wp) :: f0000,f1000,f0100,f1100,f0010,f1010,f0110,f1110, & 
                     f0001,f1001,f0101,f1101,f0011,f1011,f0111,f1111, &
                     f000,f100,f010,f110,f001,f101,f011,f111,         &
                     f00,f10,f01,f11,f0,f1
    INTEGER       :: ix0,iy0,iz0,it0,iw0,ix1,iy1,iz1,it1,iw1
    ret = .false.
    res = HUGE(res)
    ix0 = locator(x,set%x) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,set%y) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,set%z) ; IF (iz0 == 0) RETURN
    it0 = locator(t,set%t) ; IF (it0 == 0) RETURN
    iw0 = locator(w,set%w) ; IF (iw0 == 0) RETURN
    
    ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1 
    xd = (x-set%x(ix0))/(set%x(ix1)-set%x(ix0))
    yd = (y-set%y(iy0))/(set%y(iy1)-set%y(iy0))
    zd = (z-set%z(iz0))/(set%z(iz1)-set%z(iz0))
    td = (t-set%t(it0))/(set%t(it1)-set%t(it0))
    wd = (w-set%w(it0))/(set%w(it1)-set%w(it0))

    f0000 = set%data(ix0,iy0,iz0,it0,iw0)*(1d0-xd)+set%data(ix1,iy0,iz0,it0,iw0)*xd
    f1000 = set%data(ix0,iy1,iz0,it0,iw0)*(1d0-xd)+set%data(ix0,iy1,iz0,it0,iw0)*xd
    f0100 = set%data(ix0,iy0,iz1,it0,iw0)*(1d0-xd)+set%data(ix1,iy0,iz1,it0,iw0)*xd
    f1100 = set%data(ix0,iy1,iz1,it0,iw0)*(1d0-xd)+set%data(ix1,iy1,iz1,it0,iw0)*xd
    f0010 = set%data(ix0,iy0,iz0,it1,iw0)*(1d0-xd)+set%data(ix1,iy0,iz0,it1,iw0)*xd
    f1010 = set%data(ix0,iy1,iz0,it1,iw0)*(1d0-xd)+set%data(ix1,iy1,iz0,it1,iw0)*xd
    f0110 = set%data(ix0,iy0,iz1,it1,iw0)*(1d0-xd)+set%data(ix1,iy0,iz1,it1,iw0)*xd
    f1110 = set%data(ix0,iy1,iz1,it1,iw0)*(1d0-xd)+set%data(ix1,iy1,iz1,it1,iw0)*xd
    f0001 = set%data(ix0,iy0,iz0,it0,iw1)*(1d0-xd)+set%data(ix1,iy0,iz0,it0,iw1)*xd
    f1001 = set%data(ix0,iy1,iz0,it0,iw1)*(1d0-xd)+set%data(ix0,iy1,iz0,it0,iw1)*xd
    f0101 = set%data(ix0,iy0,iz1,it0,iw1)*(1d0-xd)+set%data(ix1,iy0,iz1,it0,iw1)*xd
    f1101 = set%data(ix0,iy1,iz1,it0,iw1)*(1d0-xd)+set%data(ix1,iy1,iz1,it0,iw1)*xd
    f0011 = set%data(ix0,iy0,iz0,it1,iw1)*(1d0-xd)+set%data(ix1,iy0,iz0,it1,iw1)*xd
    f1011 = set%data(ix0,iy1,iz0,it1,iw1)*(1d0-xd)+set%data(ix1,iy1,iz0,it1,iw1)*xd
    f0111 = set%data(ix0,iy0,iz1,it1,iw1)*(1d0-xd)+set%data(ix1,iy0,iz1,it1,iw1)*xd
    f1111 = set%data(ix0,iy1,iz1,it1,iw1)*(1d0-xd)+set%data(ix1,iy1,iz1,it1,iw1)*xd

    f000 = f0000*(1d0-yd) + f1000*yd
    f100 = f0100*(1d0-yd) + f1100*yd
    f010 = f0010*(1d0-yd) + f1010*yd
    f110 = f0110*(1d0-yd) + f1110*yd
    f101 = f0101*(1d0-yd) + f1101*yd
    f011 = f0011*(1d0-yd) + f1011*yd
    f111 = f0111*(1d0-yd) + f1111*yd

    f00 = f000*(1d0-zd)+f100*zd
    f10 = f010*(1d0-zd)+f110*zd 
    f01 = f001*(1d0-zd)+f101*zd 
    f11 = f011*(1d0-zd)+f111*zd 

    f0 = f00 *(1d0-td)+f10*td
    f1 = f01 *(1d0-td)+f11*td

    res = f0*(1d0-wd)+f1*wd
    ret = .true.
  END FUNCTION hdcdl5d_dset_sc

  FUNCTION hdcdl1d_dset_ve(x,set,locator,res) RESULT(ret)
    !! Hard-coded for 1D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    TYPE(dset1d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd
    INTEGER       :: nv,j,ix0,ix1
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x)
      IF (ix0==0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      ix1 = ix0+1
      xd = (x(j)-set%x(ix0))/(set%x(ix1)-set%x(ix0))
      res(j) = set%data(ix0)*(1d0-xd)+set%data(ix1)*xd
    ENDDO
    ret = .true.
  END FUNCTION hdcdl1d_dset_ve

  FUNCTION hdcdl2d_dset_ve(x,y,set,locator,res) RESULT(ret)
    !! Hard-coded for 2D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    TYPE(dset2d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd
    REAL(kind=wp) :: f0,f1
    INTEGER      :: nv,j,ix0,iy0,ix1,iy1
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x)
      iy0 = locator(y(j),set%y)
      IF (ix0==0 .OR. iy0==0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      ix1 = ix0 + 1 ; iy1 = iy0 + 1
      xd  = (x(j)-set%x(ix0))/(set%x(ix1)-set%x(ix0))
      yd  = (y(j)-set%y(iy0))/(set%y(iy1)-set%y(iy0))
      f0  = set%data(ix0,iy0)*(1d0-xd)+set%data(ix1,iy0)*xd
      f1  = set%data(ix0,iy1)*(1d0-xd)+set%data(ix1,iy1)*xd
      res(j) = f0*(1d0-yd)+f1*yd
    ENDDO
    ret = .true.
  END FUNCTION hdcdl2d_dset_ve 

  FUNCTION hdcdl3d_dset_ve(x,y,z,set,locator,res) RESULT(ret)
    !! Hard-coded for 3D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    TYPE(dset3d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd
    REAL(kind=wp) :: f00,f10,f01,f11,f0,f1
    INTEGER       :: nv,j,ix0,iy0,iz0,ix1,iy1,iz1
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x)
      iy0 = locator(y(j),set%y)
      iz0 = locator(z(j),set%z)

      IF (ix0==0 .OR. iy0==0 .OR. iz0==0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF

      ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1
      xd = (x(j)-set%x(ix0))/(set%x(ix1)-set%x(ix0))
      yd = (y(j)-set%y(iy0))/(set%y(iy1)-set%y(iy0))
      zd = (z(j)-set%z(iz0))/(set%z(iz1)-set%z(iz0))

      f00 = set%data(ix0,iy0,iz0)*(1d0-xd)+set%data(ix1,iy0,iz0)*xd
      f10 = set%data(ix0,iy1,iz0)*(1d0-xd)+set%data(ix1,iy1,iz0)*xd
      f01 = set%data(ix0,iy0,iz1)*(1d0-xd)+set%data(ix1,iy0,iz1)*xd
      f11 = set%data(ix0,iy1,iz1)*(1d0-xd)+set%data(ix1,iy1,iz1)*xd

      f0 = f00 *(1d0-yd)+f10*yd
      f1 = f00 *(1d0-yd)+f10*yd

      res(j) = f0*(1d0-zd)+f1*zd
    ENDDO
    ret = .true.
  END FUNCTION hdcdl3d_dset_ve

  FUNCTION hdcdl4d_dset_ve(x,y,z,t,set,locator,res) RESULT(ret)
    !! Hard-coded for 4D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: t       !! T coordinate of the points to interpolate
    TYPE(dset4d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd,td
    REAL(kind=wp) :: f000,f100,f010,f110,f001,f101,f011,f111, &
                     f00,f10,f01,f11,f0,f1
    INTEGER       :: nv,j,ix0,iy0,iz0,it0,ix1,iy1,iz1,it1
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x)
      iy0 = locator(y(j),set%y)
      iz0 = locator(z(j),set%z)
      it0 = locator(t(j),set%t)
    
      IF (ix0==0 .OR. iy0==0 .OR. iz0==0 .OR. it0==0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF

      ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1 ; it1=it0+1
      xd = (x(j)-set%x(ix0))/(set%x(ix1)-set%x(ix0))
      yd = (y(j)-set%y(iy0))/(set%y(iy1)-set%y(iy0))
      zd = (z(j)-set%z(iz0))/(set%z(iz1)-set%z(iz0))
      td = (t(j)-set%t(it0))/(set%t(it1)-set%t(it0))

      f000 = set%data(ix0,iy0,iz0,it0)*(1d0-xd)+set%data(ix1,iy0,iz0,it0)*xd
      f100 = set%data(ix0,iy1,iz0,it0)*(1d0-xd)+set%data(ix1,iy1,iz0,it0)*xd
      f010 = set%data(ix0,iy0,iz1,it0)*(1d0-xd)+set%data(ix1,iy0,iz1,it0)*xd
      f110 = set%data(ix0,iy1,iz1,it0)*(1d0-xd)+set%data(ix1,iy1,iz1,it0)*xd
      f001 = set%data(ix0,iy0,iz0,it1)*(1d0-xd)+set%data(ix1,iy0,iz0,it1)*xd
      f101 = set%data(ix0,iy1,iz0,it1)*(1d0-xd)+set%data(ix1,iy1,iz0,it1)*xd
      f011 = set%data(ix0,iy0,iz1,it1)*(1d0-xd)+set%data(ix1,iy0,iz1,it1)*xd
      f111 = set%data(ix0,iy1,iz1,it1)*(1d0-xd)+set%data(ix1,iy1,iz1,it1)*xd

      f00 = f000*(1d0-yd)+f100*yd
      f10 = f010*(1d0-yd)+f110*yd 
      f01 = f001*(1d0-yd)+f101*yd 
      f11 = f011*(1d0-yd)+f111*yd 

      f0 = f00 *(1d0-zd)+f10*zd
      f1 = f01 *(1d0-zd)+f11*zd

      res(j) = f0*(1d0-td)+f1*td
    ENDDO
    ret = .true.
  END FUNCTION hdcdl4d_dset_ve

  FUNCTION hdcdl5d_dset_ve(x,y,z,t,w,set,locator,res) RESULT(ret)
    !! Hard-coded for 5D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: t       !! T coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: w       !! W coordinate of the points to interpolate
    TYPE(dset5d), INTENT(in)                              :: set     !! Dataset with the tabulated function
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd,td,wd
    REAL(kind=wp) :: f0000,f1000,f0100,f1100,f0010,f1010,f0110,f1110, & 
                     f0001,f1001,f0101,f1101,f0011,f1011,f0111,f1111, &
                     f000,f100,f010,f110,f001,f101,f011,f111,         &
                     f00,f10,f01,f11,f0,f1
    INTEGER       :: nv,j,ix0,iy0,iz0,it0,iw0,ix1,iy1,iz1,it1,iw1
    ret = .false.
    nv = SIZE(x) ; ALLOCATE(res(nv))
    DO j=1,nv
      ix0 = locator(x(j),set%x)
      iy0 = locator(y(j),set%y)
      iz0 = locator(z(j),set%z)
      it0 = locator(t(j),set%t)
      iw0 = locator(w(j),set%w)

      IF (ix0==0 .OR. iy0==0 .OR. iz0==0 .OR. it0==0 .OR. iw0==0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
    
      ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1 
      xd = (x(j)-set%x(ix0))/(set%x(ix1)-set%x(ix0))
      yd = (y(j)-set%y(iy0))/(set%y(iy1)-set%y(iy0))
      zd = (z(j)-set%z(iz0))/(set%z(iz1)-set%z(iz0))
      td = (t(j)-set%t(it0))/(set%t(it1)-set%t(it0))
      wd = (w(j)-set%w(it0))/(set%w(it1)-set%w(it0))

      f0000 = set%data(ix0,iy0,iz0,it0,iw0)*(1d0-xd)+set%data(ix1,iy0,iz0,it0,iw0)*xd
      f1000 = set%data(ix0,iy1,iz0,it0,iw0)*(1d0-xd)+set%data(ix0,iy1,iz0,it0,iw0)*xd
      f0100 = set%data(ix0,iy0,iz1,it0,iw0)*(1d0-xd)+set%data(ix1,iy0,iz1,it0,iw0)*xd
      f1100 = set%data(ix0,iy1,iz1,it0,iw0)*(1d0-xd)+set%data(ix1,iy1,iz1,it0,iw0)*xd
      f0010 = set%data(ix0,iy0,iz0,it1,iw0)*(1d0-xd)+set%data(ix1,iy0,iz0,it1,iw0)*xd
      f1010 = set%data(ix0,iy1,iz0,it1,iw0)*(1d0-xd)+set%data(ix1,iy1,iz0,it1,iw0)*xd
      f0110 = set%data(ix0,iy0,iz1,it1,iw0)*(1d0-xd)+set%data(ix1,iy0,iz1,it1,iw0)*xd
      f1110 = set%data(ix0,iy1,iz1,it1,iw0)*(1d0-xd)+set%data(ix1,iy1,iz1,it1,iw0)*xd
      f0001 = set%data(ix0,iy0,iz0,it0,iw1)*(1d0-xd)+set%data(ix1,iy0,iz0,it0,iw1)*xd
      f1001 = set%data(ix0,iy1,iz0,it0,iw1)*(1d0-xd)+set%data(ix0,iy1,iz0,it0,iw1)*xd
      f0101 = set%data(ix0,iy0,iz1,it0,iw1)*(1d0-xd)+set%data(ix1,iy0,iz1,it0,iw1)*xd
      f1101 = set%data(ix0,iy1,iz1,it0,iw1)*(1d0-xd)+set%data(ix1,iy1,iz1,it0,iw1)*xd
      f0011 = set%data(ix0,iy0,iz0,it1,iw1)*(1d0-xd)+set%data(ix1,iy0,iz0,it1,iw1)*xd
      f1011 = set%data(ix0,iy1,iz0,it1,iw1)*(1d0-xd)+set%data(ix1,iy1,iz0,it1,iw1)*xd
      f0111 = set%data(ix0,iy0,iz1,it1,iw1)*(1d0-xd)+set%data(ix1,iy0,iz1,it1,iw1)*xd
      f1111 = set%data(ix0,iy1,iz1,it1,iw1)*(1d0-xd)+set%data(ix1,iy1,iz1,it1,iw1)*xd

      f000 = f0000*(1d0-yd) + f1000*yd
      f100 = f0100*(1d0-yd) + f1100*yd
      f010 = f0010*(1d0-yd) + f1010*yd
      f110 = f0110*(1d0-yd) + f1110*yd
      f101 = f0101*(1d0-yd) + f1101*yd
      f011 = f0011*(1d0-yd) + f1011*yd
      f111 = f0111*(1d0-yd) + f1111*yd

      f00 = f000*(1d0-zd)+f100*zd
      f10 = f010*(1d0-zd)+f110*zd 
      f01 = f001*(1d0-zd)+f101*zd 
      f11 = f011*(1d0-zd)+f111*zd 

      f0 = f00 *(1d0-td)+f10*td
      f1 = f01 *(1d0-td)+f11*td

      res(j) = f0*(1d0-wd)+f1*wd
    ENDDO
    ret = .true.
  END FUNCTION hdcdl5d_dset_ve

END MODULE LINTDSET

