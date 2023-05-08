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

!! file: lintgen.f90
!! summary: Linear interpolation generic interfaces
!! author: burgalat
!! date: 2010-2014 

MODULE LINTGEN
  !! Generic linear interpolation module
  USE LINT_PREC
  USE LINTCORE
  IMPLICIT NONE

  PUBLIC  :: lint_gen, hdcd_lint_gen
  PRIVATE 

  !> Interface to multi-variate linear interpolation
  !! 
  !! For __n__, the number of function's variables, user must provide: 
  !!
  !! - __n__ scalar(s) with the coordinate(s) of the point to interpolate.
  !! - __n__ vector(s) with the tabulated values of each coordinate. 
  !! - a single __n__-D array with the value of tabulated function at each
  !!   tabulated coordinate.
  !! - A function, satisfying [[lintcore(module):locate(interface)]] interface, 
  !!   that locates the neighbourhood of the point to interpolate.
  !! - A scalar that stores the output value.
  INTERFACE lint_gen
    MODULE PROCEDURE l1d_gen_sc,l2d_gen_sc,l3d_gen_sc,l4d_gen_sc,l5d_gen_sc
    MODULE PROCEDURE l1d_gen_ve,l2d_gen_ve,l3d_gen_ve,l4d_gen_ve,l5d_gen_ve
  END INTERFACE 

  !> Interface to multi-variate hard-coded linear interpolation
  !!
  !! Same remarks as in [[lintgen(module):lint_gen(interface)]] apply here.
  INTERFACE hdcd_lint_gen
    MODULE PROCEDURE hdcdl1d_gen_sc, hdcdl2d_gen_sc, hdcdl3d_gen_sc, &
                     hdcdl4d_gen_sc, hdcdl5d_gen_sc
    MODULE PROCEDURE hdcdl1d_gen_ve, hdcdl2d_gen_ve, hdcdl3d_gen_ve, &
                     hdcdl4d_gen_ve, hdcdl5d_gen_ve
  END INTERFACE

  CONTAINS

  ! GENERIC versions
  ! ----------------

  FUNCTION l1d_gen_sc(x,cx,cv, locator,res) RESULT(ret)
    !! Wrapper interface for 1D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)               :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:) :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)              :: res     !! Interpolated value
    PROCEDURE(locate)                       :: locator !! Locator function
    LOGICAL :: ret                                     !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(2,2) :: g1d 
    INTEGER                       :: i,ix
    ret = .false. 
    ix = locator(x,cx) ; IF (ix == 0) RETURN
    DO i=0,1
      g1d(i+1,1) = cx(ix+i)
      g1d(i+1,2) = cv(ix+i)
    ENDDO
    res = lintc_((/x/),g1d) 
    ret = .true.
  END FUNCTION l1d_gen_sc

  FUNCTION l2d_gen_sc(x,y,cx,cy,cv,locator,res) RESULT(ret)
    !! Wrapper interface for 2D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)                 :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                 :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)   :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)   :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)                :: res     !! Interpolated value
    PROCEDURE(locate)                         :: locator !! Locator function
    LOGICAL :: ret                                       !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(4,3) :: g2d 
    INTEGER                       :: ix0,iy0,i,a,b
    ret = .false.
    ix0 = locator(x,cx) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,cy) ; IF (iy0 == 0) RETURN
    DO i=1,4
      a=ix0+MOD((i-1),2)   ; g2d(i,1) = cx(a)
      b=iy0+MOD((i-1)/2,2) ; g2d(i,2) = cy(b) 
      g2d(i,3) = cv(a,b)
    ENDDO
    res = lintc_((/x,y/),g2d) 
    ret = .true.
  END FUNCTION l2d_gen_sc

  FUNCTION l3d_gen_sc(x,y,z,cx,cy,cz,cv,locator,res) RESULT(ret)
    !! Wrapper interface for 3D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)                   :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                   :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                   :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)     :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)     :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)     :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)                  :: res     !! Interpolated value
    PROCEDURE(locate)                           :: locator !! Locator function
    LOGICAL :: ret                                         !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(8,4) :: g3d 
    INTEGER                      :: ix0,iy0,iz0,i,a,b,c
    ret = .false.
    ix0 = locator(x,cx) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,cy) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,cz) ; IF (iz0 == 0) RETURN
    DO i=1,8
      a=ix0+MOD((i-1),2)   ; g3d(i,1) = cx(a)
      b=iy0+MOD((i-1)/2,2) ; g3d(i,2) = cy(b) 
      c=iz0+MOD((i-1)/4,2) ; g3d(i,3) = cz(c)
      g3d(i,4) = cv(a,b,c)
    ENDDO
    res = lintc_((/x,y,z/),g3d) 
    ret = .true.
  END FUNCTION l3d_gen_sc

  FUNCTION l4d_gen_sc(x,y,z,t,cx,cy,cz,ct,cv,locator,res) RESULT(ret)
    !! Wrapper interface for 4D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)                     :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                     :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                     :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                     :: t       !! T coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)       :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)       :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)       :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)       :: ct      !! Tabulated function T coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:,:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)                    :: res     !! Interpolated value
    PROCEDURE(locate)                             :: locator !! Locator function
    LOGICAL :: ret                                           !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(16,5) :: g4d 
    INTEGER                        :: ix0,iy0,iz0,it0
    INTEGER                        :: i,a,b,c,d
    ret = .false.
    ix0 = locator(x,cx) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,cy) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,cz) ; IF (iz0 == 0) RETURN
    it0 = locator(t,ct) ; IF (it0 == 0) RETURN
    DO i=1,16
      a=ix0+MOD((i-1),2)   ; g4d(i,1) = cx(a)
      b=iy0+MOD((i-1)/2,2) ; g4d(i,2) = cy(b) 
      c=iz0+MOD((i-1)/4,2) ; g4d(i,3) = cz(c)
      d=it0+MOD((i-1)/8,2) ; g4d(i,4) = ct(d)
      g4d(i,5) = cv(a,b,c,d)
    ENDDO
    res = lintc_((/x,y,z,t/),g4d) 
    ret = .true.
  END FUNCTION l4d_gen_sc

  FUNCTION l5d_gen_sc(x,y,z,t,w,cx,cy,cz,ct,cw,cv,locator,res) RESULT(ret)
    !! Wrapper interface for 5D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)                       :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                       :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                       :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                       :: t       !! T coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                       :: w       !! W coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: ct      !! Tabulated function T coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: cw      !! Tabulated function W coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:,:,:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)                      :: res     !! Interpolated value
    PROCEDURE(locate)                               :: locator !! Locator function
    LOGICAL :: ret                                             !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(32,6) :: g5d 
    INTEGER                        :: ix0,iy0,iz0,it0,iw0
    INTEGER                        :: i,a,b,c,d,e
    ret = .false.
    ix0 = locator(x,cx) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,cy) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,cz) ; IF (iz0 == 0) RETURN
    it0 = locator(t,ct) ; IF (it0 == 0) RETURN
    iw0 = locator(w,cw) ; IF (iw0 == 0) RETURN
    DO i=1,32
      a=ix0+MOD((i-1),2)    ; g5d(i,1) = cx(a)
      b=iy0+MOD((i-1)/2,2)  ; g5d(i,2) = cy(b) 
      c=iz0+MOD((i-1)/4,2)  ; g5d(i,3) = cz(c)
      d=it0+MOD((i-1)/8,2)  ; g5d(i,4) = ct(d)
      e=iw0+MOD((i-1)/16,2) ; g5d(i,5) = cw(e)
      g5d(i,6) = cv(a,b,c,d,e)
    ENDDO
    res = lintc_((/x,y,z,t,w/),g5d) 
    ret = .true.
  END FUNCTION l5d_gen_sc

  ! HARD-CODED versions
  ! -------------------

  FUNCTION hdcdl1d_gen_sc(x,cx,cv,locator,res) RESULT(ret)
    !! Hard-coded 1D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)               :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:) :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)              :: res     !! Interpolated value
    PROCEDURE(locate)                       :: locator !! Locator function
    LOGICAL :: ret                                     !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd
    INTEGER       :: ix0,ix1
    ret = .false. 
    ix0 = locator(x,cx) ; IF (ix0 == 0) RETURN
    ix1 = ix0+1
    xd = (x-cx(ix0))/(cx(ix1)-cx(ix0))
    res = cv(ix0)*(1d0-xd)+cv(ix1)*xd
    ret = .true.
  END FUNCTION hdcdl1d_gen_sc

  FUNCTION hdcdl2d_gen_sc(x,y,cx,cy,cv,locator,res) RESULT(ret)
    !! Hard-coded 2D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)                 :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                 :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)   :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)   :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)                :: res     !! Interpolated value
    PROCEDURE(locate)                         :: locator !! Locator function
    LOGICAL :: ret                                       !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd
    REAL(kind=wp) :: f0,f1
    INTEGER       :: ix0,iy0,ix1,iy1
    ret = .false.
    ix0 = locator(x,cx) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,cy) ; IF (iy0 == 0) RETURN
    ix1 = ix0 + 1 ; iy1 = iy0 + 1
    xd  = (x-cx(ix0))/(cx(ix1)-cx(ix0))
    yd  = (y-cy(iy0))/(cy(iy1)-cy(iy0))
    f0  = cv(ix0,iy0)*(1d0-xd)+cv(ix1,iy0)*xd
    f1  = cv(ix0,iy1)*(1d0-xd)+cv(ix1,iy1)*xd
    res = f0*(1d0-yd)+f1*yd
    ret = .true.
  END FUNCTION hdcdl2d_gen_sc 

  FUNCTION hdcdl3d_gen_sc(x,y,z,cx,cy,cz,cv,locator,res) RESULT(ret)
    !! Hard-coded 3D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)                   :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                   :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                   :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)     :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)     :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)     :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)                  :: res     !! Interpolated value
    PROCEDURE(locate)                           :: locator !! Locator function
    LOGICAL :: ret                                         !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd
    REAL(kind=wp) :: f00,f10,f01,f11,f0,f1
    INTEGER       :: ix0,iy0,iz0,ix1,iy1,iz1
    ret = .false.
    ix0 = locator(x,cx) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,cy) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,cz) ; IF (iz0 == 0) RETURN
    ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1
    xd = (x-cx(ix0))/(cx(ix1)-cx(ix0))
    yd = (y-cy(iy0))/(cy(iy1)-cy(iy0))
    zd = (z-cz(iz0))/(cz(iz1)-cz(iz0))

    f00 = cv(ix0,iy0,iz0)*(1d0-xd)+cv(ix1,iy0,iz0)*xd
    f10 = cv(ix0,iy1,iz0)*(1d0-xd)+cv(ix1,iy1,iz0)*xd
    f01 = cv(ix0,iy0,iz1)*(1d0-xd)+cv(ix1,iy0,iz1)*xd
    f11 = cv(ix0,iy1,iz1)*(1d0-xd)+cv(ix1,iy1,iz1)*xd

    f0 = f00 *(1d0-yd)+f10*yd
    f1 = f00 *(1d0-yd)+f10*yd

    res = f0*(1d0-zd)+f1*zd
    ret = .true.
  END FUNCTION hdcdl3d_gen_sc

  FUNCTION hdcdl4d_gen_sc(x,y,z,t,cx,cy,cz,ct,cv,locator,res) RESULT(ret)
    !! Hard-coded 4D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)                     :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                     :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                     :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                     :: t       !! T coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)       :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)       :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)       :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)       :: ct      !! Tabulated function T coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:,:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)                    :: res     !! Interpolated value
    PROCEDURE(locate)                             :: locator !! Locator function
    LOGICAL :: ret                                           !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd,td
    REAL(kind=wp) :: f000,f100,f010,f110,f001,f101,f011,f111, &
                     f00,f10,f01,f11,f0,f1
    INTEGER       :: ix0,iy0,iz0,it0,ix1,iy1,iz1,it1
    ret = .false.
    ix0 = locator(x,cx) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,cy) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,cz) ; IF (iz0 == 0) RETURN
    it0 = locator(t,ct) ; IF (it0 == 0) RETURN
    
    ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1 ; it1=it0+1
    xd = (x-cx(ix0))/(cx(ix1)-cx(ix0))
    yd = (y-cy(iy0))/(cy(iy1)-cy(iy0))
    zd = (z-cz(iz0))/(cz(iz1)-cz(iz0))
    td = (t-ct(it0))/(ct(it1)-ct(it0))

    f000 = cv(ix0,iy0,iz0,it0)*(1d0-xd)+cv(ix1,iy0,iz0,it0)*xd
    f100 = cv(ix0,iy1,iz0,it0)*(1d0-xd)+cv(ix1,iy1,iz0,it0)*xd
    f010 = cv(ix0,iy0,iz1,it0)*(1d0-xd)+cv(ix1,iy0,iz1,it0)*xd
    f110 = cv(ix0,iy1,iz1,it0)*(1d0-xd)+cv(ix1,iy1,iz1,it0)*xd
    f001 = cv(ix0,iy0,iz0,it1)*(1d0-xd)+cv(ix1,iy0,iz0,it1)*xd
    f101 = cv(ix0,iy1,iz0,it1)*(1d0-xd)+cv(ix1,iy1,iz0,it1)*xd
    f011 = cv(ix0,iy0,iz1,it1)*(1d0-xd)+cv(ix1,iy0,iz1,it1)*xd
    f111 = cv(ix0,iy1,iz1,it1)*(1d0-xd)+cv(ix1,iy1,iz1,it1)*xd

    f00 = f000*(1d0-yd)+f100*yd
    f10 = f010*(1d0-yd)+f110*yd 
    f01 = f001*(1d0-yd)+f101*yd 
    f11 = f011*(1d0-yd)+f111*yd 

    f0 = f00 *(1d0-zd)+f10*zd
    f1 = f01 *(1d0-zd)+f11*zd

    res = f0*(1d0-td)+f1*td
    ret = .true.
  END FUNCTION hdcdl4d_gen_sc

  FUNCTION hdcdl5d_gen_sc(x,y,z,t,w,cx,cy,cz,ct,cw,cv,locator,res) RESULT(ret)
    !! Hard-coded 5D linear interpolation (scalar)
    !!
    !! @warning
    !! On error, __res__ output value is undefined.
    REAL(kind=wp), INTENT(in)                       :: x       !! X coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                       :: y       !! Y coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                       :: z       !! Z coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                       :: t       !! T coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in)                       :: w       !! W coordinate of the point to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: ct      !! Tabulated function T coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)         :: cw      !! Tabulated function W coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:,:,:) :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out)                      :: res     !! Interpolated value
    PROCEDURE(locate)                               :: locator !! Locator function
    LOGICAL :: ret                                             !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd,td,wd
    REAL(kind=wp) :: f0000,f1000,f0100,f1100,f0010,f1010,f0110,f1110, & 
                     f0001,f1001,f0101,f1101,f0011,f1011,f0111,f1111, &
                     f000,f100,f010,f110,f001,f101,f011,f111,         &
                     f00,f10,f01,f11,f0,f1
    INTEGER       :: ix0,iy0,iz0,it0,iw0,ix1,iy1,iz1,it1,iw1
    ret = .false.
    ix0 = locator(x,cx) ; IF (ix0 == 0) RETURN
    iy0 = locator(y,cy) ; IF (iy0 == 0) RETURN
    iz0 = locator(z,cz) ; IF (iz0 == 0) RETURN
    it0 = locator(t,ct) ; IF (it0 == 0) RETURN
    iw0 = locator(w,cw) ; IF (iw0 == 0) RETURN
    
    ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1 ; it1=it0+1 ; iw1=iw0+1
    xd = (x-cx(ix0))/(cx(ix1)-cx(ix0))
    yd = (y-cy(iy0))/(cy(iy1)-cy(iy0))
    zd = (z-cz(iz0))/(cz(iz1)-cz(iz0))
    td = (t-ct(it0))/(ct(it1)-ct(it0))
    wd = (w-cw(iw0))/(cw(iw1)-cw(iw0))

    f0000 = cv(ix0,iy0,iz0,it0,iw0)*(1d0-xd)+cv(ix1,iy0,iz0,it0,iw0)*xd
    f1000 = cv(ix0,iy1,iz0,it0,iw0)*(1d0-xd)+cv(ix0,iy1,iz0,it0,iw0)*xd
    f0100 = cv(ix0,iy0,iz1,it0,iw0)*(1d0-xd)+cv(ix1,iy0,iz1,it0,iw0)*xd
    f1100 = cv(ix0,iy1,iz1,it0,iw0)*(1d0-xd)+cv(ix1,iy1,iz1,it0,iw0)*xd
    f0010 = cv(ix0,iy0,iz0,it1,iw0)*(1d0-xd)+cv(ix1,iy0,iz0,it1,iw0)*xd
    f1010 = cv(ix0,iy1,iz0,it1,iw0)*(1d0-xd)+cv(ix1,iy1,iz0,it1,iw0)*xd
    f0110 = cv(ix0,iy0,iz1,it1,iw0)*(1d0-xd)+cv(ix1,iy0,iz1,it1,iw0)*xd
    f1110 = cv(ix0,iy1,iz1,it1,iw0)*(1d0-xd)+cv(ix1,iy1,iz1,it1,iw0)*xd
    f0001 = cv(ix0,iy0,iz0,it0,iw1)*(1d0-xd)+cv(ix1,iy0,iz0,it0,iw1)*xd
    f1001 = cv(ix0,iy1,iz0,it0,iw1)*(1d0-xd)+cv(ix0,iy1,iz0,it0,iw1)*xd
    f0101 = cv(ix0,iy0,iz1,it0,iw1)*(1d0-xd)+cv(ix1,iy0,iz1,it0,iw1)*xd
    f1101 = cv(ix0,iy1,iz1,it0,iw1)*(1d0-xd)+cv(ix1,iy1,iz1,it0,iw1)*xd
    f0011 = cv(ix0,iy0,iz0,it1,iw1)*(1d0-xd)+cv(ix1,iy0,iz0,it1,iw1)*xd
    f1011 = cv(ix0,iy1,iz0,it1,iw1)*(1d0-xd)+cv(ix1,iy1,iz0,it1,iw1)*xd
    f0111 = cv(ix0,iy0,iz1,it1,iw1)*(1d0-xd)+cv(ix1,iy0,iz1,it1,iw1)*xd
    f1111 = cv(ix0,iy1,iz1,it1,iw1)*(1d0-xd)+cv(ix1,iy1,iz1,it1,iw1)*xd

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
  END FUNCTION hdcdl5d_gen_sc

  !--------------------
  ! VECTORIZED VERSIONS
  !--------------------

  ! GENERIC versions
  ! ----------------

  FUNCTION l1d_gen_ve(x,cx,cv,locator,res) RESULT(ret)
    !! Wrapper interface for 1D linear interpolation (vector) 
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(2,2) :: g1d 
    INTEGER                       :: n,j,i,ix
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix = locator(x(j),cx) 
      IF (ix == 0) THEN
        DEALLOCATE(res) ; RETURN 
      ENDIF
      DO i=0,1
        g1d(i+1,1) = cx(ix+i)
        g1d(i+1,2) = cv(ix+i)
      ENDDO
      res(j) = lintc_((/x(j)/),g1d) 
    ENDDO
    ret = .true.
  END FUNCTION l1d_gen_ve

  FUNCTION l2d_gen_ve(x,y,cx,cy,cv,locator,res) RESULT(ret)
    !! Wrapper interface for 2D linear interpolation (vector) 
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:)             :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(4,3) :: g2d 
    INTEGER                       :: n,j,ix0,iy0,i,a,b
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix0 = locator(x(j),cx)
      iy0 = locator(y(j),cy) 
      IF (ix0 == 0 .OR. iy0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      DO i=1,4
        a=ix0+MOD((i-1),2)   ; g2d(i,1) = cx(a)
        b=iy0+MOD((i-1)/2,2) ; g2d(i,2) = cy(b) 
        g2d(i,3) = cv(a,b)
      ENDDO
      res(j) = lintc_((/x(j),y(j)/),g2d) 
   ENDDO
    ret = .true.
  END FUNCTION l2d_gen_ve

  FUNCTION l3d_gen_ve(x,y,z,cx,cy,cz,cv,locator,res) RESULT(ret)
    !! Wrapper interface for 3D linear interpolation (vector) 
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:)           :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(8,4) :: g3d 
    INTEGER                       :: ix0,iy0,iz0
    INTEGER                       :: n,j,i,a,b,c
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix0 = locator(x(j),cx)
      iy0 = locator(y(j),cy)
      iz0 = locator(z(j),cz)
      IF (ix0 == 0 .OR. iy0 == 0 .OR. iz0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      DO i=1,8
        a=ix0+MOD((i-1),2)   ; g3d(i,1) = cx(a)
        b=iy0+MOD((i-1)/2,2) ; g3d(i,2) = cy(b) 
        c=iz0+MOD((i-1)/4,2) ; g3d(i,3) = cz(c)
        g3d(i,4) = cv(a,b,c)
      ENDDO
      res(j) = lintc_((/x(j),y(j),z(j)/),g3d) 
    ENDDO
    ret = .true.
  END FUNCTION l3d_gen_ve

  FUNCTION l4d_gen_ve(x,y,z,t,cx,cy,cz,ct,cv,locator,res) RESULT(ret)
    !! Wrapper interface for 4D linear interpolation (vector) 
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: t       !! T coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: ct      !! Tabulated function T coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:,:)         :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(16,5) :: g4d 
    INTEGER                        :: ix0,iy0,iz0,it0
    INTEGER                        :: n,j,i,a,b,c,d
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix0 = locator(x(j),cx)
      iy0 = locator(y(j),cy)
      iz0 = locator(z(j),cz)
      it0 = locator(t(j),ct)
      IF (ix0 == 0 .OR. iy0 == 0 .OR. iz0 == 0 .OR. it0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      DO i=1,16
        a=ix0+MOD((i-1),2)   ; g4d(i,1) = cx(a)
        b=iy0+MOD((i-1)/2,2) ; g4d(i,2) = cy(b) 
        c=iz0+MOD((i-1)/4,2) ; g4d(i,3) = cz(c)
        d=it0+MOD((i-1)/8,2) ; g4d(i,4) = ct(d)
        g4d(i,5) = cv(a,b,c,d)
      ENDDO
      res(j) = lintc_((/x(j),y(j),z(j),t(j)/),g4d) 
    ENDDO
    ret = .true.
  END FUNCTION l4d_gen_ve

  FUNCTION l5d_gen_ve(x,y,z,t,w,cx,cy,cz,ct,cw,cv,locator,res) RESULT(ret)
    !! Wrapper interface for 5D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: t       !! T coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: w       !! W coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: ct      !! Tabulated function T coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cw      !! Tabulated function W coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:,:,:)       :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp), DIMENSION(32,6) :: g5d 
    INTEGER                        :: ix0,iy0,iz0,it0,iw0
    INTEGER                        :: n,j,i,a,b,c,d,e
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix0 = locator(x(j),cx)
      iy0 = locator(y(j),cy)
      iz0 = locator(z(j),cz)
      it0 = locator(t(j),ct)
      iw0 = locator(w(j),cw)
      IF (ix0 == 0 .OR. iy0 == 0 .OR. iz0 == 0 .OR. it0 == 0 .OR. iw0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      DO i=1,32
        a=ix0+MOD((i-1),2)    ; g5d(i,1) = cx(a)
        b=iy0+MOD((i-1)/2,2)  ; g5d(i,2) = cy(b) 
        c=iz0+MOD((i-1)/4,2)  ; g5d(i,3) = cz(c)
        d=it0+MOD((i-1)/8,2)  ; g5d(i,4) = ct(d)
        e=iw0+MOD((i-1)/16,2) ; g5d(i,5) = cw(e)
        g5d(i,6) = cv(a,b,c,d,e)
      ENDDO
      res(j) = lintc_((/x,y,z,t,w/),g5d) 
    ENDDO
    ret = .true.
  END FUNCTION l5d_gen_ve

  ! HARD-CODED versions
  ! -------------------

  FUNCTION hdcdl1d_gen_ve(x,cx,cv,locator,res) RESULT(ret)
    !! Hard-coded 1D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd
    INTEGER       :: ix0,ix1
    INTEGER       :: n,j
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix0 = locator(x(j),cx)
      IF (ix0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      ix1 = ix0+1
      xd = (x(j)-cx(ix0))/(cx(ix1)-cx(ix0))
      res(j) = cv(ix0)*(1d0-xd)+cv(ix1)*xd
    ENDDO
    ret = .true.
  END FUNCTION hdcdl1d_gen_ve

  FUNCTION hdcdl2d_gen_ve(x,y,cx,cy,cv,locator,res) RESULT(ret)
    !! Hard-coded 2D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:)             :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd
    REAL(kind=wp) :: f0,f1
    INTEGER       :: ix0,iy0,ix1,iy1
    INTEGER       :: n,j
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix0 = locator(x(j),cx)
      iy0 = locator(y(j),cy)
      IF (ix0 == 0 .OR. iy0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      ix1 = ix0 + 1 ; iy1 = iy0 + 1
      xd  = (x(j)-cx(ix0))/(cx(ix1)-cx(ix0))
      yd  = (y(j)-cy(iy0))/(cy(iy1)-cy(iy0))
      f0  = cv(ix0,iy0)*(1d0-xd)+cv(ix1,iy0)*xd
      f1  = cv(ix0,iy1)*(1d0-xd)+cv(ix1,iy1)*xd
      res(j) = f0*(1d0-yd)+f1*yd
    ENDDO
    ret = .true.
  END FUNCTION hdcdl2d_gen_ve 

  FUNCTION hdcdl3d_gen_ve(x,y,z,cx,cy,cz,cv,locator,res) RESULT(ret)
    !! Hard-coded 3D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:)           :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd
    REAL(kind=wp) :: f00,f10,f01,f11,f0,f1
    INTEGER       :: ix0,iy0,iz0,ix1,iy1,iz1
    INTEGER       :: n,j
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix0 = locator(x(j),cx)
      iy0 = locator(y(j),cy)
      iz0 = locator(z(j),cz)
      IF (ix0 == 0 .OR. iy0 == 0 .OR. iz0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1
      xd = (x(j)-cx(ix0))/(cx(ix1)-cx(ix0))
      yd = (y(j)-cy(iy0))/(cy(iy1)-cy(iy0))
      zd = (z(j)-cz(iz0))/(cz(iz1)-cz(iz0))

      f00 = cv(ix0,iy0,iz0)*(1d0-xd)+cv(ix1,iy0,iz0)*xd
      f10 = cv(ix0,iy1,iz0)*(1d0-xd)+cv(ix1,iy1,iz0)*xd
      f01 = cv(ix0,iy0,iz1)*(1d0-xd)+cv(ix1,iy0,iz1)*xd
      f11 = cv(ix0,iy1,iz1)*(1d0-xd)+cv(ix1,iy1,iz1)*xd

      f0 = f00 *(1d0-yd)+f10*yd
      f1 = f00 *(1d0-yd)+f10*yd

      res(j) = f0*(1d0-zd)+f1*zd
    ENDDO
    ret = .true.
  END FUNCTION hdcdl3d_gen_ve

  FUNCTION hdcdl4d_gen_ve(x,y,z,t,cx,cy,cz,ct,cv,locator,res) RESULT(ret)
    !! Hard-coded 4D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: t       !! T coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: ct      !! Tabulated function T coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:,:)         :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd,td
    REAL(kind=wp) :: f000,f100,f010,f110,f001,f101,f011,f111, &
                     f00,f10,f01,f11,f0,f1
    INTEGER       :: ix0,iy0,iz0,it0,ix1,iy1,iz1,it1
    INTEGER       :: n,j
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix0 = locator(x(j),cx)
      iy0 = locator(y(j),cy)
      iz0 = locator(z(j),cz)
      it0 = locator(t(j),ct)
      IF (ix0 == 0 .OR. iy0 == 0 .OR. iz0 == 0 .OR. it0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1 ; it1=it0+1
      xd = (x(j)-cx(ix0))/(cx(ix1)-cx(ix0))
      yd = (y(j)-cy(iy0))/(cy(iy1)-cy(iy0))
      zd = (z(j)-cz(iz0))/(cz(iz1)-cz(iz0))
      td = (t(j)-ct(it0))/(ct(it1)-ct(it0))

      f000 = cv(ix0,iy0,iz0,it0)*(1d0-xd)+cv(ix1,iy0,iz0,it0)*xd
      f100 = cv(ix0,iy1,iz0,it0)*(1d0-xd)+cv(ix1,iy1,iz0,it0)*xd
      f010 = cv(ix0,iy0,iz1,it0)*(1d0-xd)+cv(ix1,iy0,iz1,it0)*xd
      f110 = cv(ix0,iy1,iz1,it0)*(1d0-xd)+cv(ix1,iy1,iz1,it0)*xd
      f001 = cv(ix0,iy0,iz0,it1)*(1d0-xd)+cv(ix1,iy0,iz0,it1)*xd
      f101 = cv(ix0,iy1,iz0,it1)*(1d0-xd)+cv(ix1,iy1,iz0,it1)*xd
      f011 = cv(ix0,iy0,iz1,it1)*(1d0-xd)+cv(ix1,iy0,iz1,it1)*xd
      f111 = cv(ix0,iy1,iz1,it1)*(1d0-xd)+cv(ix1,iy1,iz1,it1)*xd

      f00 = f000*(1d0-yd)+f100*yd
      f10 = f010*(1d0-yd)+f110*yd 
      f01 = f001*(1d0-yd)+f101*yd 
      f11 = f011*(1d0-yd)+f111*yd 

      f0 = f00 *(1d0-zd)+f10*zd
      f1 = f01 *(1d0-zd)+f11*zd
  
      res(j) = f0*(1d0-td)+f1*td
    ENDDO
    ret = .true.
  END FUNCTION hdcdl4d_gen_ve

  FUNCTION hdcdl5d_gen_ve(x,y,z,t,w,cx,cy,cz,ct,cw,cv,locator,res) RESULT(ret)
    !! Hard-coded 5D linear interpolation (vector)
    !!
    !! @warning
    !! On error, __res__ output vector is undefined (i.e. not allocated).
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: x       !! X coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: y       !! Y coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: z       !! Z coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: t       !! T coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: w       !! W coordinate of the points to interpolate
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cx      !! Tabulated function X coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cy      !! Tabulated function Y coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cz      !! Tabulated function Z coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: ct      !! Tabulated function T coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:)               :: cw      !! Tabulated function W coordinates
    REAL(kind=wp), INTENT(in), DIMENSION(:,:,:,:,:)       :: cv      !! Tabulated function values
    REAL(kind=wp), INTENT(out), DIMENSION(:), ALLOCATABLE :: res     !! Interpolated values
    PROCEDURE(locate)                                     :: locator !! Locator function
    LOGICAL :: ret                                                   !! .true. on success, .false. otherwise
    REAL(kind=wp) :: xd,yd,zd,td,wd
    REAL(kind=wp) :: f0000,f1000,f0100,f1100,f0010,f1010,f0110,f1110, & 
                     f0001,f1001,f0101,f1101,f0011,f1011,f0111,f1111, &
                     f000,f100,f010,f110,f001,f101,f011,f111,         &
                     f00,f10,f01,f11,f0,f1
    INTEGER       :: ix0,iy0,iz0,it0,iw0,ix1,iy1,iz1,it1,iw1
    INTEGER       :: n,j
    ret = .false.
    n = SIZE(x) ; ALLOCATE(res(n))
    DO j=1,n 
      ix0 = locator(x(j),cx)
      iy0 = locator(y(j),cy)
      iz0 = locator(z(j),cz)
      it0 = locator(t(j),ct)
      iw0 = locator(w(j),cw)
      IF (ix0 == 0 .OR. iy0 == 0 .OR. iz0 == 0 .OR. it0 == 0 .OR. iw0 == 0) THEN
        DEALLOCATE(res) ; RETURN
      ENDIF
      ix1=ix0+1 ; iy1=iy0+1 ; iz1=iz0+1 ; 
      xd = (x(j)-cx(ix0))/(cx(ix1)-cx(ix0))
      yd = (y(j)-cy(iy0))/(cy(iy1)-cy(iy0))
      zd = (z(j)-cz(iz0))/(cz(iz1)-cz(iz0))
      td = (t(j)-ct(it0))/(ct(it1)-ct(it0))
      wd = (w(j)-cw(it0))/(cw(it1)-cw(it0))

      f0000 = cv(ix0,iy0,iz0,it0,iw0)*(1d0-xd)+cv(ix1,iy0,iz0,it0,iw0)*xd
      f1000 = cv(ix0,iy1,iz0,it0,iw0)*(1d0-xd)+cv(ix0,iy1,iz0,it0,iw0)*xd
      f0100 = cv(ix0,iy0,iz1,it0,iw0)*(1d0-xd)+cv(ix1,iy0,iz1,it0,iw0)*xd
      f1100 = cv(ix0,iy1,iz1,it0,iw0)*(1d0-xd)+cv(ix1,iy1,iz1,it0,iw0)*xd
      f0010 = cv(ix0,iy0,iz0,it1,iw0)*(1d0-xd)+cv(ix1,iy0,iz0,it1,iw0)*xd
      f1010 = cv(ix0,iy1,iz0,it1,iw0)*(1d0-xd)+cv(ix1,iy1,iz0,it1,iw0)*xd
      f0110 = cv(ix0,iy0,iz1,it1,iw0)*(1d0-xd)+cv(ix1,iy0,iz1,it1,iw0)*xd
      f1110 = cv(ix0,iy1,iz1,it1,iw0)*(1d0-xd)+cv(ix1,iy1,iz1,it1,iw0)*xd
      f0001 = cv(ix0,iy0,iz0,it0,iw1)*(1d0-xd)+cv(ix1,iy0,iz0,it0,iw1)*xd
      f1001 = cv(ix0,iy1,iz0,it0,iw1)*(1d0-xd)+cv(ix0,iy1,iz0,it0,iw1)*xd
      f0101 = cv(ix0,iy0,iz1,it0,iw1)*(1d0-xd)+cv(ix1,iy0,iz1,it0,iw1)*xd
      f1101 = cv(ix0,iy1,iz1,it0,iw1)*(1d0-xd)+cv(ix1,iy1,iz1,it0,iw1)*xd
      f0011 = cv(ix0,iy0,iz0,it1,iw1)*(1d0-xd)+cv(ix1,iy0,iz0,it1,iw1)*xd
      f1011 = cv(ix0,iy1,iz0,it1,iw1)*(1d0-xd)+cv(ix1,iy1,iz0,it1,iw1)*xd
      f0111 = cv(ix0,iy0,iz1,it1,iw1)*(1d0-xd)+cv(ix1,iy0,iz1,it1,iw1)*xd
      f1111 = cv(ix0,iy1,iz1,it1,iw1)*(1d0-xd)+cv(ix1,iy1,iz1,it1,iw1)*xd

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
  END FUNCTION hdcdl5d_gen_ve

END MODULE LINTGEN
