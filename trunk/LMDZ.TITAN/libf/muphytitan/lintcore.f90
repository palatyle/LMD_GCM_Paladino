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

!! file: lintcore.f90
!! summary: linear interpolation core function file
!! author: burgalat
!! date: 2010-2014

MODULE LINTCORE
  !! Core module of the library.
  !! 
  !! This module contains a single function that performs the linear 
  !! interpolation of a single _N_-D point between \(2^{N}\) adjacents 
  !! points.
  USE LINT_PREC
  IMPLICIT NONE

  PRIVATE :: wp ! from LINT_PREC

  INTERFACE
    FUNCTION locate(value,vector) RESULT(idx)
      !! Locate the nearest default value in vector
      !!
      !! the method should search the subscript of the nearest value by default in
      !! in the input vector.
      IMPORT wp
      REAL(kind=wp), INTENT(in)               :: value   !! value to search 
      REAL(kind=wp), INTENT(in), DIMENSION(:) :: vector  !! Vector to search in
      INTEGER :: idx                                     !! Subscript of the nearest value in vector
    END FUNCTION locate
  END INTERFACE


  CONTAINS 

  FUNCTION lintc_(point,grid) RESULT(res)
    !! Multivariate linear interpolation core function
    !! 
    !! The method computes multivariate linear interpolation at the given __point__ using its 
    !! neighbours given in __grid__.
    !!
    !! @warning
    !! In order to get a correct result, __grid__ must be ordered so first dimensions vary first 
    !! (see [Generic method](page/index.html/#generic-method) section of main documentation).
    !!
    !! @warning 
    !! The method in its current version does not check array boundaries. This operation should be 
    !! performed in wrappers of the function !
    INTEGER, PARAMETER :: np = 2 
    REAL(kind=wp), INTENT(in), DIMENSION(:)   :: point 
      !! Coordinates of the point to compute.
      !!
      !! For __N__-D interpolation, ut should be a vector of __N__ points.
    REAL(kind=wp), INTENT(in), DIMENSION(:,:) :: grid
      !! Grid of values used for interpolation. 
      !!
      !! For __N__-D interpolation, it should be a 2D-array of \(2^{N}\) rows and \(N+1\) columns. 
      !! Each row corresponds to a point with N coordinates, the last column is reserved for the
      !! value of the point.
    REAL(kind=wp) :: res
      !! Interpolated value
    REAL(kind=wp), DIMENSION(:),ALLOCATABLE :: val
    REAL(kind=wp)                           :: cd
    INTEGER                                 :: nv,mi,ngp,cp,i,j,k
    nv = SIZE(point) ; mi = np**nv
    ALLOCATE(val(2*mi-1)) 
    val(1:mi) = grid(:,nv+1) ; val(mi+1:2*mi-1) = 0._wp
    ! Computes the QnD linear interpolation
    cp = 0
    DO i=1,nv
      cd = (point(i)-grid(1,i))/(grid(mi,i)-grid(1,i))
      k = 1 ; ngp = np**(nv-i+1) ; cp = cp + ngp
      DO j=1,ngp,np 
        val(cp+k) = val(j+cp-ngp) * (1._wp - cd) + val(j+cp-ngp+1)*cd
        k = k + 1
      ENDDO
    ENDDO
    res = val(cp+k-1) 
    DEALLOCATE(val) ! useless normally
    RETURN
  END FUNCTION lintc_

END MODULE LINTCORE
