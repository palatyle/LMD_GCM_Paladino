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

!! file: locators.f90
!! summary: Locator functions definition file
!! author: burgalat 
!! date: 2010-2014


MODULE LOCATORS
  !! Locator functions definition module.
  !! 
  !! This module defines some locator functions which search value in ordered vectors (with no 
  !! duplicates). Two algorithms are provided :
  !!
  !! - The binary search algorithm.
  !! - A simple algorithm for direct access assuming searched vector is regularly
  !!   spaced.
  !!
  !! All these functions satisfy the interface required for the __locator__ argument of linear 
  !! interpolation functions.
  USE LINT_PREC
  IMPLICIT NONE

  PRIVATE :: wp ! from LINT_PREC

  CONTAINS

  FUNCTION locate(value,vector) RESULT(res)
    !! Basic binary search algorithm
    REAL(kind=wp), INTENT(in)               :: value  !! Value to search
    REAL(kind=wp), INTENT(in), DIMENSION(:) :: vector !! Input vector to search in
    INTEGER :: res
      !! Lowest subscript of the nearest value or __0__ if value is out of range.
    REAL(kind=wp) :: l,u
    INTEGER       :: jl,jm,ju, nv
    res = 0 ; nv = SIZE(vector) ; l = vector(1) ; u = vector(nv)
    ! Check for out of range value
    IF ((value>l.AND.value>u).OR.(value<l.AND.value<u)) RETURN
    ! Search in the array
    jl=0 ; ju=nv+1
    DO WHILE (ju-jl > 1)
      res=(ju+jl)/2 
      IF (res == 0) RETURN   ! should never happen
      IF((u>=l).EQV.(value >= vector(res))) THEN
        jl=res
      ELSE
        ju=res
      ENDIF
    ENDDO
    res = jl
    RETURN
  END FUNCTION locate

  FUNCTION locate_ext(value,vector) RESULT(res)
    !! Basic binary search algorithm with extrapolation
    !! 
    !! The function performs the same computation than [[locators(module):locate(function)]] except
    !! that if __value__ is out of range, __1__ or __SIZE(vector)-1__ is returned with respect 
    !! to the nearest __vector__'s extremum.
    REAL(kind=wp), INTENT(in)               :: value  !! Value to search
    REAL(kind=wp), INTENT(in), DIMENSION(:) :: vector !! Input vector to search in
    INTEGER :: res                                    !! Lowest subscript of the nearest value 
    REAL(kind=wp) :: l,u
    INTEGER       :: jl,jm,ju, nv
    nv = SIZE(vector) ; l = vector(1) ; u= vector(nv)
    ! Check for out of range value
    IF ((value>l.AND.value>u).OR.(value<l.AND.value<u)) THEN 
      res=1 ; IF (ABS(l-value) > ABS(u-value)) res=nv-1
      RETURN
    ENDIF
    ! Search in the array
    jl=0 ; ju=nv+1
    DO WHILE (ju-jl > 1)
      res=(ju+jl)/2
      IF (res == 0) RETURN   ! should never happen
      IF((u>=l).EQV.(value >= vector(res))) THEN
        jl=res
      ELSE
        ju=res
      ENDIF
    ENDDO
    res = jl
    RETURN
  END FUNCTION locate_ext

  FUNCTION locate_reg(value,vector) RESULT(res)
    !! Direct subscript access locator method
    !! 
    !! The function assumes __vector__ is regularly spaced and computes directly
    !! the lowest subscript using __vector__ step increment.
    REAL(kind=wp), INTENT(in)               :: value  !! Value to search
    REAL(kind=wp), INTENT(in), DIMENSION(:) :: vector !! Input vector to search in
    INTEGER :: res                                    
      !! Lowest subscript of the nearest value or __0__ if value is out of range.
    INTEGER       :: nv
    REAL(kind=wp) :: step,l,u
    res = 0 
    nv = SIZE(vector) 
    l = vector(1) ; u= vector(nv)
    IF ((value>l.AND.value>u).OR.(value<l.AND.value<u)) RETURN 
    step = (vector(nv)-vector(1))/(nv-1.)
    res = MIN(1+FLOOR((value-l)/step),nv-1)
    RETURN
  END FUNCTION locate_reg 

  FUNCTION locate_reg_ext(value,vector) RESULT(res)
    !! Direct subscript access locator method with extrapolation
    !!  
    !! The function performs the same computation than [[locators(module):locate_reg(function)]] 
    !! except that if __value__ is out of range, __1__ or __SIZE(vector)-1__ is returned 
    !! with respect to the nearest __vector__'s extremum.
    REAL(kind=wp), INTENT(in)               :: value  !! Value to search
    REAL(kind=wp), INTENT(in), DIMENSION(:) :: vector !! Input vector to search in
    INTEGER :: res                                    !! Lowest subscript of the nearest value 
    INTEGER       :: nv
    REAL(kind=wp) :: step,l,u
    res = 0 
    nv = SIZE(vector) 
    l = vector(1) ; u= vector(nv)
    IF ((value>l.AND.value>u).OR.(value<l.AND.value<u)) THEN 
      res=1 ; IF (ABS(l-value) > ABS(u-value)) res = nv -1
      RETURN
    ENDIF 
    step = (vector(nv)-vector(1))/(nv-1.)
    res = MIN(1+FLOOR((value-l)/step),nv-1)
    RETURN
  END FUNCTION locate_reg_ext 

END MODULE LOCATORS
