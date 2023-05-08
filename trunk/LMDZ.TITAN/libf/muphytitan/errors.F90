! Copyright Jérémie Burgalat (2010-2015,2017)
! 
! jeremie.burgalat@univ-reims.fr
! 
! This software is a computer program whose purpose is to provide configuration 
! file and command line arguments parsing features to Fortran programs.
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

!! file: errors.F90
!! summary: Errors handling source file.
!! author: J. Burgalat
!! date: 2013-2015,2017

#include "defined.h"

MODULE ERRORS
  !! Error handler module
  !!
  !! This module provides a single derived type, [[error(type)]] which is used in all
  !! other parts of the library in order to handle errors.
  USE, INTRINSIC :: ISO_C_BINDING
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : stdout=>OUTPUT_UNIT, stderr=>ERROR_UNIT

  IMPLICIT NONE

  PUBLIC

  PRIVATE :: error_equals,error_equals_int,error_differs,error_differs_int, &
             msg_length


  INTEGER, PARAMETER :: msg_length = 250 !! Length of error message.

  TYPE, PUBLIC :: error
    !! Define an error
    !!
    !! The following derived type represents in the simplest way (I believe) an error which 
    !! stores:
    !!
    !! - An integer to numerically identify the error
    !! - A string (250 chars max) with an appropriate error message
    !! - A bounded procedure to get a string representation of the error (if bounded 
    !!   procedures are supported by the library).
    !! - internal subroutines for derived type IO WRITE statement (if Derived IO 
    !!   subroutines are supported by the library).
    !!
    !! error type comes also with two operators ("==", "/=") to compare error type with
    !! another one or an integer.
    !! If an error is not initialized explicitly, then it is set to [[errors(module):noerror(variable)]].
    CHARACTER(len=msg_length) :: msg = "No error"
      !! Message associated to the error
      !! @note
      !! The message should be short (250 characters maximum) and explicit.
    INTEGER :: id = 0
      !! Numerical identifier of the error
      !! @note 
      !! The error identifier is used to test the equality/inequality of two error objects. 
#if HAVE_FTNPROC
    CONTAINS
      PROCEDURE, PUBLIC :: to_string => error_to_string
        !! Get a string representation of the error
#endif
  END TYPE error

  INTERFACE 
    !! Clean subroutine interface
    SUBROUTINE clean_callback(err)
      !! A subroutine that may perform cleaning computation(s) before exit
      IMPORT error
      IMPLICIT NONE
      TYPE(error), INTENT(in) :: err
        !! An error object with the input error 
    END SUBROUTINE clean_callback
  END INTERFACE

  INTERFACE 
    subroutine abort_() bind(C, name="abort")
    end subroutine
  END INTERFACE

  INTERFACE assert
    !! _Raise_ an assertion.
    !!
    !! An assertion can be understood as a development error that should be raised in production mode.
    MODULE PROCEDURE :: assert_r,assert_w
  END INTERFACE assert

  !> error equality operator 
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE error_equals, error_equals_int
  END INTERFACE

  !> error inequality operator 
  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE error_differs, error_differs_int
  END INTERFACE

  !> The no error error !
  TYPE(error), PUBLIC, PARAMETER :: noerror = error("No error",0)

  CONTAINS

!===============================================================================
! error TYPE RELATED METHODS
!===============================================================================

  FUNCTION error_equals(this, other) RESULT(res)
    !! Check if two error objects are equivalent
    TYPE(error), INTENT(in) :: this, & !! The first error object to compare 
                               other   !! The second error object to compare 
    LOGICAL :: res                     !! .true. if __this__ and __other__ identifiers are the same, .false. otherwise
    res = (this%id == other%id)
    RETURN
  END FUNCTION error_equals

  FUNCTION error_equals_int(this, id) RESULT(res)
    !! Check if an error id is equal to a given integer
    TYPE(error), INTENT(in) :: this !! An error object reference
    INTEGER, INTENT(in)     :: id   !! An integer to compare to __this__ identifier
    LOGICAL :: res                  !! .true. if __this__ identifier and __id__ have the same value, .false. otherwise
    res = (this%id == id)
    RETURN
  END FUNCTION error_equals_int

  FUNCTION error_differs(this, other) RESULT(res)
    !! Check if two error objects are different
    TYPE(error), INTENT(in) :: this, & !! The first error object to compare 
                               other   !! The second error object to compare 
    LOGICAL :: res                     !! .false. if __this__ and __other__ identifiers are the same, .true. otherwise
    res = (this%id /= other%id)
    RETURN
  END FUNCTION error_differs

  FUNCTION error_differs_int(this, id) RESULT(res)
    !! Check if an error id is different from a given integer
    TYPE(error), INTENT(in) :: this !! An error object reference
    INTEGER, INTENT(in)     :: id   !! An integer to compare to __this__ identifier
    LOGICAL :: res                  !! .false. if __this__ identifier and __id__ have the same value, .true. otherwise
    res = (this%id /= id)
    RETURN
  END FUNCTION error_differs_int

  FUNCTION error_to_string(this,progname,as_warning) RESULT(str)
    !! (simple) String representation of the error
    !!
    !! The function returns a very simple formatted string with the error.
    OBJECT(error), INTENT(in)              :: this
      !! An error object reference
    CHARACTER(len=*), INTENT(in), OPTIONAL :: progname
      !! An optional string with the name of the program
    LOGICAL, INTENT(in), OPTIONAL          :: as_warning 
      !! An optional boolean flag to print the message as warning rather than as error (default to .false.).
    CHARACTER(len=:), ALLOCATABLE :: str
      !! An allocatable string with the string representation of the error
    CHARACTER(len=:), ALLOCATABLE :: pref 
    pref = "error: "
    IF (PRESENT(as_warning)) THEN ; IF (as_warning) pref = "warning: " ; ENDIF
    IF (PRESENT(progname)) THEN
      IF (LEN_TRIM(progname) /=0) THEN
        str = TRIM(progname)//': '//pref//TRIM(this%msg)
      ELSE
        str = pref//TRIM(this%msg)
      ENDIF
    ELSE
      str = pref//TRIM(this%msg)
    ENDIF
    RETURN 
  END FUNCTION error_to_string

  SUBROUTINE aborting(err)
    !! Abort the program with specific exit code
    !!
    !! The method prints the message of the given error object and
    !! stops the program using exit() subroutine.
    TYPE(error), INTENT(in) :: err
      !! An error object
    IF (err /= 0) THEN
      WRITE(*,'(a)') error_to_string(err)
      CALL EXIT(err%id)
    ENDIF
  END SUBROUTINE aborting

  SUBROUTINE assert_r(test,reason)
    !! _Raise_ an assertion.
    !! 
    !! The method raises an assertion and stops the execution if __test__ is .false.
    !! 
    !! @note
    !! If ISO_C_BINDING module is available, the method calls the method abort from the C standard library. Doing so,
    !! developer is able to debug the source code by getting the backtrace of the execution.
    !! In other situation, the method simply uses the Fortran STOP statement which makes its usage... useless. 
   LOGICAL, INTENT(in)          :: test
     !! Expression to test.
   CHARACTER(len=*), INTENT(in) :: reason
     !! Optional assertion reason.
   IF (.NOT.test) THEN
     WRITE(stderr,'(a)') "assertion: "//reason
     call abort_()
   ENDIF
  END SUBROUTINE assert_r

  SUBROUTINE assert_w(test,where,reason)
    !! _Raise_ an assertion.
    !! 
    !! The method raises an assertion and stops the execution if __test__ is .false.
    !! 
    !! See [[errors(module):assert_r(subroutine)]] remark.
    LOGICAL, INTENT(in)         :: test
     !! Expression to test.
   CHARACTER(len=*), INTENT(in) :: where
     !! Optional _location_ of the assertion.
   CHARACTER(len=*), INTENT(in) :: reason
     !! Optional assertion reason.
   IF (.NOT.test) THEN
     WRITE(stderr,'(a)') "assertion in "//where//": "//reason
     call abort_()
   ENDIF
  END SUBROUTINE assert_w

  FUNCTION free_lun() RESULT(lu)
    !> Get the first free logical unit
    !!
    !! The function loops from 7 to 9999 and returns the first free logical unit.
    !! @note
    !! According to Fortran standard, the maximum value for a lun is processor
    !! dependent. I just assume that [7,9999] is a valid range and I believe that 
    !! 9992 files to be opened is far enough for any program !
    !! @note
    !! If you intend to use loggers object from this library, you should keep in
    !! mind that loggers open files with the first free logical unit. Consequently
    !! if you need to perform I/O operations you should use this function to get a
    !! free lun instead of just randomly set a lun ! 
    INTEGER :: lu
      !! First free logical unit in the range [7,9999]  or -1 if no lun is available
    INTEGER, PARAMETER :: mxlu = 9999
    LOGICAL :: notfree
    lu = 6 ; notfree = .true.
    DO WHILE(notfree.AND.lu<=mxlu)
      lu=lu+1 ; INQUIRE(unit=lu,OPENED=notfree)
    ENDDO
    IF (lu >= mxlu) lu = -1
  END FUNCTION free_lun



END MODULE ERRORS

