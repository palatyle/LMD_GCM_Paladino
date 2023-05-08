!*****************************************************************************
! Copyright(C) 2002-2011 Intel Corporation. All Rights Reserved.
!
! The source code, information  and  material ("Material") contained herein is
! owned  by Intel Corporation or its suppliers or licensors, and title to such
! Material remains  with Intel Corporation  or its suppliers or licensors. The
! Material  contains proprietary information  of  Intel or  its  suppliers and
! licensors. The  Material is protected by worldwide copyright laws and treaty
! provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
! modified, published, uploaded, posted, transmitted, distributed or disclosed
! in any way  without Intel's  prior  express written  permission. No  license
! under  any patent, copyright  or  other intellectual property rights  in the
! Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
! implication, inducement,  estoppel or  otherwise.  Any  license  under  such
! intellectual  property  rights must  be express  and  approved  by  Intel in
! writing.
!
! *Third Party trademarks are the property of their respective owners.
!
! Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
! this  notice or  any other notice embedded  in Materials by Intel or Intel's
! suppliers or licensors in any way.
!
!*****************************************************************************
! Content:
!    Intel(R) Math Kernel Library (MKL)
!    Discrete Fourier Transform Interface (DFTI)
!*****************************************************************************

MODULE MKL_DFT_TYPE

  TYPE, PUBLIC :: DFTI_DESCRIPTOR
     PRIVATE
     INTEGER :: dontuse
     ! Structure of this type is not used in Fortran code
     ! the pointer to this type is used only
  END TYPE DFTI_DESCRIPTOR

  !======================================================================
  ! These real type kind parameters are not for direct use
  !======================================================================

  INTEGER, PARAMETER :: DFTI_SPKP = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: DFTI_DPKP = SELECTED_REAL_KIND(15,307)

  !======================================================================
  ! Descriptor configuration parameters [default values in brackets]
  !======================================================================

  ! Domain for forward transform. No default value
  INTEGER, PARAMETER :: DFTI_FORWARD_DOMAIN = 0

  ! Dimensionality, or rank. No default value
  INTEGER, PARAMETER :: DFTI_DIMENSION = 1

  ! Length(s) of transform. No default value
  INTEGER, PARAMETER :: DFTI_LENGTHS = 2

  ! Floating point precision. No default value
  INTEGER, PARAMETER :: DFTI_PRECISION = 3

  ! Scale factor for forward transform [1.0]
  INTEGER, PARAMETER :: DFTI_FORWARD_SCALE = 4

  ! Scale factor for backward transform [1.0]
  INTEGER, PARAMETER :: DFTI_BACKWARD_SCALE = 5

  ! Exponent sign for forward transform [DFTI_NEGATIVE]
  ! INTEGER, PARAMETER :: DFTI_FORWARD_SIGN = 6 ! NOT IMPLEMENTED

  ! Number of data sets to be transformed [1]
  INTEGER, PARAMETER :: DFTI_NUMBER_OF_TRANSFORMS = 7

  ! Storage of finite complex-valued sequences in complex domain
  ! [DFTI_COMPLEX_COMPLEX]
  INTEGER, PARAMETER :: DFTI_COMPLEX_STORAGE = 8

  ! Storage of finite real-valued sequences in real domain
  ! [DFTI_REAL_REAL]
  INTEGER, PARAMETER :: DFTI_REAL_STORAGE = 9

  ! Storage of finite complex-valued sequences in conjugate-even
  ! domain [DFTI_COMPLEX_REAL]
  INTEGER, PARAMETER :: DFTI_CONJUGATE_EVEN_STORAGE = 10

  ! Placement of result [DFTI_INPLACE]
  INTEGER, PARAMETER :: DFTI_PLACEMENT = 11

  ! Generalized strides for input data layout
  ! [tigth, col-major for Fortran]
  INTEGER, PARAMETER :: DFTI_INPUT_STRIDES = 12

  ! Generalized strides for output data layout
  ! [tigth, col-major for Fortran]
  INTEGER, PARAMETER :: DFTI_OUTPUT_STRIDES = 13

  ! Distance between first input elements for multiple transforms [0]
  INTEGER, PARAMETER :: DFTI_INPUT_DISTANCE = 14

  ! Distance between first output elements for multiple transforms [0]
  INTEGER, PARAMETER :: DFTI_OUTPUT_DISTANCE = 15

  ! Effort spent in initialization [DFTI_MEDIUM]
  ! INTEGER, PARAMETER :: DFTI_INITIALIZATION_EFFORT = 16 ! NOT IMPLEMENTED

  ! Use of workspace during computation [DFTI_ALLOW]
  INTEGER, PARAMETER :: DFTI_WORKSPACE = 17

  ! Ordering of the result [DFTI_ORDERED]
  INTEGER, PARAMETER :: DFTI_ORDERING = 18

  ! Possible transposition of result [DFTI_NONE]
  INTEGER, PARAMETER :: DFTI_TRANSPOSE = 19

  ! User-settable descriptor name [""]
  INTEGER, PARAMETER :: DFTI_DESCRIPTOR_NAME = 20

  ! Packing format for DFTI_COMPLEX_REAL storage of finite
  ! conjugate-even sequences [DFTI_CCS_FORMAT]
  INTEGER, PARAMETER :: DFTI_PACKED_FORMAT = 21

  ! Commit status of the descriptor. Read-only parameter
  INTEGER, PARAMETER :: DFTI_COMMIT_STATUS = 22

  ! Version string for this DFTI implementation. Read-only parameter
  INTEGER, PARAMETER :: DFTI_VERSION = 23

  ! Ordering of the forward transform. Read-only parameter
  ! INTEGER, PARAMETER :: DFTI_FORWARD_ORDERING = 24 ! NOT IMPLEMENTED

  ! Ordering of the backward transform. Read-only parameter
  ! INTEGER, PARAMETER :: DFTI_BACKWARD_ORDERING = 25 ! NOT IMPLEMENTED

  ! Number of user threads that share the descriptor [1]
  INTEGER, PARAMETER :: DFTI_NUMBER_OF_USER_THREADS = 26

  !======================================================================
  ! Values of the descriptor configuration parameters
  !======================================================================

  ! DFTI_COMMIT_STATUS
  INTEGER, PARAMETER :: DFTI_COMMITTED = 30
  INTEGER, PARAMETER :: DFTI_UNCOMMITTED = 31

  ! DFTI_FORWARD_DOMAIN
  INTEGER, PARAMETER :: DFTI_COMPLEX = 32
  INTEGER, PARAMETER :: DFTI_REAL = 33
  ! INTEGER, PARAMETER :: DFTI_CONJUGATE_EVEN = 34 ! NOT IMPLEMENTED

  ! DFTI_PRECISION
  INTEGER, PARAMETER :: DFTI_SINGLE = 35
  INTEGER, PARAMETER :: DFTI_DOUBLE = 36

  ! DFTI_PRECISION for reduced size of statically linked application.
  ! Recommended use: modify statement 'USE MKL_DFTI' in your program,
  ! so that it reads as either of:
  ! USE MKL_DFTI, FORGET=>DFTI_SINGLE, DFTI_SINGLE=>DFTI_SINGLE_R
  ! USE MKL_DFTI, FORGET=>DFTI_DOUBLE, DFTI_DOUBLE=>DFTI_DOUBLE_R
  ! where word 'FORGET' can be any name not used in the program.
  REAL(DFTI_SPKP), PARAMETER :: DFTI_SINGLE_R = 35
  REAL(DFTI_DPKP), PARAMETER :: DFTI_DOUBLE_R = 36

  ! DFTI_FORWARD_SIGN
  ! INTEGER, PARAMETER :: DFTI_NEGATIVE = 37 ! NOT IMPLEMENTED
  ! INTEGER, PARAMETER :: DFTI_POSITIVE = 38 ! NOT IMPLEMENTED

  ! DFTI_COMPLEX_STORAGE and DFTI_CONJUGATE_EVEN_STORAGE
  INTEGER, PARAMETER :: DFTI_COMPLEX_COMPLEX = 39
  INTEGER, PARAMETER :: DFTI_COMPLEX_REAL = 40

  ! DFTI_REAL_STORAGE
  INTEGER, PARAMETER :: DFTI_REAL_COMPLEX = 41
  INTEGER, PARAMETER :: DFTI_REAL_REAL = 42

  ! DFTI_PLACEMENT
  INTEGER, PARAMETER :: DFTI_INPLACE = 43 ! Result overwrites input
  INTEGER, PARAMETER :: DFTI_NOT_INPLACE  = 44 ! Have another place for result

  ! DFTI_INITIALIZATION_EFFORT
  ! INTEGER, PARAMETER :: DFTI_LOW = 45 ! NOT IMPLEMENTED
  ! INTEGER, PARAMETER :: DFTI_MEDIUM = 46 ! NOT IMPLEMENTED
  ! INTEGER, PARAMETER :: DFTI_HIGH = 47 ! NOT IMPLEMENTED

  ! DFTI_ORDERING
  INTEGER, PARAMETER :: DFTI_ORDERED = 48
  INTEGER, PARAMETER :: DFTI_BACKWARD_SCRAMBLED = 49
  ! INTEGER, PARAMETER :: DFTI_FORWARD_SCRAMBLED  = 50 ! NOT IMPLEMENTED

  ! Allow/avoid certain usages
  INTEGER, PARAMETER :: DFTI_ALLOW = 51 ! Allow transposition or workspace
  INTEGER, PARAMETER :: DFTI_AVOID = 52 ! Avoid auxiliary storage
  INTEGER, PARAMETER :: DFTI_NONE = 53

  ! DFTI_PACKED_FORMAT
  ! (for storing congugate-even finite sequence in real array)
  INTEGER, PARAMETER :: DFTI_CCS_FORMAT = 54  ! Complex conjugate-symmetric
  INTEGER, PARAMETER :: DFTI_PACK_FORMAT = 55 ! Pack format for real DFT
  INTEGER, PARAMETER :: DFTI_PERM_FORMAT = 56 ! Perm format for real DFT
  INTEGER, PARAMETER :: DFTI_CCE_FORMAT = 57  ! Complex conjugate-even

  !======================================================================
  ! Error classes
  !======================================================================
  INTEGER, PARAMETER :: DFTI_NO_ERROR = 0
  INTEGER, PARAMETER :: DFTI_MEMORY_ERROR = 1
  INTEGER, PARAMETER :: DFTI_INVALID_CONFIGURATION = 2
  INTEGER, PARAMETER :: DFTI_INCONSISTENT_CONFIGURATION = 3
  INTEGER, PARAMETER :: DFTI_MULTITHREADED_ERROR = 4
  INTEGER, PARAMETER :: DFTI_BAD_DESCRIPTOR = 5
  INTEGER, PARAMETER :: DFTI_UNIMPLEMENTED = 6
  INTEGER, PARAMETER :: DFTI_MKL_INTERNAL_ERROR = 7
  INTEGER, PARAMETER :: DFTI_NUMBER_OF_THREADS_ERROR = 8
  INTEGER, PARAMETER :: DFTI_1D_LENGTH_EXCEEDS_INT32 = 9

  ! Maximum length of error string
  INTEGER, PARAMETER :: DFTI_MAX_MESSAGE_LENGTH = 80

  ! Maximum length of user-settable descriptor name
  INTEGER, PARAMETER :: DFTI_MAX_NAME_LENGTH = 10

  ! Maximum length of MKL version string
  INTEGER, PARAMETER :: DFTI_VERSION_LENGTH = 198

  ! (deprecated parameter)
  INTEGER, PARAMETER :: DFTI_ERROR_CLASS = 60

END MODULE MKL_DFT_TYPE

