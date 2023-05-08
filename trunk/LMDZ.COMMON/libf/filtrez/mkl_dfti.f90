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


MODULE MKL_DFTI

  USE MKL_DFT_TYPE

  INTERFACE DftiCreateDescriptor

     FUNCTION dfti_create_descriptor_1d(desc, precision, domain, dim, length)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_create_descriptor_1d
       !DEC$ ATTRIBUTES REFERENCE :: dfti_create_descriptor_1d
       INTEGER dfti_create_descriptor_1d
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       INTEGER, INTENT(IN) :: precision
       INTEGER, INTENT(IN) :: domain
       INTEGER, INTENT(IN) :: dim, length
     END FUNCTION dfti_create_descriptor_1d

     FUNCTION dfti_create_descriptor_highd(desc, precision, domain, dim,length)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_create_descriptor_highd
       !DEC$ ATTRIBUTES REFERENCE :: dfti_create_descriptor_highd
       INTEGER dfti_create_descriptor_highd
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       INTEGER, INTENT(IN) :: precision
       INTEGER, INTENT(IN) :: domain
       INTEGER, INTENT(IN) :: dim
       INTEGER, INTENT(IN), DIMENSION(*) :: length
     END FUNCTION dfti_create_descriptor_highd

     FUNCTION dfti_create_descriptor_s_1d(desc, s, dom, one, dim)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_create_descriptor_s_1d
       !DEC$ ATTRIBUTES REFERENCE :: dfti_create_descriptor_s_1d
       INTEGER dfti_create_descriptor_s_1d
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_SPKP), INTENT(IN) :: s
       INTEGER, INTENT(IN) :: dom
       INTEGER, INTENT(IN) :: one
       INTEGER, INTENT(IN) :: dim
     END FUNCTION dfti_create_descriptor_s_1d

     FUNCTION dfti_create_descriptor_s_md(desc, s, dom, many, dims)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_create_descriptor_s_md
       !DEC$ ATTRIBUTES REFERENCE :: dfti_create_descriptor_s_md
       INTEGER dfti_create_descriptor_s_md
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_SPKP), INTENT(IN) :: s
       INTEGER, INTENT(IN) :: dom
       INTEGER, INTENT(IN) :: many
       INTEGER, INTENT(IN), DIMENSION(*) :: dims
     END FUNCTION dfti_create_descriptor_s_md

     FUNCTION dfti_create_descriptor_d_1d(desc, d, dom, one, dim)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_create_descriptor_d_1d
       !DEC$ ATTRIBUTES REFERENCE :: dfti_create_descriptor_d_1d
       INTEGER dfti_create_descriptor_d_1d
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_DPKP), INTENT(IN) :: d
       INTEGER, INTENT(IN) :: dom
       INTEGER, INTENT(IN) :: one
       INTEGER, INTENT(IN) :: dim
     END FUNCTION dfti_create_descriptor_d_1d

     FUNCTION dfti_create_descriptor_d_md(desc, d, dom, many, dims)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_create_descriptor_d_md
       !DEC$ ATTRIBUTES REFERENCE :: dfti_create_descriptor_d_md
       INTEGER dfti_create_descriptor_d_md
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_DPKP), INTENT(IN) :: d
       INTEGER, INTENT(IN) :: dom
       INTEGER, INTENT(IN) :: many
       INTEGER, INTENT(IN), DIMENSION(*) :: dims
     END FUNCTION dfti_create_descriptor_d_md

  END INTERFACE

  INTERFACE DftiCopyDescriptor

     FUNCTION dfti_copy_descriptor_external(desc, new_desc)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_copy_descriptor_external
       !DEC$ ATTRIBUTES REFERENCE :: dfti_copy_descriptor_external
       INTEGER dfti_copy_descriptor_external
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       TYPE(DFTI_DESCRIPTOR), POINTER :: new_desc
     END FUNCTION dfti_copy_descriptor_external

  END INTERFACE

  INTERFACE DftiCommitDescriptor

     FUNCTION dfti_commit_descriptor_external(desc)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_commit_descriptor_external
       !DEC$ ATTRIBUTES REFERENCE :: dfti_commit_descriptor_external
       INTEGER dfti_commit_descriptor_external
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_commit_descriptor_external

  END INTERFACE

  INTERFACE DftiSetValue

     FUNCTION dfti_set_value_intval(desc, OptName, IntVal)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_set_value_intval
       !DEC$ ATTRIBUTES REFERENCE :: dfti_set_value_intval
       INTEGER dfti_set_value_intval
       INTEGER, INTENT(IN) :: OptName
       INTEGER, INTENT(IN) :: IntVal
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_set_value_intval

     FUNCTION dfti_set_value_sglval(desc, OptName, sglval)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_set_value_sglval
       !DEC$ ATTRIBUTES REFERENCE :: dfti_set_value_sglval
       INTEGER dfti_set_value_sglval
       INTEGER, INTENT(IN) :: OptName
       REAL(DFTI_SPKP), INTENT(IN) :: sglval
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_set_value_sglval

     FUNCTION dfti_set_value_dblval(desc, OptName, DblVal)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_set_value_dblval
       !DEC$ ATTRIBUTES REFERENCE :: dfti_set_value_dblval
       INTEGER dfti_set_value_dblval
       INTEGER, INTENT(IN) :: OptName
       REAL(DFTI_DPKP), INTENT(IN) :: DblVal
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_set_value_dblval

     FUNCTION dfti_set_value_intvec(desc, OptName, IntVec)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_set_value_intvec
       !DEC$ ATTRIBUTES REFERENCE :: dfti_set_value_intvec
       INTEGER dfti_set_value_intvec
       INTEGER, INTENT(IN) :: OptName
       INTEGER, INTENT(IN), DIMENSION(*) :: IntVec
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_set_value_intvec

     FUNCTION dfti_set_value_chars(desc, OptName, Chars)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_set_value_chars
       !DEC$ ATTRIBUTES REFERENCE :: dfti_set_value_chars
       INTEGER dfti_set_value_chars
       INTEGER, INTENT(IN) :: OptName
       CHARACTER(*), INTENT(IN) :: Chars
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_set_value_chars

  END INTERFACE

  INTERFACE DftiGetValue

     FUNCTION dfti_get_value_intval(desc, OptName, IntVal)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_get_value_intval
       !DEC$ ATTRIBUTES REFERENCE :: dfti_get_value_intval
       INTEGER dfti_get_value_intval
       INTEGER, INTENT(IN) :: OptName
       INTEGER, INTENT(OUT) :: IntVal
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_get_value_intval

     FUNCTION dfti_get_value_sglval(desc, OptName, sglval)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_get_value_sglval
       !DEC$ ATTRIBUTES REFERENCE :: dfti_get_value_sglval
       INTEGER dfti_get_value_sglval
       INTEGER, INTENT(IN) :: OptName
       REAL(DFTI_SPKP), INTENT(OUT) :: sglval
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_get_value_sglval

     FUNCTION dfti_get_value_dblval(desc, OptName, DblVal)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_get_value_dblval
       !DEC$ ATTRIBUTES REFERENCE :: dfti_get_value_dblval
       INTEGER dfti_get_value_dblval
       INTEGER, INTENT(IN) :: OptName
       REAL(DFTI_DPKP), INTENT(OUT) :: DblVal
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_get_value_dblval

     FUNCTION dfti_get_value_intvec(desc, OptName, IntVec)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_get_value_intvec
       !DEC$ ATTRIBUTES REFERENCE :: dfti_get_value_intvec
       INTEGER dfti_get_value_intvec
       INTEGER, INTENT(IN) :: OptName
       INTEGER, INTENT(OUT), DIMENSION(*) :: IntVec
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_get_value_intvec

     FUNCTION dfti_get_value_chars(desc, OptName, Chars)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_get_value_chars
       !DEC$ ATTRIBUTES REFERENCE :: dfti_get_value_chars
       INTEGER dfti_get_value_chars
       INTEGER, INTENT(IN) :: OptName
       CHARACTER(*), INTENT(OUT) :: Chars
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_get_value_chars

  END INTERFACE

  INTERFACE DftiComputeForward

     FUNCTION dfti_compute_forward_s(desc,sSrcDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_s
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_s
       INTEGER dfti_compute_forward_s
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_SPKP), INTENT(INOUT), DIMENSION(*) :: sSrcDst
     END FUNCTION dfti_compute_forward_s

     FUNCTION dfti_compute_forward_c(desc,cSrcDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_c
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_c
       INTEGER dfti_compute_forward_c
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_SPKP), INTENT(INOUT), DIMENSION(*) :: cSrcDst
     END FUNCTION dfti_compute_forward_c

     FUNCTION dfti_compute_forward_ss(desc,sSrcDstRe,sSrcDstIm)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_ss
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_ss
       INTEGER dfti_compute_forward_ss
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_SPKP), DIMENSION(*) :: sSrcDstRe
       REAL(DFTI_SPKP), DIMENSION(*) :: sSrcDstIm
     END FUNCTION dfti_compute_forward_ss

     FUNCTION dfti_compute_forward_sc(desc,sSrc,cDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_sc
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_sc
       INTEGER dfti_compute_forward_sc
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_SPKP), INTENT(IN), DIMENSION(*) :: sSrc
       COMPLEX(DFTI_SPKP), INTENT(OUT), DIMENSION(*) :: cDst
     END FUNCTION dfti_compute_forward_sc

     FUNCTION dfti_compute_forward_cc(desc,cSrc,cDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_cc
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_cc
       INTEGER dfti_compute_forward_cc
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_SPKP), INTENT(IN), DIMENSION(*) :: cSrc
       COMPLEX(DFTI_SPKP), INTENT(OUT), DIMENSION(*) :: cDst
     END FUNCTION dfti_compute_forward_cc

     FUNCTION dfti_compute_forward_ssss(desc,sSrcRe,sSrcIm,sDstRe,sDstIm)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_ssss
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_ssss
       INTEGER dfti_compute_forward_ssss
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_SPKP), INTENT(IN), DIMENSION(*) :: sSrcRe
       REAL(DFTI_SPKP), INTENT(IN), DIMENSION(*) :: sSrcIm
       REAL(DFTI_SPKP), INTENT(OUT), DIMENSION(*) :: sDstRe
       REAL(DFTI_SPKP), INTENT(OUT), DIMENSION(*) :: sDstIm
     END FUNCTION dfti_compute_forward_ssss

     FUNCTION dfti_compute_forward_d(desc,dSrcDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_d
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_d
       INTEGER dfti_compute_forward_d
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_DPKP), INTENT(INOUT), DIMENSION(*) :: dSrcDst
     END FUNCTION dfti_compute_forward_d

     FUNCTION dfti_compute_forward_z(desc,zSrcDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_z
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_z
       INTEGER dfti_compute_forward_z
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_DPKP), INTENT(INOUT), DIMENSION(*) :: zSrcDst
     END FUNCTION dfti_compute_forward_z

     FUNCTION dfti_compute_forward_dd(desc,dSrcDstRe,dSrcDstIm)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_dd
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_dd
       INTEGER dfti_compute_forward_dd
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_DPKP), DIMENSION(*) :: dSrcDstRe
       REAL(DFTI_DPKP), DIMENSION(*) :: dSrcDstIm
     END FUNCTION dfti_compute_forward_dd

     FUNCTION dfti_compute_forward_dz(desc,dSrc,zDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_dz
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_dz
       INTEGER dfti_compute_forward_dz
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_DPKP), INTENT(IN), DIMENSION(*) :: dSrc
       COMPLEX(DFTI_DPKP), INTENT(OUT), DIMENSION(*) :: zDst
     END FUNCTION dfti_compute_forward_dz

     FUNCTION dfti_compute_forward_zz(desc,zSrc,zDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_zz
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_zz
       INTEGER dfti_compute_forward_zz
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_DPKP), INTENT(IN), DIMENSION(*) :: zSrc
       COMPLEX(DFTI_DPKP), INTENT(OUT), DIMENSION(*) :: zDst
     END FUNCTION dfti_compute_forward_zz

     FUNCTION dfti_compute_forward_dddd(desc,dSrcRe,dSrcIm,dDstRe,dDstIm)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_forward_dddd
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_forward_dddd
       INTEGER dfti_compute_forward_dddd
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_DPKP), INTENT(IN), DIMENSION(*) :: dSrcRe
       REAL(DFTI_DPKP), INTENT(IN), DIMENSION(*) :: dSrcIm
       REAL(DFTI_DPKP), INTENT(OUT), DIMENSION(*) :: dDstRe
       REAL(DFTI_DPKP), INTENT(OUT), DIMENSION(*) :: dDstIm
     END FUNCTION dfti_compute_forward_dddd

  END INTERFACE DftiComputeForward

  INTERFACE DftiComputeBackward

     FUNCTION dfti_compute_backward_s(desc,sSrcDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_s
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_s
       INTEGER dfti_compute_backward_s
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_SPKP), INTENT(INOUT), DIMENSION(*) :: sSrcDst
     END FUNCTION dfti_compute_backward_s

     FUNCTION dfti_compute_backward_c(desc,cSrcDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_c
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_c
       INTEGER dfti_compute_backward_c
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_SPKP), INTENT(INOUT), DIMENSION(*) :: cSrcDst
     END FUNCTION dfti_compute_backward_c

     FUNCTION dfti_compute_backward_ss(desc,sSrcDstRe,sSrcDstIm)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_ss
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_ss
       INTEGER dfti_compute_backward_ss
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_SPKP), DIMENSION(*) :: sSrcDstRe
       REAL(DFTI_SPKP), DIMENSION(*) :: sSrcDstIm
     END FUNCTION dfti_compute_backward_ss

     FUNCTION dfti_compute_backward_cs(desc,cSrc,sDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_cs
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_cs
       INTEGER dfti_compute_backward_cs
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_SPKP), INTENT(IN), DIMENSION(*) :: cSrc
       REAL(DFTI_SPKP), INTENT(OUT), DIMENSION(*) :: sDst
     END FUNCTION dfti_compute_backward_cs

     FUNCTION dfti_compute_backward_cc(desc,cSrc,cDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_cc
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_cc
       INTEGER dfti_compute_backward_cc
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_SPKP), INTENT(IN), DIMENSION(*) :: cSrc
       COMPLEX(DFTI_SPKP), INTENT(OUT), DIMENSION(*) :: cDst
     END FUNCTION dfti_compute_backward_cc

     FUNCTION dfti_compute_backward_ssss(desc,sSrcRe,sSrcIm,sDstRe,sDstIm)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_ssss
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_ssss
       INTEGER dfti_compute_backward_ssss
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_SPKP), INTENT(IN), DIMENSION(*) :: sSrcRe
       REAL(DFTI_SPKP), INTENT(IN), DIMENSION(*) :: sSrcIm
       REAL(DFTI_SPKP), INTENT(OUT), DIMENSION(*) :: sDstRe
       REAL(DFTI_SPKP), INTENT(OUT), DIMENSION(*) :: sDstIm
     END FUNCTION dfti_compute_backward_ssss

     FUNCTION dfti_compute_backward_d(desc,dSrcDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_d
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_d
       INTEGER dfti_compute_backward_d
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_DPKP), INTENT(INOUT), DIMENSION(*) :: dSrcDst
     END FUNCTION dfti_compute_backward_d

     FUNCTION dfti_compute_backward_z(desc,zSrcDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_z
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_z
       INTEGER dfti_compute_backward_z
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_DPKP), INTENT(INOUT), DIMENSION(*) :: zSrcDst
     END FUNCTION dfti_compute_backward_z

     FUNCTION dfti_compute_backward_dd(desc,dSrcDstRe,dSrcDstIm)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_dd
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_dd
       INTEGER dfti_compute_backward_dd
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_DPKP), DIMENSION(*) :: dSrcDstRe
       REAL(DFTI_DPKP), DIMENSION(*) :: dSrcDstIm
     END FUNCTION dfti_compute_backward_dd

     FUNCTION dfti_compute_backward_zd(desc,zSrc,dDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_zd
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_zd
       INTEGER dfti_compute_backward_zd
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_DPKP), INTENT(IN), DIMENSION(*) :: zSrc
       REAL(DFTI_DPKP), INTENT(OUT), DIMENSION(*) :: dDst
     END FUNCTION dfti_compute_backward_zd

     FUNCTION dfti_compute_backward_zz(desc,zSrc,zDst)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_zz
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_zz
       INTEGER dfti_compute_backward_zz
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       COMPLEX(DFTI_DPKP), INTENT(IN), DIMENSION(*) :: zSrc
       COMPLEX(DFTI_DPKP), INTENT(OUT), DIMENSION(*) :: zDst
     END FUNCTION dfti_compute_backward_zz

     FUNCTION dfti_compute_backward_dddd(desc,dSrcRe,dSrcIm,dDstRe,dDstIm)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_compute_backward_dddd
       !DEC$ ATTRIBUTES REFERENCE :: dfti_compute_backward_dddd
       INTEGER dfti_compute_backward_dddd
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
       REAL(DFTI_DPKP), INTENT(IN), DIMENSION(*) :: dSrcRe
       REAL(DFTI_DPKP), INTENT(IN), DIMENSION(*) :: dSrcIm
       REAL(DFTI_DPKP), INTENT(OUT), DIMENSION(*) :: dDstRe
       REAL(DFTI_DPKP), INTENT(OUT), DIMENSION(*) :: dDstIm
     END FUNCTION dfti_compute_backward_dddd

  END INTERFACE DftiComputeBackward

  INTERFACE DftiFreeDescriptor

     FUNCTION dfti_free_descriptor_external(desc)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_free_descriptor_external
       !DEC$ ATTRIBUTES REFERENCE :: dfti_free_descriptor_external
       INTEGER dfti_free_descriptor_external
       TYPE(DFTI_DESCRIPTOR), POINTER :: desc
     END FUNCTION dfti_free_descriptor_external

  END INTERFACE

  INTERFACE DftiErrorClass

     FUNCTION dfti_error_class_external(Status, ErrorClass)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_error_class_external
       !DEC$ ATTRIBUTES REFERENCE :: dfti_error_class_external
       LOGICAL dfti_error_class_external
       INTEGER, INTENT(IN) :: Status
       INTEGER, INTENT(IN) :: ErrorClass
     END FUNCTION dfti_error_class_external

  END INTERFACE

  INTERFACE DftiErrorMessage

     FUNCTION dfti_error_message_external(Status)
       USE MKL_DFT_TYPE
       !DEC$ ATTRIBUTES C :: dfti_error_message_external
       !DEC$ ATTRIBUTES REFERENCE :: dfti_error_message_external
       CHARACTER(LEN=DFTI_MAX_MESSAGE_LENGTH) :: dfti_error_message_external
       INTEGER, INTENT(IN) :: Status
     END FUNCTION dfti_error_message_external

  END INTERFACE

END MODULE MKL_DFTI

