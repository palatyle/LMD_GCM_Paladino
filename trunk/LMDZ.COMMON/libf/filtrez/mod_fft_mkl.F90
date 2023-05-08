MODULE mod_fft_mkl
#ifdef FFT_MKL

  USE MKL_DFTI
  
  REAL,SAVE                :: scale_factor
  INTEGER,SAVE             :: vsize
  INTEGER,PARAMETER        :: inc=1
 
!  TYPE FFT_HANDLE
!    TYPE(DFTI_DESCRIPTOR), POINTER :: Pt
!    LOGICAL :: IsAllocated
!  END TYPE FFT_HANDLE
  
!  TYPE(FFT_HANDLE),SAVE,ALLOCATABLE :: Forward_Handle(:)
!  TYPE(FFT_HANDLE),SAVE,ALLOCATABLE :: Backward_Handle(:)
!!$OMP THREADPRIVATE(Forward_Handle,Backward_Handle)  
  
CONTAINS
  
  SUBROUTINE Init_fft(iim,nb_vect_max)
    IMPLICIT NONE
    INTEGER :: iim
    INTEGER :: nb_vect_max
    REAL    :: rtmp=1.
    COMPLEX :: ctmp
    INTEGER :: itmp=1
    INTEGER :: isign=0
    INTEGER :: ierr
    
    vsize=iim
    scale_factor=1./SQRT(1.*vsize)
!    ALLOCATE(Forward_Handle(nb_vect_max))
!    ALLOCATE(Backward_Handle(nb_vect_max))
    
!    Forward_Handle(:)%IsAllocated=.FALSE.
!    Backward_Handle(:)%IsAllocated=.FALSE.
    
!    ALLOCATE(Table_forward(2*vsize+64))
!    ALLOCATE(Table_backward(2*vsize+64))
!    
!    CALL DZFFTM(isign,vsize,itmp,scale_factor,rtmp,vsize+inc,ctmp,vsize/2+1,table_forward,rtmp,ierr)
!    
!    CALL ZDFFTM(isign,vsize,itmp,scale_factor,ctmp,vsize/2+1,rtmp,vsize+inc,table_backward,rtmp,ierr)

!    ierr = DftiCreateDescriptor( FFT_Handle, DFTI_DOUBLE, DFTI_REAL, 1, vsize )
!    ierr = DftiSetValue(FFT_Handle,DFTI_NUMBER_OF_TRANSFORMS,nb_vect)
!    ierr = DftiSetValue(FFT_Handle,DFTI_FORWARD_SCALE,scale_factor)
!    ierr = DftiSetValue(FFT_Handle,DFTI_BACKWARD_SCALE,scale_factor)
!    ierr = DftiSetValue(FFT_Handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
!    ierr = DftiSetValue(Desc_Handle, DFTI_INPUT_DISTANCE, vsize+inc)
!    ierr = DftiSetValue(Desc_Handle, DFTI_OUTPUT_DISTANCE, vsize)
!    ierr = DftiCommitDescriptor( FFT_HANDLE )

  END SUBROUTINE Init_fft
  
  
  SUBROUTINE fft_forward(vect,TF_vect,nb_vect)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: nb_vect
    REAL,INTENT(IN)     :: vect((vsize+inc)*nb_vect)
    COMPLEX,INTENT(OUT) :: TF_vect((vsize/2+1)*nb_vect)
    REAL                :: work(4*vsize*nb_vect)
    INTEGER             :: ierr
    INTEGER, PARAMETER :: isign=-1
    REAL               :: vect_out((vsize+inc)*nb_vect)
    TYPE(DFTI_DESCRIPTOR), POINTER :: FFT_Handle
    
!    IF ( .NOT. Forward_handle(nb_vect)%IsAllocated) THEN
      ierr = DftiCreateDescriptor( FFT_Handle, DFTI_DOUBLE, DFTI_REAL, 1, vsize )
      ierr = DftiSetValue(FFT_Handle,DFTI_NUMBER_OF_TRANSFORMS,nb_vect)
      ierr = DftiSetValue(FFT_Handle,DFTI_FORWARD_SCALE,scale_factor)
      ierr = DftiSetValue(FFT_Handle,DFTI_BACKWARD_SCALE,scale_factor)
      ierr = DftiSetValue(FFT_Handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
      ierr = DftiSetValue(FFT_Handle, DFTI_INPUT_DISTANCE, vsize+inc)
      ierr = DftiSetValue(FFT_Handle, DFTI_OUTPUT_DISTANCE, (vsize/2+1)*2)
      ierr = DftiCommitDescriptor( FFT_Handle )
!      Forward_handle(nb_vect)%IsAllocated=.TRUE.
!    ENDIF
    
    ierr = DftiComputeForward( FFT_Handle, vect, TF_vect )

    ierr = DftiFreeDescriptor( FFT_Handle )

!    ierr = DftiCreateDescriptor( FFT_Handle, DFTI_DOUBLE, DFTI_REAL, 1, vsize )
!    ierr = DftiSetValue(FFT_Handle,DFTI_NUMBER_OF_TRANSFORMS,nb_vect)
!    ierr = DftiSetValue(FFT_Handle,DFTI_FORWARD_SCALE,scale_factor)
!    ierr = DftiSetValue(FFT_Handle,DFTI_BACKWARD_SCALE,scale_factor)
!    ierr = DftiSetValue(FFT_Handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
!    ierr = DftiSetValue(FFT_HANDLE, DFTI_INPUT_DISTANCE, vsize/2+1)
!    ierr = DftiSetValue(FFT_HANDLE, DFTI_OUTPUT_DISTANCE, vsize+inc)
!    ierr = DftiCommitDescriptor( FFT_HANDLE )
!    ierr = DftiComputeBackward( FFT_HANDLE, TF_vect, vect_out )
!    ierr = DftiFreeDescriptor( FFT_HANDLE )

!    CALL DZFFTM(isign,vsize,nb_vect,scale_factor,vect,vsize+inc,TF_vect,vsize/2+1,table_forward,work,ierr)
  
  END SUBROUTINE fft_forward
  
  SUBROUTINE fft_backward(TF_vect,vect,nb_vect)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: nb_vect
    REAL,INTENT(OUT)    :: vect((vsize+inc)*nb_vect)
    COMPLEX,INTENT(IN ) :: TF_vect((vsize/2+1)*nb_vect)
    REAL                :: work(4*vsize*nb_vect)
    INTEGER             :: ierr
    INTEGER, PARAMETER :: isign=1
    TYPE(DFTI_DESCRIPTOR),POINTER :: FFT_Handle
!    CALL ZDFFTM(isign,vsize,nb_vect,scale_factor,TF_vect,vsize/2+1,vect,vsize+inc,table_backward,work,ierr)
!    IF ( .NOT. Backward_handle(nb_vect)%IsAllocated) THEN
      ierr = DftiCreateDescriptor( FFT_Handle, DFTI_DOUBLE, DFTI_REAL, 1, vsize )
      ierr = DftiSetValue(FFT_Handle,DFTI_NUMBER_OF_TRANSFORMS,nb_vect)
      ierr = DftiSetValue(FFT_Handle,DFTI_FORWARD_SCALE,scale_factor)
      ierr = DftiSetValue(FFT_Handle,DFTI_BACKWARD_SCALE,scale_factor)
      ierr = DftiSetValue(FFT_Handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
      ierr = DftiSetValue(FFT_Handle, DFTI_INPUT_DISTANCE,  (vsize/2+1)*2)
      ierr = DftiSetValue(FFT_Handle, DFTI_OUTPUT_DISTANCE, vsize+inc)
      ierr = DftiCommitDescriptor( FFT_Handle )
!      Backward_handle(nb_vect)%IsAllocated=.TRUE.
!    ENDIF
    ierr = DftiComputeBackward( FFT_Handle, TF_vect, vect )
    ierr = DftiFreeDescriptor( FFT_Handle)
  
  END SUBROUTINE fft_backward

#endif
  
END MODULE mod_fft_mkl
