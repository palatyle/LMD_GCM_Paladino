MODULE mod_fft_mathkeisan
#ifdef FFT_MATHKEISAN

  REAL,SAVE,ALLOCATABLE    :: Table_forward(:)
  REAL,SAVE,ALLOCATABLE    :: Table_backward(:)
  REAL,SAVE                :: scale_factor
  INTEGER,SAVE             :: vsize
  INTEGER,PARAMETER        :: inc=2

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
    ALLOCATE(Table_forward(2*vsize+64))
    ALLOCATE(Table_backward(2*vsize+64))
    
    CALL DZFFTM(isign,vsize,itmp,scale_factor,rtmp,vsize+inc,ctmp,vsize/2+1,table_forward,rtmp,ierr)
    
    CALL ZDFFTM(isign,vsize,itmp,scale_factor,ctmp,vsize/2+1,rtmp,vsize+inc,table_backward,rtmp,ierr)

    
  END SUBROUTINE Init_fft
  
  
  SUBROUTINE fft_forward(vect,TF_vect,nb_vect)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: nb_vect
    REAL,INTENT(IN)     :: vect(vsize+inc,nb_vect)
    COMPLEX,INTENT(OUT) :: TF_vect(vsize/2+1,nb_vect)
    REAL                :: work(4*vsize*nb_vect)
    INTEGER             :: ierr
    INTEGER, PARAMETER :: isign=-1
    
    work=0
    CALL DZFFTM(isign,vsize,nb_vect,scale_factor,vect,vsize+inc,TF_vect,vsize/2+1,table_forward,work,ierr)
  
  END SUBROUTINE fft_forward
  
  SUBROUTINE fft_backward(TF_vect,vect,nb_vect)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: nb_vect
    REAL,INTENT(OUT)    :: vect(vsize+inc,nb_vect)
    COMPLEX,INTENT(IN ) :: TF_vect(vsize/2+1,nb_vect)
    REAL                :: work(4*vsize*nb_vect)
    INTEGER             :: ierr
    INTEGER, PARAMETER :: isign=1
    
    work(:)=0
    CALL ZDFFTM(isign,vsize,nb_vect,scale_factor,TF_vect,vsize/2+1,vect,vsize+inc,table_backward,work,ierr)
  
  END SUBROUTINE fft_backward

#endif
  
END MODULE mod_fft_mathkeisan


