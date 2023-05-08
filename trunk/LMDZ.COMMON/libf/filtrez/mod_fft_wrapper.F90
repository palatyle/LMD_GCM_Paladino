MODULE mod_fft_wrapper

  INTEGER,SAVE             :: vsize
  INTEGER,PARAMETER        :: inc=1

CONTAINS
  
  SUBROUTINE Init_fft(iim,nb)
  IMPLICIT NONE
    INTEGER :: iim
    INTEGER :: nb
    
    STOP "wrapper fft : une FFT doit etre specifiee a l'aide d'une clee CPP, sinon utiliser le filtre classique"
  END SUBROUTINE Init_fft
  
  
  SUBROUTINE fft_forward(vect,TF_vect,nb_vect)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: nb_vect
    REAL,INTENT(IN)     :: vect(vsize+inc,nb_vect)
    COMPLEX,INTENT(INOUT) :: TF_vect(vsize/2+1,nb_vect)
    
    STOP "wrapper fft : une FFT doit etre specifiee a l'aide d'une clee CPP, sinon utiliser le filtre classique"
    
  END SUBROUTINE fft_forward
  
  SUBROUTINE fft_backward(TF_vect,vect,nb_vect)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: nb_vect
    REAL,INTENT(INOUT)    :: vect(vsize+inc,nb_vect)
    COMPLEX,INTENT(IN ) :: TF_vect(vsize/2+1,nb_vect)
  
    STOP "wrapper fft : une FFT doit etre specifiee a l'aide d'une clee CPP, sinon utiliser le filtre classique"
    
  END SUBROUTINE fft_backward
  
END MODULE mod_fft_wrapper
