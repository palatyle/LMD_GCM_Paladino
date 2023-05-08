!
! $Id: mod_fft_fftw.F90 1403 2010-07-01 09:02:53Z fairhead $
!

MODULE mod_fft_fftw

#ifdef FFT_FFTW

  REAL, SAVE                   :: scale_factor
  INTEGER, SAVE                :: vsize
  INTEGER, PARAMETER           :: inc=1
  
  INTEGER*8, ALLOCATABLE, DIMENSION(:), SAVE :: plan_forward
  INTEGER*8, ALLOCATABLE, DIMENSION(:), SAVE :: plan_backward
  
CONTAINS
  
  SUBROUTINE Init_fft(iim,nvectmax)
  IMPLICIT NONE
#include <fftw3.f>
    INTEGER :: iim
    INTEGER :: nvectmax

    INTEGER :: itmp

    INTEGER               :: rank
    INTEGER               :: howmany
    INTEGER               :: istride, idist
    INTEGER               :: ostride, odist
    INTEGER, DIMENSION(1) :: n_array, inembed, onembed

    REAL,    DIMENSION(iim+1,nvectmax) :: dbidon
    COMPLEX, DIMENSION(iim/2+1,nvectmax) :: cbidon

    vsize = iim
    scale_factor = 1./SQRT(1.*vsize)

    dbidon = 0
    cbidon = 0

    ALLOCATE(plan_forward(nvectmax))
    ALLOCATE(plan_backward(nvectmax))
    
    WRITE(*,*)"!---------------------!"
    WRITE(*,*)"!                     !"
    WRITE(*,*)"! INITIALISATION FFTW !"
    WRITE(*,*)"!                     !"
    WRITE(*,*)"!---------------------!"
    
! On initialise tous les plans 
    DO itmp = 1, nvectmax
       rank       = 1
       n_array(1) = iim
       howmany    = itmp
       inembed(1) = iim + 1 ; onembed(1) = iim/2 + 1
       istride    = 1       ; ostride    = 1
       idist      = iim + 1 ; odist      = iim/2 + 1

       CALL dfftw_plan_many_dft_r2c(plan_forward(itmp), rank, n_array, howmany, &
            & dbidon, inembed, istride, idist, &
            & cbidon, onembed, ostride, odist, &
            & FFTW_ESTIMATE)

       rank       = 1
       n_array(1) = iim
       howmany    = itmp
       inembed(1) = iim/2 + 1 ; onembed(1) = iim + 1
       istride    = 1         ; ostride    = 1
       idist      = iim/2 + 1 ; odist      = iim + 1
       CALL dfftw_plan_many_dft_c2r(plan_backward(itmp), rank, n_array, howmany, &
            & cbidon, inembed, istride, idist, &
            & dbidon, onembed, ostride, odist, &
            & FFTW_ESTIMATE)

    ENDDO

    WRITE(*,*)"!-------------------------!"
    WRITE(*,*)"!                         !"
    WRITE(*,*)"! FIN INITIALISATION FFTW !"
    WRITE(*,*)"!                         !"
    WRITE(*,*)"!-------------------------!"

  END SUBROUTINE Init_fft
  
  
  SUBROUTINE fft_forward(vect,TF_vect,nb_vect)
    IMPLICIT NONE
#include <fftw3.f>
    INTEGER,INTENT(IN)     :: nb_vect
    REAL,INTENT(IN)        :: vect(vsize+inc,nb_vect)
    COMPLEX,INTENT(OUT) :: TF_vect(vsize/2+1,nb_vect)

    CALL dfftw_execute_dft_r2c(plan_forward(nb_vect),vect,TF_vect)

    TF_vect = scale_factor * TF_vect

  END SUBROUTINE fft_forward
  
  SUBROUTINE fft_backward(TF_vect,vect,nb_vect)
    IMPLICIT NONE
#include <fftw3.f>
    INTEGER,INTENT(IN)     :: nb_vect
    REAL,INTENT(OUT)       :: vect(vsize+inc,nb_vect)
    COMPLEX,INTENT(IN ) :: TF_vect(vsize/2+1,nb_vect)

    CALL dfftw_execute_dft_c2r(plan_backward(nb_vect),TF_vect,vect)

    vect = scale_factor * vect

  END SUBROUTINE fft_backward

#endif
  
END MODULE mod_fft_fftw
