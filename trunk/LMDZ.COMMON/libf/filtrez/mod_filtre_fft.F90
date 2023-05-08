!
! $Id: mod_filtre_fft.F90 1403 2010-07-01 09:02:53Z fairhead $
!

MODULE mod_filtre_fft

  LOGICAL,SAVE :: use_filtre_fft
  REAL,SAVE,ALLOCATABLE :: Filtre_u(:,:)
  REAL,SAVE,ALLOCATABLE :: Filtre_v(:,:)
  REAL,SAVE,ALLOCATABLE :: Filtre_inv(:,:)

CONTAINS
  
  SUBROUTINE Init_filtre_fft(coeffu,modfrstu,jfiltnu,jfiltsu,coeffv,modfrstv,jfiltnv,jfiltsv)
    USE mod_fft
    IMPLICIT NONE
    include 'dimensions.h'
    REAL,   INTENT(IN) :: coeffu(iim,jjm)
    INTEGER,INTENT(IN) :: modfrstu(jjm)
    INTEGER,INTENT(IN) :: jfiltnu
    INTEGER,INTENT(IN) :: jfiltsu
    REAL,   INTENT(IN) :: coeffv(iim,jjm)
    INTEGER,INTENT(IN) :: modfrstv(jjm)
    INTEGER,INTENT(IN) :: jfiltnv
    INTEGER,INTENT(IN) :: jfiltsv
    
    INTEGER            :: index_vp(iim)
    INTEGER            :: i,j
    INTEGER            :: l,ll_nb

    index_vp(1)=1
    DO i=1,iim/2
      index_vp(i+1)=i*2
    ENDDO
    
    DO i=1,iim/2-1
      index_vp(iim/2+i+1)=iim-2*i+1
    ENDDO
    
    ALLOCATE(Filtre_u(iim,jjm))
    ALLOCATE(Filtre_v(iim,jjm))
    ALLOCATE(Filtre_inv(iim,jjm))
  
    
    DO j=2,jfiltnu
      DO i=1,iim
        IF (index_vp(i) < modfrstu(j)) THEN
          Filtre_u(i,j)=0
        ELSE
          Filtre_u(i,j)=coeffu(index_vp(i),j)
        ENDIF
      ENDDO
    ENDDO
    
    DO j=jfiltsu,jjm
      DO i=1,iim
        IF (index_vp(i) < modfrstu(j)) THEN
          Filtre_u(i,j)=0
        ELSE
          Filtre_u(i,j)=coeffu(index_vp(i),j)
        ENDIF
      ENDDO
    ENDDO
 
    DO j=1,jfiltnv
      DO i=1,iim
        IF (index_vp(i) < modfrstv(j)) THEN
          Filtre_v(i,j)=0
        ELSE
          Filtre_v(i,j)=coeffv(index_vp(i),j)
        ENDIF
      ENDDO
    ENDDO
   
    DO j=jfiltsv,jjm
      DO i=1,iim
        IF (index_vp(i) < modfrstv(j)) THEN
          Filtre_v(i,j)=0
        ELSE
          Filtre_v(i,j)=coeffv(index_vp(i),j)
        ENDIF
      ENDDO
    ENDDO
         
    DO j=2,jfiltnu
      DO i=1,iim
        IF (index_vp(i) < modfrstu(j)) THEN
          Filtre_inv(i,j)=0
        ELSE
          Filtre_inv(i,j)=coeffu(index_vp(i),j)/(1.+coeffu(index_vp(i),j))
        ENDIF
      ENDDO
    ENDDO

    DO j=jfiltsu,jjm
      DO i=1,iim
        IF (index_vp(i) < modfrstu(j)) THEN
          Filtre_inv(i,j)=0
        ELSE
          Filtre_inv(i,j)=coeffu(index_vp(i),j)/(1.+coeffu(index_vp(i),j))
        ENDIF
      ENDDO
    ENDDO
    
#ifdef FFT_FFTW

    WRITE (*,*)"COTH jfiltnu,jfiltsu,jfiltnv,jjm-jfiltsv"
    WRITE (*,*)jfiltnu,jfiltsu,jfiltnv,jjm-jfiltsv
    WRITE (*,*)MAX(jfiltnu-2,jjm-jfiltsu,jfiltnv-2,jjm-jfiltsv)+1
    CALL Init_FFT(iim,(llm+1)*(MAX(jfiltnu-2,jjm-jfiltsu,jfiltnv-2,jjm-jfiltsv)+1))
#else    
    CALL Init_FFT(iim,(jjm+1)*(llm+1))
#endif        
    
  END SUBROUTINE Init_filtre_fft
  
  SUBROUTINE Filtre_u_fft(vect_inout,nlat,jj_begin,jj_end,nbniv)
    USE mod_fft
#ifdef CPP_PARA
    USE parallel_lmdz,ONLY : OMP_CHUNK
#endif
    IMPLICIT NONE
    include 'dimensions.h'
    INTEGER,INTENT(IN) :: nlat
    INTEGER,INTENT(IN) :: jj_begin
    INTEGER,INTENT(IN) :: jj_end
    INTEGER,INTENT(IN) :: nbniv
    REAL,INTENT(INOUT) :: vect_inout(iim+1,nlat,nbniv)

    REAL               :: vect(iim+inc,jj_end-jj_begin+1,nbniv)
    COMPLEX         :: TF_vect(iim/2+1,jj_end-jj_begin+1,nbniv)
    INTEGER            :: nb_vect
    INTEGER :: i,j,l
    INTEGER :: ll_nb
    
    ll_nb=0
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO l=1,nbniv
      ll_nb=ll_nb+1
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim+1
          vect(i,j,ll_nb)=vect_inout(i,j+jj_begin-1,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

    nb_vect=(jj_end-jj_begin+1)*ll_nb

    CALL FFT_forward(vect,TF_vect,nb_vect)

    DO l=1,ll_nb
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim/2+1
          TF_vect(i,j,l)=TF_vect(i,j,l)*Filtre_u(i,jj_begin+j-1)
        ENDDO
      ENDDO
    ENDDO
  
    CALL FFT_backward(TF_vect,vect,nb_vect)
      
      
    ll_nb=0
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO l=1,nbniv
      ll_nb=ll_nb+1
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim+1
          vect_inout(i,j+jj_begin-1,l)=vect(i,j,ll_nb)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

  END SUBROUTINE Filtre_u_fft
  

  SUBROUTINE Filtre_v_fft(vect_inout,nlat,jj_begin,jj_end,nbniv)
    USE mod_fft
#ifdef CPP_PARA
    USE parallel_lmdz,ONLY : OMP_CHUNK
#endif
    IMPLICIT NONE
    INCLUDE 'dimensions.h'
    INTEGER,INTENT(IN) :: nlat
    INTEGER,INTENT(IN) :: jj_begin
    INTEGER,INTENT(IN) :: jj_end
    INTEGER,INTENT(IN) :: nbniv
    REAL,INTENT(INOUT) :: vect_inout(iim+1,nlat,nbniv)

    REAL               :: vect(iim+inc,jj_end-jj_begin+1,nbniv)
    COMPLEX            :: TF_vect(iim/2+1,jj_end-jj_begin+1,nbniv)
    INTEGER            :: nb_vect
    INTEGER :: i,j,l
    INTEGER :: ll_nb
    
    ll_nb=0
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO l=1,nbniv
      ll_nb=ll_nb+1
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim+1
          vect(i,j,ll_nb)=vect_inout(i,j+jj_begin-1,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

    
    nb_vect=(jj_end-jj_begin+1)*ll_nb

    CALL FFT_forward(vect,TF_vect,nb_vect)
  
    DO l=1,ll_nb
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim/2+1
          TF_vect(i,j,l)=TF_vect(i,j,l)*Filtre_v(i,jj_begin+j-1)
        ENDDO
      ENDDO
    ENDDO
  
    CALL FFT_backward(TF_vect,vect,nb_vect)
    
    
    ll_nb=0
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO l=1,nbniv
      ll_nb=ll_nb+1
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim+1
          vect_inout(i,j+jj_begin-1,l)=vect(i,j,ll_nb)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  
  END SUBROUTINE Filtre_v_fft


  SUBROUTINE Filtre_inv_fft(vect_inout,nlat,jj_begin,jj_end,nbniv)
    USE mod_fft
#ifdef CPP_PARA
    USE parallel_lmdz,ONLY : OMP_CHUNK
#endif
    IMPLICIT NONE
    INCLUDE 'dimensions.h'
    INTEGER,INTENT(IN) :: nlat
    INTEGER,INTENT(IN) :: jj_begin
    INTEGER,INTENT(IN) :: jj_end
    INTEGER,INTENT(IN) :: nbniv
    REAL,INTENT(INOUT) :: vect_inout(iim+1,nlat,nbniv)

     REAL               :: vect(iim+inc,jj_end-jj_begin+1,nbniv)
    COMPLEX            :: TF_vect(iim/2+1,jj_end-jj_begin+1,nbniv)
    INTEGER            :: nb_vect
    INTEGER :: i,j,l
    INTEGER :: ll_nb
    
    ll_nb=0
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO l=1,nbniv
      ll_nb=ll_nb+1
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim+1
          vect(i,j,ll_nb)=vect_inout(i,j+jj_begin-1,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

    nb_vect=(jj_end-jj_begin+1)*ll_nb

    CALL FFT_forward(vect,TF_vect,nb_vect)
  
    DO l=1,ll_nb
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim/2+1
          TF_vect(i,j,l)=TF_vect(i,j,l)*Filtre_inv(i,jj_begin+j-1)
        ENDDO
      ENDDO
    ENDDO
  
    CALL FFT_backward(TF_vect,vect,nb_vect)

    ll_nb=0
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO l=1,nbniv
      ll_nb=ll_nb+1
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim+1
          vect_inout(i,j+jj_begin-1,l)=vect(i,j,ll_nb)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

  END SUBROUTINE Filtre_inv_fft  
   
END MODULE mod_filtre_fft
 
