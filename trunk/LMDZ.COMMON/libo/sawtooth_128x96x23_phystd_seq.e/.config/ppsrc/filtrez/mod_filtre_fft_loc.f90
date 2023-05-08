










MODULE mod_filtre_fft_loc

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
    
    
!    CALL Init_FFT(iim,(jjm+1)*(llm+1))
        
    
  END SUBROUTINE Init_filtre_fft
  
  SUBROUTINE Filtre_u_fft(vect_inout,jjb,jje,jj_begin,jj_end,nbniv)
    USE mod_fft
    IMPLICIT NONE
    include 'dimensions.h'
    INTEGER,INTENT(IN) :: jjb
    INTEGER,INTENT(IN) :: jje
    INTEGER,INTENT(IN) :: jj_begin
    INTEGER,INTENT(IN) :: jj_end
    INTEGER,INTENT(IN) :: nbniv
    REAL,INTENT(INOUT) :: vect_inout(iim+1,jjb:jje,nbniv)

    REAL               :: vect(iim+inc,jj_end-jj_begin+1,nbniv)
    COMPLEX            :: TF_vect(iim/2+1,jj_end-jj_begin+1,nbniv)
    INTEGER            :: nb_vect
    INTEGER :: i,j,l
    INTEGER :: ll_nb
!    REAL               :: vect_tmp(iim+inc,jj_end-jj_begin+1,nbniv)
    
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

!    vect_tmp=vect

    CALL FFT_forward(vect,TF_vect,nb_vect)

!    CALL FFT_forward(vect,TF_vect_test,nb_vect)
!      PRINT *,"XXXXXXXXXXXXX Filtre_u_FFT xxxxxxxxxxxx"
!      DO j=1,jj_end-jj_begin+1
!      DO i=1,iim/2+1
!         PRINT *,"====",i,j,"----->",TF_vect_test(i,j,1)
!       ENDDO
!      ENDDO

    DO l=1,ll_nb
      DO j=1,jj_end-jj_begin+1
        DO i=1,iim/2+1
          TF_vect(i,j,l)=TF_vect(i,j,l)*Filtre_u(i,jj_begin+j-1)
        ENDDO
      ENDDO
    ENDDO
       
    CALL FFT_backward(TF_vect,vect,nb_vect)
!    CALL FFT_backward(TF_vect_test,vect_test,nb_vect)
          
!      PRINT *,"XXXXXXXXXXXXX Filtre_u_FFT xxxxxxxxxxxx"
!      DO j=1,jj_end-jj_begin+1
!         DO i=1,iim
!           PRINT *,"====",i,j,"----->",vect_test(i,j,1)
!         ENDDO
!      ENDDO

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
  

  SUBROUTINE Filtre_v_fft(vect_inout,jjb,jje,jj_begin,jj_end,nbniv)
    USE mod_fft
    IMPLICIT NONE
    INCLUDE 'dimensions.h'
    INTEGER,INTENT(IN) :: jjb
    INTEGER,INTENT(IN) :: jje
    INTEGER,INTENT(IN) :: jj_begin
    INTEGER,INTENT(IN) :: jj_end
    INTEGER,INTENT(IN) :: nbniv
    REAL,INTENT(INOUT) :: vect_inout(iim+1,jjb:jje,nbniv)

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


  SUBROUTINE Filtre_inv_fft(vect_inout,jjb,jje,jj_begin,jj_end,nbniv)
    USE mod_fft
    IMPLICIT NONE
    INCLUDE 'dimensions.h'
    INTEGER,INTENT(IN) :: jjb
    INTEGER,INTENT(IN) :: jje
    INTEGER,INTENT(IN) :: jj_begin
    INTEGER,INTENT(IN) :: jj_end
    INTEGER,INTENT(IN) :: nbniv
    REAL,INTENT(INOUT) :: vect_inout(iim+1,jjb:jje,nbniv)

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
  
  
!  SUBROUTINE get_ll_index(nbniv,ll_index,ll_nb)
!  IMPLICIT NONE
!    INTEGER,INTENT(IN)  :: nbniv
!    INTEGER,INTENT(OUT) :: ll_index(nbniv)
!    INTEGER,INTENT(OUT) :: ll_nb
!
!    INTEGER :: l,ll_begin, ll_end
!   INTEGER :: omp_rank,omp_size
!   INTEGER :: OMP_GET_NUM_THREADS
!   INTEGER :: omp_chunk
!   EXTERNAL OMP_GET_NUM_THREADS
!   INTEGER :: OMP_GET_THREAD_NUM
!   EXTERNAL OMP_GET_THREAD_NUM
!
!   
!   omp_size=OMP_GET_NUM_THREADS()
!   omp_rank=OMP_GET_THREAD_NUM()    
!   omp_chunk=nbniv/omp_size+min(1,MOD(nbniv,omp_size))
!   
!   ll_begin=omp_rank*OMP_CHUNK+1
!   ll_nb=0
!   DO WHILE (ll_begin<=nbniv)
!     ll_end=min(ll_begin+OMP_CHUNK-1,nbniv)
!     DO l=ll_begin,ll_end
!       ll_nb=ll_nb+1
!       ll_index(ll_nb)=l
!     ENDDO
!     ll_begin=ll_begin+omp_size*OMP_CHUNK
!   ENDDO
!  
!  END SUBROUTINE get_ll_index
   
END MODULE mod_filtre_fft_loc
 
