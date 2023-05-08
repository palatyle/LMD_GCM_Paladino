! 
! $Id: parallel.F90 1810 2013-07-24 08:06:39Z emillour $
!
  MODULE parallel_lmdz
  USE mod_const_mpi
#ifdef CPP_IOIPSL
      use IOIPSL, only: getin
#else
! if not using IOIPSL, we still need to use (a local version of) getin
      use ioipsl_getincom, only: getin
#endif    
    
    LOGICAL,SAVE :: using_mpi=.TRUE.
    LOGICAL,SAVE :: using_omp
    
    integer, save :: mpi_size
    integer, save :: mpi_rank
    integer, save :: jj_begin
    integer, save :: jj_end
    integer, save :: jj_nb
    integer, save :: ij_begin
    integer, save :: ij_end
    logical, save :: pole_nord
    logical, save :: pole_sud
    
    integer, allocatable, save, dimension(:) :: jj_begin_para
    integer, allocatable, save, dimension(:) :: jj_end_para
    integer, allocatable, save, dimension(:) :: jj_nb_para
    integer, save :: OMP_CHUNK
    integer, save :: omp_rank
    integer, save :: omp_size  
!$OMP THREADPRIVATE(omp_rank)

! Ehouarn: add "dummy variables" (which are in dyn3d_mem/parallel_lmdz.F90)
! so that calfis_loc compiles even if using dyn3dpar
    integer,save  :: jjb_u
    integer,save  :: jje_u
    integer,save  :: jjnb_u
    integer,save  :: jjb_v
    integer,save  :: jje_v
    integer,save  :: jjnb_v    

    integer,save  :: ijb_u
    integer,save  :: ije_u
    integer,save  :: ijnb_u    
    
    integer,save  :: ijb_v
    integer,save  :: ije_v
    integer,save  :: ijnb_v    

 contains
 
    subroutine init_parallel
    USE vampir
    implicit none
#ifdef CPP_MPI
      include 'mpif.h'
#endif
#include "dimensions.h"
#include "paramet.h"
#include "iniprint.h"

      integer :: ierr
      integer :: i,j
      integer :: type_size
      integer, dimension(3) :: blocklen,type
      integer :: comp_id
      character(len=4)  :: num
      character(len=20) :: filename
 
#ifdef CPP_OMP    
      INTEGER :: OMP_GET_NUM_THREADS
      EXTERNAL OMP_GET_NUM_THREADS
      INTEGER :: OMP_GET_THREAD_NUM
      EXTERNAL OMP_GET_THREAD_NUM
#endif  

#ifdef CPP_MPI
       using_mpi=.TRUE.
#else
       using_mpi=.FALSE.
#endif
      

#ifdef CPP_OMP
       using_OMP=.TRUE.
#else
       using_OMP=.FALSE.
#endif
      
      call InitVampir
      
      IF (using_mpi) THEN
#ifdef CPP_MPI
        call MPI_COMM_SIZE(COMM_LMDZ,mpi_size,ierr)
        call MPI_COMM_RANK(COMM_LMDZ,mpi_rank,ierr)
#endif
      ELSE
        mpi_size=1
        mpi_rank=0
      ENDIF


! Open text output file with mpi_rank in suffix of file name 
      IF (lunout /= 5 .and. lunout /= 6) THEN
         WRITE(num,'(I4.4)') mpi_rank
         filename='lmdz.out_'//num
         IF (mpi_rank .NE. 0) THEN
            OPEN(UNIT=lunout,FILE=TRIM(filename),ACTION='write', &
               STATUS='unknown',FORM='formatted',IOSTAT=ierr) 
         ENDIF
      ENDIF

      
      allocate(jj_begin_para(0:mpi_size-1))
      allocate(jj_end_para(0:mpi_size-1))
      allocate(jj_nb_para(0:mpi_size-1))
      
      do i=0,mpi_size-1
        jj_nb_para(i)=(jjm+1)/mpi_size
        if ( i < MOD((jjm+1),mpi_size) ) jj_nb_para(i)=jj_nb_para(i)+1
        
        if (jj_nb_para(i) <= 1 ) then
          
         write(lunout,*)"Arret : le nombre de bande de lattitude par process est trop faible (<2)."
         write(lunout,*)" ---> diminuez le nombre de CPU ou augmentez la taille en lattitude"
          
#ifdef CPP_MPI
          IF (using_mpi) call MPI_ABORT(COMM_LMDZ,-1, ierr)
#endif          
        endif
        
      enddo
      
!      jj_nb_para(0)=11
!      jj_nb_para(1)=25
!      jj_nb_para(2)=25
!      jj_nb_para(3)=12      

      j=1
      
      do i=0,mpi_size-1 
        
        jj_begin_para(i)=j
        jj_end_para(i)=j+jj_Nb_para(i)-1
        j=j+jj_Nb_para(i)
      
      enddo
      
      jj_begin = jj_begin_para(mpi_rank)
      jj_end   = jj_end_para(mpi_rank)
      jj_nb    = jj_nb_para(mpi_rank)
      
      ij_begin=(jj_begin-1)*iip1+1
      ij_end=jj_end*iip1
      
      if (mpi_rank.eq.0) then
        pole_nord=.TRUE.
      else 
        pole_nord=.FALSE.
      endif
      
      if (mpi_rank.eq.mpi_size-1) then
        pole_sud=.TRUE.
      else 
        pole_sud=.FALSE.
      endif
        
      write(lunout,*)"init_parallel: jj_begin",jj_begin
      write(lunout,*)"init_parallel: jj_end",jj_end
      write(lunout,*)"init_parallel: ij_begin",ij_begin
      write(lunout,*)"init_parallel: ij_end",ij_end

!$OMP PARALLEL

#ifdef CPP_OMP
!$OMP MASTER
        omp_size=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP BARRIER
        omp_rank=OMP_GET_THREAD_NUM()    

!Config  Key  = omp_chunk
!Config  Desc = taille des blocs openmp
!Config  Def  = 1
!Config  Help = defini la taille des packets d'it�ration openmp
!Config         distribue a chaque tache lors de l'entree dans une
!Config         boucle parallelisee

!$OMP MASTER
      omp_chunk=(llm+1)/omp_size
      IF (MOD(llm+1,omp_size)/=0) omp_chunk=omp_chunk+1
      CALL getin('omp_chunk',omp_chunk)
!$OMP END MASTER
!$OMP BARRIER        
#else    
        omp_size=1
        omp_rank=0
#endif
!$OMP END PARALLEL         
    
    end subroutine init_parallel

    
    subroutine SetDistrib(jj_Nb_New)
    implicit none

#include "dimensions.h"
#include "paramet.h"

      INTEGER,dimension(0:MPI_Size-1) :: jj_Nb_New
      INTEGER :: i  
  
      jj_Nb_Para=jj_Nb_New
      
      jj_begin_para(0)=1
      jj_end_para(0)=jj_Nb_Para(0)
      
      do i=1,mpi_size-1 
        
        jj_begin_para(i)=jj_end_para(i-1)+1
        jj_end_para(i)=jj_begin_para(i)+jj_Nb_para(i)-1
      
      enddo
      
      jj_begin = jj_begin_para(mpi_rank)
      jj_end   = jj_end_para(mpi_rank)
      jj_nb    = jj_nb_para(mpi_rank)
      
      ij_begin=(jj_begin-1)*iip1+1
      ij_end=jj_end*iip1

    end subroutine SetDistrib



    
    subroutine Finalize_parallel
#ifdef CPP_XIOS
    ! ug Pour les sorties XIOS
	USE wxios
#endif
#ifdef CPP_COUPLE
! Use of Oasis-MCT coupler 
#if defined CPP_OMCT
    use mod_prism
#else
    use mod_prism_proto
#endif
! Ehouarn: surface_data module is in 'phylmd' ...
      use surface_data, only : type_ocean
      implicit none
#else
      implicit none
! without the surface_data module, we declare (and set) a dummy 'type_ocean'
      character(len=6),parameter :: type_ocean="dummy"
#endif
! #endif of #ifdef CPP_EARTH

      include "dimensions.h"
      include "paramet.h"
#ifdef CPP_MPI
      include 'mpif.h'
#endif      

      integer :: ierr
      integer :: i

      if (allocated(jj_begin_para)) deallocate(jj_begin_para)
      if (allocated(jj_end_para))   deallocate(jj_end_para)
      if (allocated(jj_nb_para))    deallocate(jj_nb_para)

      if (type_ocean == 'couple') then
#ifdef CPP_XIOS
    !Fermeture propre de XIOS
      CALL wxios_close()
#else
#ifdef CPP_COUPLE
         call prism_terminate_proto(ierr)
         IF (ierr .ne. PRISM_Ok) THEN
            call abort_gcm('Finalize_parallel',' Probleme dans prism_terminate_proto ',1)
         endif
#endif 
#endif
      else
#ifdef CPP_XIOS
    !Fermeture propre de XIOS
      CALL wxios_close()
#endif
#ifdef CPP_MPI
         IF (using_mpi) call MPI_FINALIZE(ierr)
#endif
      end if
      
    end subroutine Finalize_parallel
        
    subroutine Pack_Data(Field,ij,ll,row,Buffer)
    implicit none

#include "dimensions.h"
#include "paramet.h"

      integer, intent(in) :: ij,ll,row
      real,dimension(ij,ll),intent(in) ::Field
      real,dimension(ll*iip1*row), intent(out) :: Buffer 
            
      integer :: Pos
      integer :: i,l
      
      Pos=0
      do l=1,ll
        do i=1,row*iip1
          Pos=Pos+1
          Buffer(Pos)=Field(i,l)
        enddo
      enddo
      
    end subroutine Pack_data 
     
    subroutine Unpack_Data(Field,ij,ll,row,Buffer)
    implicit none

#include "dimensions.h"
#include "paramet.h"

      integer, intent(in) :: ij,ll,row
      real,dimension(ij,ll),intent(out) ::Field
      real,dimension(ll*iip1*row), intent(in) :: Buffer 
            
      integer :: Pos
      integer :: i,l
      
      Pos=0
      
      do l=1,ll
        do i=1,row*iip1
          Pos=Pos+1
          Field(i,l)=Buffer(Pos)
        enddo
      enddo
      
    end subroutine UnPack_data

    
    SUBROUTINE barrier
    IMPLICIT NONE
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif
    INTEGER :: ierr
    
!$OMP CRITICAL (MPI)      
#ifdef CPP_MPI
      IF (using_mpi) CALL MPI_Barrier(COMM_LMDZ,ierr)
#endif
!$OMP END CRITICAL (MPI)
    
    END SUBROUTINE barrier
       
      
    subroutine exchange_hallo(Field,ij,ll,up,down)
    USE Vampir
    implicit none
#include "dimensions.h"
#include "paramet.h"    
#ifdef CPP_MPI
    include 'mpif.h'
#endif    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: up,down
      
      INTEGER :: ierr
      LOGICAL :: SendUp,SendDown
      LOGICAL :: RecvUp,RecvDown
      INTEGER, DIMENSION(4) :: Request
#ifdef CPP_MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE,4) :: Status
#else
      INTEGER, DIMENSION(1,4) :: Status
#endif
      INTEGER :: NbRequest
      REAL, dimension(:),allocatable :: Buffer_Send_up,Buffer_Send_down
      REAL, dimension(:),allocatable :: Buffer_Recv_up,Buffer_Recv_down
      INTEGER :: Buffer_size      

      IF (using_mpi) THEN

        CALL barrier
      
        call VTb(VThallo)
      
        SendUp=.TRUE.
        SendDown=.TRUE.
        RecvUp=.TRUE.
        RecvDown=.TRUE.
          
        IF (pole_nord) THEN
          SendUp=.FALSE.
          RecvUp=.FALSE.
        ENDIF
    
        IF (pole_sud) THEN
          SendDown=.FALSE.
          RecvDown=.FALSE.
        ENDIF
        
        if (up.eq.0) then
          SendDown=.FALSE.
          RecvUp=.FALSE.
        endif
      
        if (down.eq.0) then
          SendUp=.FALSE.
          RecvDown=.FALSE.
        endif
      
        NbRequest=0
  
        IF (SendUp) THEN
          NbRequest=NbRequest+1
          buffer_size=down*iip1*ll
          allocate(Buffer_Send_up(Buffer_size))
          call PACK_Data(Field(ij_begin,1),ij,ll,down,Buffer_Send_up)
!$OMP CRITICAL (MPI)
#ifdef CPP_MPI
          call MPI_ISSEND(Buffer_send_up,Buffer_Size,MPI_REAL8,MPI_Rank-1,1,     &
                          COMM_LMDZ,Request(NbRequest),ierr)
#endif
!$OMP END CRITICAL (MPI)
        ENDIF
  
        IF (SendDown) THEN
          NbRequest=NbRequest+1
           
          buffer_size=up*iip1*ll
          allocate(Buffer_Send_down(Buffer_size))
          call PACK_Data(Field(ij_end+1-up*iip1,1),ij,ll,up,Buffer_send_down)
        
!$OMP CRITICAL (MPI)
#ifdef CPP_MPI
          call MPI_ISSEND(Buffer_send_down,Buffer_Size,MPI_REAL8,MPI_Rank+1,1,     &
                          COMM_LMDZ,Request(NbRequest),ierr)
#endif
!$OMP END CRITICAL (MPI)
        ENDIF
    
  
        IF (RecvUp) THEN
          NbRequest=NbRequest+1
          buffer_size=up*iip1*ll
          allocate(Buffer_recv_up(Buffer_size))
              
!$OMP CRITICAL (MPI)
#ifdef CPP_MPI
          call MPI_IRECV(Buffer_recv_up,Buffer_size,MPI_REAL8,MPI_Rank-1,1,  &
                          COMM_LMDZ,Request(NbRequest),ierr)
#endif
!$OMP END CRITICAL (MPI)
     
       
        ENDIF
  
        IF (RecvDown) THEN
          NbRequest=NbRequest+1
          buffer_size=down*iip1*ll
          allocate(Buffer_recv_down(Buffer_size))
        
!$OMP CRITICAL (MPI)
#ifdef CPP_MPI
          call MPI_IRECV(Buffer_recv_down,Buffer_size,MPI_REAL8,MPI_Rank+1,1,     &
                          COMM_LMDZ,Request(NbRequest),ierr)
#endif
!$OMP END CRITICAL (MPI)
        
        ENDIF
  
#ifdef CPP_MPI
        if (NbRequest > 0) call MPI_WAITALL(NbRequest,Request,Status,ierr)
#endif
        IF (RecvUp)  call Unpack_Data(Field(ij_begin-up*iip1,1),ij,ll,up,Buffer_Recv_up)
        IF (RecvDown) call Unpack_Data(Field(ij_end+1,1),ij,ll,down,Buffer_Recv_down)  

        call VTe(VThallo)
        call barrier
      
      ENDIF  ! using_mpi
      
      RETURN
      
    end subroutine exchange_Hallo
    

    subroutine Gather_Field(Field,ij,ll,rank)
    implicit none
#include "dimensions.h"
#include "paramet.h" 
#include "iniprint.h"
#ifdef CPP_MPI
    include 'mpif.h'
#endif    
      INTEGER :: ij,ll,rank
      REAL, dimension(ij,ll) :: Field
      REAL, dimension(:),allocatable :: Buffer_send   
      REAL, dimension(:),allocatable :: Buffer_Recv
      INTEGER, dimension(0:MPI_Size-1) :: Recv_count, displ
      INTEGER :: ierr
      INTEGER ::i
      
      IF (using_mpi) THEN

        if (ij==ip1jmp1) then 
           allocate(Buffer_send(iip1*ll*(jj_end-jj_begin+1)))
           call Pack_Data(Field(ij_begin,1),ij,ll,jj_end-jj_begin+1,Buffer_send)
        else if (ij==ip1jm) then
           allocate(Buffer_send(iip1*ll*(min(jj_end,jjm)-jj_begin+1)))
           call Pack_Data(Field(ij_begin,1),ij,ll,min(jj_end,jjm)-jj_begin+1,Buffer_send)
        else
           write(lunout,*)ij  
        stop 'erreur dans Gather_Field'
        endif
        
        if (MPI_Rank==rank) then
          allocate(Buffer_Recv(ij*ll))

!CDIR NOVECTOR
          do i=0,MPI_Size-1
             
            if (ij==ip1jmp1) then 
              Recv_count(i)=(jj_end_para(i)-jj_begin_para(i)+1)*ll*iip1
            else if (ij==ip1jm) then
              Recv_count(i)=(min(jj_end_para(i),jjm)-jj_begin_para(i)+1)*ll*iip1
            else
              stop 'erreur dans Gather_Field'
            endif
                   
            if (i==0) then 
              displ(i)=0 
            else
              displ(i)=displ(i-1)+Recv_count(i-1)
            endif
            
          enddo
          
        else
          ! Ehouarn: When in debug mode, ifort complains (for call MPI_GATHERV
          !          below) about Buffer_Recv() being not allocated.
          !          So make a dummy allocation.
          allocate(Buffer_Recv(1))
        endif ! of if (MPI_Rank==rank)
  
!$OMP CRITICAL (MPI)
#ifdef CPP_MPI
        call MPI_GATHERV(Buffer_send,(min(ij_end,ij)-ij_begin+1)*ll,MPI_REAL8,   &
                          Buffer_Recv,Recv_count,displ,MPI_REAL8,rank,COMM_LMDZ,ierr)
#endif
!$OMP END CRITICAL (MPI)
      
        if (MPI_Rank==rank) then                  
      
          if (ij==ip1jmp1) then 
            do i=0,MPI_Size-1
              call Unpack_Data(Field((jj_begin_para(i)-1)*iip1+1,1),ij,ll,                 &
                               jj_end_para(i)-jj_begin_para(i)+1,Buffer_Recv(displ(i)+1))
            enddo
          else if (ij==ip1jm) then
            do i=0,MPI_Size-1
               call Unpack_Data(Field((jj_begin_para(i)-1)*iip1+1,1),ij,ll,                       &
                               min(jj_end_para(i),jjm)-jj_begin_para(i)+1,Buffer_Recv(displ(i)+1))
            enddo
          endif
        endif 
      ENDIF ! using_mpi
      
    end subroutine Gather_Field


    subroutine AllGather_Field(Field,ij,ll)
    implicit none
#include "dimensions.h"
#include "paramet.h"    
#ifdef CPP_MPI
    include 'mpif.h'
#endif    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: ierr
      
      IF (using_mpi) THEN
        call Gather_Field(Field,ij,ll,0)
!$OMP CRITICAL (MPI)
#ifdef CPP_MPI
      call MPI_BCAST(Field,ij*ll,MPI_REAL8,0,COMM_LMDZ,ierr)
#endif
!$OMP END CRITICAL (MPI)
      ENDIF
      
    end subroutine AllGather_Field
    
   subroutine Broadcast_Field(Field,ij,ll,rank)
    implicit none
#include "dimensions.h"
#include "paramet.h"    
#ifdef CPP_MPI
    include 'mpif.h'
#endif    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: rank
      INTEGER :: ierr
      
      IF (using_mpi) THEN
      
!$OMP CRITICAL (MPI)
#ifdef CPP_MPI
      call MPI_BCAST(Field,ij*ll,MPI_REAL8,rank,COMM_LMDZ,ierr)
#endif
!$OMP END CRITICAL (MPI)
      
      ENDIF
    end subroutine Broadcast_Field
        
   
!  Subroutine verif_hallo(Field,ij,ll,up,down)
!    implicit none
!#include "dimensions.h"
!#include "paramet.h"    
!    include 'mpif.h'
!    
!      INTEGER :: ij,ll
!      REAL, dimension(ij,ll) :: Field
!      INTEGER :: up,down 
!      
!      REAL,dimension(ij,ll): NewField
!      
!      NewField=0
!      
!      ijb=ij_begin
!      ije=ij_end
!      if (pole_nord) 
!      NewField(ij_be       

  end MODULE parallel_lmdz
