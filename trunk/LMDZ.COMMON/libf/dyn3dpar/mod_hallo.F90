module mod_Hallo
USE mod_const_mpi, ONLY: COMM_LMDZ,MPI_REAL_LMDZ
USE parallel_lmdz, ONLY: using_mpi, mpi_size, mpi_rank, omp_chunk, omp_rank, &
                         pole_nord, pole_sud, jj_begin, jj_end, &
                         jj_begin_para, jj_end_para
implicit none
  logical,save :: use_mpi_alloc
  integer, parameter :: MaxRequest=200
  integer, parameter :: MaxProc=80
  integer, parameter :: MaxBufferSize=1024*1024*16
  integer, parameter :: ListSize=1000
  
  integer,save       :: MaxBufferSize_Used
!$OMP THREADPRIVATE( MaxBufferSize_Used)

   real,save,pointer,dimension(:) :: Buffer
!$OMP THREADPRIVATE(Buffer)

   integer,save,dimension(Listsize) :: Buffer_Pos
   integer,save :: Index_Pos
!$OMP THREADPRIVATE(Buffer_Pos,Index_pos)
   
  type Hallo
    real, dimension(:,:),pointer :: Field
    integer :: offset
    integer :: size
    integer :: NbLevel
    integer :: Stride
  end type Hallo
  
  type request_SR
    integer :: NbRequest=0
    integer :: Pos
    integer :: Index 
    type(Hallo),dimension(MaxRequest) :: Hallo
    integer :: MSG_Request
  end type request_SR

  type request
    type(request_SR),dimension(0:MaxProc-1) :: RequestSend
    type(request_SR),dimension(0:MaxProc-1) :: RequestRecv
    integer :: tag=1
  end type request
  
    
  contains

  subroutine Init_mod_hallo
    implicit none

    Index_Pos=1
    Buffer_Pos(Index_Pos)=1
    MaxBufferSize_Used=0

    IF (use_mpi_alloc .AND. using_mpi) THEN
      CALL create_global_mpi_buffer
    ELSE 
      CALL create_standard_mpi_buffer
    ENDIF
     
  end subroutine init_mod_hallo

  SUBROUTINE create_standard_mpi_buffer
  IMPLICIT NONE
    
    ALLOCATE(Buffer(MaxBufferSize))
    
  END SUBROUTINE create_standard_mpi_buffer
  
  SUBROUTINE create_global_mpi_buffer
  IMPLICIT NONE
#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif  
    POINTER (Pbuffer,MPI_Buffer(MaxBufferSize))
    REAL :: MPI_Buffer
#ifdef CPP_MPI
    INTEGER(KIND=MPI_ADDRESS_KIND) :: BS 
#else
    INTEGER(KIND=8) :: BS
#endif
    INTEGER :: i,ierr

!  Allocation du buffer MPI
      Bs=8*MaxBufferSize
!$OMP CRITICAL (MPI)
#ifdef CPP_MPI
      CALL MPI_ALLOC_MEM(BS,MPI_INFO_NULL,Pbuffer,ierr)
#endif
!$OMP END CRITICAL (MPI)
      DO i=1,MaxBufferSize
	MPI_Buffer(i)=i
      ENDDO
     
      CALL  Associate_buffer(MPI_Buffer)
      
  CONTAINS
     
     SUBROUTINE Associate_buffer(MPI_Buffer)
     IMPLICIT NONE
       REAL,DIMENSION(:),target :: MPI_Buffer  

         Buffer=>MPI_Buffer
 
      END SUBROUTINE  Associate_buffer
                                      
  END SUBROUTINE create_global_mpi_buffer
 
      
  subroutine allocate_buffer(Size,Index,Pos)
  implicit none
    integer :: Size
    integer :: Index
    integer :: Pos

    if (Buffer_pos(Index_pos)+Size>MaxBufferSize_Used) MaxBufferSize_Used=Buffer_pos(Index_pos)+Size  
    if (Buffer_pos(Index_pos)+Size>MaxBufferSize) then
      print *,'STOP :: La taille de MaxBufferSize dans mod_hallo.F90 est trop petite !!!!'
      stop
    endif
    
    if (Index_pos>=ListSize) then
      print *,'STOP :: La taille de ListSize dans mod_hallo.F90 est trop petite !!!!'
      stop
    endif
     
    Pos=Buffer_Pos(Index_Pos)
    Buffer_Pos(Index_pos+1)=Buffer_Pos(Index_Pos)+Size
    Index_Pos=Index_Pos+1
    Index=Index_Pos
    
  end subroutine allocate_buffer
     
  subroutine deallocate_buffer(Index)
  implicit none
    integer :: Index
    
    Buffer_Pos(Index)=-1
    
    do while (Buffer_Pos(Index_Pos)==-1 .and. Index_Pos>1)
      Index_Pos=Index_Pos-1
    end do

  end subroutine deallocate_buffer  
  
  subroutine SetTag(a_request,tag)
  implicit none
    type(request):: a_request
    integer :: tag
    
    a_request%tag=tag
  end subroutine SetTag
  
  
  subroutine Init_Hallo(Field,Stride,NbLevel,offset,size,NewHallo)
    integer :: Stride
    integer :: NbLevel
    integer :: size
    integer :: offset
    real, dimension(Stride,NbLevel),target :: Field
    type(Hallo) :: NewHallo
    
    NewHallo%Field=>Field
    NewHallo%Stride=Stride
    NewHallo%NbLevel=NbLevel
    NewHallo%size=size
    NewHallo%offset=offset
    
    
  end subroutine Init_Hallo
  
  subroutine Register_SendField(Field,ij,ll,offset,size,target,a_request)
  implicit none

#include "dimensions.h"
#include "paramet.h"    
    
      INTEGER :: ij,ll,offset,size,target
      REAL, dimension(ij,ll) :: Field
      type(request),target :: a_request
      type(request_SR),pointer :: Ptr_request

      Ptr_Request=>a_request%RequestSend(target)
      Ptr_Request%NbRequest=Ptr_Request%NbRequest+1
      if (Ptr_Request%NbRequest>=MaxRequest) then
        print *,'STOP :: La taille de MaxRequest dans mod_hallo.F90 est trop petite !!!!'
        stop
      endif      
      call Init_Hallo(Field,ij,ll,offset,size,Ptr_request%Hallo(Ptr_Request%NbRequest))
      
   end subroutine Register_SendField      
      
  subroutine Register_RecvField(Field,ij,ll,offset,size,target,a_request)
  implicit none

#include "dimensions.h"
#include "paramet.h"    
    
      INTEGER :: ij,ll,offset,size,target
      REAL, dimension(ij,ll) :: Field
      type(request),target :: a_request
      type(request_SR),pointer :: Ptr_request

      Ptr_Request=>a_request%RequestRecv(target)
      Ptr_Request%NbRequest=Ptr_Request%NbRequest+1
      
      if (Ptr_Request%NbRequest>=MaxRequest) then
        print *,'STOP :: La taille de MaxRequest dans mod_hallo.F90 est trop petite !!!!'
        stop
      endif   
            
      call Init_Hallo(Field,ij,ll,offset,size,Ptr_request%Hallo(Ptr_Request%NbRequest))

      
   end subroutine Register_RecvField      
  
  subroutine Register_SwapField(FieldS,FieldR,ij,ll,jj_Nb_New,a_request)
  
      implicit none
#include "dimensions.h"
#include "paramet.h"    
    
    INTEGER :: ij,ll
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    type(request) :: a_request
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    integer ::i,jje,jjb
    
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo
    
    do i=0,MPI_Size-1
      if (i /= MPI_Rank) then
        jjb=max(jj_begin_new(i),jj_begin)
        jje=min(jj_end_new(i),jj_end)
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_SendField(FieldS,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
        jjb=max(jj_begin_new(MPI_Rank),jj_begin_Para(i))
        jje=min(jj_end_new(MPI_Rank),jj_end_Para(i))
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_RecvField(FieldR,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
      endif
    enddo
    
  end subroutine Register_SwapField    
  
  
    subroutine Register_SwapFieldHallo(FieldS,FieldR,ij,ll,jj_Nb_New,Up,Down,a_request)
  
      implicit none
#include "dimensions.h"
#include "paramet.h"    
    
    INTEGER :: ij,ll,Up,Down
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    type(request) :: a_request
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    integer ::i,jje,jjb
    
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo
    
    do i=0,MPI_Size-1
      jj_begin_New(i)=max(1,jj_begin_New(i)-Up)
      jj_end_New(i)=min(jjp1,jj_end_new(i)+Down)
    enddo
   
    do i=0,MPI_Size-1
      if (i /= MPI_Rank) then
        jjb=max(jj_begin_new(i),jj_begin)
        jje=min(jj_end_new(i),jj_end)
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_SendField(FieldS,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
        jjb=max(jj_begin_new(MPI_Rank),jj_begin_Para(i))
        jje=min(jj_end_new(MPI_Rank),jj_end_Para(i))
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_RecvField(FieldR,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
      endif
    enddo
    
  end subroutine Register_SwapFieldHallo
  
  subroutine Register_Hallo(Field,ij,ll,RUp,Rdown,SUp,SDown,a_request)
  
      implicit none
#include "dimensions.h"
#include "paramet.h"    
#ifdef CPP_MPI
    include 'mpif.h'
#endif    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: Sup,Sdown,rup,rdown
      type(request) :: a_request
      type(Hallo),pointer :: PtrHallo
      LOGICAL :: SendUp,SendDown
      LOGICAL :: RecvUp,RecvDown
   
 
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
      
      if (Sup.eq.0) then
        SendUp=.FALSE.
       endif
      
      if (Sdown.eq.0) then
        SendDown=.FALSE.
      endif

      if (Rup.eq.0) then
        RecvUp=.FALSE.
      endif
      
      if (Rdown.eq.0) then
        RecvDown=.FALSE.
      endif
      
      IF (SendUp) THEN
        call Register_SendField(Field,ij,ll,jj_begin,SUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (SendDown) THEN
        call Register_SendField(Field,ij,ll,jj_end-SDown+1,SDown,MPI_Rank+1,a_request)
      ENDIF
    
  
      IF (RecvUp) THEN
        call Register_RecvField(Field,ij,ll,jj_begin-Rup,RUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (RecvDown) THEN
        call Register_RecvField(Field,ij,ll,jj_end+1,RDown,MPI_Rank+1,a_request)
      ENDIF
  
    end subroutine Register_Hallo
    
    subroutine SendRequest(a_Request)
      implicit none

#include "dimensions.h"
#include "paramet.h"
#ifdef CPP_MPI
      include 'mpif.h'
#endif

      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer :: SizeBuffer
      integer :: i,rank,l,ij,Pos,ierr
      integer :: offset
      real,dimension(:,:),pointer :: Field
      integer :: Nb
       
      do rank=0,MPI_SIZE-1
      
        Req=>a_Request%RequestSend(rank)
        
        SizeBuffer=0
        do i=1,Req%NbRequest
          PtrHallo=>Req%Hallo(i)
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
          DO l=1,PtrHallo%NbLevel
            SizeBuffer=SizeBuffer+PtrHallo%size*iip1
          ENDDO
!$OMP ENDDO NOWAIT          
        enddo
      
        if (SizeBuffer>0) then
       
          call allocate_buffer(SizeBuffer,Req%Index,Req%pos)

          Pos=Req%Pos
          do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            offset=(PtrHallo%offset-1)*iip1+1
            Nb=iip1*PtrHallo%size-1
            Field=>PtrHallo%Field

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)           
            do l=1,PtrHallo%NbLevel
!cdir NODEP
              do ij=0,Nb
	        Buffer(Pos+ij)=Field(Offset+ij,l)
	      enddo
              
              Pos=Pos+Nb+1
            enddo
!$OMP END DO NOWAIT            
          enddo
    
!$OMP CRITICAL (MPI)
         
#ifdef CPP_MPI
         call MPI_ISSEND(Buffer(req%Pos),SizeBuffer,MPI_REAL_LMDZ,rank,a_request%tag+1000*omp_rank,     &
                         COMM_LMDZ,Req%MSG_Request,ierr)
#endif
         IF (.NOT.using_mpi) THEN
           PRINT *,'Erreur, echange MPI en mode sequentiel !!!'
           STOP
         ENDIF
!         PRINT *,"-------------------------------------------------------------------"
!         PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->"
!         PRINT *,"Requete envoye au proc :",rank,"tag :",a_request%tag+1000*omp_rank
!         PRINT *,"Taille du message :",SizeBuffer,"requete no :",Req%MSG_Request
!         PRINT *,"-------------------------------------------------------------------"
!$OMP END CRITICAL (MPI)
        endif

    enddo
   
           
      do rank=0,MPI_SIZE-1
         
          Req=>a_Request%RequestRecv(rank)
          SizeBuffer=0
          
	  do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
            DO l=1,PtrHallo%NbLevel
              SizeBuffer=SizeBuffer+PtrHallo%size*iip1
            ENDDO
!$OMP ENDDO NOWAIT          
          enddo
        
          if (SizeBuffer>0) then

             call allocate_buffer(SizeBuffer,Req%Index,Req%Pos)
!$OMP CRITICAL (MPI)

#ifdef CPP_MPI
             call MPI_IRECV(Buffer(Req%Pos),SizeBuffer,MPI_REAL_LMDZ,rank,a_request%tag+1000*omp_rank,     &
                           COMM_LMDZ,Req%MSG_Request,ierr)
#endif             
             IF (.NOT.using_mpi) THEN
               PRINT *,'Erreur, echange MPI en mode sequentiel !!!'
               STOP
             ENDIF

!         PRINT *,"-------------------------------------------------------------------"
!         PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->"
!         PRINT *,"Requete en attente du proc :",rank,"tag :",a_request%tag+1000*omp_rank
!         PRINT *,"Taille du message :",SizeBuffer,"requete no :",Req%MSG_Request
!         PRINT *,"-------------------------------------------------------------------"

!$OMP END CRITICAL (MPI)
          endif
      
      enddo
                        
   end subroutine SendRequest 
   
   subroutine WaitRequest(a_Request)
   implicit none
   
#include "dimensions.h"
#include "paramet.h"
#ifdef CPP_MPI
      include 'mpif.h'   
#endif
      
      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer, dimension(2*mpi_size) :: TabRequest
#ifdef CPP_MPI
      integer, dimension(MPI_STATUS_SIZE,2*mpi_size) :: TabStatus
#else
      integer, dimension(1,2*mpi_size) :: TabStatus
#endif
      integer :: NbRequest
      integer :: i,rank,pos,ij,l,ierr
      integer :: offset
      integer :: Nb
      
      
      NbRequest=0
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
     
      if (NbRequest>0) then
!$OMP CRITICAL (MPI)
!        PRINT *,"-------------------------------------------------------------------"
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"en attente"
!        PRINT *,"No des requetes :",TabRequest(1:NbRequest)
#ifdef CPP_MPI
        call MPI_WAITALL(NbRequest,TabRequest,TabStatus,ierr)
#endif
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"complete"
!        PRINT *,"-------------------------------------------------------------------"
!$OMP END CRITICAL (MPI)
      endif
      do rank=0,MPI_Size-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          Pos=Req%Pos
          do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            offset=(PtrHallo%offset-1)*iip1+1
	    Nb=iip1*PtrHallo%size-1

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
	    do l=1,PtrHallo%NbLevel
!cdir NODEP
              do ij=0,Nb
	        PtrHallo%Field(offset+ij,l)=Buffer(Pos+ij)
	      enddo

              Pos=Pos+Nb+1
	    enddo
!$OMP ENDDO NOWAIT	    
          enddo
        endif
      enddo
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
              
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
     
      a_request%tag=1
    end subroutine WaitRequest
     
   subroutine WaitSendRequest(a_Request)
   implicit none
   
#include "dimensions.h"
#include "paramet.h"
#ifdef CPP_MPI
      include 'mpif.h'   
#endif      
      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer, dimension(mpi_size) :: TabRequest
#ifdef CPP_MPI
      integer, dimension(MPI_STATUS_SIZE,mpi_size) :: TabStatus
#else
      integer, dimension(1,mpi_size) :: TabStatus
#endif
      integer :: NbRequest
      integer :: i,rank,pos,ij,l,ierr
      integer :: offset
      
      
      NbRequest=0
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
      

      if (NbRequest>0) THEN 
!$OMP CRITICAL (MPI)     
!        PRINT *,"-------------------------------------------------------------------"
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"en attente"
!        PRINT *,"No des requetes :",TabRequest(1:NbRequest)
#ifdef CPP_MPI
        call MPI_WAITALL(NbRequest,TabRequest,TabStatus,ierr)
#endif
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"complete"
!        PRINT *,"-------------------------------------------------------------------"

!$OMP END CRITICAL (MPI)
      endif      
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
              
      a_request%tag=1
    end subroutine WaitSendRequest
    
   subroutine WaitRecvRequest(a_Request)
   implicit none
   
#include "dimensions.h"
#include "paramet.h"
#ifdef CPP_MPI
      include 'mpif.h'   
#endif
      
      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer, dimension(mpi_size) :: TabRequest
#ifdef CPP_MPI
      integer, dimension(MPI_STATUS_SIZE,mpi_size) :: TabStatus
#else
      integer, dimension(1,mpi_size) :: TabStatus
#endif
      integer :: NbRequest
      integer :: i,rank,pos,ij,l,ierr
      integer :: offset,Nb
      
      
      NbRequest=0
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
     
      
      if (NbRequest>0) then
!$OMP CRITICAL (MPI)     
!        PRINT *,"-------------------------------------------------------------------"
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"en attente"
!        PRINT *,"No des requetes :",TabRequest(1:NbRequest)
#ifdef CPP_MPI
        call MPI_WAITALL(NbRequest,TabRequest,TabStatus,ierr)
#endif
!        PRINT *,"Process de rang",mpi_rank,"Task : ",omp_rank,"--->",NbRequest,"complete"
!        PRINT *,"-------------------------------------------------------------------"
!$OMP END CRITICAL (MPI)     
      endif
      
      do rank=0,MPI_Size-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          Pos=Req%Pos
          do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            offset=(PtrHallo%offset-1)*iip1+1
	    Nb=iip1*PtrHallo%size-1
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
	    do l=1,PtrHallo%NbLevel
!cdir NODEP
              do ij=0,Nb
	        PtrHallo%Field(offset+ij,l)=Buffer(Pos+ij)
	      enddo
                 Pos=Pos+Nb+1
            enddo
!$OMP END DO NOWAIT
          enddo
        endif
      enddo
      
           
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
     
      a_request%tag=1
    end subroutine WaitRecvRequest
    
    
    
    subroutine CopyField(FieldS,FieldR,ij,ll,jj_Nb_New)
  
      implicit none
#include "dimensions.h"
#include "paramet.h"    
    
    INTEGER :: ij,ll,l
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    integer ::i,jje,jjb,ijb,ije
    
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo
    
    jjb=max(jj_begin,jj_begin_new(MPI_Rank))
    jje=min(jj_end,jj_end_new(MPI_Rank))
    if (ij==ip1jm) jje=min(jje,jjm)

    if (jje >= jjb) then
      ijb=(jjb-1)*iip1+1
      ije=jje*iip1

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      do l=1,ll
        FieldR(ijb:ije,l)=FieldS(ijb:ije,l)
      enddo
!$OMP ENDDO NOWAIT
    endif


  end subroutine CopyField    

  subroutine CopyFieldHallo(FieldS,FieldR,ij,ll,jj_Nb_New,Up,Down)
  
      implicit none
#include "dimensions.h"
#include "paramet.h"    
    
    INTEGER :: ij,ll,Up,Down
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New

    integer ::i,jje,jjb,ijb,ije,l

     
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo

        
    jjb=max(jj_begin,jj_begin_new(MPI_Rank)-Up)
    jje=min(jj_end,jj_end_new(MPI_Rank)+Down)
    if (ij==ip1jm) jje=min(jje,jjm)
    
    
    if (jje >= jjb) then
      ijb=(jjb-1)*iip1+1
      ije=jje*iip1

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      do l=1,ll
        FieldR(ijb:ije,l)=FieldS(ijb:ije,l)
      enddo
!$OMP ENDDO NOWAIT

    endif
   end subroutine CopyFieldHallo        
          
end module mod_Hallo 
   
