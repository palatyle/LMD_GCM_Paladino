!
! $Id: bands.F90 2351 2015-08-25 15:14:59Z emillour $
!
  module Bands
  
    integer, parameter :: bands_caldyn=1
    integer, parameter :: bands_vanleer=2
    integer, parameter :: bands_dissip=3
    
    INTEGER,dimension(:),allocatable :: jj_Nb_Caldyn
    INTEGER,dimension(:),allocatable :: jj_Nb_vanleer
    INTEGER,dimension(:),allocatable :: jj_Nb_vanleer2
    INTEGER,dimension(:),allocatable :: jj_Nb_dissip
    INTEGER,dimension(:),allocatable :: jj_Nb_physic
    INTEGER,dimension(:),allocatable :: jj_Nb_physic_bis
    INTEGER,dimension(:),allocatable :: distrib_phys
  
  contains
  
  subroutine AllocateBands
    USE parallel_lmdz
    implicit none
    
    allocate(jj_Nb_Caldyn(0:MPI_Size-1))
    allocate(jj_Nb_vanleer(0:MPI_Size-1))
    allocate(jj_Nb_vanleer2(0:MPI_Size-1))
    allocate(jj_Nb_dissip(0:MPI_Size-1))
    allocate(jj_Nb_physic(0:MPI_Size-1))
    allocate(jj_Nb_physic_bis(0:MPI_Size-1))
    allocate(distrib_phys(0:MPI_Size-1))
  
  end subroutine AllocateBands
  
  subroutine Read_distrib
    USE parallel_lmdz
    implicit none

    include "dimensions.h"
      integer :: i,j
      character (len=4) :: siim,sjjm,sllm,sproc
      character (len=255) :: filename
      integer :: unit_number=10
      integer :: ierr
    
      call AllocateBands
      write(siim,'(i3)') iim
      write(sjjm,'(i3)') jjm
      write(sllm,'(i3)') llm
      write(sproc,'(i3)') mpi_size
      filename='Bands_'//TRIM(ADJUSTL(siim))//'x'//TRIM(ADJUSTL(sjjm))//'x'//TRIM(ADJUSTL(sllm))//'_'  &
                        //TRIM(ADJUSTL(sproc))//'prc.dat'    
       
      OPEN(UNIT=unit_number,FILE=trim(filename),STATUS='old',FORM='formatted',IOSTAT=ierr)
      
      if (ierr==0) then
      
         do i=0,mpi_size-1
          read (unit_number,*) j,jj_nb_caldyn(i)
        enddo
      
        do i=0,mpi_size-1
          read (unit_number,*) j,jj_nb_vanleer(i)
        enddo
      
        do i=0,mpi_size-1
          read (unit_number,*) j,jj_nb_dissip(i)
        enddo
      
        do i=0,mpi_size-1
          read (unit_number,*) j,distrib_phys(i)
        enddo
	
	CLOSE(unit_number)  
  
      else
        do i=0,mpi_size-1
          jj_nb_caldyn(i)=(jjm+1)/mpi_size
	  if (i<MOD(jjm+1,mpi_size)) jj_nb_caldyn(i)=jj_nb_caldyn(i)+1
        enddo
      
        jj_nb_vanleer(:)=jj_nb_caldyn(:)
        jj_nb_dissip(:)=jj_nb_caldyn(:)
        
	do i=0,mpi_size-1
	  distrib_phys(i)=(iim*(jjm-1)+2)/mpi_size
	  IF (i<MOD(iim*(jjm-1)+2,mpi_size)) distrib_phys(i)=distrib_phys(i)+1
	enddo
      endif
  
   end subroutine Read_distrib
   
   
   SUBROUTINE  Set_Bands 
     USE parallel_lmdz
     IMPLICIT NONE
     INCLUDE 'dimensions.h'    
     INTEGER :: i, ij
     INTEGER :: jj_para_begin(0:mpi_size-1)
     INTEGER :: jj_para_end(0:mpi_size-1)
       
      do i=0,mpi_size-1
         jj_nb_vanleer2(i)=(jjm+1)/mpi_size
	 if (i<MOD(jjm+1,mpi_size)) jj_nb_vanleer2(i)=jj_nb_vanleer2(i)+1
      enddo
          
      jj_para_begin(0)=1
      ij=distrib_phys(0)+iim-1
      jj_para_end(0)=((ij-1)/iim)+1
      
      DO i=1,mpi_Size-1
        ij=ij+1
        jj_para_begin(i)=((ij-1)/iim)+1
        ij=ij+distrib_phys(i)-1
        jj_para_end(i)=((ij-1)/iim)+1
      ENDDO
 
      do i=0,MPI_Size-1
        jj_Nb_physic(i)=jj_para_end(i)-jj_para_begin(i)+1
        if (i/=0) then
          if (jj_para_begin(i)==jj_para_end(i-1)) then
            jj_Nb_physic(i-1)=jj_Nb_physic(i-1)-1
          endif
        endif
      enddo
      
      do i=0,MPI_Size-1
        jj_Nb_physic_bis(i)=jj_para_end(i)-jj_para_begin(i)+1
        if (i/=0) then
          if (jj_para_begin(i)==jj_para_end(i-1)) then
            jj_Nb_physic_bis(i)=jj_Nb_physic_bis(i)-1
          else
	    jj_Nb_physic_bis(i-1)=jj_Nb_physic_bis(i-1)+1
	    jj_Nb_physic_bis(i)=jj_Nb_physic_bis(i)-1
	  endif
        endif
      enddo
      
    end subroutine Set_Bands


    subroutine AdjustBands_caldyn
      use times
      USE parallel_lmdz
      implicit none

      real :: minvalue,maxvalue
      integer :: min_proc,max_proc
      integer :: i,j
      real,allocatable,dimension(:) :: value
      integer,allocatable,dimension(:) :: index
      real :: tmpvalue
      integer :: tmpindex
      
      allocate(value(0:mpi_size-1))
      allocate(index(0:mpi_size-1))
        
  
      call allgather_timer_average

      do i=0,mpi_size-1
        value(i)=timer_average(jj_nb_caldyn(i),timer_caldyn,i)
	index(i)=i
      enddo
      
      do i=0,mpi_size-2
        do j=i+1,mpi_size-1
	  if (value(i)>value(j)) then
	    tmpvalue=value(i)
	    value(i)=value(j)
	    value(j)=tmpvalue
	    
	    tmpindex=index(i)
	    index(i)=index(j)
	    index(j)=tmpindex
	   endif
	 enddo
      enddo
      
      maxvalue=value(mpi_size-1)
      max_proc=index(mpi_size-1)           
           
      do i=0,mpi_size-2
        minvalue=value(i)
        min_proc=index(i)
        if (jj_nb_caldyn(max_proc)>2) then
          if (timer_iteration(jj_nb_caldyn(min_proc)+1,timer_caldyn,min_proc)<=1 ) then
             jj_nb_caldyn(min_proc)=jj_nb_caldyn(min_proc)+1
             jj_nb_caldyn(max_proc)=jj_nb_caldyn(max_proc)-1
	     exit
           else
             if (timer_average(jj_nb_caldyn(min_proc)+1,timer_caldyn,min_proc)                 &
	        -timer_delta(jj_nb_caldyn(min_proc)+1,timer_caldyn,min_proc) < maxvalue) then
               jj_nb_caldyn(min_proc)=jj_nb_caldyn(min_proc)+1
               jj_nb_caldyn(max_proc)=jj_nb_caldyn(max_proc)-1
               exit
	     endif
           endif
         endif
      enddo
      
      deallocate(value)
      deallocate(index)
         
    end subroutine AdjustBands_caldyn
    
    subroutine AdjustBands_vanleer
      use times
      USE parallel_lmdz
      implicit none

      real :: minvalue,maxvalue
      integer :: min_proc,max_proc
      integer :: i,j
      real,allocatable,dimension(:) :: value
      integer,allocatable,dimension(:) :: index
      real :: tmpvalue
      integer :: tmpindex
      
      allocate(value(0:mpi_size-1))
      allocate(index(0:mpi_size-1))
        
  
      call allgather_timer_average

      do i=0,mpi_size-1
        value(i)=timer_average(jj_nb_vanleer(i),timer_vanleer,i)
	index(i)=i
      enddo
      
      do i=0,mpi_size-2
        do j=i+1,mpi_size-1
	  if (value(i)>value(j)) then
	    tmpvalue=value(i)
	    value(i)=value(j)
	    value(j)=tmpvalue
	    
	    tmpindex=index(i)
	    index(i)=index(j)
	    index(j)=tmpindex
	   endif
	 enddo
      enddo
      
      maxvalue=value(mpi_size-1)
      max_proc=index(mpi_size-1)           
           
      do i=0,mpi_size-2
        minvalue=value(i)
        min_proc=index(i)

        if (jj_nb_vanleer(max_proc)>2) then
          if (timer_average(jj_nb_vanleer(min_proc)+1,timer_vanleer,min_proc)==0. .or. &
             timer_average(jj_nb_vanleer(max_proc)-1,timer_vanleer,max_proc)==0.) then
             jj_nb_vanleer(min_proc)=jj_nb_vanleer(min_proc)+1
             jj_nb_vanleer(max_proc)=jj_nb_vanleer(max_proc)-1
	     exit
           else
             if (timer_average(jj_nb_vanleer(min_proc)+1,timer_vanleer,min_proc) < maxvalue) then
               jj_nb_vanleer(min_proc)=jj_nb_vanleer(min_proc)+1
               jj_nb_vanleer(max_proc)=jj_nb_vanleer(max_proc)-1
               exit
	     endif
           endif
         endif
      enddo
      
      deallocate(value)
      deallocate(index)
         
    end subroutine AdjustBands_vanleer

    subroutine AdjustBands_dissip
      use times
      USE parallel_lmdz
      implicit none

      real :: minvalue,maxvalue
      integer :: min_proc,max_proc
      integer :: i,j
      real,allocatable,dimension(:) :: value
      integer,allocatable,dimension(:) :: index
      real :: tmpvalue
      integer :: tmpindex
      
      allocate(value(0:mpi_size-1))
      allocate(index(0:mpi_size-1))
        
  
      call allgather_timer_average

      do i=0,mpi_size-1
        value(i)=timer_average(jj_nb_dissip(i),timer_dissip,i)
	index(i)=i
      enddo
      
      do i=0,mpi_size-2
        do j=i+1,mpi_size-1
	  if (value(i)>value(j)) then
	    tmpvalue=value(i)
	    value(i)=value(j)
	    value(j)=tmpvalue
	    
	    tmpindex=index(i)
	    index(i)=index(j)
	    index(j)=tmpindex
	   endif
	 enddo
      enddo
      
      maxvalue=value(mpi_size-1)
      max_proc=index(mpi_size-1)           
           
      do i=0,mpi_size-2
        minvalue=value(i)
        min_proc=index(i)

        if (jj_nb_dissip(max_proc)>3) then
          if (timer_iteration(jj_nb_dissip(min_proc)+1,timer_dissip,min_proc)<=1) then
             jj_nb_dissip(min_proc)=jj_nb_dissip(min_proc)+1
             jj_nb_dissip(max_proc)=jj_nb_dissip(max_proc)-1
	     exit
           else
             if (timer_average(jj_nb_dissip(min_proc)+1,timer_dissip,min_proc)         &
	        - timer_delta(jj_nb_dissip(min_proc)+1,timer_dissip,min_proc) < maxvalue) then
               jj_nb_dissip(min_proc)=jj_nb_dissip(min_proc)+1
               jj_nb_dissip(max_proc)=jj_nb_dissip(max_proc)-1
               exit
	     endif
           endif
         endif
      enddo
      
      deallocate(value)
      deallocate(index)
         
    end subroutine AdjustBands_dissip

    subroutine AdjustBands_physic
      use times
#ifdef CPP_PHYS
! Ehouarn: what follows is only related to // physics
      USE mod_phys_lmdz_para, only : klon_mpi_para_nb
#endif
      USE parallel_lmdz
      implicit none

      integer :: i,Index
      real,allocatable,dimension(:) :: value
      integer,allocatable,dimension(:) :: Inc
      real :: medium
      integer :: NbTot,sgn
      
      allocate(value(0:mpi_size-1))
      allocate(Inc(0:mpi_size-1))
        
  
      call allgather_timer_average
      
      medium=0
      do i=0,mpi_size-1
        value(i)=timer_average(jj_nb_physic(i),timer_physic,i)
	medium=medium+value(i)
      enddo    
      
      medium=medium/mpi_size      
      NbTot=0
#ifdef CPP_PHYS
      do i=0,mpi_size-1
        Inc(i)=nint(klon_mpi_para_nb(i)*(medium-value(i))/value(i))
        NbTot=NbTot+Inc(i)  
      enddo
      
      if (NbTot>=0) then
        Sgn=1
      else
        Sgn=-1
	NbTot=-NbTot
      endif
      
      Index=0
      do i=1,NbTot
        Inc(Index)=Inc(Index)-Sgn
	Index=Index+1
	if (Index>mpi_size-1) Index=0
      enddo
      
      do i=0,mpi_size-1
        distrib_phys(i)=klon_mpi_para_nb(i)+inc(i)
      enddo
#endif     
    end subroutine AdjustBands_physic

    subroutine WriteBands
    USE parallel_lmdz
    implicit none
    include "dimensions.h"

      integer :: i,j
      character (len=4) :: siim,sjjm,sllm,sproc
      character (len=255) :: filename
      integer :: unit_number=10
      integer :: ierr
  
      write(siim,'(i3)') iim
      write(sjjm,'(i3)') jjm
      write(sllm,'(i3)') llm
      write(sproc,'(i3)') mpi_size

      filename='Bands_'//TRIM(ADJUSTL(siim))//'x'//TRIM(ADJUSTL(sjjm))//'x'//TRIM(ADJUSTL(sllm))//'_'  &
                        //TRIM(ADJUSTL(sproc))//'prc.dat'    
      
      OPEN(UNIT=unit_number,FILE=trim(filename),STATUS='replace',FORM='formatted',IOSTAT=ierr)
      
      if (ierr==0) then
        
!	write (unit_number,*) '*** Bandes caldyn ***'
	do i=0,mpi_size-1
          write (unit_number,*) i,jj_nb_caldyn(i)
        enddo
        
!	write (unit_number,*) '*** Bandes vanleer ***' 
        do i=0,mpi_size-1
          write (unit_number,*) i,jj_nb_vanleer(i)
        enddo
       
!        write (unit_number,*) '*** Bandes dissip ***'
        do i=0,mpi_size-1
          write (unit_number,*) i,jj_nb_dissip(i)
        enddo
        
	do i=0,mpi_size-1
          write (unit_number,*) i,distrib_phys(i)
        enddo
	
        CLOSE(unit_number)   
      else 
        print *,'probleme lors de l ecriture des bandes'
      endif
       
    end subroutine WriteBands
  
  end module Bands
  
  
