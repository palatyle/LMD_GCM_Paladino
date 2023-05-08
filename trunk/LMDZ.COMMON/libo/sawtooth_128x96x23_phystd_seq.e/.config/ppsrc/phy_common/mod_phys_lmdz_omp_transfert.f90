










!
!$Header$
!
MODULE mod_phys_lmdz_omp_transfert

  PRIVATE
  
  INTEGER,PARAMETER :: grow_factor=1.5
  INTEGER,PARAMETER :: size_min=1024
  
  CHARACTER(LEN=size_min),SAVE            :: buffer_c
!  INTEGER,SAVE                            :: size_c=0
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)   :: buffer_i
  INTEGER,SAVE                            :: size_i=0
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)      :: buffer_r
  INTEGER,SAVE                            :: size_r=0
  LOGICAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: buffer_l
  INTEGER,SAVE                            :: size_l=0


  
  
  INTERFACE bcast_omp
    MODULE PROCEDURE bcast_omp_c,                                                     &
                     bcast_omp_i,bcast_omp_i1,bcast_omp_i2,bcast_omp_i3,bcast_omp_i4, &
                     bcast_omp_r,bcast_omp_r1,bcast_omp_r2,bcast_omp_r3,bcast_omp_r4, &
                     bcast_omp_l,bcast_omp_l1,bcast_omp_l2,bcast_omp_l3,bcast_omp_l4
  END INTERFACE

  INTERFACE scatter_omp
    MODULE PROCEDURE scatter_omp_i,scatter_omp_i1,scatter_omp_i2,scatter_omp_i3, &
                     scatter_omp_r,scatter_omp_r1,scatter_omp_r2,scatter_omp_r3, &
                     scatter_omp_l,scatter_omp_l1,scatter_omp_l2,scatter_omp_l3
  END INTERFACE

  
  INTERFACE gather_omp
    MODULE PROCEDURE gather_omp_i,gather_omp_i1,gather_omp_i2,gather_omp_i3, &
                     gather_omp_r,gather_omp_r1,gather_omp_r2,gather_omp_r3, &
                     gather_omp_l,gather_omp_l1,gather_omp_l2,gather_omp_l3  
  END INTERFACE
  
  
  INTERFACE reduce_sum_omp
    MODULE PROCEDURE reduce_sum_omp_i,reduce_sum_omp_i1,reduce_sum_omp_i2,reduce_sum_omp_i3,reduce_sum_omp_i4, &
                     reduce_sum_omp_r,reduce_sum_omp_r1,reduce_sum_omp_r2,reduce_sum_omp_r3,reduce_sum_omp_r4
  END INTERFACE 


  PUBLIC bcast_omp,scatter_omp,gather_omp,reduce_sum_omp, omp_barrier

CONTAINS

  SUBROUTINE omp_barrier
  IMPLICIT NONE

!$OMP BARRIER

  END SUBROUTINE omp_barrier
  
  SUBROUTINE check_buffer_i(buff_size)
  IMPLICIT NONE
  INTEGER :: buff_size

!$OMP BARRIER
!$OMP MASTER
    IF (buff_size>size_i) THEN
      IF (ALLOCATED(buffer_i)) DEALLOCATE(buffer_i)
      size_i=MAX(size_min,INT(grow_factor*buff_size))
      ALLOCATE(buffer_i(size_i))
    ENDIF
!$OMP END MASTER
!$OMP BARRIER
  
  END SUBROUTINE check_buffer_i
  
  SUBROUTINE check_buffer_r(buff_size)
  IMPLICIT NONE
  INTEGER :: buff_size

!$OMP BARRIER
!$OMP MASTER
    IF (buff_size>size_r) THEN
      IF (ALLOCATED(buffer_r)) DEALLOCATE(buffer_r)
      size_r=MAX(size_min,INT(grow_factor*buff_size))
      ALLOCATE(buffer_r(size_r))
    ENDIF
!$OMP END MASTER
!$OMP BARRIER
  
  END SUBROUTINE check_buffer_r
  
  SUBROUTINE check_buffer_l(buff_size)
  IMPLICIT NONE
  INTEGER :: buff_size

!$OMP BARRIER
!$OMP MASTER
    IF (buff_size>size_l) THEN
      IF (ALLOCATED(buffer_l)) DEALLOCATE(buffer_l)
      size_l=MAX(size_min,INT(grow_factor*buff_size))
      ALLOCATE(buffer_l(size_l))
    ENDIF
!$OMP END MASTER
!$OMP BARRIER
  
  END SUBROUTINE check_buffer_l
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Broadcast --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! -- Les chaine de charactère -- !!

  SUBROUTINE bcast_omp_c(var)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(INOUT) :: Var
    
    CALL bcast_omp_cgen(Var,len(Var),buffer_c)
    
  END SUBROUTINE bcast_omp_c

!! -- Les entiers -- !!
  
  SUBROUTINE bcast_omp_i(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var
    INTEGER :: Var_tmp(1)
    
    Var_tmp(1)=Var
    CALL check_buffer_i(1)
    CALL bcast_omp_igen(Var_tmp,1,buffer_i)
    Var=Var_tmp(1)

  END SUBROUTINE bcast_omp_i


  SUBROUTINE bcast_omp_i1(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:)
   
    CALL check_buffer_i(size(Var))
    CALL bcast_omp_igen(Var,size(Var),buffer_i)

  END SUBROUTINE bcast_omp_i1


  SUBROUTINE bcast_omp_i2(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:)
   
    CALL check_buffer_i(size(Var))
    CALL bcast_omp_igen(Var,size(Var),buffer_i)

  END SUBROUTINE bcast_omp_i2


  SUBROUTINE bcast_omp_i3(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:)

    CALL check_buffer_i(size(Var))
    CALL bcast_omp_igen(Var,size(Var),buffer_i)

  END SUBROUTINE bcast_omp_i3


  SUBROUTINE bcast_omp_i4(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:,:)
   
    CALL check_buffer_i(size(Var))
    CALL bcast_omp_igen(Var,size(Var),buffer_i)

  END SUBROUTINE bcast_omp_i4


!! -- Les reels -- !!

  SUBROUTINE bcast_omp_r(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var
    REAL :: Var_tmp(1)
    
    Var_tmp(1)=Var
    CALL check_buffer_r(1)
    CALL bcast_omp_rgen(Var_tmp,1,buffer_r)
    Var=Var_tmp(1)

  END SUBROUTINE bcast_omp_r


  SUBROUTINE bcast_omp_r1(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:)
   
    CALL check_buffer_r(size(Var))
    CALL bcast_omp_rgen(Var,size(Var),buffer_r)

  END SUBROUTINE bcast_omp_r1


  SUBROUTINE bcast_omp_r2(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:)
   
    CALL check_buffer_r(size(Var))
    CALL bcast_omp_rgen(Var,size(Var),buffer_r)

  END SUBROUTINE bcast_omp_r2


  SUBROUTINE bcast_omp_r3(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:)

    CALL check_buffer_r(size(Var))
    CALL bcast_omp_rgen(Var,size(Var),buffer_r)

  END SUBROUTINE bcast_omp_r3


  SUBROUTINE bcast_omp_r4(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:,:)
   
    CALL check_buffer_r(size(Var))
    CALL bcast_omp_rgen(Var,size(Var),buffer_r)

  END SUBROUTINE bcast_omp_r4

  
!! -- Les booleans -- !!

  SUBROUTINE bcast_omp_l(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var
    LOGICAL :: Var_tmp(1)
    
    Var_tmp(1)=Var
    CALL check_buffer_l(1)
    CALL bcast_omp_lgen(Var_tmp,1,buffer_l)
    Var=Var_tmp(1)

  END SUBROUTINE bcast_omp_l


  SUBROUTINE bcast_omp_l1(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:)
   
    CALL check_buffer_l(size(Var))
    CALL bcast_omp_lgen(Var,size(Var),buffer_l)

  END SUBROUTINE bcast_omp_l1


  SUBROUTINE bcast_omp_l2(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:)
   
    CALL check_buffer_l(size(Var))
    CALL bcast_omp_lgen(Var,size(Var),buffer_l)

  END SUBROUTINE bcast_omp_l2


  SUBROUTINE bcast_omp_l3(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:)

    CALL check_buffer_l(size(Var))
    CALL bcast_omp_lgen(Var,size(Var),buffer_l)

  END SUBROUTINE bcast_omp_l3


  SUBROUTINE bcast_omp_l4(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:,:)
   
    CALL check_buffer_l(size(Var))
    CALL bcast_omp_lgen(Var,size(Var),buffer_l)

  END SUBROUTINE bcast_omp_l4



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Scatter   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scatter_omp_i(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut

    CALL Check_buffer_i(size(VarIn))   
    CALL scatter_omp_igen(VarIn,Varout,1,buffer_i)
    
  END SUBROUTINE scatter_omp_i


  SUBROUTINE scatter_omp_i1(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL Check_buffer_i(size(VarIn))   
    CALL scatter_omp_igen(VarIn,Varout,Size(VarOut,2),buffer_i)
    
  END SUBROUTINE scatter_omp_i1
  
  
  SUBROUTINE scatter_omp_i2(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL Check_buffer_i(size(VarIn))   
    CALL scatter_omp_igen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3),buffer_i)

  END SUBROUTINE scatter_omp_i2


  SUBROUTINE scatter_omp_i3(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL Check_buffer_i(size(VarIn))   
    CALL scatter_omp_igen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4),buffer_i)
  
  END SUBROUTINE scatter_omp_i3




  SUBROUTINE scatter_omp_r(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut

    CALL Check_buffer_r(size(VarIn))   
    CALL scatter_omp_rgen(VarIn,Varout,1,buffer_r)
    
  END SUBROUTINE scatter_omp_r


  SUBROUTINE scatter_omp_r1(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL Check_buffer_r(size(VarIn))   
    CALL scatter_omp_rgen(VarIn,Varout,Size(VarOut,2),buffer_r)
        
  END SUBROUTINE scatter_omp_r1
  
  
  SUBROUTINE scatter_omp_r2(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL Check_buffer_r(size(VarIn))   
    CALL scatter_omp_rgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3),buffer_r)

  END SUBROUTINE scatter_omp_r2


  SUBROUTINE scatter_omp_r3(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL Check_buffer_r(size(VarIn))   
    CALL scatter_omp_rgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4),buffer_r)
  
  END SUBROUTINE scatter_omp_r3
  


  SUBROUTINE scatter_omp_l(VarIn, VarOut)
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut

    CALL Check_buffer_l(size(VarIn))   
    CALL scatter_omp_lgen(VarIn,Varout,1,buffer_l)
    
  END SUBROUTINE scatter_omp_l


  SUBROUTINE scatter_omp_l1(VarIn, VarOut)
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL Check_buffer_l(size(VarIn))   
    CALL scatter_omp_lgen(VarIn,Varout,Size(VarOut,2),buffer_l)
    
  END SUBROUTINE scatter_omp_l1
  
  
  SUBROUTINE scatter_omp_l2(VarIn, VarOut)
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL Check_buffer_l(size(VarIn))   
    CALL scatter_omp_lgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3),buffer_l)

  END SUBROUTINE scatter_omp_l2


  SUBROUTINE scatter_omp_l3(VarIn, VarOut)
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL Check_buffer_l(size(VarIn))   
    CALL scatter_omp_lgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4),buffer_l)
  
  END SUBROUTINE scatter_omp_l3  
  

  SUBROUTINE gather_omp_i(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut

    CALL Check_buffer_i(size(VarOut))   
    CALL gather_omp_igen(VarIn,Varout,1,buffer_i)
    
  END SUBROUTINE gather_omp_i


  SUBROUTINE gather_omp_i1(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL Check_buffer_i(size(VarOut))   
    CALL gather_omp_igen(VarIn,Varout,Size(VarIn,2),buffer_i)
    
  END SUBROUTINE gather_omp_i1


  SUBROUTINE gather_omp_i2(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL Check_buffer_i(size(VarOut))   
    CALL gather_omp_igen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3),buffer_i)
          
  END SUBROUTINE gather_omp_i2
  

  SUBROUTINE gather_omp_i3(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL Check_buffer_i(size(VarOut))   
    CALL gather_omp_igen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4),buffer_i)
    
  END SUBROUTINE gather_omp_i3



  SUBROUTINE gather_omp_r(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut

    CALL Check_buffer_r(size(VarOut))   
    CALL gather_omp_rgen(VarIn,Varout,1,buffer_r)
    
  END SUBROUTINE gather_omp_r


  SUBROUTINE gather_omp_r1(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL Check_buffer_r(size(VarOut))   
    CALL gather_omp_rgen(VarIn,Varout,Size(VarIn,2),buffer_r)
        
  END SUBROUTINE gather_omp_r1


  SUBROUTINE gather_omp_r2(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL Check_buffer_r(size(VarOut))   
    CALL gather_omp_rgen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3),buffer_r)
    
  END SUBROUTINE gather_omp_r2
  

  SUBROUTINE gather_omp_r3(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    CALL Check_buffer_r(size(VarOut))       
    CALL gather_omp_rgen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4),buffer_r)
    
  END SUBROUTINE gather_omp_r3


  SUBROUTINE gather_omp_l(VarIn, VarOut)
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut

    CALL Check_buffer_l(size(VarOut))   
    CALL gather_omp_lgen(VarIn,Varout,1,buffer_l)
    
  END SUBROUTINE gather_omp_l


  SUBROUTINE gather_omp_l1(VarIn, VarOut)
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL Check_buffer_l(size(VarOut))   
    CALL gather_omp_lgen(VarIn,Varout,Size(VarIn,2),buffer_l)
    
  END SUBROUTINE gather_omp_l1


  SUBROUTINE gather_omp_l2(VarIn, VarOut)
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL Check_buffer_l(size(VarOut))   
    CALL gather_omp_lgen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3),buffer_l)
    
  END SUBROUTINE gather_omp_l2
  

  SUBROUTINE gather_omp_l3(VarIn, VarOut)
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL Check_buffer_l(size(VarOut))   
    CALL gather_omp_lgen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4),buffer_l)
    
  END SUBROUTINE gather_omp_l3




  SUBROUTINE reduce_sum_omp_i(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN)  :: VarIn
    INTEGER,INTENT(OUT) :: VarOut
    INTEGER             :: VarIn_tmp(1)
    INTEGER             :: VarOut_tmp(1)
    
    VarIn_tmp(1)=VarIn
    CALL Check_buffer_i(1)   
    CALL reduce_sum_omp_igen(VarIn_tmp,Varout_tmp,1,buffer_i)
    VarOut=VarOut_tmp(1)
    
  END SUBROUTINE reduce_sum_omp_i

  SUBROUTINE reduce_sum_omp_i1(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
    CALL Check_buffer_i(size(VarIn))   
    CALL reduce_sum_omp_igen(VarIn,Varout,Size(VarIn),buffer_i)
   
  END SUBROUTINE reduce_sum_omp_i1
  
  
  SUBROUTINE reduce_sum_omp_i2(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL Check_buffer_i(size(VarIn))   
    CALL reduce_sum_omp_igen(VarIn,Varout,Size(VarIn),buffer_i)
  
  END SUBROUTINE reduce_sum_omp_i2


  SUBROUTINE reduce_sum_omp_i3(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL Check_buffer_i(size(VarIn))   
    CALL reduce_sum_omp_igen(VarIn,Varout,Size(VarIn),buffer_i)
  
  END SUBROUTINE reduce_sum_omp_i3


  SUBROUTINE reduce_sum_omp_i4(VarIn, VarOut)
    IMPLICIT NONE

    INTEGER,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
  
    CALL Check_buffer_i(size(VarIn))   
    CALL reduce_sum_omp_igen(VarIn,Varout,Size(VarIn),buffer_i)
  
  END SUBROUTINE reduce_sum_omp_i4


  SUBROUTINE reduce_sum_omp_r(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN)  :: VarIn
    REAL,INTENT(OUT) :: VarOut
    REAL             :: VarIn_tmp(1)
    REAL             :: VarOut_tmp(1)
    
    VarIn_tmp(1)=VarIn
    CALL Check_buffer_r(1)   
    CALL reduce_sum_omp_rgen(VarIn_tmp,Varout_tmp,1,buffer_r)
    VarOut=VarOut_tmp(1)
  
  END SUBROUTINE reduce_sum_omp_r

  SUBROUTINE reduce_sum_omp_r1(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    CALL Check_buffer_r(size(VarIn))   
    CALL reduce_sum_omp_rgen(VarIn,Varout,Size(VarIn),buffer_r)
   
  END SUBROUTINE reduce_sum_omp_r1
  
  
  SUBROUTINE reduce_sum_omp_r2(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL Check_buffer_r(size(VarIn))   
    CALL reduce_sum_omp_rgen(VarIn,Varout,Size(VarIn),buffer_r)
  
  END SUBROUTINE reduce_sum_omp_r2


  SUBROUTINE reduce_sum_omp_r3(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL Check_buffer_r(size(VarIn))   
    CALL reduce_sum_omp_rgen(VarIn,Varout,Size(VarIn),buffer_r)
  
  END SUBROUTINE reduce_sum_omp_r3


  SUBROUTINE reduce_sum_omp_r4(VarIn, VarOut)
    IMPLICIT NONE

    REAL,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
  
    CALL Check_buffer_r(size(VarIn))   
    CALL reduce_sum_omp_rgen(VarIn,Varout,Size(VarIn),buffer_r)
  
  END SUBROUTINE reduce_sum_omp_r4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    LES ROUTINES GENERIQUES    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE bcast_omp_cgen(Var,Nb,Buff)
  IMPLICIT NONE
    
    CHARACTER(LEN=*),INTENT(INOUT) :: Var
    CHARACTER(LEN=*),INTENT(INOUT) :: Buff
    INTEGER,INTENT(IN) :: Nb
    
    INTEGER :: i
  
  !$OMP MASTER
      Buff=Var
  !$OMP END MASTER
  !$OMP BARRIER

    DO i=1,Nb
      Var=Buff
    ENDDO
  !$OMP BARRIER      
  
  END SUBROUTINE bcast_omp_cgen


      
  SUBROUTINE bcast_omp_igen(Var,Nb,Buff)
  IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: Nb
    INTEGER,DIMENSION(Nb),INTENT(INOUT) :: Var
    INTEGER,DIMENSION(Nb),INTENT(INOUT) :: Buff

    INTEGER :: i
    
  !$OMP MASTER
    DO i=1,Nb
      Buff(i)=Var(i)
    ENDDO
  !$OMP END MASTER
  !$OMP BARRIER

    DO i=1,Nb
      Var(i)=Buff(i)
    ENDDO
  !$OMP BARRIER        

  END SUBROUTINE bcast_omp_igen


  SUBROUTINE bcast_omp_rgen(Var,Nb,Buff)
  IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: Nb
    REAL,DIMENSION(Nb),INTENT(INOUT) :: Var
    REAL,DIMENSION(Nb),INTENT(INOUT) :: Buff

    INTEGER :: i
    
  !$OMP MASTER
    DO i=1,Nb
      Buff(i)=Var(i)
    ENDDO
  !$OMP END MASTER
  !$OMP BARRIER

    DO i=1,Nb
      Var(i)=Buff(i)
    ENDDO
  !$OMP BARRIER        

  END SUBROUTINE bcast_omp_rgen

  SUBROUTINE bcast_omp_lgen(Var,Nb,Buff)
  IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: Nb
    LOGICAL,DIMENSION(Nb),INTENT(INOUT) :: Var
    LOGICAL,DIMENSION(Nb),INTENT(INOUT) :: Buff
  
    INTEGER :: i
    
  !$OMP MASTER
    DO i=1,Nb
      Buff(i)=Var(i)
    ENDDO
  !$OMP END MASTER
  !$OMP BARRIER

    DO i=1,Nb
      Var(i)=Buff(i)
    ENDDO
  !$OMP BARRIER        

  END SUBROUTINE bcast_omp_lgen


  SUBROUTINE scatter_omp_igen(VarIn,VarOut,dimsize,Buff)
    USE mod_phys_lmdz_omp_data
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN),DIMENSION(klon_mpi,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(klon_omp,dimsize) :: VarOut
    INTEGER,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

    INTEGER :: i,ij
    
  !$OMP MASTER
    DO i=1,dimsize
      DO ij=1,klon_mpi
        Buff(ij,i)=VarIn(ij,i)
      ENDDO
    ENDDO  
  !$OMP END MASTER
  !$OMP BARRIER
 
    DO i=1,dimsize
      DO ij=1,klon_omp
        VarOut(ij,i)=Buff(klon_omp_begin-1+ij,i)
      ENDDO
    ENDDO
  !$OMP BARRIER  
 
  END SUBROUTINE scatter_omp_igen


  SUBROUTINE scatter_omp_rgen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN),DIMENSION(klon_mpi,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(klon_omp,dimsize) :: VarOut
    REAL,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

    INTEGER :: i,ij
    
  !$OMP MASTER
    DO i=1,dimsize
      DO ij=1,klon_mpi
        Buff(ij,i)=VarIn(ij,i)
      ENDDO
    ENDDO  
  !$OMP END MASTER
  !$OMP BARRIER

    DO i=1,dimsize
      DO ij=1,klon_omp
        VarOut(ij,i)=Buff(klon_omp_begin-1+ij,i)
      ENDDO
    ENDDO
  !$OMP BARRIER  

  END SUBROUTINE scatter_omp_rgen


  SUBROUTINE scatter_omp_lgen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN),DIMENSION(klon_mpi,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(klon_omp,dimsize) :: VarOut
    LOGICAL,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

    INTEGER :: i,ij
    
 !$OMP MASTER 
    DO i=1,dimsize
      DO ij=1,klon_mpi
        Buff(ij,i)=VarIn(ij,i)
      ENDDO
    ENDDO  
  !$OMP END MASTER
  !$OMP BARRIER

    DO i=1,dimsize
      DO ij=1,klon_omp
        VarOut(ij,i)=Buff(klon_omp_begin-1+ij,i)
      ENDDO
    ENDDO
  !$OMP BARRIER  

  END SUBROUTINE scatter_omp_lgen





  SUBROUTINE gather_omp_igen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN),DIMENSION(klon_omp,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(klon_mpi,dimsize) :: VarOut
    INTEGER,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

    INTEGER :: i,ij
    
    DO i=1,dimsize
      DO ij=1,klon_omp
        Buff(klon_omp_begin-1+ij,i)=VarIn(ij,i)
      ENDDO
    ENDDO
  !$OMP BARRIER  
  
  
  !$OMP MASTER
    DO i=1,dimsize
      DO ij=1,klon_mpi
        VarOut(ij,i)=Buff(ij,i)
      ENDDO
    ENDDO  
  !$OMP END MASTER
  !$OMP BARRIER

  END SUBROUTINE gather_omp_igen


  SUBROUTINE gather_omp_rgen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN),DIMENSION(klon_omp,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(klon_mpi,dimsize) :: VarOut
    REAL,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

    INTEGER :: i,ij
    
    DO i=1,dimsize
      DO ij=1,klon_omp
        Buff(klon_omp_begin-1+ij,i)=VarIn(ij,i)
      ENDDO
    ENDDO
  !$OMP BARRIER  


  !$OMP MASTER
    DO i=1,dimsize
      DO ij=1,klon_mpi
        VarOut(ij,i)=Buff(ij,i)
      ENDDO
    ENDDO  
  !$OMP END MASTER
  !$OMP BARRIER

  END SUBROUTINE gather_omp_rgen


  SUBROUTINE gather_omp_lgen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN),DIMENSION(klon_omp,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(klon_mpi,dimsize) :: VarOut
    LOGICAL,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

    INTEGER :: i,ij
    
    DO i=1,dimsize
      DO ij=1,klon_omp
        Buff(klon_omp_begin-1+ij,i)=VarIn(ij,i)
      ENDDO
    ENDDO
  !$OMP BARRIER  


  !$OMP MASTER
    DO i=1,dimsize
      DO ij=1,klon_mpi
        VarOut(ij,i)=Buff(ij,i)
      ENDDO
    ENDDO  
  !$OMP END MASTER
  !$OMP BARRIER

  END SUBROUTINE gather_omp_lgen


  SUBROUTINE reduce_sum_omp_igen(VarIn,VarOut,dimsize,Buff)
  IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN),DIMENSION(dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(dimsize) :: VarOut
    INTEGER,INTENT(INOUT),DIMENSION(dimsize) :: Buff

    INTEGER :: i

  !$OMP MASTER
    Buff(:)=0
  !$OMP END MASTER
  !$OMP BARRIER
  
  !$OMP CRITICAL     
    DO i=1,dimsize
      Buff(i)=Buff(i)+VarIn(i)
    ENDDO
  !$OMP END CRITICAL
  !$OMP BARRIER  
  
  !$OMP MASTER
    DO i=1,dimsize
      VarOut(i)=Buff(i)
    ENDDO
  !$OMP END MASTER
  !$OMP BARRIER
  
  END SUBROUTINE reduce_sum_omp_igen

  SUBROUTINE reduce_sum_omp_rgen(VarIn,VarOut,dimsize,Buff)
  IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN),DIMENSION(dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(dimsize) :: VarOut
    REAL,INTENT(INOUT),DIMENSION(dimsize) :: Buff

    INTEGER :: i

  !$OMP MASTER
    Buff(:)=0
  !$OMP END MASTER
  !$OMP BARRIER
  
  !$OMP CRITICAL     
    DO i=1,dimsize
      Buff(i)=Buff(i)+VarIn(i)
    ENDDO
  !$OMP END CRITICAL
  !$OMP BARRIER  
  
  !$OMP MASTER
    DO i=1,dimsize
      VarOut(i)=Buff(i)
    ENDDO
  !$OMP END MASTER
  !$OMP BARRIER
  
  END SUBROUTINE reduce_sum_omp_rgen

END MODULE mod_phys_lmdz_omp_transfert
