MODULE buffer_mod

PRIVATE
  REAL,PARAMETER :: grow_factor=1.5

  INTEGER, POINTER, SAVE :: buffer_i(:)
  INTEGER,SAVE :: size_buffer_i = 0
!$OMP THREADPRIVATE(buffer_i,size_buffer_i)

  REAL,POINTER,SAVE      :: buffer_r(:)
  INTEGER,SAVE :: size_buffer_r = 0 
!$OMP THREADPRIVATE(buffer_r,size_buffer_r)
  
  LOGICAL,POINTER,SAVE   :: buffer_l(:)
  INTEGER,SAVE :: size_buffer_l = 0
!$OMP THREADPRIVATE(buffer_l,size_buffer_l)

  CHARACTER,POINTER,SAVE :: buffer_c(:)
  INTEGER,SAVE :: size_buffer_c = 0
!$OMP THREADPRIVATE(buffer_c,size_buffer_c)

INTERFACE get_buffer
  MODULE PROCEDURE get_buffer_i, get_buffer_r, get_buffer_l, get_buffer_c
END INTERFACE
  
PUBLIC :: get_buffer

CONTAINS

  SUBROUTINE get_buffer_i(buff,buff_size)
  IMPLICIT NONE
    INTEGER,POINTER    :: buff(:)
    INTEGER,INTENT(IN) :: buff_size 

    IF (buff_size>size_buffer_i) THEN
      DEALLOCATE(buffer_i)
      size_buffer_i=MAX(2,INT(size_buffer_i*grow_factor))
      ALLOCATE(buffer_i(size_buffer_i))
    ENDIF
    
    buff=>buffer_i
  END SUBROUTINE get_buffer_i

  SUBROUTINE get_buffer_r(buff,buff_size)
  IMPLICIT NONE
    REAL,POINTER       :: buff(:)
    INTEGER,INTENT(IN) :: buff_size 

    IF (buff_size>size_buffer_r) THEN
      DEALLOCATE(buffer_r)
      size_buffer_r=MAX(2,INT(size_buffer_r*grow_factor))
      ALLOCATE(buffer_r(size_buffer_r))
    ENDIF
    
    buff=>buffer_r
  END SUBROUTINE get_buffer_r

  SUBROUTINE get_buffer_l(buff,buff_size)
  IMPLICIT NONE
    LOGICAL,POINTER    :: buff(:)
    INTEGER,INTENT(IN) :: buff_size 

    IF (buff_size>size_buffer_l) THEN
      DEALLOCATE(buffer_l)
      size_buffer_l=MAX(2,INT(size_buffer_l*grow_factor))
      ALLOCATE(buffer_l(size_buffer_l))
    ENDIF
    
    buff=>buffer_l
  END SUBROUTINE get_buffer_l
  
  SUBROUTINE get_buffer_c(buff,buff_size)
  IMPLICIT NONE
    CHARACTER,POINTER  :: buff(:)
    INTEGER,INTENT(IN) :: buff_size 

    IF (buff_size>size_buffer_c) THEN
      DEALLOCATE(buffer_c)
      size_buffer_c=MAX(2,INT(size_buffer_c*grow_factor))
      ALLOCATE(buffer_c(size_buffer_c))
    ENDIF
    
    buff=>buffer_c
  END SUBROUTINE get_buffer_c
  
END MODULE buffer_mod
