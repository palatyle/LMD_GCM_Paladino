










MODULE sort_mod

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Sorting algorithms
  ! ~~~~~~~~~~~~~~~~~~~
  ! Contains quicksort (qsort) and insertion sort (isort)
  ! -> Sorts an array (real only for now) by increasing values
  ! -> If wanted, output indexes of permutations if you want to apply the same sorting to other arrays
  !          ( this is quicker than having other arrays going through the whole process )
  ! 
  ! TODO : Extend interface to integers - Enable decreasing order
  ! 
  ! Author : J. Vatant d'Ollone - 2018
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  INTERFACE qsort
     MODULE PROCEDURE qsort_r
     MODULE PROCEDURE qsort_outp_r
  END INTERFACE qsort

  INTERFACE isort
     MODULE PROCEDURE isort_r
     MODULE PROCEDURE isort_outp_r
  END INTERFACE isort

CONTAINS

  ! 1-pivot quicksort
  RECURSIVE SUBROUTINE qsort_r(A,n)

    IMPLICIT NONE

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! I/O
    !logical,               intent(in)    :: o !! Order T=increasing/F=decreasing
    integer,               intent(in)    :: n !! Size of array
    real,    dimension(n), intent(inout) :: A !! Array to sort

    ! Local variables
    integer :: left, right
    real    :: random
    real    :: pivot
    real    :: temp
    integer :: marker
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (n > 1) then

       call random_number(random)
       pivot = A(int(random*real(n-1))+1) ! random pivot (not best performance, but avoids worst-case)
       left = 0
       right = n + 1

       do while (left < right)
          right = right - 1
          do while (A(right) > pivot)
             right = right - 1
          enddo
          left = left + 1
          do while (A(left) < pivot)
             left = left + 1
          enddo
          if (left < right) then
             temp     = A(left)
             A(left)  = A(right)
             A(right) = temp
          endif
       enddo

       if (left == right) then
          marker = left + 1
       else
          marker = left
       endif

       call qsort(A(:marker-1), marker-1)
       call qsort(A(marker:),   n-marker+1)

    endif

  END SUBROUTINE qsort_r


  ! 1-pivot quicksort + output array of permutations
  RECURSIVE SUBROUTINE qsort_outp_r(A,n,P)

    IMPLICIT NONE

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! I/O
    !logical,               intent(in)    :: o !! Order T=increasing/F=decreasing
    integer,               intent(in)    :: n !! Size of array
    integer, dimension(n), intent(inout) :: P !! Output array of permutations
    real,    dimension(n), intent(inout) :: A !! Array to sort

    ! Local variables
    integer :: i
    integer :: left, right
    real    :: random
    real    :: pivot
    real    :: tempA
    integer :: tempP
    integer :: marker
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (n > 1) then

       call random_number(random)
       pivot = A(int(random*real(n-1))+1) ! random pivot (not best performance, but avoids worst-case)
       left = 0
       right = n + 1

       do while (left < right)
          right = right - 1
          do while (A(right) > pivot)
             right = right - 1
          enddo
          left = left + 1
          do while (A(left) < pivot)
             left = left + 1
          enddo
          if (left < right) then
             tempA    = A(left)
             tempP    = P(left)
             A(left)  = A(right)
             P(left)  = P(right)
             A(right) = tempA
             P(right) = tempP
          endif
       enddo

       if (left == right) then
          marker = left + 1
       else
          marker = left
       endif

       call qsort(A(:marker-1), marker-1, P(:marker-1))
       call qsort(A(marker:),   n-marker+1, P(marker:))

    endif

  END SUBROUTINE qsort_outp_r


  ! Insertion sort
  SUBROUTINE isort_r(A,n)
    IMPLICIT NONE

    INTEGER,               INTENT(in)    :: n !! Size of array
    REAL,    DIMENSION(n), INTENT(inout) :: A !! Array to sort

    INTEGER ::  i, j
    REAL    ::  x

    DO i = 2, n
       x = a(i)
       j = i - 1
       DO WHILE (j >= 1)
          IF (a(j) <= x) EXIT
          a(j + 1) = a(j)
          j = j - 1
       ENDDO
       a(j + 1) = x
    ENDDO
  END SUBROUTINE isort_r


  ! Insertion sort + output array of permutations
  SUBROUTINE isort_outp_r(A,n,P)
    IMPLICIT NONE

    INTEGER,               INTENT(in)    :: n !! Size of array
    INTEGER, DIMENSION(n), INTENT(inout) :: P !! Output array of permutations
    REAL,    DIMENSION(n), INTENT(inout) :: A !! Array to sort

    INTEGER ::  i, j, px
    REAL    ::  x

    DO i = 2, n
       x = a(i)
       px = p(i)
       j = i - 1
       DO WHILE (j >= 1)
          IF (a(j) <= x) EXIT
          a(j + 1) = a(j)
          p(j + 1) = p(j)
          j = j - 1
       ENDDO
       a(j + 1) = x
       p(j + 1) = px
    ENDDO
  END SUBROUTINE isort_outp_r


END MODULE sort_mod

