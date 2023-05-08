










! $Id$
module interpolation

  ! From Press et al., 1996, version 2.10a
  ! B3 Interpolation and Extrapolation

  IMPLICIT NONE 

contains

  pure FUNCTION locate(xx,x)

    REAL, DIMENSION(:), INTENT(IN) :: xx
    REAL, INTENT(IN) :: x
    INTEGER  locate

    ! Given an array xx(1:N), and given a value x, returns a value j,
    ! between 0 and N, such that x is between xx(j) and xx(j + 1). xx
    ! must be monotonic, either increasing or decreasing. j = 0 or j =
    ! N is returned to indicate that x is out of range. This
    ! procedure should not be called with a zero-sized array argument.
    ! See notes.

    INTEGER  n,jl,jm,ju
    LOGICAL  ascnd

    !----------------------------

    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    ! (True if ascending order of table, false otherwise.)
    ! Initialize lower and upper limits:
    jl=0
    ju=n+1
    do while (ju-jl > 1)
       jm=(ju+jl)/2 ! Compute a midpoint,
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm ! and replace either the lower limit
       else
          ju=jm ! or the upper limit, as appropriate.
       end if
    end do
    ! {ju == jl + 1}

    ! {(ascnd .and. xx(jl) <= x < xx(jl+1)) 
    !  .neqv. 
    !  (.not. ascnd .and. xx(jl+1) <= x < xx(jl))}

    ! Then set the output, being careful with the endpoints:
    if (x == xx(1)) then
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if

  END FUNCTION locate

  !***************************

  pure SUBROUTINE hunt(xx,x,jlo)

    ! Given an array xx(1:N ), and given a value x, returns a value
    ! jlo such that x is between xx(jlo) and xx(jlo+1). xx must be
    ! monotonic, either increasing or decreasing. jlo = 0 or jlo = N is
    ! returned to indicate that x is out of range. jlo on input is taken as
    ! the initial guess for jlo on output.
    ! Modified so that it uses the information "jlo = 0" on input.

    INTEGER, INTENT(INOUT) :: jlo
    REAL, INTENT(IN) :: x
    REAL, DIMENSION(:), INTENT(IN) :: xx
    INTEGER  n,inc,jhi,jm
    LOGICAL  ascnd, hunt_up

    !-----------------------------------------------------

    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    ! (True if ascending order of table, false otherwise.)
    if (jlo < 0 .or. jlo > n) then
       ! Input guess not useful. Go immediately to bisection.
       jlo=0
       jhi=n+1
    else
       inc=1 ! Set the hunting increment.
       if (jlo == 0) then
          hunt_up = .true.
       else
          hunt_up = x >= xx(jlo) .eqv. ascnd
       end if
       if (hunt_up) then ! Hunt up:
          do
             jhi=jlo+inc
             if (jhi > n) then ! Done hunting, since off end of table.
                jhi=n+1
                exit
             else
                if (x < xx(jhi) .eqv. ascnd) exit
                jlo=jhi ! Not done hunting,
                inc=inc+inc ! so double the increment
             end if
          end do ! and try again.
       else ! Hunt down:
          jhi=jlo
          do
             jlo=jhi-inc
             if (jlo < 1) then ! Done hunting, since off end of table.
                jlo=0
                exit
             else
                if (x >= xx(jlo) .eqv. ascnd) exit
                jhi=jlo ! Not done hunting,
                inc=inc+inc ! so double the increment
             end if
          end do ! and try again.
       end if
    end if ! Done hunting, value bracketed.

    do ! Hunt is done, so begin the final bisection phase:
       if (jhi-jlo <= 1) then
          if (x == xx(n)) jlo=n-1
          if (x == xx(1)) jlo=1
          exit
       else
          jm=(jhi+jlo)/2
          if (x >= xx(jm) .eqv. ascnd) then
             jlo=jm
          else
             jhi=jm
          end if
       end if
    end do

  END SUBROUTINE hunt

end module interpolation
