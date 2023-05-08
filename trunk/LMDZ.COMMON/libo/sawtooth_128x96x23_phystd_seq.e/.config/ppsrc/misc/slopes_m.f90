










module slopes_m

  ! Author: Lionel GUEZ

  implicit none

  interface slopes
     ! This generic function computes second order slopes with Van
     ! Leer slope-limiting, given cell averages. Reference: Dukowicz,
     ! 1987, SIAM Journal on Scientific and Statistical Computing, 8,
     ! 305.

     ! The only difference between the specific functions is the rank
     ! of the first argument and the equal rank of the result.

     ! real, intent(in), rank >= 1:: f ! (n, ...) cell averages, n must be >= 1
     ! real, intent(in):: x(:) ! (n + 1) cell edges
     ! real slopes, same shape as f ! (n, ...)

     module procedure slopes1, slopes2, slopes3, slopes4
  end interface

  private
  public slopes

contains

  pure function slopes1(f, x)

    real, intent(in):: f(:)
    real, intent(in):: x(:)
    real slopes1(size(f))

    ! Local:
    integer n, i
    real xc(size(f)) ! (n) cell centers
    real h

    !------------------------------------------------------

    n = size(f)
    forall (i = 1:n) xc(i) = (x(i) + x(i + 1)) / 2.
    slopes1(1) = 0.
    slopes1(n) = 0.

    do i = 2, n - 1
       if (f(i) >= max(f(i - 1), f(i + 1)) .or. f(i) &
            <= min(f(i - 1), f(i + 1))) then
          ! Local extremum
          slopes1(i) = 0.
       else
          ! (f(i - 1), f(i), f(i + 1)) strictly monotonous

          ! Second order slope:
          slopes1(i) = (f(i + 1) - f(i - 1)) / (xc(i + 1) - xc(i - 1))

          ! Slope limitation:
          h = abs(x(i + 1) - xc(i))
          slopes1(i) = sign(min(abs(slopes1(i)), abs(f(i + 1) - f(i)) / h, &
               abs(f(i) - f(i - 1)) / h), slopes1(i))
       end if
    end do

  end function slopes1

  !*************************************************************

  pure function slopes2(f, x)

    real, intent(in):: f(:, :)
    real, intent(in):: x(:)
    real slopes2(size(f, 1), size(f, 2))

    ! Local:
    integer n, i, j
    real xc(size(f, 1)) ! (n) cell centers
    real h(2:size(f, 1) - 1), delta_xc(2:size(f, 1) - 1) ! (2:n - 1)

    !------------------------------------------------------

    n = size(f, 1)
    forall (i = 1:n) xc(i) = (x(i) + x(i + 1)) / 2.

    forall (i = 2:n - 1) 
       h(i) = abs(x(i + 1) - xc(i))
       delta_xc(i) = xc(i + 1) - xc(i - 1)
    end forall

    do j = 1, size(f, 2)
       slopes2(1, j) = 0.

       do i = 2, n - 1
          if (f(i, j) >= max(f(i - 1, j), f(i + 1, j)) .or. &
               f(i, j) <= min(f(i - 1, j), f(i + 1, j))) then
             ! Local extremum
             slopes2(i, j) = 0.
          else
             ! (f(i - 1, j), f(i, j), f(i + 1, j))
             ! strictly monotonous

             ! Second order slope:
             slopes2(i, j) = (f(i + 1, j) - f(i - 1, j)) / delta_xc(i)

             ! Slope limitation:
             slopes2(i, j) = sign(min(abs(slopes2(i, j)), &
                  abs(f(i + 1, j) - f(i, j)) / h(i), &
                  abs(f(i, j) - f(i - 1, j)) / h(i)), slopes2(i, j))
          end if
       end do

       slopes2(n, j) = 0.
    end do

  end function slopes2

  !*************************************************************

  pure function slopes3(f, x)

    real, intent(in):: f(:, :, :)
    real, intent(in):: x(:)
    real slopes3(size(f, 1), size(f, 2), size(f, 3))

    ! Local:
    integer n, i, j, k
    real xc(size(f, 1)) ! (n) cell centers
    real h(2:size(f, 1) - 1), delta_xc(2:size(f, 1) - 1) ! (2:n - 1)

    !------------------------------------------------------

    n = size(f, 1)
    forall (i = 1:n) xc(i) = (x(i) + x(i + 1)) / 2.

    forall (i = 2:n - 1) 
       h(i) = abs(x(i + 1) - xc(i))
       delta_xc(i) = xc(i + 1) - xc(i - 1)
    end forall

    do k = 1, size(f, 3)
       do j = 1, size(f, 2)
          slopes3(1, j, k) = 0.

          do i = 2, n - 1
             if (f(i, j, k) >= max(f(i - 1, j, k), f(i + 1, j, k)) .or. &
                  f(i, j, k) <= min(f(i - 1, j, k), f(i + 1, j, k))) then
                ! Local extremum
                slopes3(i, j, k) = 0.
             else
                ! (f(i - 1, j, k), f(i, j, k), f(i + 1, j, k))
                ! strictly monotonous

                ! Second order slope:
                slopes3(i, j, k) = (f(i + 1, j, k) - f(i - 1, j, k)) &
                     / delta_xc(i)

                ! Slope limitation:
                slopes3(i, j, k) = sign(min(abs(slopes3(i, j, k)), &
                     abs(f(i + 1, j, k) - f(i, j, k)) / h(i), &
                     abs(f(i, j, k) - f(i - 1, j, k)) / h(i)), slopes3(i, j, k))
             end if
          end do

          slopes3(n, j, k) = 0.
       end do
    end do

  end function slopes3

  !*************************************************************

  pure function slopes4(f, x)

    real, intent(in):: f(:, :, :, :)
    real, intent(in):: x(:)
    real slopes4(size(f, 1), size(f, 2), size(f, 3), size(f, 4))

    ! Local:
    integer n, i, j, k, l
    real xc(size(f, 1)) ! (n) cell centers
    real h(2:size(f, 1) - 1), delta_xc(2:size(f, 1) - 1) ! (2:n - 1)

    !------------------------------------------------------

    n = size(f, 1)
    forall (i = 1:n) xc(i) = (x(i) + x(i + 1)) / 2.

    forall (i = 2:n - 1) 
       h(i) = abs(x(i + 1) - xc(i))
       delta_xc(i) = xc(i + 1) - xc(i - 1)
    end forall

    do l = 1, size(f, 4)
       do k = 1, size(f, 3)
          do j = 1, size(f, 2)
             slopes4(1, j, k, l) = 0.

             do i = 2, n - 1
                if (f(i, j, k, l) >= max(f(i - 1, j, k, l), f(i + 1, j, k, l)) &
                     .or. f(i, j, k, l) &
                     <= min(f(i - 1, j, k, l), f(i + 1, j, k, l))) then
                   ! Local extremum
                   slopes4(i, j, k, l) = 0.
                else
                   ! (f(i - 1, j, k, l), f(i, j, k, l), f(i + 1, j, k, l))
                   ! strictly monotonous

                   ! Second order slope:
                   slopes4(i, j, k, l) = (f(i + 1, j, k, l) &
                        - f(i - 1, j, k, l)) / delta_xc(i)

                   ! Slope limitation:
                   slopes4(i, j, k, l) = sign(min(abs(slopes4(i, j, k, l)), &
                        abs(f(i + 1, j, k, l) - f(i, j, k, l)) / h(i), &
                        abs(f(i, j, k, l) - f(i - 1, j, k, l)) / h(i)), &
                        slopes4(i, j, k, l))
                end if
             end do

             slopes4(n, j, k, l) = 0.
          end do
       end do
    end do

  end function slopes4

end module slopes_m
