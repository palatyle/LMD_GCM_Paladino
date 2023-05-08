










module regr1_conserv_m

  ! Author: Lionel GUEZ

  use assert_eq_m, only: assert_eq
  use assert_m, only: assert
  use interpolation, only: locate

  implicit none

  interface regr1_conserv

     ! This generic procedure regrids a piecewise linear function (not
     ! necessarily continuous) by averaging it. This is a conservative
     ! regridding. The regridding operation is done on the first
     ! dimension of the input array. Input are positions of cell
     ! edges. The target grid should be included in the source grid:
     ! no extrapolation is allowed.

     ! The only difference between the procedures is the rank of the
     ! first argument.

     ! real, intent(in), rank >= 1:: vs ! (ns, ...)
     ! averages of cells of the source grid
     ! vs(is, ...) for [xs(is), xs(is + 1)]

     ! real, intent(in):: xs(:) ! (ns + 1)
     ! edges of cells of the source grid, in strictly ascending order

     ! real, intent(in):: xt(:) ! (nt + 1)
     ! edges of cells of the target grid, in strictly ascending order

     ! real, intent(in), optional, rank >= 1:: slope ! (ns, ...)
     ! same rank as vs
     ! slopes inside cells of the source grid
     ! slope(is, ...) for [xs(is), xs(is + 1)]
     ! If not present, slopes are taken equal to 0. The regridding is
     ! then first order.

     ! real, intent(out), rank >= 1:: vt(nt, ...) 
     ! has the same rank as vs and slope
     ! averages of cells of the target grid
     ! vt(it, ...) for  [xt(it), xt(it + 1)]

     ! ns and nt must be >= 1.

     ! See notes for explanations on the algorithm and justification
     ! of algorithmic choices.

     module procedure regr11_conserv, regr12_conserv, regr13_conserv, &
          regr14_conserv
  end interface regr1_conserv

  private
  public regr1_conserv

contains

  subroutine regr11_conserv(vs, xs, xt, vt, slope)

    ! vs and slope have rank 1.

    real, intent(in):: vs(:)
    real, intent(in):: xs(:)
    real, intent(in):: xt(:)
    real, intent(out):: vt(:)
    real, intent(in), optional:: slope(:)

    ! Local:
    integer is, it, ns, nt
    logical slope_present

    !---------------------------------------------

    ns = assert_eq(size(vs), size(xs) - 1, "regr11_conserv ns")
    nt = assert_eq(size(xt) - 1, size(vt), "regr11_conserv nt")

    ! Quick check on sort order:
    call assert(xs(1) < xs(2), "regr11_conserv xs bad order")
    call assert(xt(1) < xt(2), "regr11_conserv xt bad order")

    call assert(xs(1) <= xt(1) .and. xt(nt + 1) <= xs(ns + 1), &
         "regr11_conserv extrapolation")
    slope_present = present(slope)

    is = locate(xs, xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
    do it = 1, nt
       ! 1 <= is <= ns
       ! xs(is) <= xt(it) < xs(is + 1)
       if (xt(it + 1) <= xs(is + 1)) then
          vt(it) = mean_lin(xt(it), xt(it + 1))
       else
          vt(it) = mean_lin(xt(it), xs(is + 1)) * (xs(is + 1) - xt(it))
          is = is + 1
          do while (xs(is + 1) < xt(it + 1))
             ! 1 <= is <= ns - 1
             vt(it) = vt(it) + (xs(is + 1) - xs(is)) * vs(is)
             is = is + 1
          end do
          ! 1 <= is <= ns
          vt(it) = (vt(it) + mean_lin(xs(is), xt(it + 1)) * (xt(it + 1) &
               - xs(is))) / (xt(it + 1) - xt(it))
       end if

       if (xs(is + 1) == xt(it + 1)) is = is + 1
       ! 1 <= is <= ns .or. it == nt
    end do

  contains

    real function mean_lin(a, b)

      ! mean in [a, b] of the linear function in [xs(is), xs(is + 1)]

      real, intent(in):: a, b

      !---------------------------------------------

      if (slope_present) then
         mean_lin = slope(is) / 2. * (a + b - xs(is) - xs(is + 1)) + vs(is)
      else
         mean_lin = vs(is)
      end if

    end function mean_lin

  end subroutine regr11_conserv

  !********************************************

  subroutine regr12_conserv(vs, xs, xt, vt, slope)

    ! vs and slope have rank 2.

    real, intent(in):: vs(:, :)
    real, intent(in):: xs(:)
    real, intent(in):: xt(:)
    real, intent(out):: vt(:, :)
    real, intent(in), optional:: slope(:, :)

    ! Local:
    integer is, it, ns, nt, n2
    logical slope_present

    !---------------------------------------------

    ns = assert_eq(size(vs, 1), size(xs) - 1, "regr12_conserv ns")
    nt = assert_eq(size(xt) - 1, size(vt, 1), "regr12_conserv nt")
    n2 = assert_eq(size(vs, 2), size(vt, 2), "regr12_conserv n2")

    ! Quick check on sort order:
    call assert(xs(1) < xs(2), "regr12_conserv xs bad order")
    call assert(xt(1) < xt(2), "regr12_conserv xt bad order")

    call assert(xs(1) <= xt(1) .and. xt(nt + 1) <= xs(ns + 1), &
         "regr12_conserv extrapolation")
    slope_present = present(slope)

    is = locate(xs, xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
    do it = 1, nt
       ! 1 <= is <= ns
       ! xs(is) <= xt(it) < xs(is + 1)
       if (xt(it + 1) <= xs(is + 1)) then
          vt(it, :) = mean_lin(xt(it), xt(it + 1))
       else
          vt(it, :) = mean_lin(xt(it), xs(is + 1)) * (xs(is + 1) - xt(it))
          is = is + 1
          do while (xs(is + 1) < xt(it + 1))
             ! 1 <= is <= ns - 1
             vt(it, :) = vt(it, :) + (xs(is + 1) - xs(is)) * vs(is, :)
             is = is + 1
          end do
          ! 1 <= is <= ns
          vt(it, :) = (vt(it, :) + mean_lin(xs(is), xt(it + 1)) * (xt(it + 1) &
               - xs(is))) / (xt(it + 1) - xt(it))
       end if

       if (xs(is + 1) == xt(it + 1)) is = is + 1
       ! 1 <= is <= ns .or. it == nt
    end do

  contains

    function mean_lin(a, b)

      ! mean in [a, b] of the linear function in [xs(is), xs(is + 1)]

      real, intent(in):: a, b
      real mean_lin(n2)

      !---------------------------------------------

      if (slope_present) then
         mean_lin = slope(is, :) / 2. * (a + b - xs(is) - xs(is + 1)) &
              + vs(is, :)
      else
         mean_lin = vs(is, :)
      end if

    end function mean_lin

  end subroutine regr12_conserv

  !********************************************

  subroutine regr13_conserv(vs, xs, xt, vt, slope)

    ! vs and slope have rank 3.

    real, intent(in):: vs(:, :, :)
    real, intent(in):: xs(:)
    real, intent(in):: xt(:)
    real, intent(out):: vt(:, :, :)
    real, intent(in), optional:: slope(:, :, :)

    ! Local:
    integer is, it, ns, nt, n2, n3
    logical slope_present

    !---------------------------------------------

    ns = assert_eq(size(vs, 1), size(xs) - 1, "regr13_conserv ns")
    nt = assert_eq(size(xt) - 1, size(vt, 1), "regr13_conserv nt")
    n2 = assert_eq(size(vs, 2), size(vt, 2), "regr13_conserv n2")
    n3 = assert_eq(size(vs, 3), size(vt, 3), "regr13_conserv n3")

    ! Quick check on sort order:
    call assert(xs(1) < xs(2), "regr13_conserv xs bad order")
    call assert(xt(1) < xt(2), "regr13_conserv xt bad order")

    call assert(xs(1) <= xt(1) .and. xt(nt + 1) <= xs(ns + 1), &
         "regr13_conserv extrapolation")
    slope_present = present(slope)

    is = locate(xs, xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
    do it = 1, nt
       ! 1 <= is <= ns
       ! xs(is) <= xt(it) < xs(is + 1)
       if (xt(it + 1) <= xs(is + 1)) then
          vt(it, :, :) = mean_lin(xt(it), xt(it + 1))
       else
          vt(it, :, :) = mean_lin(xt(it), xs(is + 1)) * (xs(is + 1) - xt(it))
          is = is + 1
          do while (xs(is + 1) < xt(it + 1))
             ! 1 <= is <= ns - 1
             vt(it, :, :) = vt(it, :, :) + (xs(is + 1) - xs(is)) * vs(is, :, :)
             is = is + 1
          end do
          ! 1 <= is <= ns
          vt(it, :, :) = (vt(it, :, :) + mean_lin(xs(is), xt(it + 1)) &
               * (xt(it + 1) - xs(is))) / (xt(it + 1) - xt(it))
       end if

       if (xs(is + 1) == xt(it + 1)) is = is + 1
       ! 1 <= is <= ns .or. it == nt
    end do

  contains

    function mean_lin(a, b)

      ! mean in [a, b] of the linear function in [xs(is), xs(is + 1)]

      real, intent(in):: a, b
      real mean_lin(n2, n3)

      !---------------------------------------------

      if (slope_present) then
         mean_lin = slope(is, :, :) / 2. * (a + b - xs(is) - xs(is + 1)) &
              + vs(is, :, :)
      else
         mean_lin = vs(is, :, :)
      end if

    end function mean_lin

  end subroutine regr13_conserv

  !********************************************

  subroutine regr14_conserv(vs, xs, xt, vt, slope)

    ! vs and slope have rank 4.

    real, intent(in):: vs(:, :, :, :)
    real, intent(in):: xs(:)
    real, intent(in):: xt(:)
    real, intent(out):: vt(:, :, :, :)
    real, intent(in), optional:: slope(:, :, :, :)

    ! Local:
    integer is, it, ns, nt, n2, n3, n4
    logical slope_present

    !---------------------------------------------

    ns = assert_eq(size(vs, 1), size(xs) - 1, "regr14_conserv ns")
    nt = assert_eq(size(xt) - 1, size(vt, 1), "regr14_conserv nt")
    n2 = assert_eq(size(vs, 2), size(vt, 2), "regr14_conserv n2")
    n3 = assert_eq(size(vs, 3), size(vt, 3), "regr14_conserv n3")
    n4 = assert_eq(size(vs, 4), size(vt, 4), "regr14_conserv n4")

    ! Quick check on sort order:
    call assert(xs(1) < xs(2), "regr14_conserv xs bad order")
    call assert(xt(1) < xt(2), "regr14_conserv xt bad order")

    call assert(xs(1) <= xt(1) .and. xt(nt + 1) <= xs(ns + 1), &
         "regr14_conserv extrapolation")
    slope_present = present(slope)

    is = locate(xs, xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
    do it = 1, nt
       ! 1 <= is <= ns
       ! xs(is) <= xt(it) < xs(is + 1)
       if (xt(it + 1) <= xs(is + 1)) then
          vt(it, :, :, :) = mean_lin(xt(it), xt(it + 1))
       else
          vt(it, :, :, :) = mean_lin(xt(it), xs(is + 1)) * (xs(is + 1) - xt(it))
          is = is + 1
          do while (xs(is + 1) < xt(it + 1))
             ! 1 <= is <= ns - 1
             vt(it, :, :, :) = vt(it, :, :, :) + (xs(is + 1) - xs(is)) &
                  * vs(is, :, :, :)
             is = is + 1
          end do
          ! 1 <= is <= ns
          vt(it, :, :, :) = (vt(it, :, :, :) + mean_lin(xs(is), xt(it + 1)) &
               * (xt(it + 1) - xs(is))) / (xt(it + 1) - xt(it))
       end if

       if (xs(is + 1) == xt(it + 1)) is = is + 1
       ! 1 <= is <= ns .or. it == nt
    end do

  contains

    function mean_lin(a, b)

      ! mean in [a, b] of the linear function in [xs(is), xs(is + 1)]

      real, intent(in):: a, b
      real mean_lin(n2, n3, n4)

      !---------------------------------------------

      if (slope_present) then
         mean_lin = slope(is, :, :, :) / 2. * (a + b - xs(is) - xs(is + 1)) &
              + vs(is, :, :, :)
      else
         mean_lin = vs(is, :, :, :)
      end if

    end function mean_lin

  end subroutine regr14_conserv

end module regr1_conserv_m
