










module PCHFE_95_m

  implicit none

contains

  SUBROUTINE PCHFE_95(X, F, D, SKIP, XE, FE, IERR)

    ! PURPOSE  Evaluate a piecewise cubic Hermite function at an array of
    !            points.  May be used by itself for Hermite interpolation,
    !            or as an evaluator for PCHIM or PCHIC.
    ! CATEGORY  E3
    ! KEYWORDS  CUBIC HERMITE EVALUATION, HERMITE INTERPOLATION, PCHIP,
    !             PIECEWISE CUBIC EVALUATION

    !          PCHFE:  Piecewise Cubic Hermite Function Evaluator
    ! Evaluates the cubic Hermite function defined by  X, F, D  at
    ! the points  XE.

    use assert_eq_m, only: assert_eq

    REAL, intent(in):: X(:) ! real array of independent variable values
    ! The elements of X must be strictly increasing.

    REAL, intent(in):: F(:) ! real array of function values
    ! F(I) is the value corresponding to X(I).

    REAL, intent(in):: D(:) ! real array of derivative values
    ! D(I) is the value corresponding to X(I).

    LOGICAL, intent(inout):: SKIP 
    ! request to skip checks for validity of "x"
    ! If "skip" is false then "pchfe" will check that size(x) >= 2 and
    ! "x" is in strictly ascending order.
    ! Setting "skip" to true will save time in case these checks have
    ! already been performed (say, in "PCHIM" or "PCHIC").
    ! "SKIP" will be set to TRUE on normal return.

    real, intent(in):: XE(:) ! points at which the function is to be evaluated
    ! NOTES:
    ! 1. The evaluation will be most efficient if the elements of XE
    ! are increasing relative to X.
    ! That is,   XE(J) .GE. X(I)
    ! implies    XE(K) .GE. X(I),  all K.GE.J
    ! 2. If any of the XE are outside the interval [X(1),X(N)], values
    ! are extrapolated from the nearest extreme cubic, and a warning
    ! error is returned.

    real, intent(out):: FE(:) ! values of the cubic Hermite function
    ! defined by X, F, D at the points XE

    integer, intent(out):: IERR ! error flag
    ! Normal return:
    ! IERR = 0  no error
    ! Warning error:
    ! IERR > 0  means that extrapolation was performed at IERR points
    ! "Recoverable" errors:
    !              IERR = -1  if N < 2
    !              IERR = -3  if the X-array is not strictly increasing
    !              IERR = -4  if NE < 1
    ! NOTE: The above errors are checked in the order listed, and
    ! following arguments have **NOT** been validated.

    ! Variables local to the procedure:

    INTEGER  N, NE

    !---------------------------------------

    n = assert_eq(size(x), size(f), size(d), "PCHFE_95 n")
    ne = assert_eq(size(xe), size(fe), "PCHFE_95 ne")
    call PCHFE(N, X, F, D, 1, SKIP, NE, XE, FE, IERR)

  end SUBROUTINE PCHFE_95

end module PCHFE_95_m
