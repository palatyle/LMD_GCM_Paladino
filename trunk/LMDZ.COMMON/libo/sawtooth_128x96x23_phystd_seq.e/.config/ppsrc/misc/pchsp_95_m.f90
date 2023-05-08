










module pchsp_95_m

  implicit none

contains

  function pchsp_95(x, f, ibeg, iend, vc_beg, vc_end)

    ! PURPOSE: Set derivatives needed to determine the Hermite
    ! representation of the cubic spline interpolant to given data,
    ! with specified boundary conditions.

    ! Part of the "pchip" package.

    ! CATEGORY: E1A

    ! KEYWORDS: cubic hermite interpolation, piecewise cubic
    ! interpolation, spline interpolation

    ! DESCRIPTION: "pchsp" stands for "Piecewise Cubic Hermite Spline"
    ! Computes the Hermite representation of the cubic spline
    ! interpolant to the data given in X and F satisfying the boundary
    ! conditions specified by Ibeg, iend, vc_beg and VC_end.

    ! The resulting piecewise cubic Hermite function may be evaluated
    ! by "pchfe" or "pchfd".

    ! NOTE: This is a modified version of C. de Boor's cubic spline
    ! routine "cubspl".

    ! REFERENCE: Carl de Boor, A Practical Guide to Splines, Springer,
    ! 2001, pages 43-47

    use assert_eq_m, only: assert_eq

    real, intent(in):: x(:)
    ! independent variable values
    ! The elements of X must be strictly increasing:
    !                X(I-1) < X(I),  I = 2...N.
    !           (Error return if not.)
    ! (error if size(x) < 2)

    real, intent(in):: f(:)
    !     dependent variable values to be interpolated
    !  F(I) is value corresponding to X(I).

    INTEGER, intent(in):: ibeg
    !     desired boundary condition at beginning of data

    !        IBEG = 0  to set pchsp_95(1) so that the third derivative is con-
    !              tinuous at X(2).  This is the "not a knot" condition
    !              provided by de Boor's cubic spline routine CUBSPL.
    !              This is the default boundary condition.
    !        IBEG = 1  if first derivative at X(1) is given in VC_BEG.
    !        IBEG = 2  if second derivative at X(1) is given in VC_BEG.
    !        IBEG = 3  to use the 3-point difference formula for pchsp_95(1).
    !              (Reverts to the default boundary condition if size(x) < 3 .)
    !        IBEG = 4  to use the 4-point difference formula for pchsp_95(1).
    !              (Reverts to the default boundary condition if size(x) < 4 .)

    !          NOTES:
    !           1. An error return is taken if IBEG is out of range.
    !           2. For the "natural" boundary condition, use IBEG=2 and
    !              VC_BEG=0.

    INTEGER, intent(in):: iend
    !           IC(2) = IEND, desired condition at end of data.
    !  IEND may take on the same values as IBEG, but applied to
    !  derivative at X(N). In case IEND = 1 or 2, The value is given in VC_END.

    !          NOTES:
    !           1. An error return is taken if IEND is out of range.
    !           2. For the "natural" boundary condition, use IEND=2 and
    !              VC_END=0.

    REAL, intent(in), optional:: vc_beg
    ! desired boundary value, as indicated above.
    !           VC_BEG need be set only if IBEG = 1 or 2 .

    REAL, intent(in), optional:: vc_end
    ! desired boundary value, as indicated above.
    !           VC_END need be set only if Iend = 1 or 2 .

    real pchsp_95(size(x))
    ! derivative values at the data points
    !           These values will determine the cubic spline interpolant
    !           with the requested boundary conditions.
    !           The value corresponding to X(I) is stored in
    !                PCHSP_95(I),  I=1...N.

    ! LOCAL VARIABLES:
    real wk(2, size(x)) ! real array of working storage
    INTEGER n ! number of data points
    integer ierr, ic(2)
    real vc(2)

    !-------------------------------------------------------------------

    n = assert_eq(size(x), size(f), "pchsp_95 n")
    if ((ibeg == 1 .or. ibeg == 2) .and. .not. present(vc_beg)) then
       print *, "vc_beg required for IBEG = 1 or 2"
       stop 1
    end if
    if ((iend == 1 .or. iend == 2) .and. .not. present(vc_end)) then
       print *, "vc_end required for IEND = 1 or 2"
       stop 1
    end if
    ic = (/ibeg, iend/)
    if (present(vc_beg)) vc(1) = vc_beg
    if (present(vc_end)) vc(2) = vc_end
    call PCHSP(IC, VC, N, X, F, pchsp_95, 1, WK, size(WK), IERR)
    if (ierr /= 0) stop 1

  END function pchsp_95

end module pchsp_95_m
