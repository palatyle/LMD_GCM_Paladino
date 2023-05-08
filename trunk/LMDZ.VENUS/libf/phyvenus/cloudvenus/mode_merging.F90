
!  SUBROUTINE   MERGING   Mode merging calculations
!  FUNCTION     U         For the donor mode
!  FUNCTION     ALPHA
!  FUNCTION     BETA
!  SUBROUTINE   GetN1R1
!  SUBROUTINE   GetN2R2


!*****************************************************************************
SUBROUTINE MERGING(rini,ni,sigmai,sigmaf,dM0i,dM3i,dM0f,dM3f)

  ! Example: mode 1 merge to mode 2
  ! dMi and dMf are the fraction of moment need to be 
  ! additioned to M_2 and replace for M_1

  use donnees
  IMPLICIT NONE

  ! Inputs
  real, intent(in) :: rini, ni, sigmai, sigmaf
  ! Outputs
  real, intent(out) :: dM0i,dM3i,dM0f,dM3f
  ! Local varialbes
  real :: Nt1,Rt1,Nt2,Rt2
  ! Functions
  real :: moment

  ! A mode is free to move anywhere in size space, it has a preferred position,
  ! called home base rhb1 for mode 1 and rhb2 for mode 2

  ! Calculation of the new N1,R1,N2 and R2 after mode-merging
  CALL GetN1R1(ni,rini,Nt1,Rt1,sigmai)
  CALL GetN2R2(ni,rini,Nt2,Rt2,sigmai,sigmaf)

  ! BE CAREFUL: here the moment are not normalized
  dM0i = moment(0,Nt1,Rt1,sigmai)
  dM3i = moment(3,Nt1,Rt1,sigmai)

  dM0f = moment(0,Nt2,Rt2,sigmaf)
  dM3f = moment(3,Nt2,Rt2,sigmaf)

  RETURN
END SUBROUTINE MERGING


!*****************************************************************************
FUNCTION U(rini,k,sigmai)

  !*     Computes [ln(Redge)-ln(Ri) - k*ln(sigmai)**2]/ [SQRT(2)*ln(sigmai)]
  !*     for the donor mode
  !*
  !*     INPUTS
  !*           rini     mean radius of donor mode before mode-merging
  !*           k        moment order used for computation
  !*
  !*     OUTPUTS
  !*           U        a function of the log-N distrib parameters
  !*
  !*     See also: J.B. notes about mode-merging

  use donnees
  IMPLICIT NONE

  real, intent(in) :: rini, k, sigmai
  real :: U

  U = (log(redge)-log(rini)-k*log(sigmai)**2)/(sqrt(2.)*log(sigmai))

  RETURN
END FUNCTION U


!*****************************************************************************
FUNCTION ALPHA(rini,k,sigmai)

  !*     Computes: 1 + erf(u(param,k)) for the donor mode
  !*
  !*     INPUTS
  !*           rini     Mean radius of donor mode before mode-merging
  !*           k        Moment order used for computation
  !*
  !*     OUTPUTS
  !*           ALPHA    error of the donor mode
  !*
  !*     See also: J.B. notes about mode-merging

  IMPLICIT NONE

  ! Imputs
  real, intent(in) :: rini, k, sigmai
  ! Function
  REAL ::  U, ALPHA

  ALPHA = 1.0D0 + erf(U(rini,k,sigmai))

  RETURN
END FUNCTION ALPHA


!*****************************************************************************
FUNCTION BETA(rini,k,sigmai)

  !*     Computes: 1 - erf(u(param,k)) for the donor mode
  !*
  !*     INPUTS
  !*           rini     Mean radius of donor mode before mode-merging
  !*           k        Moment order used for computation
  !*
  !*     OUTPUTS
  !*           BETA     Error of the donor mode
  !*     
  !*     See also: J.B. notes about mode-merging

  IMPLICIT NONE

  ! Imputs
  real, intent(in) :: rini, k, sigmai
  ! Function
  REAL :: U, BETA

  BETA = 1.D0 - erf(U(rini,k,sigmai))

  RETURN
END FUNCTION BETA


!*****************************************************************************
SUBROUTINE GetN1R1(ni,rini,Nt1,Rt1,sigmai)

  !*     Computes the final Mean radius and total number particles for donor mode
  !*
  !*     INPUTS
  !*           rini     Mean radius of donor mode before mode-merging
  !*           ni       Total number of particles at initial state
  !*
  !*     OUTPUTS
  !*     A tuple (N1,R1) with the total number of particles and Mean radius of mode 1
  !*           N1       Total number of particles of donor mode after mode-merging
  !*           R1       Mean radius of donor mode after mode-merging
  !*
  !*    See also : J.B. notes about mode-merging

  IMPLICIT NONE

  ! Imputs
  real, intent(in) :: ni,rini,sigmai
  ! Outputs
  real, intent(out) :: Nt1,Rt1
  ! Function
  REAL :: ALPHA
  ! Local variables
  REAL :: cn1, cr1, p, q, a, b

  p = 3.0D0
  q = 6.0D0

  a = ALPHA(rini,p,sigmai)
  b = ALPHA(rini,q,sigmai)

  cn1 = ni/2.0D0
  cr1 = rini

  Nt1 = cn1 * a / (a/b)**(p/(p-q))
  Rt1 = cr1 * (a / b)**(1.D0/(p-q))

  RETURN
END SUBROUTINE GetN1R1


!*****************************************************************************
SUBROUTINE GetN2R2(ni,rini,N2,R2,sigmai,sigmaf)

  !*     Computes the final Mean radius and total number of particles for 
  !*     the receptor mode
  !*
  !*     FUNCTION
  !*           BETA
  !*     
  !*     INPUTS
  !*           rini     Mean radius of donor mode before mode-merging
  !*           ni       Total number of particles at initial state
  !*
  !*     OUTPUTS
  !*     A tuple (N2,R2) with the total number of particles and Mean radius of mode 2
  !*           N2       Total number of particles of receptor mode after mode-merging
  !*           R2       Mean radius of receptor mode after mode-merging
  !*
  !*     See also : J.B. notes about mode-merging

  IMPLICIT NONE

  ! Imputs
  REAL :: sigmai, sigmaf
  real :: ni, rini, N2, R2
  ! Local variables
  real :: a, b, cn2, cr2, p, q
  ! Function
  REAL :: BETA

  p = 3.0D0
  q = 6.0D0

  a = BETA(rini,p,sigmai)
  b = BETA(rini,q,sigmai)

  cn2 = ni/2.D0 *exp(-0.5D0*p*q * (log(sigmai)**2 - log(sigmaf)**2))
  cr2 = rini * exp(0.5D0* (p+q) * (log(sigmai)**2 - log(sigmaf)**2))

  N2 = cn2 * a/(a/b)**(p/(p-q))
  R2 = cr2 * (a/b)**(1.D0/(p-q))

  RETURN
END SUBROUTINE GetN2R2
