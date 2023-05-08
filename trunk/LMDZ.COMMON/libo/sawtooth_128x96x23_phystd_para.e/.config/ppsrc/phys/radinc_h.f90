










module radinc_h

  implicit none

  include "bands.h"
  include "scatterers.h"

!======================================================================
!
!     RADINC.H 
!
!     Includes for the radiation code; RADIATION LAYERS, LEVELS,
!     number of spectral intervals. . .
! 
!======================================================================

!     RADIATION parameters

!     In radiation code, layer 1 corresponds to the stratosphere.  Level
!     1 is the top of the stratosphere.  The dummy layer is at the same
!     temperature as the (vertically isothermal) stratosphere, and
!     any time it is explicitly needed, the appropriate quantities will
!     be dealt with (aka "top". . .)

!     L_NLEVRAD corresponds to the surface - i.e., the GCM Level that
!     is at the surface.  PLEV(L_NLEVRAD) = P(J,I)+PTROP, 
!     PLEV(2) = PTROP, PLEV(1) = ptop

!     L_NLAYRAD is the number of radiation code layers
!     L_NLEVRAD is the number of radiation code levels.  Level N is the
!               top of layer N. 
!
!     L_NSPECTI is the number of IR spectral intervals
!     L_NSPECTV is the number of Visual(or Solar) spectral intervals
!     L_NGAUSS  is the number of Gauss points for K-coefficients
!               GAUSS POINT 17 (aka the last one) is the special case
!
!     L_NPREF   is the number of reference pressures that the 
!               k-coefficients are calculated on
!     L_PINT    is the number of Lagrange interpolated reference
!               pressures for the gas k-coefficients - now for a
!		smaller p-grid than before
!     L_NTREF   is the number of reference temperatures for the
!               k-coefficients
!     L_TAUMAX  is the largest optical depth - larger ones are set
!               to this value
!
!     L_REFVAR  The number of different mixing ratio values for
!               the k-coefficients. Variable component of the mixture
!		can in princple be anything: currently it's H2O.
!
!     NAERKIND  The number of radiatively active aerosol types
!
!     NSIZEMAX  The maximum number of aerosol particle sizes
!
!----------------------------------------------------------------------

      integer,save :: L_NLAYRAD  ! = nbp_lev ! set by ini_radinc_h
      integer,save :: L_LEVELS   ! = 2*(nbp_lev-1)+3 ! set by ini_radinc_h
      integer,save :: L_NLEVRAD  ! = nbp_lev+1 ! set by ini_radinc_h
!$OMP THREADPRIVATE(L_NLAYRAD,L_LEVELS,L_NLEVRAD)

      ! These are set in sugas_corrk
      ! [uses allocatable arrays] -- AS 12/2011
      integer :: L_NPREF, L_NTREF, L_REFVAR, L_PINT, L_NGAUSS  !L_NPREF, L_NTREF, L_REFVAR, L_PINT, L_NGAUSS read by master in sugas_corrk

      integer, parameter :: L_NSPECTI = NBinfrared
      integer, parameter :: L_NSPECTV = NBvisible

!      integer, parameter :: NAERKIND  = 2 ! set in scatterers.h
      real,    parameter :: L_TAUMAX  = 35

      ! For Planck function integration: 
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Integration boundary temperatures are NTstart/NTfac and Ntstop/NTfac
      ! -- JVO 20 : Now read boundary T and integration dT as inputs in callphys.def
      !             NTstart, Nstop and NTfac then set by ini_radinc_h
      !             Smart user can adjust values depending he's running hot or cold atm
      !             Default is wide range : 30K-1500K, with 0.1K step
      !             ->  NTstart=300, Nstop=15000, NTfac=10
      integer :: NTstart, NTstop
      real*8  :: NTfac

      ! Maximum number of grain size classes for aerosol convolution:
      ! This must correspond to size of largest dataset used for aerosol 
      ! optical properties in datagcm folder.
      integer, parameter :: nsizemax = 60

      character(len=100),save :: corrkdir
!$OMP THREADPRIVATE(corrkdir)

      character(len=100),save :: banddir
!$OMP THREADPRIVATE(banddir)

contains

  subroutine ini_radinc_h(nbp_lev,tplanckmin,tplanckmax,dtplanck)
  ! Initialize module variables
  implicit none
  integer,intent(in) :: nbp_lev
  real*8, intent(in) :: tplanckmin
  real*8, intent(in) :: tplanckmax
  real*8, intent(in) :: dtplanck
  
  L_NLAYRAD = nbp_lev
  L_LEVELS  = 2*(nbp_lev-1)+3
  L_NLEVRAD = nbp_lev+1

  NTfac   = 1.D0 / dtplanck
  NTstart = int(tplanckmin * NTfac)
  NTstop  = int(tplanckmax * NTfac)
 
  end subroutine

end module radinc_h
