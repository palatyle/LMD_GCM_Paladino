module yomlw_h
! Coefficients for the longwave radiation subroutines
use dimradmars_mod, only: nir, nabsmx, npademx
implicit none

  real,save :: at(2,2)
  real,save :: bt(2,2)
  real,save :: tref= 200.0 ! tref: temperature dependence of the absorption
  real,save :: xp(6,nir)
  real,save :: tstand = 200.0 ! reference temperature for the Planck function
  real,save :: ga(npademx,nabsmx)
  real,save :: gb(npademx,nabsmx)
  real,save :: cst_voigt(2,nabsmx)
  real,save :: gcp ! = g/cpp (set in callradite)

! Number of layers on which LTE calculations (in lw and sw) are performed 
! (Computed in nlthermeq) :
  integer,save :: nlaylte

  real,save,allocatable :: xi(:,:,:,:)
  real,save,allocatable :: xi_ground(:,:)
  real,save,allocatable :: xi_emis(:,:,:)

contains

! Allocate array (subroutine ini_yomlw_h is called by phys_state_var_init)
  subroutine ini_yomlw_h(ngrid)
    
    use dimradmars_mod, only: nuco2,nflev
    implicit none
    
    integer,intent(in) :: ngrid ! number of atmospheric columns
    
    allocate(xi(ngrid,nuco2,0:nflev+1,0:nflev+1))
    allocate(xi_ground(ngrid,nuco2))
    allocate(xi_emis(ngrid,nuco2,nflev-1))
   
    xi (:,:,:,:)=0. ! initialisation previously done in lwmain at firstcall
 
    ! initialize xp(),at(),bt(),ga(),gb(),cst_voigt()

    ! What follows was taken from old "sulw.F" routine (by Jean-Jacques 
    ! Morcrette, ECMWF, further simplified by F. Forget 01/2000),
    ! adapted to F90
    
    ! ROOTS AND WEIGHTS FOR THE 2-POINT GAUSSIAN QUADRATURE
!     DATA (RT1(IG1),IG1=1,2) / -0.577350269, +0.577350269 /
!     DATA (WG1(IG1),IG1=1,2) /  1.0        ,  1.0         /
    ! COEFFICIENTS OF THE POLYNOMIALS GIVING THE PLANCK FUNCTIONS
    xp = reshape ( (/ &
           0.63849788E+01, 0.30969419E+02, 0.44790835E+02, &
           0.52651048E+01,-0.18799237E+02, 0.92836181E+01, &
           0.26166790E+02, 0.12348011E+03, 0.17868306E+03, &
           0.33657659E+02,-0.66869343E+02, 0.21017507E+02, &
           0.11101254E+02, 0.86037325E+02, 0.25892695E+03, &
           0.35582991E+03, 0.16958020E+03,-0.41311413E+02, &
           0.47045285E+02, 0.12234377E+03, 0.61873275E+02, &
           -0.31971883E+02, 0.59168472E+01, 0.91927407E+01 &
                    /) , (/ 6,nir /) )

    ! temperature dependency of absorber amounts:
    at = reshape( (/ &
           0.694E-03, 0.272E-02, 0.275E-02, 0.178E-01  &
                   /) , (/ 2,2 /) )
    bt = reshape( (/ &
           0.328E-05, 0.298E-05,-0.705E-04,-0.163E-04  &
                   /) , (/ 2,2 /) )

    ga(1:4,1) = (/ 0.288231E-04, 0.170794E-01,-0.339714E-01, 0.000000E+00 /)
    gb(1:4,1) = (/ 0.288231E-04, 0.145426E-01, 0.543812E+00, 0.100000E+01 /)
    
    ga(1:4,2) = (/ 0.289299E-01, 0.190634E+01, 0.384061E+01, 0.000000E+00 /)
    gb(1:4,2) = (/ 0.289299E-01, 0.189485E+01, 0.600363E+01, 0.100000E+01 /)
    
    cst_voigt = reshape( (/ &
                  0.500E-02, 0.100E-01, 0.150E-01, 0.100E+00  &
                         /) , (/ 2,2 /) )
  
  end subroutine ini_yomlw_h
   
  subroutine end_yomlw_h

    implicit none

    if (allocated(xi)) deallocate(xi)
    if (allocated(xi_ground)) deallocate(xi_ground)
    if (allocated(xi_emis)) deallocate(xi_emis)

  end subroutine end_yomlw_h
 
end module yomlw_h
