MODULE donnees

  implicit none

  real, parameter :: ZERO = 2.220446049250313E-016
  real, parameter :: PI = 3.1415926536D0
  real, save :: E          ! Constant for microphys. processes


  ! THERMO
  real, parameter :: kbz = 1.3806488D-23 !m2.kg.s-2.K-2 Boltzmann constant
  real, parameter :: RGAS = 8.31441D0    !J/(mole K)    Gas constant
  real, parameter :: LSA = 6.322D5       !J/kg     Latent heat of evap of SA 
  real, parameter :: MAIR = 43.45D-3     !kg/mole  Molecular mass of venus air
  real, parameter :: MWV = 18.0153D-3    !kg/mole  Molecular mass of water
  real, parameter :: CSA = 0.4430057096D0!=MAIR/MSA   

  real, save :: PPSA, PPWV ! Partial pressure of SA and water
  real, save :: RH     ! Relative Humidity
  real, save :: KAIR   ! Thermal conductivity of air (J/m sec K)
  real, save :: KN     ! Knudsen number
  real, save :: D      ! Molecular diffusivity of the air (m2/s)
  real, save :: akn    ! Approx slip-flow correction coef
!  real, parameter :: akn = 1.591D0 ! Approx slip-flow correction coeff 


  ! SEDIMENTATION and PRODUCTION
  !> Fractal dimension of fractal aerosols.
  REAL, SAVE :: df = 3.D0  ! 3 = sphere
  !> Monomer radius (m).
  REAL, SAVE :: rmono = 6.66D-8    
  REAL, SAVE :: rpla = 6051800.D0 !m Venus radius
  !> Planet acceleration due to gravity constant (ground) (\(m.s^{-2}\)).
  REAL, SAVE :: g0 = 6.673D-11*(4.867D24/(6051800.D0)**2) 
!!$  !> Air molecules mean radius (m).   !SG   Y|°O°|Y 
  REAL, SAVE :: air_rad  = 1.75D-10


  ! MODES  (H2SO4-H2O)
  real, save :: redge  ! Virtual limits between the two modes
  real, save :: r1, r2 ! Mean raidus of the distribution (m)
  real, save :: N1, N2 ! Total number concentration (#/m3)
  real, save :: NTOT   ! Total number of droplets (mode 1 and 2)


  ! DROPLET AND GAS
  real, save :: WSA, WSAEQ, SH2SO4, RHOSA, RHOsasat

  ! Heterogeneous nucleation parameter (Franck's thesis, p.34)
  real, parameter :: desorpco2 = 6.D-20  !J   Absorp./desorp. of 1 h2o molecule at the CCN surface
  real, parameter :: surfdifco2 = 6.D-21 !J   Activation energy for surface diffusion
  real, parameter :: nusco2 = 1.D+13     !s-1 Jump frequency

  ! CONSERVATION
  real, save :: CHECKSA
  real, save :: CHECKWV

  ! H2SO4
  ! vo1h2so4 = masse 1 moléc. AS en kg / densité AS en kg/m3
  ! vo1h2so4 = m_point/rho = 1.628656574D-25/1.8302D3
  real, parameter :: vo1h2so4 = 8.898790154D-29 !m3
  real, parameter :: m0h2so4 = 1.628640074D-25  !kg       Wettability
  real, parameter :: MSA = 98.08D-3             !kg/mole  Molecular mass of SA
  real, save :: ST         ! Surface tension of sulfuric acid solution/vapor (N/m)

  real, save :: sigh2so4 ! tension superficielle de H2SO4
  real, save :: h2so4_m3 ! Number concentration of H2SO4 (m-3)

  ! RADII GRID (m)
  real, save, dimension(:), allocatable :: rad_cld,ri,rs,vol
  real, save ::  vratio




  ! FIADERO'S CORRECTION

  !! This flag enables/disables the __Fiadero__ correction alogrithm for fractal mode settling velocity
  !! computation.
  !!
  !! @bug
  !! Currently, the Fiadero correction creates instatibilities on the vertical structure. It seems to be
  !! related to the coupling between the two moments. In order to reduce the instabilities, settling
  !! velocity of moments are forced to be the same, see [[mm_globals(module):mm_wsed_m0(variable)]] and
  !! [[mm_globals(module):mm_wsed_m3(variable)]]).
  LOGICAL, SAVE          :: no_fiadero_w = .false.
  !> Minimum ratio for __Fiadero__ correction.
  !!
  !! When [[mm_globals(module):mm_no_fiadero_w(variable)]] is disabled, this variable defines the minimum
  !! value of the moment's ratio between two adjacents vertical cells to be used within the correction.
  REAL, SAVE :: fiadero_min  = 0.1D0
  !> Maximum ratio for __Fiadero__ correction.
  !!
  !! When [[mm_globals(module):mm_no_fiadero_w(variable)]] is disabled, this variable defines the maximum
  !! value of the moment's ratio between two adjacents vertical cells to be used within the correction.
  REAL, SAVE :: fiadero_max  = 10.D0

  LOGICAL, SAVE :: wsed_m0 = .false. !! Force all aerosols moments to fall at M0 settling velocity.
  LOGICAL, SAVE :: wsed_m3 = .false. !! Force all aerosols moments to fall at M3 settling velocity.



contains
  !  subroutine   build_radius_grid
  !  SUBROUTINE   logdist   Calculates the lognormal size distribution function with bins


!******************************************************************************
  SUBROUTINE build_radius_grid(nbins,rmi,rma,vratio)
!******************************************************************************
    INTEGER, INTENT(in)                                :: nbins
    REAL, INTENT(in)                     :: rmi,rma
    REAL, INTENT(out)                    :: vratio
    INTEGER :: i

    if (allocated(rad_cld)) deallocate(rad_cld)
    if (allocated(ri)) deallocate(ri)
    if (allocated(rs)) deallocate(rs)
    if (allocated(vol)) deallocate(vol)
    ALLOCATE(rad_cld(nbins),ri(nbins),rs(nbins),vol(nbins))
  
    vratio=(rma/rmi)**(3.D0/(nbins-1))
    vol(1)=4.D0/3.D0*pi*rmi**3

    DO i=1,nbins
      rad_cld(i)=rmi*vratio**(i/3.D0)
      IF (i < nbins) vol(i+1)=vol(1)*vratio**i
      ri(i) = (2.D0/(vratio+1.))**(1.D0/3.D0)*rad_cld(i)
      rs(i) = (2.D0*vratio/(vratio+1.D0))**(1.D0/3.D0)*rad_cld(i)
    ENDDO

  END SUBROUTINE build_radius_grid


!******************************************************************************
  SUBROUTINE logdist(sigma,ntotal,r0,rad,n_aer)
!******************************************************************************

    USE free_param
    implicit none
    
    ! Inputs / Outputs
    real,intent(in) :: sigma, ntotal, r0
    real,intent(out) :: n_aer(nbin)
    real,intent(in) :: rad(nbin)

    ! Local variables
    real :: lognr,dr,n_aer_erf(nbin),sqr2sig,logsig
    integer :: i

    dr = ((2.D0*vratio)/(vratio +1.D0))**(1.D0/3.D0) - (2.D0/(vratio + 1.D0))**(1.D0/3.D0)

    logsig = log(sigma)
    sqr2sig = logsig*sqrt(2.0)

    OPEN(666,FILE="logn.dat")

    DO i=1,nbin
       lognr = 1.D0 / (sqrt(2.D0*PI)*logsig*rad(i))               !m-1
       lognr = lognr * exp(-0.5D0*(log(rad(i)/r_aer)/logsig)**2)  !m-1
       n_aer(i) = lognr*dr*rad(i) !#
       n_aer_erf(i) =  0.5 * (erf(log(rs(i)/r_aer)/sqr2sig) - erf(log(ri(i)/r_aer)/sqr2sig))
       WRITE(666,'(6(ES15.7,2X))') rad(i),lognr,n_aer(i), n_aer_erf(i), rad(i)*dr, rs(i)-ri(i)
    END DO

    CLOSE(666)

!    write(*,'(a,2(ES15.7,1X))') "logdist: nAER(norm), ERF method = ", SUM(n_aer(:)),SUM(n_aer_erf)
    n_aer(:) = n_aer(:)*ntotal
!    write(*,'((a),2(ES15.7,1X))') '  -->      total', SUM(n_aer(:)), ntotal
!    write(*,'((a),ES15.7)')       '  --> rel. error', abs(SUM(n_aer(:))-ntotal)/ntotal
!    write(*,'(a)') 'Using erf METHOD (more accurate)'
    n_aer(:) = n_aer_erf(:) * ntotal

  END SUBROUTINE logdist

END MODULE donnees
