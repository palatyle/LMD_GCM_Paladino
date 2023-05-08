MODULE muphy_diag
  !! Microphysics diagnostcs vairables module.
  !! 
  !! The module contains all the variables related to output diagnostics of the microphysics.
  !! Such variables are not (and should not be) used as input in the GCM except for output data writing.
  !!
  !! The module also contains two methods:
  !! 
  !! - ini_diag_arrays(ngrid,nlayer,nices)
  !! - free_diag_arrays()

  USE tracer_h
  IMPLICIT NONE
  REAL(kind=8), ALLOCATABLE, DIMENSION(:)     :: mmd_aer_prec   !! Aerosols precipitations (both modes) (m).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:)     :: mmd_ccn_prec   !! CCN precipitations (m).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:)   :: mmd_aer_s_flux !! Spherical aerosol mass flux (\(kg.m^{-2}.s^{-1}\)).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:)   :: mmd_aer_f_flux !! Fractal aerosol mass flux (\(kg.m^{-2}.s^{-1}\)).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:)   :: mmd_ccn_flux   !! CCN mass flux (\(kg.m^{-2}.s^{-1}\)).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: mmd_ice_fluxes !! Ice sedimentation fluxes (\(kg.m^{-2}.s^{-1}\)).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: mmd_gazs_sat   !! Condensible gaz saturation ratios (--).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:)   :: mmd_ice_prec   !! Ice precipitations (m).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:)   :: mmd_rc_sph     !! Spherical mode characteristic radius (m).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:)   :: mmd_rc_fra     !! Fractal mode characteristic radius (m).
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:)   :: mmd_rc_cld     !! Cloud drop radius (m).

  !$OMP TRHEADPRIVATE(mmd_aer_prec,mmd_ccnprec,mmd_aer_s_flux,mmd_aer_f_flux,mmd_ccn_flux,mmd_ice_fluxes)
  !$OMP TRHEADPRIVATE(mmd_gazs_sat,mmd_ice_prec,mmd_rc_sph,mmd_rc_fra,mmd_rc_cld)

  CONTAINS

  SUBROUTINE ini_diag_arrays(ngrid,nlayer,nices)
    !! Initialize the variables associated to microphysics diagnostics.
    INTEGER, INTENT(in) :: ngrid  !! Number of points of the horizontal grid.
    INTEGER, INTENT(in) :: nlayer !! Number of points of the vertical grid (layers).
    INTEGER, INTENT(in) :: nices  !! Number of condensible species (cloud microphysics).
    ALLOCATE(mmd_aer_prec(ngrid))
    ALLOCATE(mmd_ccn_prec(ngrid))
    ALLOCATE(mmd_aer_s_flux(ngrid,nlayer))
    ALLOCATE(mmd_aer_f_flux(ngrid,nlayer))
    ALLOCATE(mmd_ccn_flux(ngrid,nlayer))
    ALLOCATE(mmd_ice_fluxes(ngrid,nlayer,nices))
    ALLOCATE(mmd_gazs_sat(ngrid,nlayer,nices))
    ALLOCATE(mmd_ice_prec(ngrid,nices))
    ALLOCATE(mmd_rc_sph(ngrid,nlayer))
    ALLOCATE(mmd_rc_fra(ngrid,nlayer))
    ALLOCATE(mmd_rc_cld(ngrid,nlayer))

    mmd_aer_prec(:)       = 0d0
    mmd_ccn_prec(:)       = 0d0
    mmd_aer_s_flux(:,:)   = 0d0
    mmd_aer_f_flux(:,:)   = 0d0
    mmd_ccn_flux(:,:)     = 0d0
    mmd_ice_fluxes(:,:,:) = 0d0
    mmd_gazs_sat(:,:,:)   = 0d0
    mmd_ice_prec(:,:)     = 0d0
    mmd_rc_sph(:,:)       = 0d0
    mmd_rc_fra(:,:)       = 0d0
    mmd_rc_cld(:,:)       = 0d0
      
  END SUBROUTINE ini_diag_arrays

  SUBROUTINE free_diag_arrays()
    !! Free memory of the variables associated to microphysics diagnostics.
    IF (ALLOCATED(mmd_aer_prec))   DEALLOCATE(mmd_aer_prec)
    IF (ALLOCATED(mmd_ccn_prec))   DEALLOCATE(mmd_ccn_prec)
    IF (ALLOCATED(mmd_aer_s_flux)) DEALLOCATE(mmd_aer_s_flux)
    IF (ALLOCATED(mmd_aer_f_flux)) DEALLOCATE(mmd_aer_f_flux)
    IF (ALLOCATED(mmd_ccn_flux))   DEALLOCATE(mmd_ccn_flux)
    IF (ALLOCATED(mmd_ice_fluxes)) DEALLOCATE(mmd_ice_fluxes)
    IF (ALLOCATED(mmd_gazs_sat))   DEALLOCATE(mmd_gazs_sat)
    IF (ALLOCATED(mmd_ice_prec))   DEALLOCATE(mmd_ice_prec)
    IF (ALLOCATED(mmd_rc_sph))     DEALLOCATE(mmd_rc_sph)
    IF (ALLOCATED(mmd_rc_fra))     DEALLOCATE(mmd_rc_fra)
    IF (ALLOCATED(mmd_rc_cld))     DEALLOCATE(mmd_rc_cld)
  END SUBROUTINE free_diag_arrays
END MODULE muphy_diag
