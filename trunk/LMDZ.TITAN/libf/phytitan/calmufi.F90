

SUBROUTINE calmufi(dt, plev, zlev, play, zlay, g3d, temp, pq, zdqfi, zdq)
  !! Interface subroutine to YAMMS model for Titan LMDZ GCM.
  !!
  !! The subroutine computes the microphysics processes for a single vertical column.
  !!
  !! - All input vectors are assumed to be defined from GROUND to TOP of the atmosphere.
  !! - All output vectors are defined from GROUND to TOP of the atmosphere.
  !! - Only tendencies are returned.
  !!
  !! @important
  !! The method assumes global initialization of YAMMS model (and extras) has been already
  !! done elsewhere.
  !! 
  !! Authors : J.Burgalat, J.Vatant d'Ollone - 2017
  !!
  USE MMP_GCM
  USE tracer_h
  USE callkeys_mod, only : callclouds
  USE muphy_diag
  IMPLICIT NONE

  REAL(kind=8), INTENT(IN) :: dt  !! Physics timestep (s).
  
  REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: plev  !! Pressure levels (Pa).
  REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: zlev  !! Altitude levels (m).
  REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: play  !! Pressure layers (Pa).
  REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: zlay  !! Altitude at the center of each layer (m).
  REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: g3d   !! Latitude-Altitude depending gravitational acceleration (m.s-2).
  REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: temp  !! Temperature at the center of each layer (K).

  REAL(kind=8), DIMENSION(:,:,:), INTENT(IN)  :: pq    !! Tracers (\(X.kg^{-1}}\)).
  REAL(kind=8), DIMENSION(:,:,:), INTENT(IN)  :: zdqfi !! Tendency from former processes for tracers (\(X.kg^{-1}}\)).
  REAL(kind=8), DIMENSION(:,:,:), INTENT(OUT) :: zdq   !! Microphysical tendency for tracers (\(X.kg^{-1}}\)).
  
  REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: zq !! Local tracers updated from former processes (\(X.kg^{-1}}\)).
  
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: m0as !! 0th order moment of the spherical mode (\(m^{-2}\)).
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: m3as !! 3rd order moment of the spherical mode (\(m^{3}.m^{-2}\)).
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: m0af !! 0th order moment of the fractal mode (\(m^{-2}\)).
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: m3af !! 3rd order moment of the fractal mode (\(m^{3}.m^{-2}\)).

  REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: m0n  !! 0th order moment of the CCN distribution (\(m^{-2}\)).
  REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: m3n  !! 3rd order moment of the CCN distribution (\(m^{3}.m^{-2}\)).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: m3i  !! 3rd order moments of the ice components (\(m^{3}.m^{-2}\)).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: gazs !! Condensible species gazs molar fraction (\(mol.mol^{-1}\)).

  REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: dm0as !! Tendency of the 0th order moment of the spherical mode distribution (\(m^{-2}\)).
  REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: dm3as !! Tendency of the 3rd order moment of the spherical mode distribution (\(m^{3}.m^{-2}\)). 
  REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: dm0af !! Tendency of the 0th order moment of the fractal mode distribution (\(m^{-2}\)).
  REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: dm3af !! Tendency of the 3rd order moment of the fractal mode distribution (\(m^{3}.m^{-2}\)).
  REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: dm0n  !! Tendency of the 0th order moment of the _CCN_ distribution (\(m^{-2}\)).
  REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: dm3n  !! Tendency of the 3rd order moment of the _CCN_ distribution (\(m^{3}.m^{-2}\)).
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: dm3i  !! Tendencies of the 3rd order moments of each ice components (\(m^{3}.m^{-2}\)). 
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: dgazs !! Tendencies of each condensible gaz species !(\(mol.mol^{-1}\)). 

  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE ::  int2ext !! (\(m^{-2}\)).
  TYPE(error) :: err

  INTEGER :: ilon, i,nices
  INTEGER :: nq,nlon,nlay

  ! Read size of arrays
  nq    = size(pq,DIM=3)
  nlon  = size(play,DIM=1)
  nlay  = size(play,DIM=2)
  nices = size(ices_indx)
  ! Conversion intensive to extensive
  ALLOCATE( int2ext(nlon,nlay) )  

  ! Allocate arrays
  ALLOCATE( zq(nlon,nlay,nq) )

  ALLOCATE( m0as(nlay) )
  ALLOCATE( m3as(nlay) )
  ALLOCATE( m0af(nlay) )
  ALLOCATE( m3af(nlay) )
  ALLOCATE( m0n(nlay) )
  ALLOCATE( m3n(nlay) )
  ALLOCATE( m3i(nlay,nices) )
  ALLOCATE( gazs(nlay,nices) )

  ALLOCATE( dm0as(nlay) )
  ALLOCATE( dm3as(nlay) )
  ALLOCATE( dm0af(nlay) )
  ALLOCATE( dm3af(nlay) )
  ALLOCATE( dm0n(nlay) )
  ALLOCATE( dm3n(nlay) )
  ALLOCATE( dm3i(nlay,nices) )
  ALLOCATE( dgazs(nlay,nices) )

  ! Initialization of zdq here since intent=out and no action performed on every tracers
  zdq(:,:,:) = 0.D0

  ! Initialize tracers updated with former processes from physics
  zq(:,:,:) = pq(:,:,:) + zdqfi(:,:,:)*dt

  ! Loop on horizontal grid points
  DO ilon = 1, nlon
    ! Convert tracers to extensive ( except for gazs where we work with molar mass ratio )
    ! We suppose a given order of tracers !
    int2ext(ilon,:) = ( plev(ilon,1:nlay)-plev(ilon,2:nlay+1) ) / g3d(ilon,1:nlay)
    
    ! Check because of the threshold of small tracers values in the dynamics
    ! WARNING : With this patch it enables to handles the small values required
    ! by YAMMS, but still it might still leads to some unphysical values of
    ! radii inside YAMMS, harmless for now, but who might be a problem the day
    ! you'll want to compute optics from the radii.
    WHERE (pq(ilon,:,1) > 2.D-200 .AND. pq(ilon,:,2) > 2.D-200) ! Test on both moments to avoid divergences if one hit the threshold but not the other
      m0as(:) = zq(ilon,:,1) * int2ext(ilon,:)                  ! It can still be a pb if both m0 and m3 has been set to epsilon at the beginning of dynamics
      m3as(:) = zq(ilon,:,2) * int2ext(ilon,:)                  ! then mixed, even though both are above the threshold, their ratio can be nonsense
    ELSEWHERE
      m0as(:)=0.D0
      m3as(:)=0.D0
    ENDWHERE
    WHERE (pq(ilon,:,3) > 2.D-200 .AND. pq(ilon,:,4) > 2.D-200)
      m0af(:) = zq(ilon,:,3) * int2ext(ilon,:)
      m3af(:) = zq(ilon,:,4) * int2ext(ilon,:)
    ELSEWHERE
      m0af(:)=0.D0
      m3af(:)=0.D0
    ENDWHERE
    
    if (callclouds) then ! if call clouds
      m0n(:) = zq(ilon,:,5) * int2ext(ilon,:)
      m3n(:) = zq(ilon,:,6) * int2ext(ilon,:)
      do i=1,nices
        m3i(:,nices) = zq(ilon,:,6+i) * int2ext(ilon,:)
        gazs(:,i)    = zq(ilon,:,ices_indx(i)) / rat_mmol(ices_indx(i)) ! For gazs we work on the full tracer array !!
        ! We use the molar mass ratio from GCM in case there is discrepancy with the mm one
      enddo
    endif

    ! Initialize YAMMS atmospheric column
    err = mm_column_init(plev(ilon,:),zlev(ilon,:),play(ilon,:),zlay(ilon,:),temp(ilon,:)) ; IF (err /= 0) call abort_program(err)
    ! Initialize YAMMS aerosols moments column
    err = mm_aerosols_init(m0as,m3as,m0af,m3af) ; IF (err /= 0) call abort_program(err)
    IF (callclouds) THEN ! call clouds
      err = mm_clouds_init(m0n,m3n,m3i,gazs) ; IF (err /= 0) call abort_program(err)
    ENDIF

    ! Check on size (???)

    ! initializes tendencies:
    !dm0as(:) = 0._mm_wp ; dm3as(:) = 0._mm_wp ; dm0af(:) = 0._mm_wp ; dm3af(:) = 0._mm_wp
    !dm0n(:) = 0._mm_wp ; dm3n(:) = 0._mm_wp ; dm3i(:,:) = 0._mm_wp ; dgazs(:,:) = 0._mm_wp

    dm0as(:) = 0.D0 ; dm3as(:) = 0.D0 ; dm0af(:) = 0.D0 ; dm3af(:) = 0.D0
    dm0n(:) = 0.D0 ; dm3n(:) = 0.D0 ; dm3i(:,:) = 0.D0 ; dgazs(:,:) = 0.D0
    ! call microphysics

    IF (callclouds) THEN ! call clouds
      IF(.NOT.mm_muphys(dm0as,dm3as,dm0af,dm3af,dm0n,dm3n,dm3i,dgazs)) &
        call abort_program(error("mm_muphys aborted -> initialization not done !",-1))
    ELSE
      IF (.NOT.mm_muphys(dm0as,dm3as,dm0af,dm3af)) & 
        call abort_program(error("mm_muphys aborted -> initialization not done !",-1))
    ENDIF
    ! save diags (if no clouds, relevant arrays will be set to 0 !)
    call mm_diagnostics(mmd_aer_prec(ilon),mmd_aer_s_flux(ilon,:),mmd_aer_f_flux(ilon,:),  &
                        mmd_ccn_prec(ilon),mmd_ccn_flux(ilon,:), mmd_ice_prec(ilon,:),   &
                        mmd_ice_fluxes(ilon,:,:),mmd_gazs_sat(ilon,:,:))
    call mm_get_radii(mmd_rc_sph(ilon,:),mmd_rc_fra(ilon,:),mmd_rc_cld(ilon,:))

    ! Convert tracers back to intensives ( except for gazs where we work with molar mass ratio )
    ! We suppose a given order of tracers !
 
    zdq(ilon,:,1) = dm0as(:) / int2ext(ilon,:)
    zdq(ilon,:,2) = dm3as(:) / int2ext(ilon,:)
    zdq(ilon,:,3) = dm0af(:) / int2ext(ilon,:)
    zdq(ilon,:,4) = dm3af(:) / int2ext(ilon,:)
    
    if (callclouds) then ! if call clouds
      zdq(ilon,:,5) = dm0n(:) / int2ext(ilon,:)
      zdq(ilon,:,6) = dm3n(:) / int2ext(ilon,:)
      do i=1,nices
        zdq(ilon,:,6+i) = dm3i(:,nices) / int2ext(ilon,:)
        zdq(ilon,:,ices_indx(i)) = dgazs(:,i) * rat_mmol(ices_indx(i)) ! For gazs we work on the full tracer array !!
        ! We use the molar mass ratio from GCM in case there is discrepancy with the mm one
      enddo
    endif
    
  END DO ! loop on ilon

  ! YAMMS gives a tendency which is integrated for all the timestep but in the GCM 
  ! we want to have routines spitting tendencies in s-1 -> let's divide !
  zdq(:,:,:) = zdq(:,:,:) / dt

END SUBROUTINE calmufi
