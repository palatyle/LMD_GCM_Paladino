!IDEAL:MODEL_LAYER:INITIALIZATION
!

!  This MODULE holds the routines which are used to perform various initializations
!  for the individual domains.  

!  This MODULE CONTAINS the following routines:

!  initialize_field_test - 1. Set different fields to different constant
!                             values.  This is only a test.  If the correct
!                             domain is not found (based upon the "id")
!                             then a fatal error is issued.               

!-----------------------------------------------------------------------

MODULE module_initialize

   USE module_domain
   USE module_io_domain
   USE module_state_description
   USE module_model_constants
   USE module_bc
   USE module_timing
   USE module_configure
   USE module_init_utilities
#ifdef DM_PARALLEL
   USE module_dm
#endif


CONTAINS


!-------------------------------------------------------------------
! this is a wrapper for the solver-specific init_domain routines.
! Also dereferences the grid variables and passes them down as arguments.
! This is crucial, since the lower level routines may do message passing
! and this will get fouled up on machines that insist on passing down
! copies of assumed-shape arrays (by passing down as arguments, the 
! data are treated as assumed-size -- ie. f77 -- arrays and the copying
! business is avoided).  Fie on the F90 designers.  Fie and a pox.

   SUBROUTINE init_domain ( grid )

   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid 
   !  Local data.
   INTEGER                :: dyn_opt 
   INTEGER :: idum1, idum2

   CALL nl_get_dyn_opt( 1,dyn_opt )
   
   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )

   IF (      dyn_opt .eq. 1 &
        .or. dyn_opt .eq. 2 &
        .or. dyn_opt .eq. 3 &
                                       ) THEN
     CALL init_domain_rk( grid &
!
#include <em_actual_new_args.inc>
!
                        )

   ELSE
     WRITE(0,*)' init_domain: unknown or unimplemented dyn_opt = ',dyn_opt
     CALL wrf_error_fatal ( ' init_domain: unknown or unimplemented dyn_opt ' )
   ENDIF

   END SUBROUTINE init_domain

!-------------------------------------------------------------------

   SUBROUTINE init_domain_rk ( grid &
!
# include <em_dummy_new_args.inc>
!
)
   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid

# include <em_dummy_new_decl.inc>

   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k

   ! Local data
!****Mars
    REAL    :: x_param,y_param,rho_param,dilat
    REAL    :: mulu, mulv, addu, addv
!****Mars
   INTEGER, PARAMETER :: nl_max = 1000
   REAL, DIMENSION(nl_max) :: zk, p_in, theta, tk, rho, u, v, qv, pd_in
   INTEGER :: nl_in


   INTEGER :: icm,jcm, ii, im1, jj, jm1, loop, error, fid, nxc, nyc
   REAL    :: u_mean,v_mean, f0, p_surf, p_level, qvf, z_at_v, z_at_u
   REAL    :: xrad, yrad, zrad, rad, delt, cof1, cof2
!   REAL, EXTERNAL :: interp_0
   REAL    :: hm, xa
   REAL    :: pi

!  stuff from original initialization that has been dropped from the Registry 
   REAL    :: vnu, xnu, xnus, dinit0, cbh, p0_temp, t0_temp, zd, zt
   REAL    :: qvf1, qvf2, pd_surf
   INTEGER :: it
   real :: thtmp, ptmp, temp(3)

   LOGICAL :: moisture_init
   LOGICAL :: stretch_grid, dry_sounding

  INTEGER :: xs , xe , ys , ye
  REAL :: mtn_ht
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor

!!MARS
 REAL :: lon_input, lat_input, alt_input, tsurf_input
 INTEGER :: ierr
!!MARS
 REAL, DIMENSION(nl_max) :: profdustq,profdustn
 REAL, DIMENSION(nl_max) :: prescribed_sw,prescribed_lw

      REAL :: pfu, pfd, phm
      INTEGER :: hypsometric_opt = 1 ! classic
      !INTEGER :: hypsometric_opt = 2 ! Wee et al. 2012 correction

#ifdef DM_PARALLEL
#    include <em_data_calls.inc>
#endif

   call init_module_model_constants

   SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
   kds = grid%sd31 ; kde = grid%ed31 ;
   ids = grid%sd32 ; ide = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   kms = grid%sm31 ; kme = grid%em31 ;
   ims = grid%sm32 ; ime = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   kts = grid%sp31 ; kte = grid%ep31 ;   ! note that tile is entire patch
   its = grid%sp32 ; ite = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XYZ )
   ids = grid%sd31 ; ide = grid%ed31 ;
   jds = grid%sd32 ; jde = grid%ed32 ;
   kds = grid%sd33 ; kde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   jms = grid%sm32 ; jme = grid%em32 ;
   kms = grid%sm33 ; kme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   jts = grid%sp32 ; jte = grid%ep32 ;   ! note that tile is entire patch
   kts = grid%sp33 ; kte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XZY )
   ids = grid%sd31 ; ide = grid%ed31 ;
   kds = grid%sd32 ; kde = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   kms = grid%sm32 ; kme = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   kts = grid%sp32 ; kte = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch

   END SELECT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!MARS : NOFILE      === no mountain
!!MARS : FILE xa=0.  === a linear slope with elevation hm
!!MARS : FILE        === mountain height hm, width xa
  open(unit=22,file='ze_hill',form='formatted',status='old',iostat=ierr)
  IF (ierr .eq. 0) THEN
    rewind(22) 
    read(22,*) hm, xa
    write(6,*) 'hm, xa ', hm, xa
    close(22)
  ENDIF
!!MARS
!!MARS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   icm = ide/2
!****Mars
   jcm = jde/2
!
!   xa1  = 5000./500.
!   xal1 = 4000./500.
!   pii  = 2.*asin(1.0)
!   hm1  = 250.
!!   hm1  = 1000.
!****Mars

   delt = 3.
!   delt = 10.

!****Mars
   stretch_grid = .true.
!   stretch_grid = .false.
!****Mars
!   z_scale = .50
   z_scale = .40
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

!!!MARS
!!!MARS
!  open(unit=16,file='input_vert',form='formatted',status='old')
!  rewind(16) 
!  read(16,*) delt, z_scale
!  write(6,*) 'delt, z_scale are  ', delt, z_scale
!  close(16)
!!!MARS
!!!MARS

   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

! here we check to see if the boundary conditions are set properly

   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

   moisture_init = .true.

    grid%itimestep=0

#ifdef DM_PARALLEL
   CALL wrf_dm_bcast_bytes( icm , IWORDSIZE )
   CALL wrf_dm_bcast_bytes( jcm , IWORDSIZE )
#endif

    CALL nl_set_mminlu(1,'    ')
    CALL nl_set_iswater(1,0)
    CALL nl_set_cen_lat(1,40.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_map_proj(1,0)


!  here we initialize data we currently is not initialized 
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msft(i,j)     = 1.
         grid%msfu(i,j)     = 1.
         grid%msfv(i,j)     = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
         grid%f(i,j)        = 0.

      END DO
   END DO

    DO j = jts, jte
    DO k = kts, kte
      DO i = its, ite
         grid%em_ww(i,k,j)     = 0.
      END DO
   END DO
   END DO

   grid%step_number = 0

!! set up the grid
!
!   IF (stretch_grid) THEN ! exponential stretch for eta (nearly constant dz)
!     DO k=1, kde
!      grid%em_znw(k) = (exp(-(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
!                                (1.-exp(-1./z_scale))
!     ENDDO
!   ELSE
!     DO k=1, kde
!      grid%em_znw(k) = 1. - float(k-1)/float(kde-1)
!     ENDDO
!   ENDIF

!!MARS
!!MARS
  open(unit=12,file='levels',form='formatted',status='old')
  rewind(12)
  DO k=1, kde
  read(12,*) grid%em_znw(k)
  write(6,*) 'read level ', k,grid%em_znw(k)
  ENDDO
  close(12)
!!MARS
!!MARS
   

   DO k=1, kde-1
    grid%em_dnw(k) = grid%em_znw(k+1) - grid%em_znw(k)
    grid%em_rdnw(k) = 1./grid%em_dnw(k)
    grid%em_znu(k) = 0.5*(grid%em_znw(k+1)+grid%em_znw(k))
   ENDDO
   DO k=2, kde-1
    grid%em_dn(k) = 0.5*(grid%em_dnw(k)+grid%em_dnw(k-1))
    grid%em_rdn(k) = 1./grid%em_dn(k)
    grid%em_fnp(k) = .5* grid%em_dnw(k  )/grid%em_dn(k)
    grid%em_fnm(k) = .5* grid%em_dnw(k-1)/grid%em_dn(k)
   ENDDO

   cof1 = (2.*grid%em_dn(2)+grid%em_dn(3))/(grid%em_dn(2)+grid%em_dn(3))*grid%em_dnw(1)/grid%em_dn(2) 
   cof2 =     grid%em_dn(2)        /(grid%em_dn(2)+grid%em_dn(3))*grid%em_dnw(1)/grid%em_dn(3) 
   grid%cf1  = grid%em_fnp(2) + cof1
   grid%cf2  = grid%em_fnm(2) - cof1 - cof2
   grid%cf3  = cof2       

   grid%cfn  = (.5*grid%em_dnw(kde-1)+grid%em_dn(kde-1))/grid%em_dn(kde-1)
   grid%cfn1 = -.5*grid%em_dnw(kde-1)/grid%em_dn(kde-1)
   grid%rdx = 1./config_flags%dx
   grid%rdy = 1./config_flags%dy

!  get the sounding from the ascii sounding file, first get dry sounding and 
!  calculate base state

     !!!!
     !!!! user-modified wind speed
     !!!!
     mulu = 1.  !! default 
     mulv = 1.  !! default
     addu = 0.  !! default
     addv = 0.  !! default
     IF (config_flags%init_MU .ne. 0.) mulu = config_flags%init_MU
     IF (config_flags%init_MV .ne. 0.) mulv = config_flags%init_MV
     IF (config_flags%init_U  .ne. 0.) addu = config_flags%init_U
     IF (config_flags%init_V  .ne. 0.) addv = config_flags%init_V
     write(6,*) ' coeff for winds: ', mulu, mulv, addu, addv

  dry_sounding = .true.
  IF ( wrf_dm_on_monitor() ) THEN
  write(6,*) ' getting dry sounding for base state '

  CALL get_sounding( zk, p_in, pd_in, theta, tk, rho, u, v, qv, dry_sounding, nl_max, nl_in, &
                        mulu, mulv, addu, addv )
  ENDIF
  CALL wrf_dm_bcast_real( zk , nl_max )
  CALL wrf_dm_bcast_real( p_in , nl_max )
  CALL wrf_dm_bcast_real( pd_in , nl_max )
  CALL wrf_dm_bcast_real( theta , nl_max )
  CALL wrf_dm_bcast_real( tk , nl_max )
  CALL wrf_dm_bcast_real( rho , nl_max )
  CALL wrf_dm_bcast_real( u , nl_max )
  CALL wrf_dm_bcast_real( v , nl_max )
  CALL wrf_dm_bcast_real( qv , nl_max )
  CALL wrf_dm_bcast_integer ( nl_in , 1 ) 

  write(6,*) ' returned from reading sounding, nl_in is ',nl_in

!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure

  grid%p_top = interp_0( p_in, zk, config_flags%ztop, nl_in )

!!MARS
!!MARS
  open(unit=14,file='input_coord',form='formatted',status='old')
  rewind(14)
  read(14,*) lon_input
  read(14,*) lat_input
  close(14)
  write(6,*) ' lon is ',lon_input
  write(6,*) ' lat is ',lat_input
!!MARS
!!MARS

!!MARS
!!MARS
  open(unit=18,file='input_more',form='formatted',status='old')
  rewind(18)
  read(18,*) alt_input, tsurf_input
  close(18)
  write(6,*) ' alt is ',alt_input
  write(6,*) ' tsurf is ',tsurf_input
!!MARS
!!MARS

  DO j=jts,jte
  DO i=its,ite
!!MARS
IF (ierr .eq. 0) THEN
  write(6,*) ' IDEALIZED TOPOGRAPHY '
  IF (xa .ne. 0.) THEN
      !!!2D hill
      !grid%ht(i,j) = alt_input + hm/(1.+(float(i-icm)/xa)**2)
      !    grid%ht(i,j) = hm1*exp(-(( float(i-icm)/xa1)**2))   &
      !               *( (cos(pii*float(i-icm)/xal1))**2 )
    IF (hm .gt. 0.) THEN
      write(6,*) '3D hill. height, width: ',hm,xa
      write(6,*) 'input sounding is out of the mountain'
      grid%ht(i,j) = alt_input + hm/(1.+(float(i-icm)/xa)**2+(float(j-jcm)/xa)**2)
    ELSE IF (hm .lt. 0.) THEN
      write(6,*) '3D crater. height, width: ',hm,xa
      write(6,*) 'input sounding is at the bottom of crater'
      grid%ht(i,j) = (alt_input - hm) + hm/(1.+(float(i-icm)/xa)**2+(float(j-jcm)/xa)**2)
      !! AS: cannot use same formula as hill because it would force ideal.exe to extrapolate
      !!       which is not possible given how interp_0 is written
    ELSE
      write(6,*) 'Nothing. Height is 0. Flat topography'
      grid%ht(i,j) = alt_input
    ENDIF
  ELSE
    write(6,*) 'linear slope '
    write(6,*) 'height ',hm
    IF (hm .gt. 0.) THEN  
      grid%ht(i,j) = alt_input + hm * float(i)  
    ELSE IF (hm .lt. 0.) THEN
      !! see above, crater case
      grid%ht(i,j) = (alt_input - hm) + hm * float(i)
    ELSE
      write(6,*) 'Nothing. Height is 0. Flat topography'
      grid%ht(i,j) = alt_input
    ENDIF
  ENDIF
      !!!3D crater    
      !!        grid%ht(i,j) = hm - hm/(1.+(float(i-icm)/xa)**2+(float(j-jcm)/xa)**2)
      !!3D crater w/ rims
      !  x_param = float(i-icm)
      !  y_param = float(j-jcm)
      !  dilat = xa/2
      !  rho_param = sqrt(x_param**2 + y_param**2)
      !  ! revolution surface ; seed is a fourth order polynom
      !  grid%ht(i,j) = (rho_param+6*dilat)*(rho_param+10*dilat)
      !  grid%ht(i,j) = (rho_param-6*dilat)*(rho_param-10*dilat)*grid%ht(i,j)
      !  ! flat terrain elsewhere - smooth gradient (no abrupt fall)
      !  grid%ht(i,j) = grid%ht(i,j)*(tanh(rho_param+7*dilat)/2 - tanh(rho_param-7*dilat)/2)
      !  grid%ht(i,j) = hm - (hm*.4/1500)*grid%ht(i,j)/(dilat**4)
      !           !NONONONONON  
      !           !grid%ht(i,j) = grid%ht(i,j) + alt_input  
      !   !if (rho_param .GE. dilat*10) ht(i,j) = hm
ELSE
    write(6,*) ' FLAT SURFACE '
    grid%ht(i,j) = alt_input
    !grid%ht(i,j) = 0.
ENDIF
    grid%tsk(i,j) = tsurf_input
    grid%m_tsurf(i,j) = tsurf_input
!!MARS
    grid%xlat(i,j) = lat_input
    grid%xlong(i,j) = lon_input!+float(i)*config_flags%dx/59000.
    grid%m_emiss(i,j)=0.95
    grid%m_co2ice(i,j)=0.
    grid%m_h2oice(i,j)=0.
!! >> Used for restarts only:
    grid%m_q2(i,:,j)=0.
    grid%m_fluxrad(i,j)=0.
    grid%m_wstar(i,j)=0.
!! <<
    write(6,*) 'NOTE TO SELF. slpx and slpy set to 0 which means no slope insolation.'
    grid%slpx(i,j) = 0.
    grid%slpy(i,j) = 0.
   DO k=1,config_flags%num_soil_layers
    grid%m_tsoil(i,k,j) = 0.
   ENDDO
    !!! COMMENT THE LINES BELOW IF YOU DON'T WANT CORIOLIS TERMS
    !!! cf. doc WRF2008 page 11 for e and f expressions
    grid%e(i,j) = 2. * EOMEG * COS(pi*lat_input/180.)
    grid%f(i,j) = 2. * EOMEG * SIN(pi*lat_input/180.)
    write(6,*) 'CALCULATE CORIOLIS TERM', grid%f(i,j),grid%e(i,j)
!!MARS
  ENDDO
  ENDDO

  xs=ide/2 -3
  xs=ids   -3
  xe=xs + 6
  ys=jde/2 -3
  ye=ys + 6
  mtn_ht = 500
#ifdef MTN
  DO j=max(ys,jds),min(ye,jde-1)
  DO i=max(xs,ids),min(xe,ide-1)
     grid%ht(i,j) = mtn_ht * 0.25 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) ) * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )
  ENDDO
  ENDDO
#endif
#ifdef EW_RIDGE
  DO j=max(ys,jds),min(ye,jde-1)
  DO i=ids,ide
     grid%ht(i,j) = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )
  ENDDO
  ENDDO
#endif
#ifdef NS_RIDGE
  DO j=jds,jde
  DO i=max(xs,ids),min(xe,ide-1)
     grid%ht(i,j) = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) )
  ENDDO
  ENDDO
#endif
  DO j=jts,jte
  DO i=its,ite
    grid%em_phb(i,1,j) = g * grid%ht(i,j)
    grid%em_ph0(i,1,j) = g * grid%ht(i,j)
  ENDDO
  ENDDO

!!!dans hill_2d OK
!    grid%em_phb(i,1,j) = g*grid%ht(i,j)
!    grid%em_php(i,1,j) = 0.
!    grid%em_ph0(i,1,j) = grid%em_phb(i,1,j)

  DO J = jts, jte
  DO I = its, ite

    p_surf = interp_0( p_in, zk, grid%em_phb(i,1,j)/g, nl_in )
    grid%em_mub(i,j) = p_surf-grid%p_top

!  this is dry hydrostatic sounding (base state), so given grid%em_p (coordinate),
!  interp theta (from interp) and compute 1/rho from eqn. of state

    DO K = 1, kte-1
      p_level = grid%em_znu(k)*(p_surf - grid%p_top) + grid%p_top
      grid%em_pb(i,k,j) = p_level

! OLD METHOD
!      grid%em_t_init(i,k,j) = interp_0( theta, p_in, p_level, nl_in ) - t0
! NEW METHOD: Wee et al. 2012
! interpolate temperature. then convert to potential temperature.
      grid%em_t_init(i,k,j) = interp_0( tk, p_in, p_level, nl_in )
      !! l un ou l autre pareil
      grid%em_t_init(i,k,j) = - t0 + (grid%em_t_init(i,k,j) * ((p1000mb/p_level)**rcp))
      !grid%em_t_init(i,k,j) = - t0 + (grid%em_t_init(i,k,j) * ((610./p_level)**(1.0/3.9)))
      grid%em_alb(i,k,j) = (r_d/p1000mb)*(grid%em_t_init(i,k,j)+t0)*(grid%em_pb(i,k,j)/p1000mb)**cvpm
    ENDDO

!  calc hydrostatic balance (alternatively we could interp the geopotential from the
!  sounding, but this assures that the base state is in exact hydrostatic balance with
!  respect to the model eqns.

   IF (hypsometric_opt == 1) THEN
    DO k  = 2,kte
      grid%em_phb(i,k,j) = grid%em_phb(i,k-1,j) - grid%em_dnw(k-1)*grid%em_mub(i,j)*grid%em_alb(i,k-1,j)
    ENDDO
   ELSE IF (hypsometric_opt == 2) THEN
    DO k = 2,kte
      pfu = grid%em_mub(i,j)*grid%em_znw(k)   + grid%p_top
      pfd = grid%em_mub(i,j)*grid%em_znw(k-1)   + grid%p_top
      phm = grid%em_mub(i,j)*grid%em_znu(k-1)   + grid%p_top
      grid%em_phb(i,k,j) = grid%em_phb(i,k-1,j) + grid%em_alb(i,k-1,j)*phm*LOG(pfd/pfu)
    END DO
   END IF


  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
    write(6,*) ' ptop is ',grid%p_top
    write(6,*) ' base state grid%em_mub(1,1), p_surf is ',grid%em_mub(1,1),grid%em_mub(1,1)+grid%p_top
  ENDIF

!  calculate full state for each column - this includes moisture.

!!!!!MARS MARS
!  write(6,*) ' getting moist sounding for full state '
!  dry_sounding = .false.
  dry_sounding = .true.
  CALL get_sounding( zk, p_in, pd_in, theta, tk, rho, u, v, qv, dry_sounding, nl_max, nl_in, &
                              mulu, mulv, addu, addv )

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)

!  At this point grid%p_top is already set. find the DRY mass in the column 
!  by interpolating the DRY pressure.  

   pd_surf = interp_0( pd_in, zk, grid%em_phb(i,1,j)/g, nl_in )

!  compute the perturbation mass and the full mass

    grid%em_mu_1(i,j) = pd_surf-grid%p_top - grid%em_mub(i,j)
    grid%em_mu_2(i,j) = grid%em_mu_1(i,j)
    grid%em_mu0(i,j) = grid%em_mu_1(i,j) + grid%em_mub(i,j)

! given the dry pressure and coordinate system, interp the potential
! temperature and qv

    do k=1,kde-1

      p_level = grid%em_znu(k)*(pd_surf - grid%p_top) + grid%p_top

      moist(i,k,j,P_QV) = interp_0( qv, pd_in, p_level, nl_in )
      grid%em_t_1(i,k,j)          = interp_0( theta, pd_in, p_level, nl_in ) - t0
      grid%em_t_2(i,k,j)          = grid%em_t_1(i,k,j)
      

    enddo

!  integrate the hydrostatic equation (from the RHS of the bigstep
!  vertical momentum equation) down from the top to get grid%em_p.
!  first from the top of the model to the top pressure

    k = kte-1  ! top level

    qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k,j,P_QV))
    qvf2 = 1./(1.+qvf1)
    qvf1 = qvf1*qvf2

!    grid%em_p(i,k,j) = - 0.5*grid%em_mu_1(i,j)/grid%em_rdnw(k)
    grid%em_p(i,k,j) = - 0.5*(grid%em_mu_1(i,j)+qvf1*grid%em_mub(i,j))/grid%em_rdnw(k)/qvf2
    qvf = 1. + rvovrd*moist(i,k,j,P_QV)
    grid%em_alt(i,k,j) = (r_d/p1000mb)*(grid%em_t_1(i,k,j)+t0)*qvf* &
                (((grid%em_p(i,k,j)+grid%em_pb(i,k,j))/p1000mb)**cvpm)
    grid%em_al(i,k,j) = grid%em_alt(i,k,j) - grid%em_alb(i,k,j)

!  down the column

    do k=kte-2,1,-1
      qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k+1,j,P_QV))
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2
      grid%em_p(i,k,j) = grid%em_p(i,k+1,j) - (grid%em_mu_1(i,j) + qvf1*grid%em_mub(i,j))/qvf2/grid%em_rdn(k+1)
      qvf = 1. + rvovrd*moist(i,k,j,P_QV)
      grid%em_alt(i,k,j) = (r_d/p1000mb)*(grid%em_t_1(i,k,j)+t0)*qvf* &
                  (((grid%em_p(i,k,j)+grid%em_pb(i,k,j))/p1000mb)**cvpm)
      grid%em_al(i,k,j) = grid%em_alt(i,k,j) - grid%em_alb(i,k,j)
    enddo

!  this is the hydrostatic equation used in the model after the
!  small timesteps.  In the model, grid%em_al (inverse density)
!  is computed from the geopotential.


    grid%em_ph_1(i,1,j) = 0.
   IF (hypsometric_opt == 1) THEN
    DO k  = 2,kte
      grid%em_ph_1(i,k,j) = grid%em_ph_1(i,k-1,j) - (1./grid%em_rdnw(k-1))*(       &
                   (grid%em_mub(i,j)+grid%em_mu_1(i,j))*grid%em_al(i,k-1,j)+ &
                    grid%em_mu_1(i,j)*grid%em_alb(i,k-1,j)  )
                                                   
      grid%em_ph_2(i,k,j) = grid%em_ph_1(i,k,j) 
      grid%em_ph0(i,k,j) = grid%em_ph_1(i,k,j) + grid%em_phb(i,k,j)
    ENDDO
   ELSE IF (hypsometric_opt == 2) THEN

             ! Alternative hydrostatic eq.: dZ = -al*p*dLOG(p), where p is dry pressure.
             ! Note that al*p approximates Rd*T and dLOG(p) does z.
             ! Here T varies mostly linear with z, the first-order integration produces better result.

               grid%em_ph_2(i,1,j) = grid%em_phb(i,1,j)
               DO k = 2,kte
                  pfu = grid%em_mu0(i,j)*grid%em_znw(k)   + grid%p_top
                  pfd = grid%em_mu0(i,j)*grid%em_znw(k-1) + grid%p_top
                  phm = grid%em_mu0(i,j)*grid%em_znu(k-1) + grid%p_top
                  grid%em_ph_2(i,k,j) = grid%em_ph_2(i,k-1,j) + grid%em_alt(i,k-1,j)*phm*LOG(pfd/pfu)
               END DO

               DO k = 1,kte
                  grid%em_ph_2(i,k,j) = grid%em_ph_2(i,k,j) - grid%em_phb(i,k,j)
                  grid%em_ph_1(i,k,j) = grid%em_ph_2(i,k,j)
               END DO

   END IF



    IF ( wrf_dm_on_monitor() ) THEN
    if((i==2) .and. (j==2)) then
     write(6,*) ' grid%em_ph_1 calc ',grid%em_ph_1(2,1,2),grid%em_ph_1(2,2,2),&
                              grid%em_mu_1(2,2)+grid%em_mub(2,2),grid%em_mu_1(2,2), &
                              grid%em_alb(2,1,2),grid%em_al(1,2,1),grid%em_rdnw(1)
    endif
    ENDIF

  ENDDO
  ENDDO

!#if 0

!  thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

  DO J = jts, min(jde-1,jte)
    yrad = config_flags%dy*float(j-nyc)/10000.
yrad = 0.
    DO I = its, min(ide-1,ite)
      xrad = config_flags%dx*float(i-nxc)/10000.
xrad = 0.
      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%em_ph_1(i,k,j)+grid%em_ph_1(i,k+1,j)  &
                   +grid%em_phb(i,k,j)+grid%em_phb(i,k+1,j))/g
        zrad = (zrad-1500.)/1500.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
           grid%em_t_1(i,k,j)=grid%em_t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%em_t_2(i,k,j)=grid%em_t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%em_alt(i,k,j) = (r_d/p1000mb)*(grid%em_t_1(i,k,j)+t0)*qvf* &
                        (((grid%em_p(i,k,j)+grid%em_pb(i,k,j))/p1000mb)**cvpm)
           grid%em_al(i,k,j) = grid%em_alt(i,k,j) - grid%em_alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

   IF (hypsometric_opt == 1) THEN

      DO k  = 2,kte
        grid%em_ph_1(i,k,j) = grid%em_ph_1(i,k-1,j) - (1./grid%em_rdnw(k-1))*(       &
                     (grid%em_mub(i,j)+grid%em_mu_1(i,j))*grid%em_al(i,k-1,j)+ &
                      grid%em_mu_1(i,j)*grid%em_alb(i,k-1,j)  )
                                                   
        grid%em_ph_2(i,k,j) = grid%em_ph_1(i,k,j) 
        grid%em_ph0(i,k,j) = grid%em_ph_1(i,k,j) + grid%em_phb(i,k,j)
      ENDDO

   ELSE IF (hypsometric_opt == 2) THEN

             ! Alternative hydrostatic eq.: dZ = -al*p*dLOG(p), where p is dry pressure.
             ! Note that al*p approximates Rd*T and dLOG(p) does z.
             ! Here T varies mostly linear with z, the first-order integration produces better result.

               grid%em_ph_2(i,1,j) = grid%em_phb(i,1,j)
               DO k = 2,kte
                  pfu = grid%em_mu0(i,j)*grid%em_znw(k)   + grid%p_top
                  pfd = grid%em_mu0(i,j)*grid%em_znw(k-1) + grid%p_top
                  phm = grid%em_mu0(i,j)*grid%em_znu(k-1) + grid%p_top
                  grid%em_ph_2(i,k,j) = grid%em_ph_2(i,k-1,j) + grid%em_alt(i,k-1,j)*phm*LOG(pfd/pfu)
               END DO

               DO k = 1,kte
                  grid%em_ph_2(i,k,j) = grid%em_ph_2(i,k,j) - grid%em_phb(i,k,j)
                  grid%em_ph_1(i,k,j) = grid%em_ph_2(i,k,j)
               END DO

   END IF


    ENDDO
  ENDDO

!#endif

   IF ( wrf_dm_on_monitor() ) THEN
   write(6,*) ' grid%em_mu_1 from comp ', grid%em_mu_1(1,1)
   write(6,*) ' full state sounding from comp, ph, grid%em_p, grid%em_al, grid%em_t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%em_ph_1(1,k,1)+grid%em_phb(1,k,1), &
                                      grid%em_p(1,k,1)+grid%em_pb(1,k,1), grid%em_alt(1,k,1), &
                                      grid%em_t_1(1,k,1)+t0, moist(1,k,1,P_QV)
   enddo

   write(6,*) ' pert state sounding from comp, grid%em_ph_1, pp, alp, grid%em_t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%em_ph_1(1,k,1), &
                                      grid%em_p(1,k,1), grid%em_al(1,k,1), &
                                      grid%em_t_1(1,k,1), moist(1,k,1,P_QV)
   enddo
   ENDIF

! interp v

  DO J = jts, jte
  DO I = its, min(ide-1,ite)

    IF (j == jds) THEN
      z_at_v = grid%em_phb(i,1,j)/g
    ELSE IF (j == jde) THEN
      z_at_v = grid%em_phb(i,1,j-1)/g
    ELSE
      z_at_v = 0.5*(grid%em_phb(i,1,j)+grid%em_phb(i,1,j-1))/g
    END IF

    p_surf = interp_0( p_in, zk, z_at_v, nl_in )

    DO K = 1, kte-1
      p_level = grid%em_znu(k)*(p_surf - grid%p_top) + grid%p_top
      grid%em_v_1(i,k,j) = interp_0( v, p_in, p_level, nl_in )
      grid%em_v_2(i,k,j) = grid%em_v_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

! interp u

  DO J = jts, min(jde-1,jte)
  DO I = its, ite

    IF (i == ids) THEN
      z_at_u = grid%em_phb(i,1,j)/g
    ELSE IF (i == ide) THEN
      z_at_u = grid%em_phb(i-1,1,j)/g
    ELSE
      z_at_u = 0.5*(grid%em_phb(i,1,j)+grid%em_phb(i-1,1,j))/g
    END IF

    p_surf = interp_0( p_in, zk, z_at_u, nl_in )

    DO K = 1, kte-1
      p_level = grid%em_znu(k)*(p_surf - grid%p_top) + grid%p_top
      grid%em_u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
      grid%em_u_2(i,k,j) = grid%em_u_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

!  set w

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte
  DO I = its, min(ide-1,ite)
    grid%em_w_1(i,k,j) = 0.
    grid%em_w_2(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

!  set a few more things

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%h_diabatic(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
  DO k=1,kte-1
    grid%em_t_base(k) = grid%em_t_1(1,k,1)
    grid%qv_base(k) = moist(1,k,1,P_QV)
    grid%u_base(k) = grid%em_u_1(1,k,1)
    grid%v_base(k) = grid%em_v_1(1,k,1)
    grid%z_base(k) = 0.5*(grid%em_phb(1,k,1)+grid%em_phb(1,k+1,1)+grid%em_ph_1(1,k,1)+grid%em_ph_1(1,k+1,1))/g
  ENDDO
  ENDIF
  CALL wrf_dm_bcast_real( grid%em_t_base , kte )
  CALL wrf_dm_bcast_real( grid%qv_base , kte )
  CALL wrf_dm_bcast_real( grid%u_base , kte )
  CALL wrf_dm_bcast_real( grid%v_base , kte )
  CALL wrf_dm_bcast_real( grid%z_base , kte )

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
     thtmp   = grid%em_t_2(i,1,j)+t0
     ptmp    = grid%em_p(i,1,j)+grid%em_pb(i,1,j)
     temp(1) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%em_t_2(i,2,j)+t0
     ptmp    = grid%em_p(i,2,j)+grid%em_pb(i,2,j)
     temp(2) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%em_t_2(i,3,j)+t0
     ptmp    = grid%em_p(i,3,j)+grid%em_pb(i,3,j)
     temp(3) = thtmp * (ptmp/p1000mb)**rcp

!!MARS
!     grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)
     grid%tmn(I,J)=grid%tsk(I,J)-0.5
!!!MARS
!!TODO: passer la valeur a partir des donnees
!grid%mars_tsoil(I,:,J)=grid%tsk(I,J)
!!!MARS
  ENDDO
  ENDDO
    IF (planet.eq."prescribed") Then
      call read_hr(profdustq,profdustn,nl_in)
      open(unit=17,file="prescribed_sw.txt",action="write")
      open(unit=18,file="prescribed_lw.txt",action="write")
      DO k=1,kte!-1
        p_level = grid%em_znu(k)*(pd_surf - grid%p_top) + grid%p_top
        prescribed_sw(k) = interp_0( profdustq, pd_in, p_level, nl_in )
        prescribed_lw(k) = interp_0( profdustn, pd_in, p_level, nl_in )
        write (17,*) prescribed_sw(k)
        write (18,*) prescribed_lw(k)
      ENDDO
      close(unit=17)
      close(unit=18)
    ENDIF

    if (      ( config_flags%mars == 1  ) &
         .OR. ( config_flags%mars == 11 ) &
         .OR. ( config_flags%mars == 12 ) ) then
      print *, '**** INTERPOLATE HV **** RANK 2 in SCALAR'
      DO k=1,kte-1
         p_level = grid%em_znu(k)*(pd_surf - grid%p_top) + grid%p_top
         scalar(its:ite,k,jts:jte,2) = interp_0( qv, pd_in, p_level, nl_in )
         scalar(its:ite,k,jts:jte,3) = 0.
           !! water ice is set to 0 (was put into water vapor when building prof
           !from MCD)
      ENDDO
      print *, "WATER VAPOR",scalar(its,:,jts,2)
    endif

 END SUBROUTINE init_domain_rk

   SUBROUTINE init_module_initialize
   END SUBROUTINE init_module_initialize

!---------------------------------------------------------------------

!  test driver for get_sounding
!
!      implicit none
!      integer n
!      parameter(n = 1000)
!      real zk(n),p(n),theta(n),rho(n),u(n),v(n),qv(n),pd(n)
!      logical dry
!      integer nl,k
!
!      dry = .false.
!      dry = .true.
!      call get_sounding( zk, p, pd, theta, rho, u, v, qv, dry, n, nl )
!      write(6,*) ' input levels ',nl
!      write(6,*) ' sounding '
!      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
!      do k=1,nl
!        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), pd(k), theta(k), rho(k), u(k), v(k), qv(k)
!      enddo
!      end
!
!---------------------------------------------------------------------------

      subroutine get_sounding( zk, p, p_dry, theta, tk, rho, &
                               u, v, qv, dry, nl_max, nl_in, &
                               mulu, mulv, addu, addv )
      implicit none

      integer nl_max, nl_in
      real zk(nl_max), p(nl_max), theta(nl_max), tk(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
      logical dry

      integer n
      parameter(n=1000)
      logical debug

!      parameter( debug = .false.)
!****Mars
      parameter( debug = .true.)
      real mulu, mulv, addu, addv


! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

      !! special MARS
      real r_input(n)
      real cp_input(n)
      real cv_input(n)
      real cvpm_input(n)
      real pfile_input(n)
      real t_input(n)
      real rhofile_input(n)
      !! special MARS

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding

! local data

      !real p1000mb,cv,cp,r,cvpm,g
!****Mars
!      parameter (p1000mb = 1.e+05, r = 287, cp = 1003., cv = cp-r, cvpm = -cv/cp, g=9.81 )
!      parameter (p1000mb = 610., r = 192., cp = 844.6, cv = cp-r, cvpm = -cv/cp, g=3.72)
!      parameter (p1000mb = 610., r = 191., cp = 744.5, cv = cp-r, cvpm = -cv/cp, g=3.72)
!****Mars
      integer k, it, nl
      real qvf, qvf1, dz

!  first, read the sounding

      call read_sounding( p_surf, th_surf, qv_surf, &
                          h_input, th_input, qv_input, &
                          u_input, v_input, r_input, cp_input, &
                          pfile_input, t_input, rhofile_input, n, nl, debug )


      !! special MARS
      do k=1,nl
        cv_input(k) = cp_input(k) - r_input(k)
        cvpm_input(k) = - cv_input(k) / cp_input(k)
      enddo
      !! special MARS

      if(dry) then
       do k=1,nl
         qv_input(k) = 0.
       enddo
      endif

      if(debug) write(6,*) ' number of input levels = ',nl

        nl_in = nl
        if(nl_in .gt. nl_max ) then
          write(6,*) ' too many levels for input arrays ',nl_in,nl_max
          call wrf_error_fatal ( ' too many levels for input arrays ' )
        end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      do k=1,nl
        qv_input(k) = 0.001*qv_input(k)
      enddo

      p_surf = 100.*p_surf  ! convert to pascals
      qvf = 1. + rvovrd*qv_input(1) 
      rho_surf = 1./((r_d/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(rcp)
          !!!!!! rcp variable
          !rho_surf =  1./((r_input(1)/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm_input(1)))
          !pi_surf = (p_surf/p1000mb)**(r_input(1)/cp_input(1))


      if(debug) then
        write(6,*) ' surface density is ',rho_surf
        write(6,*) ' surface pi is      ',pi_surf
      end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

          qvf = 1. + rvovrd*qv_input(1) 
          qvf1 = 1. + qv_input(1)
          rho_input(1) = rho_surf
          dz = h_input(1)
          do it=1,10
!!MARS MARS
            pm_input(1) = p_surf !&
!                    - dz*(0.25*rho_surf+0.75*rho_input(1))*g*qvf1  !!! BEURK
!                    - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1  !! parce que couche 1 tres proche
            rho_input(1) = 1./((r_d/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
               !!!!!!! rcp variable
               !rho_input(1) = 1./((r_input(1)/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm_input(1)))
          enddo

! integrate up the column

          do k=2,nl
            rho_input(k) = rho_input(k-1)
            dz = h_input(k)-h_input(k-1)
               !!!!!!! rcp variable
               !dz = r_input(k) * t_input(k) * (- p_input(k) + p_input(k-1)) / p_input(k) / g
               !!dz = - cp_input(k) * (- t_input(k) + t_input(k-1)) / g
            qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
            qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here

print *, 'input', pfile_input(k), rhofile_input(k)  

            do it=1,10 !!ou moins??? non. !! trop de rho(k-1) donne une pression trop faible puis crash
                                          !! mais coeff ci-dessous vont varier la pression calcul√e
              pm_input(k) = pm_input(k-1) &
                      - dz*(0.75*rho_input(k)+0.25*rho_input(k-1))*g*qvf1
                      !- 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
              rho_input(k) = 1./((r_d/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
!!
!! marche pas 
!!
!pm_input(k) = pm_input(k-1) &
!        - 0.5*dz*(1./rho_input(k)+1./rho_input(k-1))*g*qvf1
!rho_input(k) = (r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm)
              !!!!!!! rcp variable
              !rho_input(k) = 1./((r_input(k)/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm_input(k)))
              !print *, p_input(k), pm_input(k),((r_input(k)/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm_input(k))),k
print *, it, pm_input(k), rho_input(k), dz
            enddo
          enddo


!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          tk(k) = t_input(k)
          rho(k) = rho_input(k)
          u(k) = mulu*u_input(k) + addu
          v(k) = mulv*v_input(k) + addv
          qv(k) = qv_input(k)

         !!!! direct input from file
         write(6,*) '*** DIRECT INPUT FROM FILE ***' 
         p(k) = pfile_input(k)
         p_dry(k) = pfile_input(k)
         rho(k) = rhofile_input(k)

        enddo

     if(debug) then
      write(6,*) ' sounding '
      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

     end if

      end subroutine get_sounding

!-------------------------------------------------------

      subroutine read_sounding( ps,ts,qvs,h,th,qv,u,v,r,cp,p,t,rho,n,nl,debug )
      implicit none
      integer n,nl
      real ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n),r(n),cp(n),p(n),t(n),rho(n) 
      logical end_of_file
      logical debug

      integer k

      open(unit=10,file='input_sounding',form='formatted',status='old')
      rewind(10)
      read(10,*) ps, ts, qvs
      if(debug) then
        write(6,*) ' input sounding surface parameters '
        write(6,*) ' surface pressure (mb) ',ps
        write(6,*) ' surface pot. temp (K) ',ts
        write(6,*) ' surface mixing ratio (g/kg) ',qvs
      end if

      end_of_file = .false.
      k = 0

      do while (.not. end_of_file)

        read(10,*,end=100) h(k+1), th(k+1), qv(k+1), u(k+1), v(k+1)
        k = k+1
        if(debug) write(6,'(1x,i3,5(1x,e10.3))') k, h(k), th(k), qv(k), u(k), v(k)
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo


      !!! special MARS
      open(unit=11,file='input_therm',form='formatted',status='old')
      rewind(11)
      end_of_file = .false.
      k = 0
      do while (.not. end_of_file)

        read(11,*,end=101) r(k+1), cp(k+1), p(k+1), rho(k+1), t(k+1) 
        write(*,*) k, r(k+1), cp(k+1), p(k+1), rho(k+1), t(k+1)
        k = k+1
        go to 112
 101    end_of_file = .true.
 112    continue
      enddo
      !!! special MARS



      nl = k

      close(unit=10,status = 'keep')

      end subroutine read_sounding

      subroutine read_hr(hr_sw,hr_lw,n)
      implicit none
      integer n
      real hr_sw(n),hr_lw(n)
      logical end_of_file

      integer k

! first element is the surface

      open(unit=11,file='input_hr',form='formatted',status='old')
      rewind(11)
      end_of_file = .false.
      k = 0
      do while (.not. end_of_file)

        read(11,*,end=102) hr_sw(k+1),hr_lw(k+1)
        write(*,*) k,hr_sw(k+1),hr_lw(k+1)
        k = k+1
        go to 113
 102    end_of_file = .true.
 113    continue
      enddo

      close(unit=11,status = 'keep')

      end subroutine read_hr

END MODULE module_initialize
