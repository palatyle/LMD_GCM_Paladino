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

MODULE module_initialize_ideal

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
   INTEGER :: idum1, idum2

   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )

     CALL init_domain_rk( grid &
!
#include <actual_new_args.inc>
!
                        )

   END SUBROUTINE init_domain

!-------------------------------------------------------------------

   SUBROUTINE init_domain_rk ( grid &
!
# include <dummy_new_args.inc>
!
)
   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid

# include <dummy_new_decl.inc>

   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k

   ! Local data

   INTEGER, PARAMETER :: nl_max = 1000
   REAL, DIMENSION(nl_max) :: zk, p_in, theta, rho, u, v, qv, pd_in
   REAL*8, DIMENSION(nl_max) :: pd_in8
   INTEGER :: nl_in

   INTEGER :: icm,jcm, ii, im1, jj, jm1, loop, error, fid, nxc, nyc
   REAL    :: u_mean,v_mean, f0, p_surf, p_level, qvf, z_at_v, z_at_u
   REAL*8 :: p_level8
   REAL    :: xrad, yrad, zrad, rad, delt, cof1, cof2
!   REAL, EXTERNAL :: interp_0
   REAL    :: hm
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
!  For LES, add randx
   real :: randx

!!MARS
 REAL :: lon_input, lat_input, alt_input, tsurf_input
 ! for mode 3
 REAL, DIMENSION(nl_max) :: profdustq,profdustn
 REAL, DIMENSION(nl_max) :: prescribed_sw,prescribed_lw,prescribed_dyn
 REAL, DIMENSION(nl_max) :: hrsw,hrlw,hrdyn
 REAL, DIMENSION(nl_max) :: lsf_dt,lsf_dq,lsfdt,lsfdq
 REAL, DIMENSION(nl_max) :: venus_hrdyn
 REAL, DIMENSION(nl_max) :: altitude
 REAL*8, DIMENSION(nl_max) :: trac
!!MARS

      REAL :: pfu, pfd, phm
      INTEGER :: hypsometric_opt = 1 ! classic
      !INTEGER :: hypsometric_opt = 2 ! Wee et al. 2012 correction

      LOGICAL :: logp = .true. ! use logp to interpolate (and not p)

#ifdef DM_PARALLEL
#    include <data_calls.inc>
#endif


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

   IF (planet == "mars") THEN
     stretch_grid = .false.
     !! FOR LES, set stretch to false
   ELSE
     stretch_grid = .true. !! VENUS
   ENDIF
   delt = 3.
!   z_scale = .50
!   z_scale = .10
!   z_scale = .25
!   z_scale = .15
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

! here we check to see if the boundary conditions are set properly

   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

   moisture_init = .true.

    grid%itimestep=0

#ifdef DM_PARALLEL
   CALL wrf_dm_bcast_bytes( icm , IWORDSIZE )
   CALL wrf_dm_bcast_bytes( jcm , IWORDSIZE )
#endif

    CALL nl_set_mminlu(1, '    ')
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
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
!  for LES, include Coriolis force
         grid%f(i,j)        = 0.  !!MARS MARS 1.e-4 
!!      grid%f(i,j)     = 2*EOMEG*SIN(grid%xlat(i,j)*degrad)
      END DO
   END DO

    DO j = jts, jte
    DO k = kts, kte
      DO i = its, ite
         grid%ww(i,k,j)     = 0.
      END DO
   END DO
   END DO

   grid%step_number = 0

! set up the grid

   IF (stretch_grid) THEN ! exponential stretch for eta (nearly constant dz)
     DO k=1, kde
      grid%znw(k) = (exp(-(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
                                (1.-exp(-1./z_scale))
     ENDDO
   ELSE
     !!MARS
     !!MARS
     IF (planet .ne. "mars") THEN
       open(unit=12,file='levels',form='formatted',status='old')
       rewind(12)
       DO k=1, kde
        read(12,*) grid%znw(k)
        write(6,*) 'read level ', k,grid%znw(k)
       ENDDO
       close(12)
     ENDIF
     !!MARS
     !!MARS
     !  !DO k=1, kde
     !  ! grid%znw(k) = 1. - float(k-1)/float(kde-1)
     !  !ENDDO
   ENDIF

   !! SPECIFIC FOR LES PBL MARS
   IF (planet == "mars") THEN
     !!!MARS
     grid%znw(1)=1.000
     grid%znw(2)=0.9995 !5m
     grid%znw(3)=0.9980 !20m
     grid%znw(4)=0.9950 !55m
     DO k=5, kde
       grid%znw(k) = grid%znw(4) * ( 1. - float(k-4)/float(kde-4) )
     ENDDO
   ENDIF


   DO k=1, kde-1
    grid%dnw(k) = grid%znw(k+1) - grid%znw(k)
    grid%rdnw(k) = 1./grid%dnw(k)
    grid%znu(k) = 0.5*(grid%znw(k+1)+grid%znw(k))
   ENDDO
   DO k=2, kde-1
    grid%dn(k) = 0.5*(grid%dnw(k)+grid%dnw(k-1))
    grid%rdn(k) = 1./grid%dn(k)
    grid%fnp(k) = .5* grid%dnw(k  )/grid%dn(k)
    grid%fnm(k) = .5* grid%dnw(k-1)/grid%dn(k)
   ENDDO

   cof1 = (2.*grid%dn(2)+grid%dn(3))/(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(2) 
   cof2 =     grid%dn(2)        /(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(3) 
   grid%cf1  = grid%fnp(2) + cof1
   grid%cf2  = grid%fnm(2) - cof1 - cof2
   grid%cf3  = cof2       

   grid%cfn  = (.5*grid%dnw(kde-1)+grid%dn(kde-1))/grid%dn(kde-1)
   grid%cfn1 = -.5*grid%dnw(kde-1)/grid%dn(kde-1)
   grid%rdx = 1./config_flags%dx
   grid%rdy = 1./config_flags%dy

!  get the sounding from the ascii sounding file, first get dry sounding and 
!  calculate base state

  dry_sounding = .true.
  IF ( wrf_dm_on_monitor() ) THEN
  write(6,*) ' getting dry sounding for base state '

  CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, nl_max, nl_in )
  ENDIF
  CALL wrf_dm_bcast_real( zk , nl_max )
  CALL wrf_dm_bcast_real( p_in , nl_max )
  CALL wrf_dm_bcast_real( pd_in , nl_max )
  CALL wrf_dm_bcast_real( theta , nl_max )
  CALL wrf_dm_bcast_real( rho , nl_max )
  CALL wrf_dm_bcast_real( u , nl_max )
  CALL wrf_dm_bcast_real( v , nl_max )
  CALL wrf_dm_bcast_real( qv , nl_max )
  CALL wrf_dm_bcast_integer ( nl_in , 1 ) 

  write(6,*) ' returned from reading sounding, nl_in is ',nl_in

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

!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure

  write(6,*) ' ztop above ground is ',config_flags%ztop
  write(6,*) ' real ztop is ',config_flags%ztop + alt_input
  grid%p_top = interp_0( p_in, zk, config_flags%ztop + alt_input, nl_in )

  DO j=jts,jte
  DO i=its,ite
!!MARS
    grid%ht(i,j) = alt_input
    grid%m_tsurf(i,j) = tsurf_input
!!MARS
    grid%xlat(i,j) = lat_input !+ float(j)*config_flags%dy/59000.
    grid%xlong(i,j) = lon_input !+ float(i)*config_flags%dx/59000.
    grid%m_emiss(i,j)=0.95
    grid%m_co2ice(i,j)=0.
    grid%m_h2oice(i,j)=0.
!! >> Used for restarts only:
    grid%m_q2(i,:,j)=0.
    grid%m_fluxrad(i,j)=0.
    grid%m_wstar(i,j)=0.
!! <<
    grid%slpx(i,j) = 0.
    grid%slpy(i,j) = 0.
   DO k=1,config_flags%num_soil_layers
    grid%m_tsoil(i,k,j) = 0.
   ENDDO
    grid%m_gw(i,1,j) = 0.
    grid%m_gw(i,2,j) = 0.
    grid%m_gw(i,3,j) = 0.
    grid%m_gw(i,4,j) = 0.
    grid%m_gw(i,5,j) = 0.
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
     grid%ht(i,j) = alt_input + mtn_ht * 0.25 * &
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
    grid%phb(i,1,j) = g * grid%ht(i,j)
    grid%ph0(i,1,j) = g * grid%ht(i,j)
  ENDDO
  ENDDO

  IF (.not.logp) THEN
    write(6,*) 'interpolate in p'
  ELSE
    write(6,*) 'interpolate in logp'
  ENDIF

  DO J = jts, jte
  DO I = its, ite

    p_surf = interp_0( p_in, zk, grid%phb(i,1,j)/g, nl_in )
    grid%mub(i,j) = p_surf-grid%p_top

!  this is dry hydrostatic sounding (base state), so given grid%p (coordinate),
!  interp theta (from interp) and compute 1/rho from eqn. of state

    DO K = 1, kte-1
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
      grid%pb(i,k,j) = p_level
     IF (.not.logp) THEN
      grid%t_init(i,k,j) = interp_0( theta, p_in, p_level, nl_in ) - t0
     ELSE
      grid%t_init(i,k,j) = interp_0_log( theta, p_in, p_level, nl_in ) - t0
     ENDIF
      grid%alb(i,k,j) = (r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm
    ENDDO

!  calc hydrostatic balance (alternatively we could interp the geopotential from the
!  sounding, but this assures that the base state is in exact hydrostatic balance with
!  respect to the model eqns.

   IF (hypsometric_opt == 1) THEN
    DO k  = 2,kte
      grid%phb(i,k,j) = grid%phb(i,k-1,j) - grid%dnw(k-1)*grid%mub(i,j)*grid%alb(i,k-1,j)
    ENDDO
   ELSE IF (hypsometric_opt == 2) THEN
    DO k = 2,kte
      pfu = grid%mub(i,j)*grid%znw(k)   + grid%p_top
      pfd = grid%mub(i,j)*grid%znw(k-1)   + grid%p_top
      phm = grid%mub(i,j)*grid%znu(k-1)   + grid%p_top
      grid%phb(i,k,j) = grid%phb(i,k-1,j) + grid%alb(i,k-1,j)*phm*LOG(pfd/pfu)
    END DO
   END IF


  ENDDO
  ENDDO
  IF ( wrf_dm_on_monitor() ) THEN
    write(6,*) ' ptop is ',grid%p_top
    write(6,*) ' base state grid%mub(1,1), p_surf is ',grid%mub(1,1),grid%mub(1,1)+grid%p_top
  ENDIF

!  calculate full state for each column - this includes moisture.

  write(6,*) ' getting moist sounding for full state '
  dry_sounding = .false.
  CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, nl_max, nl_in )

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)

!  At this point grid%p_top is already set. find the DRY mass in the column 
!  by interpolating the DRY pressure.  

   pd_surf = interp_0( pd_in, zk, grid%phb(i,1,j)/g, nl_in )

!  compute the perturbation mass and the full mass

    grid%mu_1(i,j) = pd_surf-grid%p_top - grid%mub(i,j)
    grid%mu_2(i,j) = grid%mu_1(i,j)
    grid%mu0(i,j) = grid%mu_1(i,j) + grid%mub(i,j)

! given the dry pressure and coordinate system, interp the potential
! temperature and qv

    do k=1,kde-1

      p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
     IF (.not.logp) THEN
      moist(i,k,j,P_QV) = interp_0( qv, pd_in, p_level, nl_in )
      grid%t_1(i,k,j)          = interp_0( theta, pd_in, p_level, nl_in ) - t0
     ELSE
      moist(i,k,j,P_QV) = interp_0_log( qv, pd_in, p_level, nl_in )
      grid%t_1(i,k,j)          = interp_0_log( theta, pd_in, p_level, nl_in ) - t0
     ENDIF
      grid%t_2(i,k,j)          = grid%t_1(i,k,j)
      

    enddo

!  integrate the hydrostatic equation (from the RHS of the bigstep
!  vertical momentum equation) down from the top to get grid%p.
!  first from the top of the model to the top pressure

    k = kte-1  ! top level

    qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k,j,P_QV))
    qvf2 = 1./(1.+qvf1)
    qvf1 = qvf1*qvf2

!    grid%p(i,k,j) = - 0.5*grid%mu_1(i,j)/grid%rdnw(k)
    grid%p(i,k,j) = - 0.5*(grid%mu_1(i,j)+qvf1*grid%mub(i,j))/grid%rdnw(k)/qvf2
    qvf = 1. + rvovrd*moist(i,k,j,P_QV)
    grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
    grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)

!  down the column

    do k=kte-2,1,-1
      qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k+1,j,P_QV))
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2
      grid%p(i,k,j) = grid%p(i,k+1,j) - (grid%mu_1(i,j) + qvf1*grid%mub(i,j))/qvf2/grid%rdn(k+1)
      qvf = 1. + rvovrd*moist(i,k,j,P_QV)
      grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                  (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
      grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
    enddo

!  this is the hydrostatic equation used in the model after the
!  small timesteps.  In the model, grid%al (inverse density)
!  is computed from the geopotential.


    grid%ph_1(i,1,j) = 0.
   IF (hypsometric_opt == 1) THEN
    DO k  = 2,kte
      grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (1./grid%rdnw(k-1))*(       &
                   (grid%mub(i,j)+grid%mu_1(i,j))*grid%al(i,k-1,j)+ &
                    grid%mu_1(i,j)*grid%alb(i,k-1,j)  )
                                                   
      grid%ph_2(i,k,j) = grid%ph_1(i,k,j) 
      grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
    ENDDO
   ELSE IF (hypsometric_opt == 2) THEN

             ! Alternative hydrostatic eq.: dZ = -al*p*dLOG(p), where p is dry pressure.
             ! Note that al*p approximates Rd*T and dLOG(p) does z.
             ! Here T varies mostly linear with z, the first-order integration produces better result.

               grid%ph_2(i,1,j) = grid%phb(i,1,j)
               DO k = 2,kte
                  pfu = grid%mu0(i,j)*grid%znw(k)   + grid%p_top
                  pfd = grid%mu0(i,j)*grid%znw(k-1) + grid%p_top
                  phm = grid%mu0(i,j)*grid%znu(k-1) + grid%p_top
                  grid%ph_2(i,k,j) = grid%ph_2(i,k-1,j) + grid%alt(i,k-1,j)*phm*LOG(pfd/pfu)
               END DO

               DO k = 1,kte
                  grid%ph_2(i,k,j) = grid%ph_2(i,k,j) - grid%phb(i,k,j)
                  grid%ph_1(i,k,j) = grid%ph_2(i,k,j)
               END DO

   END IF



    IF ( wrf_dm_on_monitor() ) THEN
    if((i==2) .and. (j==2)) then
     write(6,*) ' grid%ph_1 calc ',grid%ph_1(2,1,2),grid%ph_1(2,2,2),&
                              grid%mu_1(2,2)+grid%mub(2,2),grid%mu_1(2,2), &
                              grid%alb(2,1,2),grid%al(1,2,1),grid%rdnw(1)
    endif
    ENDIF

  ENDDO
  ENDDO

!#if 0

!  thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

! For LES, change the initial random perturbations
! For 2D test, call randx outside I-loop
! For 3D runs, call randx inside both I-J loops

  DO J = jts, min(jde-1,jte)
!   yrad = config_flags%dy*float(j-nyc)/10000.
    yrad = 0.
    DO I = its, min(ide-1,ite)
!     xrad = config_flags%dx*float(i-nxc)/10000.
      xrad = 0.
      call random_number (randx)
      randx = randx - 0.5
!     DO K = 1, kte-1
      DO K = 1, 4 

!  No bubbles for LES!
!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

!       zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
!                  +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
!       zrad = (zrad-1500.)/1500.
        zrad = 0.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
!          grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+ 0.1 *randx
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

   IF (hypsometric_opt == 1) THEN

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (1./grid%rdnw(k-1))*(       &
                     (grid%mub(i,j)+grid%mu_1(i,j))*grid%al(i,k-1,j)+ &
                      grid%mu_1(i,j)*grid%alb(i,k-1,j)  )
                                                   
        grid%ph_2(i,k,j) = grid%ph_1(i,k,j) 
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

   ELSE IF (hypsometric_opt == 2) THEN

             ! Alternative hydrostatic eq.: dZ = -al*p*dLOG(p), where p is dry pressure.
             ! Note that al*p approximates Rd*T and dLOG(p) does z.
             ! Here T varies mostly linear with z, the first-order integration produces better result.

               grid%ph_2(i,1,j) = grid%phb(i,1,j)
               DO k = 2,kte
                  pfu = grid%mu0(i,j)*grid%znw(k)   + grid%p_top
                  pfd = grid%mu0(i,j)*grid%znw(k-1) + grid%p_top
                  phm = grid%mu0(i,j)*grid%znu(k-1) + grid%p_top
                  grid%ph_2(i,k,j) = grid%ph_2(i,k-1,j) + grid%alt(i,k-1,j)*phm*LOG(pfd/pfu)
               END DO

               DO k = 1,kte
                  grid%ph_2(i,k,j) = grid%ph_2(i,k,j) - grid%phb(i,k,j)
                  grid%ph_1(i,k,j) = grid%ph_2(i,k,j)
               END DO

   END IF


    ENDDO
  ENDDO

!#endif

   IF ( wrf_dm_on_monitor() ) THEN
   write(6,*) ' grid%mu_1 from comp ', grid%mu_1(1,1)
   write(6,*) ' full state sounding from comp, ph/g, grid%p, grid%al, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, (grid%ph_1(1,k,1)+grid%phb(1,k,1))/g, &
                                      grid%p(1,k,1)+grid%pb(1,k,1), grid%alt(1,k,1), &
                                      grid%t_1(1,k,1)+t0, moist(1,k,1,P_QV)
   enddo

   write(6,*) ' pert state sounding from comp, grid%ph_1, pp, alp, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1), &
                                      grid%p(1,k,1), grid%al(1,k,1), &
                                      grid%t_1(1,k,1), moist(1,k,1,P_QV)
   enddo
   ENDIF

! interp v

  DO J = jts, jte
  DO I = its, min(ide-1,ite)

    IF (j == jds) THEN
      z_at_v = grid%phb(i,1,j)/g
    ELSE IF (j == jde) THEN
      z_at_v = grid%phb(i,1,j-1)/g
    ELSE
      z_at_v = 0.5*(grid%phb(i,1,j)+grid%phb(i,1,j-1))/g
    END IF

    p_surf = interp_0( p_in, zk, z_at_v, nl_in )

    DO K = 1, kte-1
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
     IF (.not.logp) THEN
      grid%v_1(i,k,j) = interp_0( v, p_in, p_level, nl_in )
     ELSE
      grid%v_1(i,k,j) = interp_0_log( v, p_in, p_level, nl_in )
     ENDIF
      grid%v_2(i,k,j) = grid%v_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

! interp u

  DO J = jts, min(jde-1,jte)
  DO I = its, ite

    IF (i == ids) THEN
      z_at_u = grid%phb(i,1,j)/g
    ELSE IF (i == ide) THEN
      z_at_u = grid%phb(i-1,1,j)/g
    ELSE
      z_at_u = 0.5*(grid%phb(i,1,j)+grid%phb(i-1,1,j))/g
    END IF

    p_surf = interp_0( p_in, zk, z_at_u, nl_in )

    DO K = 1, kte-1
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
     IF (.not.logp) THEN
      grid%u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
     ELSE
      grid%u_1(i,k,j) = interp_0_log( u, p_in, p_level, nl_in )
     ENDIF
      grid%u_2(i,k,j) = grid%u_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

!  set w

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte
  DO I = its, min(ide-1,ite)
    grid%w_1(i,k,j) = 0.
    grid%w_2(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

!!!MARS MARS
IF (config_flags%init_MU .ne. 0.) THEN
  grid%u_1 = grid%u_1*config_flags%init_MU
  grid%u_2 = grid%u_2*config_flags%init_MU
  print *, 'multiply zonal wind ', config_flags%init_MU
ENDIF
IF (config_flags%init_MV .ne. 0.) THEN
  grid%v_1 = grid%v_1*config_flags%init_MV
  grid%v_2 = grid%v_2*config_flags%init_MV
  print *, 'multiply meridional wind ', config_flags%init_MV
ENDIF
IF (config_flags%init_U .ne. 0.) THEN
  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%u_1(i,k,j) = config_flags%init_U
    grid%u_2(i,k,j) = config_flags%init_U
  ENDDO
  ENDDO
  ENDDO
  print *, 'constant zonal wind ', config_flags%init_U
  !!! ****** ou autre possibilité
  !!! >   grid%u_1 = grid%u_1*0. + config_flags%init_U
  !!! >   grid%u_2 = grid%u_2*0. + config_flags%init_U
ENDIF
IF (config_flags%init_V .ne. 0.) THEN
  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%v_1(i,k,j) = config_flags%init_V
    grid%v_2(i,k,j) = config_flags%init_V
  ENDDO
  ENDDO
  ENDDO
  print *, 'constant meridional wind ', config_flags%init_V
ENDIF

!!!MARS MARS


!  set a few more things

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%h_diabatic(i,k,j) = 0.
      !!!!! MARS NO WIND CASE
      !grid%u_1(i,k,j) = 0.
      !grid%u_2(i,k,j) = 0.
      !grid%v_1(i,k,j) = 0.
      !grid%v_2(i,k,j) = 0.
      !!!!! MARS NO WIND CASE
  ENDDO
  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
  DO k=1,kte-1
    grid%t_base(k) = grid%t_1(1,k,1)
    grid%qv_base(k) = moist(1,k,1,P_QV)
    grid%u_base(k) = grid%u_1(1,k,1)
    grid%v_base(k) = grid%v_1(1,k,1)
    grid%z_base(k) = 0.5*(grid%phb(1,k,1)+grid%phb(1,k+1,1)+grid%ph_1(1,k,1)+grid%ph_1(1,k+1,1))/g

!!!!! MARS SIMPLE LES (PURE BUOYANCY)
!!      grid%t_base(k)  = grid%t_init(its,k,jts)
!      grid%t_base(k) = 0.
!      grid%qv_base(k) = 0.
!      grid%u_base(k)  = 0.
!      grid%v_base(k)  = 0.
!      grid%z_base(k) = 0.
!!!!! MARS SIMPLE LES

  ENDDO
  ENDIF
  CALL wrf_dm_bcast_real( grid%t_base , kte )
  CALL wrf_dm_bcast_real( grid%qv_base , kte )
  CALL wrf_dm_bcast_real( grid%u_base , kte )
  CALL wrf_dm_bcast_real( grid%v_base , kte )
  CALL wrf_dm_bcast_real( grid%z_base , kte )

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
     thtmp   = grid%t_2(i,1,j)+t0
     ptmp    = grid%p(i,1,j)+grid%pb(i,1,j)
     temp(1) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,2,j)+t0
     ptmp    = grid%p(i,2,j)+grid%pb(i,2,j)
     temp(2) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,3,j)+t0
     ptmp    = grid%p(i,3,j)+grid%pb(i,3,j)
     temp(3) = thtmp * (ptmp/p1000mb)**rcp

!!    For LES-CBL, add 5 degrees to the surface temperature!
!!
!     grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)
!!     grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)+5.
     grid%tmn(I,J)=grid%tsk(I,J)-0.5

  ENDDO
  ENDDO

!!!!! MARS 

    ! interpolate water vapor
    if (      ( config_flags%mars == 1  ) &
         .OR. ( config_flags%mars == 11 ) &
         .OR. ( config_flags%mars == 12 ) ) then
      print *, '**** INTERPOLATE HV **** RANK 2 in SCALAR'
      DO k=1,kte-1
         p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
        IF (.not.logp) THEN
         scalar(its:ite,k,jts:jte,2) = interp_0( qv, pd_in, p_level, nl_in )
        ELSE
         scalar(its:ite,k,jts:jte,2) = interp_0_log( qv, pd_in, p_level, nl_in )
        ENDIF
         scalar(its:ite,k,jts:jte,3) = 0.
           !! water ice is set to 0 (was put into water vapor when building prof from MCD)
      ENDDO
      print *, "WATER VAPOR",scalar(its,:,jts,2)
    endif

    ! interpolate qdust
    if (      ( config_flags%mars == 11 ) &
         .OR. ( config_flags%mars == 12 ) ) then
      call read_dust(profdustq,profdustn,nl_in)
      print *, '**** INTERPOLATE DUSTQ **** RANK 4 in SCALAR'
      print *, '**** INTERPOLATE DUSTN **** RANK 5 in SCALAR'
      DO k=1,kte-1
         p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
        IF (.not.logp) THEN
         scalar(its:ite,k,jts:jte,4) = interp_0( profdustq, pd_in, p_level,nl_in )
         scalar(its:ite,k,jts:jte,5) = interp_0( profdustn, pd_in, p_level,nl_in )
        ELSE
         scalar(its:ite,k,jts:jte,4) = interp_0_log( profdustq, pd_in, p_level, nl_in )
         scalar(its:ite,k,jts:jte,5) = interp_0_log( profdustn, pd_in, p_level, nl_in )
        ENDIF
      ENDDO
      print *, "DUST Q", scalar(its,:,jts,4)
      print *, "DUST N", scalar(its,:,jts,5)
    endif

    if ( config_flags%mars == 12 ) then
      scalar(its:ite,1:kte-1,jts:jte,6) = 0.
      scalar(its:ite,1:kte-1,jts:jte,7) = 0.
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (planet .ne. "mars") THEN
    call read_dust(profdustq,profdustn,nl_in)
    DO k=1,kte!-1
       p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
       DO j = jts, jte
       DO i = its, ite
         !!! we use Q2 as a vehicle for heating rates! sick!
         grid%m_q2(i,k,j) = interp_0_log( profdustq, pd_in, p_level, nl_in ) &
                          + interp_0_log( profdustn, pd_in, p_level, nl_in )
       ENDDO
       ENDDO
       !print*,'grid%m_q2' 
       !print*,k,grid%m_q2(1,k,1)
    ENDDO
   ENDIF

    IF (planet.eq."prescribed") Then
      call read_hr(hrsw,hrlw,hrdyn,nl_in)
      open(unit=17,file="prescribed_sw.txt",action="write")
      open(unit=18,file="prescribed_lw.txt",action="write")
      open(unit=19,file="prescribed_dyn.txt",action="write")
      DO k=1,kte!-1
        p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
        prescribed_sw(k) = interp_0_log( hrsw, pd_in, p_level, nl_in )
        prescribed_lw(k) = interp_0_log( hrlw, pd_in, p_level, nl_in )
        prescribed_dyn(k) = interp_0_log( hrdyn, pd_in, p_level, nl_in )
        write (17,*) prescribed_sw(k)
        write (18,*) prescribed_lw(k)
        write (19,*) prescribed_dyn(k)
      ENDDO
      close(unit=19)
      close(unit=18)
      close(unit=17)
    ENDIF
    
    IF (planet.eq."venus") Then
      call read_hr(hrsw,hrlw,hrdyn,nl_in)
      open(unit=20,file="venus_hrdyn.txt",action="write")
      DO k=1,kte!-1
        p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
        venus_hrdyn(k) = interp_0_log( hrdyn, pd_in, p_level, nl_in )
        write (20,*) venus_hrdyn(k)
      ENDDO
      close(unit=20)
    ENDIF

    IF (planet.eq."generic") THEN
      call read_lsf(lsfdt,lsfdq,nl_in)
      open(unit=17,file="lsf.txt",action="write")
      DO k=1,kte!-1
        p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
        lsf_dt = interp_0_log( lsfdt, pd_in, p_level, nl_in )
        lsf_dq = interp_0_log( lsfdq, pd_in, p_level, nl_in )
        write (17,*) lsf_dt(k),lsf_dq(k)
      ENDDO
    ENDIF

    IF ((planet.eq."venus") .AND. ( config_flags%mars == 34 )) Then
       pd_in8(:)=pd_in(:)
       do i = 1,34
         call read_tracer(trac,num_scalar,i,nl_in)
         DO k=1,kte-1
           p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
           p_level8=p_level
           scalar(its:ite,k,jts:jte,i+1) =  interp_0_log2( trac, pd_in8, p_level8, nl_in )
         ENDDO
       ENDDO
       !close(unit=22)
    ENDIF

    open(unit=21,file="altitude.txt",action="write")
    DO k=1,kte!-1
      p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
      altitude(k) = interp_0_log( zk, pd_in, p_level, nl_in )
      write (21,*) altitude(k)
    ENDDO
    close(unit=21)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! MARS 


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

      subroutine get_sounding( zk, p, p_dry, theta, rho, &
                               u, v, qv, dry, nl_max, nl_in )
      implicit none

      integer nl_max, nl_in
      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
      logical dry

      integer n
      parameter(n=1000)
      logical debug
      parameter( debug = .true.)

! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! input therm data (element 0 is the ground so it's n+1 but n is 1000 anyway so...)

     real r_therm(n),cp_therm(n),p_therm(n),rho_therm(n),t_therm(n)

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding

! local data

!      real p1000mb,cv,cp,r,cvpm,g
!      parameter (p1000mb = 1.e+05, r = 287, cp = 1003., cv = cp-r, cvpm = -cv/cp, g=9.81 )
!      parameter (p1000mb = 610., r = 192., cp = 844.6, cv = cp-r, cvpm = -cv/cp, g=3.72)
!      parameter (p1000mb = 610., r = 191., cp = 744.5, cv = cp-r, cvpm = -cv/cp, g=3.72)
       !parameter (p1000mb = 92.e+05, r = 192., cp = 900.0, cv = cp-r, cvpm = -cv/cp, g=8.87)
      integer k, it, nl
      real qvf, qvf1, dz

      LOGICAL :: direct_from_file

      IF (planet == "mars") THEN
        direct_from_file = .false.
      ELSE
        direct_from_file = .true.
      ENDIF

!  first, read the sounding

      call read_sounding( p_surf, th_surf, qv_surf, &
                          h_input, th_input, qv_input, &
                          u_input, v_input,n, nl, debug )

! and the therm :

      call read_therm(r_therm,cp_therm,p_therm,rho_therm,t_therm,n)

      if(debug) write(6,*) ' number of input levels = ',nl
      nl_in = nl
      if(nl_in .gt. nl_max ) then
        write(6,*) ' too many levels for input arrays ',nl_in,nl_max
        call wrf_error_fatal ( ' too many levels for input arrays ' )
      end if

      IF (.NOT. direct_from_file) THEN

! To use r/cp as defined above, one has to recompute teta from T (default MCD computes
! teta for a variable r/cp)

      do k=1,nl
        th_input(k) = t_therm(k)*(p1000mb/p_therm(k))**(rcp)
      enddo
      th_surf = t_therm(1)*(p1000mb/p_therm(1))**(rcp)
! -----

      !if(dry) then
      ! do k=1,nl
      !   qv_input(k) = 0.
      ! enddo
      !endif

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      do k=1,nl
        !!!!!!!!!!!!!! MARS
        !! from mol/mol to kg/kg
        qv_input(k) = qv_input(k)*18./mwdry
        qv_input(k) = 0.001*qv_input(k)
      enddo

      p_surf = 100.*p_surf  ! convert to pascals
      qvf = 1. + rvovrd*qv_input(1) 
!!MARS
      qvf = 1.
      rho_surf = 1./((r_d/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(rcp)

      if(debug) then
        write(6,*) ' surface density is ',rho_surf
        write(6,*) ' surface pi is      ',pi_surf
      end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

          qvf = 1. + rvovrd*qv_input(1) 
          qvf1 = 1. + qv_input(1)
!!MARS
          qvf = 1.
          qvf1 = 1.
          rho_input(1) = rho_surf
          dz = h_input(1)
          do it=1,10
!            pm_input(1) = p_surf &
!                    - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
!!!MARS MARS MARS
            pm_input(1) = p_surf 
            rho_input(1) = 1./((r_d/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
          enddo

! integrate up the column

          do k=2,nl
            rho_input(k) = rho_input(k-1)
            dz = h_input(k)-h_input(k-1)
            qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
            qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here
!!MARS
          qvf = 1.
          qvf1 = 1.
!!MARS 

            do it=1,10
              pm_input(k) = pm_input(k-1) &
                      - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
              rho_input(k) = 1./((r_d/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
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

      ELSE !IF (.NOT. direct_from_file) THEN
        
        do k=1,nl
         !!!! direct input from file
         write(6,*) '*** DIRECT INPUT FROM FILE ***'
         pm_input(k) = p_therm(k)
         p_input(k) = p_therm(k)
         rho_input(k) = rho_therm(k)
        enddo

      ENDIF


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

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

      subroutine read_sounding( ps,ts,qvs,h,th,qv,u,v,n,nl,debug )
      implicit none
      integer n,nl
      real ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n)
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

      nl = k

      close(unit=10,status = 'keep')

      end subroutine read_sounding

      subroutine read_therm(r,cp,p,rho,t,n)
      implicit none
      integer n
      real r(n),cp(n),p(n),rho(n),t(n)
      logical end_of_file

      integer k

! first element is the surface

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

      close(unit=11,status = 'keep')

      end subroutine read_therm

!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_dust(pdustq,pdustn,n)
      implicit none
      integer n
      real pdustq(n+1),pdustn(n+1)
      logical end_of_file

      integer k

! first element is the surface

      open(unit=11,file='input_dust',form='formatted',status='old')
      rewind(11)
      end_of_file = .false.
      k = 0
      do while (.not. end_of_file)

        read(11,*,end=102) pdustq(k+1),pdustn(k+1)
        write(*,*) k,pdustq(k+1),pdustn(k+1)
        k = k+1
        go to 113
 102    end_of_file = .true.
 113    continue
      enddo

      close(unit=11,status = 'keep')

      end subroutine read_dust
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_hr(hr_sw,hr_lw,hr_dyn,n)
      implicit none
      integer n
      real hr_sw(n+1),hr_lw(n+1),hr_dyn(n+1)
      logical end_of_file

      integer k

! first element is the surface

      open(unit=11,file='input_hr',form='formatted',status='old')
      rewind(11)
      end_of_file = .false.
      k = 0
      do while (.not. end_of_file)

        read(11,*,end=103) hr_sw(k+1),hr_lw(k+1),hr_dyn(k+1)
        write(*,*) k,hr_sw(k+1),hr_lw(k+1),hr_dyn(k+1)
        k = k+1
        go to 114
 103    end_of_file = .true.
 114    continue
      enddo

      close(unit=11,status = 'keep')

      end subroutine read_hr

      subroutine read_lsf(dt,dq,n)
      implicit none
      integer n
      real dt(n+1),dq(n+1)
      logical end_of_file

      integer k

! first element is the surface

      open(unit=12,file='input_lsf',form='formatted',status='old')
      rewind(12)
      end_of_file = .false.
      k = 0
      do while (.not. end_of_file)

        read(12,*,end=103) dt(k+1),dq(k+1)
        write(*,*) k,dt(k+1),dq(k+1)
        k = k+1
        go to 114
 103    end_of_file = .true.
 114    continue
      enddo

      close(unit=12,status = 'keep')

      end subroutine read_lsf

      subroutine read_tracer(trace,nq,qn,n)
      implicit none
      integer n,qn,nq ! qn : number of the tracer
      real*8 tra(nq-1,n+1)
      real*8 trace(n+1) !output
      logical end_of_file

      integer k,j

! first element is the surface
      open(unit=14,file='input_tracer',form='formatted',status='old')
      rewind(14)
      end_of_file = .false.
        DO k=1,n
          read(14,*) tra(:,k)
          write(*,*) k,tra(qn,k)
      ENDDO

      close(14)
      trace(:)=tra(qn,:)
      end subroutine read_tracer


END MODULE module_initialize_ideal
