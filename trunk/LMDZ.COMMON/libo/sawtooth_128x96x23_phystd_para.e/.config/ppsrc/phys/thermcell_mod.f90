










MODULE thermcell_mod
      
      IMPLICIT NONE
      
      
! Flags for computations
                                                 !  default
      LOGICAL,SAVE :: dqimpl                     !  .true.   Flag for thermcell_dq version (True : implicit scheme || False : explicit scheme)
      LOGICAL,SAVE :: dvimpl                     !  .false.  Flag for specific u, v mixing (True : thermcell_dv2 || False : thermcell_dq)
      
! Physical parameters
      
      REAL,SAVE :: r_aspect_thermals             !  1.0      Aspect ratio of the thermals (width / height)
      REAL,SAVE :: tau_thermals                  !  0.       Relaxation time
      REAL,SAVE :: betalpha                      !  0.9      - between 0 (e=d) and 1 (rho*fraca=cst)
      REAL,SAVE :: afact                         !  2./3.    - buoyancy acceleration efficiency, between 0 and 1
      REAL,SAVE :: fact_epsilon                  !  2.e-3    - turbulence and friction deceleration
      REAL,SAVE :: nu                            !  0.       Geometrical contributions to entrainment and detrainment
      REAL,SAVE :: alpha_max                     !  0.7      Maximal permitted updraft fraction
      REAL,SAVE :: fomass_max                    !  0.5      Maximal permitted outgoing layer mass fraction
      REAL,SAVE :: pres_limit                    !  1.e5     -
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB : linf is used to set the lowest possible first level because we allow it
!      to begin higher than the surface. It is set to 2 in order to remove the
!      first layer for gas giant.
!      If there is a surface, it has to be set to 1.
!      If someone want to call more than once the thermal plume model in some
!      grid points, this variable may become a saved array of INTEGER with size
!      ngrid. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER,SAVE :: linf                       !     1
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB : d_temp is an artificial virtual potential temperature offset added in
!      layer linf which can be used to force convection to begin in it.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      REAL,SAVE :: d_temp                        !     0.
      
! Physical constants
      
      REAL,SAVE :: RTT
      REAL,SAVE :: RG
      REAL,SAVE :: RKAPPA
      REAL,SAVE :: RPI
      REAL,SAVE :: RD
      
!$OMP THREADPRIVATE(RTT, RG, RKAPPA, RPI, RD)
      
      
      CONTAINS
      
      SUBROUTINE init_thermcell_mod(g, rcp, r, pi, T_h2o_ice_liq, RV)
         
         USE ioipsl_getin_p_mod, ONLY: getin_p
         
         IMPLICIT NONE
         
         REAL g
         REAL rcp
         REAL r
         REAL pi
         REAL T_h2o_ice_liq
         REAL RV
         
         RTT = T_h2o_ice_liq
         RG = g
         RKAPPA = rcp
         RPI = pi
         RD = r
         
         print *, 'Implicit mixing scheme ?'
         dqimpl = .true.
         call getin_p('dqimpl', dqimpl)
         print *, 'dqimpl = ', dqimpl
         
         print *, 'Use specific horizontal momentum scheme ?'
         dvimpl = .false.
         call getin_p('dvimpl', dvimpl)
         print *, 'dvimpl = ', dvimpl
         
         print *, 'Plume aspect ratio ?'
         r_aspect_thermals = 1.
         call getin_p('r_aspect_thermals', r_aspect_thermals)
         print *, 'r_aspect_thermals = ', r_aspect_thermals
         
         print *, 'Relaxation time ?'
         tau_thermals = 0.
         call getin_p('tau_thermals', tau_thermals)
         print *, 'tau_thermals = ', tau_thermals
         
         print *, 'betalpha ?'
         betalpha = 0.9
         call getin_p('betalpha', betalpha)
         print *, 'betalpha = ', betalpha
         
         print *, 'Buoyancy acceleration efficiency ?'
         afact = 2./3.
         call getin_p('afact', afact)
         print *, 'afact = ', afact
         
         print *, 'Turbulence and friction factor ?'
         fact_epsilon = 2.e-3
         call getin_p('fact_epsilon', fact_epsilon)
         print *, 'fact_epsilon = ', fact_epsilon
         
         print *, 'Constant entrainment and detrainment term ?'
         nu = 0.
         call getin_p('nu', nu)
         print *, 'nu = ', nu
         
         print *, 'Maximal authorized updraft fraction ?'
         alpha_max = 0.7
         call getin_p('alpha_max', alpha_max)
         print *, 'alpha_max = ', alpha_max
         
         print *, 'Maximal authorized entrained mass flux ?'
         fomass_max = 0.5
         call getin_p('fomass_max', fomass_max)
         print *, 'fomass_max = ', fomass_max
         
         print *, 'Minimal pressure below which plume can no longer trigger ?'
         pres_limit = 1.e5
         call getin_p('pres_limit', pres_limit)
         print *, 'pres_limit = ', pres_limit
         
         print *, 'Deepest layer from which a plume can rise ?'
         linf = 1
         call getin_p('linf', linf)
         print *, 'linf = ', linf
         
         print *, 'd_temp ?'
         d_temp = 0.
         call getin_p('d_temp', d_temp)
         print *, 'd_temp = ', d_temp
         
         RETURN
      END SUBROUTINE init_thermcell_mod
      
END MODULE thermcell_mod
