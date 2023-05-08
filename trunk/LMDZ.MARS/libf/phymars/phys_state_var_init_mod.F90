MODULE phys_state_var_init_mod

CONTAINS

      SUBROUTINE phys_state_var_init(ngrid,nlayer,nq,tname, &
                                     day_ini,hour_ini,pdaysec,ptimestep, &
                                     prad,pg,pr,pcpp)

!=======================================================================
!
!   purpose:
!   -------
!
!   Allocate arrays in modules
!   Fill geometrical arrays
!   Fill a first set of physical constants
!   -- was done previously in inifis
!
!=======================================================================
!   
!   authors: Ehouarn Millour and Aymeric Spiga
!            14/04/2014
!
!   arguments:
!   ----------
!
!   input:
!   ------
!
!    ngrid                 Size of the horizontal grid.
!    nlayer                Number of vertical layers.
!    nq                    Number of tracers.
!
!=======================================================================

      use init_print_control_mod, only: init_print_control
      use slope_mod, only: ini_slope_mod,end_slope_mod
      use comsaison_h, only: ini_comsaison_h,end_comsaison_h
      use surfdat_h, only: ini_surfdat_h,end_surfdat_h
      use comgeomfi_h, only: ini_comgeomfi_h,end_comgeomfi_h
      use comsoil_h, only: ini_comsoil_h,end_comsoil_h
      use dimradmars_mod, only: ini_dimradmars_mod,end_dimradmars_mod
      use yomlw_h, only: ini_yomlw_h,end_yomlw_h
      use conc_mod, only: ini_conc_mod,end_conc_mod
      use turb_mod, only: ini_turb_mod,end_turb_mod
      use comcstfi_h, only: pi,rad,cpp,g,r,rcp
      use tracer_mod, only: ini_tracer_mod,end_tracer_mod
      use time_phylmdz_mod, only: init_time
      use co2cloud_mod, only: ini_co2cloud,end_co2cloud
      use compute_dtau_mod, only: ini_compute_dtau_mod, &
                                  end_compute_dtau_mod
      use rocketduststorm_mod, only: ini_rocketduststorm_mod, &
                                     end_rocketduststorm_mod
      use topmons_mod, only: ini_topmons_mod, &
                             end_topmons_mod
      use calchim_mod, only: ini_calchim_mod,end_calchim_mod
      use watercloud_mod, only: ini_watercloud_mod, &
                                end_watercloud_mod
      use nonoro_gwd_ran_mod, only: ini_nonoro_gwd_ran, &
                                    end_nonoro_gwd_ran

      IMPLICIT NONE
      
      INTEGER,INTENT(IN) :: ngrid,nlayer,nq
      CHARACTER(len=*),INTENT(IN) :: tname(nq)
      INTEGER,INTENT(IN) :: day_ini
      REAL,INTENT(IN) :: hour_ini
      REAL,INTENT(IN) :: pdaysec,ptimestep,prad,pg,pr,pcpp

      ! set dimension and allocate arrays in tracer_mod
      call end_tracer_mod
      call ini_tracer_mod(nq,tname)


      ! initialize "print_control" constants/flags ("prt_level","lunout", etc.)
      call init_print_control
      
      ! set parameters in comcstfi_h
      pi=2.*asin(1.) 
      rad=prad
      cpp=pcpp
      g=pg
      r=pr
      rcp=r/cpp

      ! Initialize some "temporal and calendar" related variables
      call init_time(day_ini,hour_ini,pdaysec,ptimestep)

      ! allocate "slope_mod" arrays
      call end_slope_mod
      call ini_slope_mod(ngrid)

      ! allocate "comsaison_h" arrays
      call end_comsaison_h
      call ini_comsaison_h(ngrid)

      ! allocate "surfdat_h" arrays
      call end_surfdat_h
      call ini_surfdat_h(ngrid,nq)

      ! allocate "comgeomfi_h" arrays
      call end_comgeomfi_h
      call ini_comgeomfi_h(ngrid)

      ! allocate "comsoil_h" arrays
      call end_comsoil_h
      call ini_comsoil_h(ngrid)

      ! set some variables in "dimradmars_mod"
      call end_dimradmars_mod
      call ini_dimradmars_mod(ngrid,nlayer)

      ! allocate arrays in "yomlw_h"
      call end_yomlw_h
      call ini_yomlw_h(ngrid)

      ! allocate arrays in "conc_mod" (aeronomars)
      call end_conc_mod
      call ini_conc_mod(ngrid,nlayer)

      ! allocate arrays in "turb_mod"
      call end_turb_mod
      call ini_turb_mod(ngrid,nlayer)
      
      ! allocate arrays in "co2cloud" : 
      ! Memory of the origin of the co2 particles      
      call end_co2cloud
      call ini_co2cloud(ngrid,nlayer)
      
      ! allocate arrays in "compute_dtau_mod":
      call end_compute_dtau_mod
      call ini_compute_dtau_mod(ngrid)

      ! allocate arrays in "rocketduststorm_mod":
      call end_rocketduststorm_mod
      call ini_rocketduststorm_mod(ngrid)

      ! allocate arrays in "topmons_mod":
      call end_topmons_mod
      call ini_topmons_mod(ngrid,nlayer)

      ! allocate arrays in "calchim_mod" (aeronomars)
      call end_calchim_mod
      call ini_calchim_mod(ngrid,nlayer,nq)

      ! allocate arrays in "watercloud_mod"
      call end_watercloud_mod
      call ini_watercloud_mod(ngrid,nlayer,nq)

      ! allocate arrays in "nonoro_gwd_ran_mod"
      call end_nonoro_gwd_ran
      call ini_nonoro_gwd_ran(ngrid,nlayer)


      END SUBROUTINE phys_state_var_init

END MODULE phys_state_var_init_mod
