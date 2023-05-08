!
! $Id: $
!
MODULE callphysiq_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE call_physiq(klon,llm,nqtot,tname,                              &
                       debut_split,lafin_split,                           &
                       jD_cur,jH_cur_split,zdt_split,                     &
                       zplev_omp,zplay_omp,                               &
                       zpk_omp,zphi_omp,zphis_omp,                        &
                       presnivs_omp,                                      &
                       zufi_omp,zvfi_omp,zrfi_omp,ztfi_omp,zqfi_omp,      &
                       flxwfi_omp,pducov,                                 &
                       zdufi_omp,zdvfi_omp,zdtfi_omp,zdqfi_omp,zdpsrf_omp)

  USE control_mod, ONLY: planet_type
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
  USE physiq_mod, ONLY: physiq
!  USE moyzon_mod, ONLY: plevmoy,tmoy ! planetary mean values to send to physics
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: klon ! (local) number of atmospheric columns
  INTEGER,INTENT(IN) :: llm  ! number of atmospheric layers
  INTEGER,INTENT(IN) :: nqtot ! number of tracers
  CHARACTER(len=*),INTENT(IN) :: tname(nqtot) ! tracer names
  LOGICAL,INTENT(IN) :: debut_split ! .true. if very first call to physics
  LOGICAL,INTENT(IN) :: lafin_split ! .true. if last call to physics
  REAL,INTENT(IN) :: JD_cur ! Julian day
  REAL,INTENT(IN) :: JH_cur_split ! Julian hour (fraction of day)
  REAL,INTENT(IN) :: zdt_split ! time step over which the physics are evaluated
  REAL,INTENT(IN) :: zplev_omp(klon,llm+1) ! interlayer pressure (Pa)
  REAL,INTENT(IN) :: zplay_omp(klon,llm) ! mid-layer pressure (Pa)
  REAL,INTENT(IN) :: zpk_omp(klon,llm)
  REAL,INTENT(IN) :: zphi_omp(klon,llm) ! geopotential at midlayer
  REAL,INTENT(IN) :: zphis_omp(klon) ! surface geopotential
  REAL,INTENT(IN) :: presnivs_omp(llm) ! approximate pressure of atm. layers
  REAL,INTENT(IN) :: zufi_omp(klon,llm) ! zonal wind (m/s)
  REAL,INTENT(IN) :: zvfi_omp(klon,llm) ! meridional wind (m/s)
  REAL,INTENT(IN) :: zrfi_omp(klon,llm) ! relative wind vorticity, in s-1
  REAL,INTENT(IN) :: ztfi_omp(klon,llm) ! temperature (K)
  REAL,INTENT(IN) :: zqfi_omp(klon,llm,nqtot) ! tracers (*/kg of air)
  REAL,INTENT(IN) :: flxwfi_omp(klon,llm) ! Vertical mass flux on lower mesh interfaces (kg/s) 
  REAL,INTENT(IN) :: pducov(nbp_lon+1,nbp_lat,llm) ! dynamical tendency on ucov
  ! tendencies (in */s) from the physics:
  REAL,INTENT(OUT) :: zdufi_omp(klon,llm) ! tendency on zonal winds
  REAL,INTENT(OUT) :: zdvfi_omp(klon,llm) ! tendency on meridional winds
  REAL,INTENT(OUT) :: zdtfi_omp(klon,llm) ! tendency on temperature
  REAL,INTENT(OUT) :: zdqfi_omp(klon,llm,nqtot) ! tendency on tracers
  REAL,INTENT(OUT) :: zdpsrf_omp(klon) ! tendency on surface pressure

  ! Local variables
  CHARACTER(len=11) :: modname="call_physiq"
  LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)

! Sanity check on physics package type
  IF (firstcall) THEN
    IF (planet_type.ne."venus") THEN
      CALL abort_gcm(modname,"wrong planet_type for this physics package",1)
    ENDIF
    firstcall=.false.
  ENDIF


! Call physics package with required inputs/outputs
  CALL physiq(klon,           &
              llm,            &
              nqtot,          &
              debut_split,    &
              lafin_split,    &
              jD_cur,         &
              jH_cur_split,   &
              zdt_split,      &
              zplev_omp,      &
              zplay_omp,      &
              zpk_omp,        &
              zphi_omp,       &
              zphis_omp,      &
              presnivs_omp,   &
              zufi_omp,       &
              zvfi_omp,       &
              ztfi_omp,       &
              zqfi_omp,       &
              flxwfi_omp,     &
!              plevmoy,        & ! planet-averaged mean pressure (Pa) at interfaces
!              tmoy,           & ! planet-averaged mean temperature (K) at mid-layers
              zdufi_omp,      &
              zdvfi_omp,      &
              zdtfi_omp,      &
              zdqfi_omp,      &
              zdpsrf_omp)


END SUBROUTINE call_physiq

END MODULE callphysiq_mod
