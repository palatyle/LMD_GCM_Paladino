!
! $Header$
!
MODULE surface_data

  REAL, PARAMETER        :: calice=1.0/(5.1444e+06*0.15)
  REAL, PARAMETER        :: tau_gl=86400.*5.
  REAL, PARAMETER        :: calsno=1./(2.3867e+06*.15)
  
  LOGICAL, SAVE          :: ok_veget      ! true for use of vegetation model ORCHIDEE
  !$OMP THREADPRIVATE(ok_veget)

  CHARACTER(len=6), SAVE :: type_ocean    ! force/slab/couple
  !$OMP THREADPRIVATE(type_ocean)

  ! if type_ocean=couple : version_ocean=opa8 ou nemo
  ! if type_ocean=slab   : version_ocean=sicOBS
  CHARACTER(len=6), SAVE :: version_ocean 
  !$OMP THREADPRIVATE(version_ocean)

END MODULE surface_data
