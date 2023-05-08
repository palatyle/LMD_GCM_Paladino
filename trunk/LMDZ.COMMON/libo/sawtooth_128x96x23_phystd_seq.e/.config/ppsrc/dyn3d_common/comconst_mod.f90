










MODULE comconst_mod

IMPLICIT NONE  

      INTEGER im,jm,lllm,imp1,jmp1,lllmm1,lllmp1,lcl
      REAL dtvr ! dynamical time step (in s)
      REAL daysec !length (in s) of a standard day
      REAL pi    ! something like 3.14159....
      REAL dtphys ! (s) time step for the physics
      REAL dtdiss ! (s) time step for the dissipation
      REAL rad ! (m) radius of the planet
      REAL r ! Reduced Gas constant r=R/mu
             ! with R=8.31.. J.K-1.mol-1, mu: mol mass of atmosphere (kg/mol)
      REAL cpp   ! Specific heat Cp (J.kg-1.K-1)
      REAL kappa ! kappa=r/Cp 
      REAL cotot
      REAL unsim ! = 1./iim
      REAL g ! (m/s2) gravity
      REAL omeg ! (rad/s) rotation rate of the planet
! Dissipation factors, for Earth model:
      REAL dissip_factz,dissip_zref !dissip_deltaz
! Dissipation factors, for other planets:
      REAL dissip_fac_mid,dissip_fac_up,dissip_deltaz,dissip_hdelta
      REAL dissip_pupstart
! top_bound sponge:
      INTEGER iflag_top_bound ! sponge type
      INTEGER ngroup ! parameter to group points (along longitude) near poles
      INTEGER mode_top_bound  ! sponge mode
      REAL tau_top_bound ! inverse of sponge characteristic time scale (Hz)
      REAL daylen ! length of solar day, in 'standard' day length
      REAL year_day ! Number of standard days in a year
      REAL molmass ! (g/mol) molar mass of the atmosphere

      REAL nu_venus,t0_venus ! coeffs needed for Cp(T), Venus atmosphere
      REAL ihf  ! (W/m2) intrinsic heat flux for giant planets


END MODULE comconst_mod
