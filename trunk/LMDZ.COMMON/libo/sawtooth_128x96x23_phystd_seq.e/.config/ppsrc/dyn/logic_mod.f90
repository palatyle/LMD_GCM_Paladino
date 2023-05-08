










MODULE logic_mod

IMPLICIT NONE  

  LOGICAL purmats ! true if time stepping is purely Matsuno scheme
                  ! false implies Matsuno-Leapfrog time stepping scheme
  LOGICAL forward ! true if during forward phase of Matsuno step
  LOGICAL leapf ! true if during a leapfrog time stepping step
  LOGICAL apphys ! true if during a time step when physics will be called
  LOGICAL statcl
  LOGICAL conser
  LOGICAL apdiss ! true if during a time step when dissipation will be called
  LOGICAL apdelq
  LOGICAL saison
  LOGICAL ecripar
  LOGICAL fxyhypb ! true if using hyperbolic function discretization
                  ! for latitudinal grid 
  LOGICAL ysinus ! true if using sine function discretiation
                 ! for latitudinal grid
  LOGICAL read_start ! true if reading a start.nc file to initialize fields
  LOGICAL ok_guide ! true if nudging
  LOGICAL ok_strato
  LOGICAL tidal  ! true if adding tidal forces (for Titan)
  LOGICAL ok_gradsfile
  LOGICAL ok_limit  ! true for boundary conditions file creation (limit.nc)
  LOGICAL ok_etat0  ! true for initial states creation (start.nc, startphy.nc)
  LOGICAL read_orop ! true for sub-cell scales orographic params read in file
  LOGICAL hybrid ! vertical coordinate is hybrid if true (sigma otherwise)
                 ! (only used if disvert_type==2)
  LOGICAL moyzon_mu,moyzon_ch ! used for zonal averages in Titan

  INTEGER iflag_phys ! type of physics to call: 0 none, 1: phy*** package,
                     ! 2: Held & Suarez, 101-200: aquaplanets & terraplanets
  INTEGER iflag_trac

END MODULE logic_mod
