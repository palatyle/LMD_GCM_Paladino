!-----------------------------------------------------------------------
      module datafile_mod
!  Address of the directory containing tables of data needed by the GCM
      implicit none

      ! Main directory: 'datadir':
      ! Default for LMD machines:
      character(len=300),save :: datadir='datagcm'
!$OMP THREADPRIVATE(datadir)
      
      ! Subdirectories of 'datadir':
      
      ! surfdir stores planetary topography, albedo, etc. (surface.nc files)
      character(len=12),parameter :: surfdir="surface_data"
      
      ! Default directories for correlated-k data
      ! Read in inifis_mod, set in physiq_mod
      character(LEN=100),save :: corrkdir = 'datagcm/corrk_data/50X50'
      character(LEN=100),save :: banddir  = 'datagcm/corrk_data/50X50/50x50'
!$OMP THREADPRIVATE(corrkdir,banddir)
      
      ! Default directory for microphysics
      ! Set in inifis_mod
      character(LEN=100),save :: config_mufi ='datagcm/microphysics/config.cfg'
!$OMP THREADPRIVATE(config_mufi)

      ! Default file for coupled microphysics optical properties
      ! Set in physiq_mod
      character(LEN=100),save :: haze_opt_file ='datagcm/HAZE_OPT_23x23.DAT'
!$OMP THREADPRIVATE(haze_opt_file)

      end module datafile_mod
!-----------------------------------------------------------------------
