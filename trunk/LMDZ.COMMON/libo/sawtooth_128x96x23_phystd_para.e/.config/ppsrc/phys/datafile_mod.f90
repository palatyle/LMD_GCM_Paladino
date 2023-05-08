










!-----------------------------------------------------------------------
      module datafile_mod
!  Address of the directory containing tables of data needed by the GCM
      implicit none

      ! Main directory: 'datadir':
      ! Default for Berserker @ UChicago:
!      character(len=300) :: datadir='/home/rwordsworth/datagcm'
      ! Default for Gnome Idataplex:
!      character(len=300) :: datadir='/san/home/rdword/gcm/datagcm'
      ! Default for LMD machines:
      character(len=300),save :: datadir='/home/palatyle/LMD_gen/trunk/datadir'

!$OMP THREADPRIVATE(datadir)
      
      ! Subdirectories of 'datadir':
      
      ! surfdir stores planetary topography, albedo, etc. (surface.nc files)
      character(len=12),parameter :: surfdir="surface_data"
      
      ! aerdir stores aerosol properties files (optprop_*dat files)
      character(LEN=18),parameter :: aerdir="aerosol_properties" 

      end module datafile_mod
!-----------------------------------------------------------------------
