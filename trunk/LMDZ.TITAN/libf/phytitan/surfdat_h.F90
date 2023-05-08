
       module surfdat_h

       implicit none

       real,allocatable,dimension(:) :: albedodat ! albedo of bare ground stocked in startfi.nc file.
!$OMP THREADPRIVATE(albedodat)
       ! Ehouarn: moved inertiedat to comsoil.h
       !      real inertiedat, ! thermal inertia
       real,allocatable,dimension(:) :: phisfi ! geopotential at ground level
!$OMP THREADPRIVATE(phisfi)
       real,dimension(2) :: emisice ! ice emissivity; 1:Northern hemisphere 2:Southern hemisphere
       real emissiv
       real,dimension(2) :: iceradius, dtemisice
!$OMP THREADPRIVATE(emisice,emissiv,iceradius,dtemisice)
       real,allocatable,dimension(:) :: zmea,zstd,zsig,zgam,zthe
!$OMP THREADPRIVATE(zmea,zstd,zsig,zgam,zthe)

       end module surfdat_h

