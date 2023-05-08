MODULE inigeomphy_mod

CONTAINS

SUBROUTINE inigeomphy(iim,jjm,nlayer, &
                     nbp, communicator, &
                     rlatu,rlatv,rlonu,rlonv,aire,cu,cv)
  USE mod_grid_phy_lmdz, ONLY: klon_glo,  & ! number of atmospheric columns (on full grid)
                               regular_lonlat, &  ! regular longitude-latitude grid type
                               nbp_lon, nbp_lat, nbp_lev
  USE mod_phys_lmdz_para, ONLY: klon_omp, & ! number of columns (on local omp grid)
                                klon_omp_begin, & ! start index of local omp subgrid
                                klon_omp_end, & ! end index of local omp subgrid
                                klon_mpi_begin ! start indes of columns (on local mpi grid)
  USE geometry_mod, ONLY : init_geometry
  USE physics_distribution_mod, ONLY : init_physics_distribution
  USE regular_lonlat_mod, ONLY : init_regular_lonlat, &
                                 east, west, north, south, &
                                 north_east, north_west, &
                                 south_west, south_east
  USE mod_interface_dyn_phys, ONLY :  init_interface_dyn_phys
  USE nrtype, ONLY: pi
  USE comvert_mod, ONLY: preff, ap, bp, aps, bps, presnivs, &
                         scaleheight, pseudoalt
  USE vertical_layers_mod, ONLY: init_vertical_layers
  IMPLICIT NONE

  ! =======================================================================
  ! Initialisation of the physical constants and some positional and
  ! geometrical arrays for the physics
  ! =======================================================================

  include "iniprint.h"

  INTEGER, INTENT (IN) :: nlayer ! number of atmospheric layers
  INTEGER, INTENT (IN) :: iim ! number of atmospheric columns along longitudes
  INTEGER, INTENT (IN) :: jjm ! number of atompsheric columns along latitudes
  INTEGER, INTENT(IN) :: nbp ! number of physics columns for this MPI process
  INTEGER, INTENT(IN) :: communicator ! MPI communicator
  REAL, INTENT (IN) :: rlatu(jjm+1) ! latitudes of the physics grid
  REAL, INTENT (IN) :: rlatv(jjm) ! latitude boundaries of the physics grid
  REAL, INTENT (IN) :: rlonv(iim+1) ! longitudes of the physics grid
  REAL, INTENT (IN) :: rlonu(iim+1) ! longitude boundaries of the physics grid
  REAL, INTENT (IN) :: aire(iim+1,jjm+1) ! area of the dynamics grid (m2)
  REAL, INTENT (IN) :: cu((iim+1)*(jjm+1)) ! cu coeff. (u_covariant = cu * u)
  REAL, INTENT (IN) :: cv((iim+1)*jjm) ! cv coeff. (v_covariant = cv * v)

  INTEGER :: ibegin, iend, offset
  INTEGER :: i,j,k
  CHARACTER (LEN=20) :: modname = 'iniphysiq'
  CHARACTER (LEN=80) :: abort_message
  REAL :: total_area_phy, total_area_dyn

  ! boundaries, on global grid
  REAL,ALLOCATABLE :: boundslon_reg(:,:)
  REAL,ALLOCATABLE :: boundslat_reg(:,:)

  ! global array, on full physics grid:
  REAL,ALLOCATABLE :: latfi_glo(:)
  REAL,ALLOCATABLE :: lonfi_glo(:)
  REAL,ALLOCATABLE :: cufi_glo(:)
  REAL,ALLOCATABLE :: cvfi_glo(:)
  REAL,ALLOCATABLE :: airefi_glo(:)
  REAL,ALLOCATABLE :: boundslonfi_glo(:,:)
  REAL,ALLOCATABLE :: boundslatfi_glo(:,:)

  ! local arrays, on given MPI/OpenMP domain:
  REAL,ALLOCATABLE,SAVE :: latfi(:)
  REAL,ALLOCATABLE,SAVE :: lonfi(:)
  REAL,ALLOCATABLE,SAVE :: cufi(:)
  REAL,ALLOCATABLE,SAVE :: cvfi(:)
  REAL,ALLOCATABLE,SAVE :: airefi(:)
  REAL,ALLOCATABLE,SAVE :: boundslonfi(:,:)
  REAL,ALLOCATABLE,SAVE :: boundslatfi(:,:)
  INTEGER,ALLOCATABLE,SAVE :: ind_cell_glo_fi(:)
!$OMP THREADPRIVATE (latfi,lonfi,cufi,cvfi,airefi,boundslonfi,boundslatfi,ind_cell_glo_fi)

  ! Initialize Physics distibution and parameters and interface with dynamics
  IF (iim*jjm>1) THEN ! general 3D case
    CALL init_physics_distribution(regular_lonlat,4, &
                                 nbp,iim,jjm+1,nlayer,communicator)
  ELSE ! For 1D model
    CALL init_physics_distribution(regular_lonlat,4, &
                                 1,1,1,nlayer,communicator)
  ENDIF
  CALL init_interface_dyn_phys
  
  ! init regular global longitude-latitude grid points and boundaries
  ALLOCATE(boundslon_reg(iim,2))
  ALLOCATE(boundslat_reg(jjm+1,2))
  
  DO i=1,iim
   boundslon_reg(i,east)=rlonu(i) 
   boundslon_reg(i,west)=rlonu(i+1) 
  ENDDO

  boundslat_reg(1,north)= PI/2 
  boundslat_reg(1,south)= rlatv(1)
  DO j=2,jjm
   boundslat_reg(j,north)=rlatv(j-1) 
   boundslat_reg(j,south)=rlatv(j) 
  ENDDO
  boundslat_reg(jjm+1,north)= rlatv(jjm) 
  boundslat_reg(jjm+1,south)= -PI/2

  ! Write values in module regular_lonlat_mod
  CALL init_regular_lonlat(iim,jjm+1, rlonv(1:iim), rlatu, &
                           boundslon_reg, boundslat_reg)

  ! Generate global arrays on full physics grid
  ALLOCATE(latfi_glo(klon_glo),lonfi_glo(klon_glo))
  ALLOCATE(cufi_glo(klon_glo),cvfi_glo(klon_glo))
  ALLOCATE(airefi_glo(klon_glo))
  ALLOCATE(boundslonfi_glo(klon_glo,4))
  ALLOCATE(boundslatfi_glo(klon_glo,4))

  IF (klon_glo>1) THEN ! general case
    ! North pole
    latfi_glo(1)=rlatu(1)
    lonfi_glo(1)=0.
    cufi_glo(1) = cu(1)
    cvfi_glo(1) = cv(1)
    boundslonfi_glo(1,north_east)=0
    boundslatfi_glo(1,north_east)=PI/2
    boundslonfi_glo(1,north_west)=2*PI
    boundslatfi_glo(1,north_west)=PI/2
    boundslonfi_glo(1,south_west)=2*PI
    boundslatfi_glo(1,south_west)=rlatv(1)
    boundslonfi_glo(1,south_east)=0
    boundslatfi_glo(1,south_east)=rlatv(1)
    DO j=2,jjm
      DO i=1,iim
        k=(j-2)*iim+1+i
        latfi_glo(k)= rlatu(j)
        lonfi_glo(k)= rlonv(i)
        cufi_glo(k) = cu((j-1)*(iim+1)+i)
        cvfi_glo(k) = cv((j-1)*(iim+1)+i)
        boundslonfi_glo(k,north_east)=rlonu(i)
        boundslatfi_glo(k,north_east)=rlatv(j-1)
        boundslonfi_glo(k,north_west)=rlonu(i+1)
        boundslatfi_glo(k,north_west)=rlatv(j-1)
        boundslonfi_glo(k,south_west)=rlonu(i+1)
        boundslatfi_glo(k,south_west)=rlatv(j)
        boundslonfi_glo(k,south_east)=rlonu(i)
        boundslatfi_glo(k,south_east)=rlatv(j)
      ENDDO
    ENDDO
    ! South pole
    latfi_glo(klon_glo)= rlatu(jjm+1)
    lonfi_glo(klon_glo)= 0.
    cufi_glo(klon_glo) = cu((iim+1)*jjm+1)
    cvfi_glo(klon_glo) = cv((iim+1)*jjm-iim)
    boundslonfi_glo(klon_glo,north_east)= 0
    boundslatfi_glo(klon_glo,north_east)= rlatv(jjm)
    boundslonfi_glo(klon_glo,north_west)= 2*PI
    boundslatfi_glo(klon_glo,north_west)= rlatv(jjm)
    boundslonfi_glo(klon_glo,south_west)= 2*PI
    boundslatfi_glo(klon_glo,south_west)= -PI/2
    boundslonfi_glo(klon_glo,south_east)= 0
    boundslatfi_glo(klon_glo,south_east)= -Pi/2

    ! build airefi(), mesh area on physics grid
    CALL gr_dyn_fi(1,iim+1,jjm+1,klon_glo,aire,airefi_glo)
    ! Poles are single points on physics grid
    airefi_glo(1)=sum(aire(1:iim,1))
    airefi_glo(klon_glo)=sum(aire(1:iim,jjm+1))

    ! Sanity check: do total planet area match between physics and dynamics?
    total_area_dyn=sum(aire(1:iim,1:jjm+1))
    total_area_phy=sum(airefi_glo(1:klon_glo))
    IF (total_area_dyn/=total_area_phy) THEN
      WRITE (lunout, *) 'iniphysiq: planet total surface discrepancy !!!'
      WRITE (lunout, *) '     in the dynamics total_area_dyn=', total_area_dyn
      WRITE (lunout, *) '  but in the physics total_area_phy=', total_area_phy
      IF (abs(total_area_dyn-total_area_phy)>0.00001*total_area_dyn) THEN
        ! stop here if the relative difference is more than 0.001%
        abort_message = 'planet total surface discrepancy'
        CALL abort_gcm(modname, abort_message, 1)
      ENDIF
    ENDIF
  ELSE ! klon_glo==1, running the 1D model
    ! just copy over input values
    latfi_glo(1)=rlatu(1)
    lonfi_glo(1)=rlonv(1)
    cufi_glo(1)=cu(1)
    cvfi_glo(1)=cv(1)
    airefi_glo(1)=aire(1,1)
    boundslonfi_glo(1,north_east)=rlonu(1)
    boundslatfi_glo(1,north_east)=PI/2
    boundslonfi_glo(1,north_west)=rlonu(2)
    boundslatfi_glo(1,north_west)=PI/2
    boundslonfi_glo(1,south_west)=rlonu(2)
    boundslatfi_glo(1,south_west)=rlatv(1)
    boundslonfi_glo(1,south_east)=rlonu(1)
    boundslatfi_glo(1,south_east)=rlatv(1)
  ENDIF ! of IF (klon_glo>1)

!$OMP PARALLEL 
  ! Now generate local lon/lat/cu/cv/area/bounds arrays
  ALLOCATE(latfi(klon_omp),lonfi(klon_omp),cufi(klon_omp),cvfi(klon_omp))
  ALLOCATE(airefi(klon_omp))
  ALLOCATE(boundslonfi(klon_omp,4))
  ALLOCATE(boundslatfi(klon_omp,4))
  ALLOCATE(ind_cell_glo_fi(klon_omp))


  offset = klon_mpi_begin - 1
  airefi(1:klon_omp) = airefi_glo(offset+klon_omp_begin:offset+klon_omp_end)
  cufi(1:klon_omp) = cufi_glo(offset+klon_omp_begin:offset+klon_omp_end)
  cvfi(1:klon_omp) = cvfi_glo(offset+klon_omp_begin:offset+klon_omp_end)
  lonfi(1:klon_omp) = lonfi_glo(offset+klon_omp_begin:offset+klon_omp_end)
  latfi(1:klon_omp) = latfi_glo(offset+klon_omp_begin:offset+klon_omp_end)
  boundslonfi(1:klon_omp,:) = boundslonfi_glo(offset+klon_omp_begin:offset+klon_omp_end,:)
  boundslatfi(1:klon_omp,:) = boundslatfi_glo(offset+klon_omp_begin:offset+klon_omp_end,:)
  ind_cell_glo_fi(1:klon_omp)=(/ (i,i=offset+klon_omp_begin,offset+klon_omp_end) /)

  ! copy over local grid longitudes and latitudes
  CALL init_geometry(klon_omp,lonfi,latfi,boundslonfi,boundslatfi, &
                     airefi,ind_cell_glo_fi,cufi,cvfi)

  ! copy over preff , ap(), bp(), etc 
  CALL init_vertical_layers(nlayer,preff,scaleheight, &
                            ap,bp,aps,bps,presnivs,pseudoalt)

!$OMP END PARALLEL


END SUBROUTINE inigeomphy

END MODULE inigeomphy_mod
