MODULE geometry_mod

! Store informations concerning the local (to MPI/OpenMP core) geometry

  REAL,SAVE,ALLOCATABLE :: longitude(:) ! longitude of the cell (rad)
!$OMP THREADPRIVATE(longitude)

  REAL,SAVE,ALLOCATABLE :: latitude(:)! latitude of the cell (rad)
!$OMP THREADPRIVATE(latitude)

  REAL,SAVE,ALLOCATABLE :: longitude_deg(:) ! longitude of the cell (degree)
!$OMP THREADPRIVATE(longitude_deg)

  REAL,SAVE,ALLOCATABLE :: latitude_deg(:)! latitude of the cell (degree)
!$OMP THREADPRIVATE(latitude_deg)

  REAL,SAVE,ALLOCATABLE :: boundslon(:,:)  ! boundaries of the cell (rad)
!$OMP THREADPRIVATE(boundslon)

  REAL,SAVE,ALLOCATABLE :: boundslat(:,:) ! boundaries of the cell (rad)
!$OMP THREADPRIVATE(boundslat)

  REAL,SAVE,ALLOCATABLE :: dx(:)      ! resolution of longitude cell (valid only for 2D grid)
!$OMP THREADPRIVATE(dx)
  
  REAL,SAVE,ALLOCATABLE :: dy(:)      ! resolution of latitude cell (valid only for 2D grid)
!$OMP THREADPRIVATE(dy)

  REAL,SAVE,ALLOCATABLE :: cell_area(:)      ! area of the cell
!$OMP THREADPRIVATE(cell_area)

  INTEGER,SAVE,ALLOCATABLE :: ind_cell_glo(:)      ! global index of a local cell
!$OMP THREADPRIVATE(ind_cell_glo)

CONTAINS

  SUBROUTINE init_geometry(klon,longitude_,latitude_, &
                           boundslon_,boundslat_, &
                           cell_area_,ind_cell_glo_,dx_,dy_)
  USE mod_grid_phy_lmdz, ONLY: nvertex
  USE nrtype, ONLY : PI
  IMPLICIT NONE
    INTEGER,INTENT(IN) :: klon ! number of columns for this MPI/OpenMP domain
    REAL,INTENT(IN) :: longitude_(klon)
    REAL,INTENT(IN) :: latitude_(klon)
    REAL,INTENT(IN) :: boundslon_(klon,nvertex)
    REAL,INTENT(IN) :: boundslat_(klon,nvertex)
    REAL,INTENT(IN) :: cell_area_(klon)
    INTEGER,OPTIONAL,INTENT(IN) :: ind_cell_glo_(klon)
    REAL,OPTIONAL,INTENT(IN) :: dx_(klon)
    REAL,OPTIONAL,INTENT(IN) :: dy_(klon)
    
    ALLOCATE(longitude(klon))
    ALLOCATE(latitude(klon))
    ALLOCATE(longitude_deg(klon))
    ALLOCATE(latitude_deg(klon))
    ALLOCATE(boundslon(klon,nvertex))
    ALLOCATE(boundslat(klon,nvertex))
    ALLOCATE(cell_area(klon))
    IF (PRESENT(ind_cell_glo_)) ALLOCATE(ind_cell_glo(klon))
    IF (PRESENT(dx_)) ALLOCATE(dx(klon))
    IF (PRESENT(dy_))ALLOCATE(dy(klon))
    
    longitude(:) = longitude_(:)
    latitude(:) = latitude_(:)
    longitude_deg(:) = longitude(:)*180./PI
    latitude_deg(:) = latitude(:)*180./PI
    boundslon(:,:) = boundslon_(:,:)
    boundslat(:,:) = boundslat_(:,:)
    cell_area(:) = cell_area_(:)
    IF (PRESENT(ind_cell_glo_)) ind_cell_glo(:) = ind_cell_glo_(:)
    IF (PRESENT(dx_)) dx(:) = dx_(:)
    IF (PRESENT(dy_)) dy(:) = dy_(:)
    
  END SUBROUTINE init_geometry


END MODULE geometry_mod

