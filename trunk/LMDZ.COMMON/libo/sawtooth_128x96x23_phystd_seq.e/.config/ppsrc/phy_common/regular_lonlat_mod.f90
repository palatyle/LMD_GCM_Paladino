










MODULE regular_lonlat_mod

  ! Store information on the global physics grid
  ! for a regular (longitude-latitude) grid
  
  INTEGER, PARAMETER :: north_east=1     ! boundaries of regular lontlat
  INTEGER, PARAMETER :: north_west=2     ! boundaries of regular lontlat
  INTEGER, PARAMETER :: south_west=3     ! boundaries of regular lontlat
  INTEGER, PARAMETER :: south_east=4     ! boundaries of regular lontlat
  INTEGER, PARAMETER :: east=1           ! boundaries of regular lontlat
  INTEGER, PARAMETER :: west=2           ! boundaries of regular lontlat
  INTEGER, PARAMETER :: north=1          ! boundaries of regular lontlat
  INTEGER, PARAMETER :: south=2          ! boundaries of regular lontlat

! Global definition, shared by all threads
! Do not set threadprivate directives

  REAL,SAVE,ALLOCATABLE :: lon_reg(:)      ! value of longitude cell (rad)

  REAL,SAVE,ALLOCATABLE :: lat_reg(:)      ! value of longitude cell (rad)

  REAL,SAVE,ALLOCATABLE :: boundslon_reg(:,:)      ! value of boundaries cell (1=>east, 2=>west)(rad)

  REAL,SAVE,ALLOCATABLE :: boundslat_reg(:,:)      ! value of longitude cell (1=>north, 2=>south)(rad)


CONTAINS
  
  SUBROUTINE init_regular_lonlat(nbp_lon, nbp_lat, lon_reg_, lat_reg_, boundslon_reg_, boundslat_reg_)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nbp_lon
    INTEGER,INTENT(IN) :: nbp_lat
    REAL,   INTENT(IN) :: lon_reg_(nbp_lon)  
    REAL,   INTENT(IN) :: lat_reg_(nbp_lat) 
    REAL,   INTENT(IN) :: boundslon_reg_(nbp_lon,2) 
    REAL,   INTENT(IN) :: boundslat_reg_(nbp_lat,2)  

    
    ALLOCATE(lon_reg(nbp_lon))
    lon_reg(:)=lon_reg_(:)

    ALLOCATE(lat_reg(nbp_lat)) 
    lat_reg(:)=lat_reg_(:)

    ALLOCATE(boundslon_reg(nbp_lon,2)) 
    boundslon_reg(:,:)=boundslon_reg_(:,:)

    ALLOCATE(boundslat_reg(nbp_lat,2))  
    boundslat_reg(:,:)=boundslat_reg_(:,:)

  
  END SUBROUTINE init_regular_lonlat    

END MODULE regular_lonlat_mod

