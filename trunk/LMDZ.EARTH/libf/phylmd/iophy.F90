!
! $Header$
!
module iophy
  
! abd  REAL,private,allocatable,dimension(:),save :: io_lat
! abd  REAL,private,allocatable,dimension(:),save :: io_lon
  REAL,allocatable,dimension(:),save :: io_lat
  REAL,allocatable,dimension(:),save :: io_lon
  INTEGER, save :: phys_domain_id
  
  INTERFACE histwrite_phy
    MODULE PROCEDURE histwrite2d_phy,histwrite3d_phy
  END INTERFACE


contains

  subroutine init_iophy_new(rlat,rlon)
  USE dimphy
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz
  USE ioipsl
  implicit none
  include 'dimensions.h'   
    real,dimension(klon),intent(in) :: rlon
    real,dimension(klon),intent(in) :: rlat

    REAL,dimension(klon_glo)        :: rlat_glo
    REAL,dimension(klon_glo)        :: rlon_glo
    
    INTEGER,DIMENSION(2) :: ddid
    INTEGER,DIMENSION(2) :: dsg
    INTEGER,DIMENSION(2) :: dsl
    INTEGER,DIMENSION(2) :: dpf
    INTEGER,DIMENSION(2) :: dpl
    INTEGER,DIMENSION(2) :: dhs
    INTEGER,DIMENSION(2) :: dhe 
    INTEGER :: i    

    CALL gather(rlat,rlat_glo)
    CALL bcast(rlat_glo)
    CALL gather(rlon,rlon_glo)
    CALL bcast(rlon_glo)
    
!$OMP MASTER  
    ALLOCATE(io_lat(jjm+1-1/iim))
    io_lat(1)=rlat_glo(1)
    io_lat(jjm+1-1/iim)=rlat_glo(klon_glo)
    IF (iim > 1) then
      DO i=2,jjm
        io_lat(i)=rlat_glo(2+(i-2)*iim)
      ENDDO
    ENDIF

    ALLOCATE(io_lon(iim))
    io_lon(:)=rlon_glo(2-1/iim:iim+1-1/iim)

    ddid=(/ 1,2 /)
    dsg=(/ iim, jjm+1-1/iim /)
    dsl=(/ iim, jj_nb /)
    dpf=(/ 1,jj_begin /)
    dpl=(/ iim, jj_end /)
    dhs=(/ ii_begin-1,0 /)
    if (mpi_rank==mpi_size-1) then
      dhe=(/0,0/)
    else
      dhe=(/ iim-ii_end,0 /)  
    endif
    
    call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, &
                      'APPLE',phys_domain_id)

!$OMP END MASTER
      
  end subroutine init_iophy_new

  subroutine init_iophy(lat,lon)
  USE dimphy
  USE mod_phys_lmdz_para
  use ioipsl
  implicit none
  include 'dimensions.h'   
    real,dimension(iim),intent(in) :: lon
    real,dimension(jjm+1-1/iim),intent(in) :: lat

    INTEGER,DIMENSION(2) :: ddid
    INTEGER,DIMENSION(2) :: dsg
    INTEGER,DIMENSION(2) :: dsl
    INTEGER,DIMENSION(2) :: dpf
    INTEGER,DIMENSION(2) :: dpl
    INTEGER,DIMENSION(2) :: dhs
    INTEGER,DIMENSION(2) :: dhe 

!$OMP MASTER  
    allocate(io_lat(jjm+1-1/iim))
    io_lat(:)=lat(:)
    allocate(io_lon(iim))
    io_lon(:)=lon(:)
   
    ddid=(/ 1,2 /)
    dsg=(/ iim, jjm+1-1/iim /)
    dsl=(/ iim, jj_nb /)
    dpf=(/ 1,jj_begin /)
    dpl=(/ iim, jj_end /)
    dhs=(/ ii_begin-1,0 /)
    if (mpi_rank==mpi_size-1) then
      dhe=(/0,0/)
    else
      dhe=(/ iim-ii_end,0 /)  
    endif
    
    call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, &
                      'APPLE',phys_domain_id)

!$OMP END MASTER
      
  end subroutine init_iophy
  
  subroutine histbeg_phy(name,itau0,zjulian,dtime,nhori,nid_day)
  USE dimphy
  USE mod_phys_lmdz_para
  use ioipsl
  use write_field
  implicit none
  include 'dimensions.h'
    
    character*(*), intent(IN) :: name
    integer, intent(in) :: itau0
    real,intent(in) :: zjulian
    real,intent(in) :: dtime
    integer,intent(out) :: nhori
    integer,intent(out) :: nid_day

!$OMP MASTER    
    if (is_sequential) then
      call histbeg(name,iim,io_lon, jj_nb,io_lat(jj_begin:jj_end), &
                   1,iim,1,jj_nb,itau0, zjulian, dtime, nhori, nid_day)
    else
      call histbeg(name,iim,io_lon, jj_nb,io_lat(jj_begin:jj_end), &
                   1,iim,1,jj_nb,itau0, zjulian, dtime, nhori, nid_day,phys_domain_id)
    endif
!$OMP END MASTER
  
  end subroutine histbeg_phy
  
  subroutine histwrite2d_phy(nid,name,itau,field)
  USE dimphy
  USE mod_phys_lmdz_para
  USE ioipsl
  implicit none
  include 'dimensions.h'
    
    integer,intent(in) :: nid
    character*(*), intent(IN) :: name
    integer, intent(in) :: itau
    real,dimension(:),intent(in) :: field
    REAL,dimension(klon_mpi) :: buffer_omp
    INTEGER :: index2d(iim*jj_nb)
    REAL :: Field2d(iim,jj_nb)

    IF (size(field)/=klon) CALL abort_gcm('iophy::histwrite2d','Field first dimension not equal to klon',1)
    
    CALL Gather_omp(field,buffer_omp)    
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,Field2d)
    CALL histwrite(nid,name,itau,Field2d,iim*jj_nb,index2d)
!$OMP END MASTER    
  end subroutine histwrite2d_phy


  
  subroutine histwrite3d_phy(nid,name,itau,field)
  USE dimphy
  USE mod_phys_lmdz_para

  use ioipsl
  implicit none
  include 'dimensions.h'
    
    integer,intent(in) :: nid
    character*(*), intent(IN) :: name
    integer, intent(in) :: itau
    real,dimension(:,:),intent(in) :: field  ! --> field(klon,:)
    REAL,dimension(klon_mpi,size(field,2)) :: buffer_omp
    INTEGER :: nlev
    INTEGER :: index3d(iim*jj_nb*size(field,2))
    REAL :: Field3d(iim,jj_nb,size(field,2))
    
    IF (size(field,1)/=klon) CALL abort_gcm('iophy::histwrite3d','Field first dimension not equal to klon',1)
    nlev=size(field,2)
    
    CALL Gather_omp(field,buffer_omp)
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,field3d)
    CALL histwrite(nid,name,itau,Field3d,iim*jj_nb*nlev,index3d)
!$OMP END MASTER    
  end subroutine histwrite3d_phy
  
  

end module iophy
