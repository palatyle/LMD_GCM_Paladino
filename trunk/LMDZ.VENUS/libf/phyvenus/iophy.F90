!
! $Header$
!
module iophy
  
! abd  REAL,private,allocatable,dimension(:),save :: io_lat
! abd  REAL,private,allocatable,dimension(:),save :: io_lon
  REAL,allocatable,dimension(:),save :: io_lat
  REAL,allocatable,dimension(:),save :: io_lon
  INTEGER, save :: phys_domain_id
  INTEGER, save :: npstn
  INTEGER, allocatable, dimension(:), save :: nptabij
  
  INTERFACE histwrite_phy
    MODULE PROCEDURE histwrite2d_phy,histwrite3d_phy
  END INTERFACE

  INTERFACE histbeg_phy_all
    MODULE PROCEDURE histbeg_phy,histbeg_phy_points
  END INTERFACE


contains

  subroutine init_iophy_new(rlat,rlon)
  USE dimphy, only: klon
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat, klon_glo
  USE ioipsl
  implicit none
    real,dimension(klon),intent(in) :: rlon ! longitudes, in degrees
    real,dimension(klon),intent(in) :: rlat ! latitudes, in degrees

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
    ALLOCATE(io_lat(nbp_lat))
    IF (klon_glo == 1) THEN
      io_lat(1)=rlat_glo(1)
    ELSE
      io_lat(1)=rlat_glo(1)
      io_lat(nbp_lat)=rlat_glo(klon_glo)
      DO i=2,nbp_lat-1
        io_lat(i)=rlat_glo(2+(i-2)*nbp_lon)
      ENDDO
    ENDIF

    ALLOCATE(io_lon(nbp_lon))
    IF (klon_glo == 1) THEN
      io_lon(1)=rlon_glo(1)
    ELSE 
      io_lon(1:nbp_lon)=rlon_glo(2:nbp_lon+1)
    ENDIF

! POUR VENUS, LA CARTE EST A L ENVERS !!
    io_lat = -1.*io_lat
    io_lon = -1.*io_lon

    ddid=(/ 1,2 /)
    dsg=(/ nbp_lon, nbp_lat /)
    dsl=(/ nbp_lon, jj_nb /)
    dpf=(/ 1,jj_begin /)
    dpl=(/ nbp_lon, jj_end /)
    dhs=(/ ii_begin-1,0 /)
    IF (mpi_rank==mpi_size-1) THEN
      dhe=(/0,0/)
    ELSE
      dhe=(/ nbp_lon-ii_end,0 /)  
    ENDIF
    
    call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, &
                      'APPLE',phys_domain_id)

!$OMP END MASTER
      
  end subroutine init_iophy_new

  subroutine init_iophy(lat,lon)
  USE dimphy
  USE mod_phys_lmdz_para
  use ioipsl, only: flio_dom_set
  USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat
  implicit none
    real,dimension(nbp_lon),intent(in) :: lon
    real,dimension(nbp_lat),intent(in) :: lat

    INTEGER,DIMENSION(2) :: ddid
    INTEGER,DIMENSION(2) :: dsg
    INTEGER,DIMENSION(2) :: dsl
    INTEGER,DIMENSION(2) :: dpf
    INTEGER,DIMENSION(2) :: dpl
    INTEGER,DIMENSION(2) :: dhs
    INTEGER,DIMENSION(2) :: dhe 

!$OMP MASTER  
    allocate(io_lat(nbp_lat))
    io_lat(:)=lat(:)
    allocate(io_lon(nbp_lon))
    io_lon(:)=lon(:)
   
    ddid=(/ 1,2 /)
    dsg=(/ nbp_lon, nbp_lat /)
    dsl=(/ nbp_lon, jj_nb /)
    dpf=(/ 1,jj_begin /)
    dpl=(/ nbp_lon, jj_end /)
    dhs=(/ ii_begin-1,0 /)
    if (mpi_rank==mpi_size-1) then
      dhe=(/0,0/)
    else
      dhe=(/ nbp_lon-ii_end,0 /)  
    endif
    
    call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, &
                      'APPLE',phys_domain_id)

!$OMP END MASTER
      
  end subroutine init_iophy
  
  subroutine histbeg_phy(name,itau0,zjulian,dtime,nhori,nid_day)
  USE mod_phys_lmdz_para, only: jj_begin, jj_end, jj_nb, is_sequential
  USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat
  use ioipsl, only: histbeg
  use write_field
  implicit none
    
    character*(*), intent(IN) :: name
    integer, intent(in) :: itau0
    real,intent(in) :: zjulian
    real,intent(in) :: dtime
    integer,intent(out) :: nhori
    integer,intent(out) :: nid_day

!$OMP MASTER    
    if (is_sequential) then
      call histbeg(name,nbp_lon,io_lon, jj_nb,io_lat(jj_begin:jj_end), &
                   1,nbp_lon,1,jj_nb,itau0, zjulian, dtime, nhori, nid_day)
    else
      call histbeg(name,nbp_lon,io_lon, jj_nb,io_lat(jj_begin:jj_end), &
                   1,nbp_lon,1,jj_nb,itau0, zjulian, dtime, nhori, nid_day,phys_domain_id)
    endif
!$OMP END MASTER
  
  end subroutine histbeg_phy

  subroutine histbeg_phy_points(rlon,rlat,pim,tabij,ipt,jpt, &
             plon,plat,plon_bounds,plat_bounds, &
             nname,itau0,zjulian,dtime,nnhori,nnid_day)
  USE dimphy, only: klon
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz, only: klon_glo, nbp_lon, nbp_lat
  use ioipsl, only: histbeg

  implicit none

    real,dimension(klon),intent(in) :: rlon
    real,dimension(klon),intent(in) :: rlat
    integer, intent(in) :: itau0
    real,intent(in) :: zjulian
    real,intent(in) :: dtime
    integer, intent(in) :: pim
    integer, intent(out) :: nnhori
    character(len=20), intent(in) :: nname
    INTEGER, intent(out) :: nnid_day
    integer :: i
    REAL,dimension(klon_glo)        :: rlat_glo
    REAL,dimension(klon_glo)        :: rlon_glo
    INTEGER, DIMENSION(pim), intent(in)  :: tabij
    REAL,dimension(pim), intent(in) :: plat, plon
    INTEGER,dimension(pim), intent(in) :: ipt, jpt
    REAL,dimension(pim,2), intent(out) :: plat_bounds, plon_bounds

    INTEGER, SAVE :: tabprocbeg, tabprocend
!$OMP THREADPRIVATE(tabprocbeg, tabprocend)
    INTEGER :: ip
    INTEGER, PARAMETER :: nip=1
    INTEGER :: npproc
    REAL, allocatable, dimension(:) :: npplat, npplon
    REAL, allocatable, dimension(:,:) :: npplat_bounds, npplon_bounds
    REAL, dimension(nbp_lon,nbp_lat) :: zx_lon, zx_lat

    CALL gather(rlat,rlat_glo)
    CALL bcast(rlat_glo)
    CALL gather(rlon,rlon_glo)
    CALL bcast(rlon_glo)

!$OMP MASTER
    DO i=1,pim

!    print*,'CFMIP_iophy i tabij lon lat',i,tabij(i),plon(i),plat(i)

     plon_bounds(i,1)=rlon_glo(tabij(i)-1)
     plon_bounds(i,2)=rlon_glo(tabij(i)+1)
     if(plon_bounds(i,2).LE.0..AND.plon_bounds(i,1).GE.0.) THEN
      if(rlon_glo(tabij(i)).GE.0.) THEN
       plon_bounds(i,2)=-1*plon_bounds(i,2)
      endif
     endif
     if(plon_bounds(i,2).GE.0..AND.plon_bounds(i,1).LE.0.) THEN
      if(rlon_glo(tabij(i)).LE.0.) THEN
       plon_bounds(i,2)=-1*plon_bounds(i,2)
      endif
     endif
!
     IF ( tabij(i).LE.nbp_lon) THEN
      plat_bounds(i,1)=rlat_glo(tabij(i))
     ELSE
      plat_bounds(i,1)=rlat_glo(tabij(i)-nbp_lon)
     ENDIF
     plat_bounds(i,2)=rlat_glo(tabij(i)+nbp_lon)
!
!    print*,'CFMIP_iophy point i lon lon_bds',i,plon_bounds(i,1),rlon_glo(tabij(i)),plon_bounds(i,2) 
!    print*,'CFMIP_iophy point i lat lat_bds',i,plat_bounds(i,1),rlat_glo(tabij(i)),plat_bounds(i,2) 
!
    ENDDO
    if (is_sequential) then

     npstn=pim
     IF(.NOT. ALLOCATED(nptabij)) THEN
      ALLOCATE(nptabij(pim))
     ENDIF 
     DO i=1,pim
      nptabij(i)=tabij(i)
     ENDDO

       CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,rlon_glo,zx_lon)
       if ((nbp_lon*nbp_lat).gt.1) then
       DO i = 1, nbp_lon
         zx_lon(i,1) = rlon_glo(i+1)
         zx_lon(i,nbp_lat) = rlon_glo(i+1)
       ENDDO
       endif
       CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,rlat_glo,zx_lat)

    DO i=1,pim
!    print*,'CFMIP_iophy i tabij lon lat',i,tabij(i),plon(i),plat(i)

     plon_bounds(i,1)=zx_lon(ipt(i)-1,jpt(i))
     plon_bounds(i,2)=zx_lon(ipt(i)+1,jpt(i))

     if (ipt(i).EQ.1) then
      plon_bounds(i,1)=zx_lon(nbp_lon,jpt(i))
      plon_bounds(i,2)=360.+zx_lon(ipt(i)+1,jpt(i))
     endif
 
     if (ipt(i).EQ.nbp_lon) then
      plon_bounds(i,2)=360.+zx_lon(1,jpt(i))
     endif

     plat_bounds(i,1)=zx_lat(ipt(i),jpt(i)-1)
     plat_bounds(i,2)=zx_lat(ipt(i),jpt(i)+1)

     if (jpt(i).EQ.1) then
      plat_bounds(i,1)=zx_lat(ipt(i),1)+0.001
      plat_bounds(i,2)=zx_lat(ipt(i),1)-0.001
     endif
 
     if (jpt(i).EQ.nbp_lat) then
      plat_bounds(i,1)=zx_lat(ipt(i),nbp_lat)+0.001
      plat_bounds(i,2)=zx_lat(ipt(i),nbp_lat)-0.001
     endif
!
!    print*,'CFMIP_iophy point i lon lon_bds',i,plon_bounds(i,1),rlon(tabij(i)),plon_bounds(i,2) 
!    print*,'CFMIP_iophy point i lat lat_bds',i,plat_bounds(i,1),rlat(tabij(i)),plat_bounds(i,2) 
!
    ENDDO
!    print*,'iophy is_sequential nname, nnhori, nnid_day=',trim(nname),nnhori,nnid_day
     call histbeg(nname,pim,plon,plon_bounds, & 
                           plat,plat_bounds, &
                           itau0, zjulian, dtime, nnhori, nnid_day)
    else
     npproc=0
     DO ip=1, pim
      tabprocbeg=klon_mpi_begin
      tabprocend=klon_mpi_end
      IF(tabij(ip).GE.tabprocbeg.AND.tabij(ip).LE.tabprocend) THEN
       npproc=npproc+1
       npstn=npproc
      ENDIF 
     ENDDO
!    print*,'CFMIP_iophy mpi_rank npstn',mpi_rank,npstn
     IF(.NOT. ALLOCATED(nptabij)) THEN
      ALLOCATE(nptabij(npstn))
      ALLOCATE(npplon(npstn), npplat(npstn))
      ALLOCATE(npplon_bounds(npstn,2), npplat_bounds(npstn,2))
     ENDIF
     npproc=0
     DO ip=1, pim
      IF(tabij(ip).GE.tabprocbeg.AND.tabij(ip).LE.tabprocend) THEN
       npproc=npproc+1
       nptabij(npproc)=tabij(ip)
!      print*,'mpi_rank npproc ip plon plat tabij=',mpi_rank,npproc,ip, &
!      plon(ip),plat(ip),tabij(ip)
       npplon(npproc)=plon(ip)
       npplat(npproc)=plat(ip)
       npplon_bounds(npproc,1)=plon_bounds(ip,1)
       npplon_bounds(npproc,2)=plon_bounds(ip,2)
       npplat_bounds(npproc,1)=plat_bounds(ip,1)
       npplat_bounds(npproc,2)=plat_bounds(ip,2)
!!!
!!! print qui sert a reordonner les points stations selon l'ordre CFMIP
!!! ne pas enlever
        print*,'iophy_mpi rank ip lon lat',mpi_rank,ip,plon(ip),plat(ip)
!!!
      ENDIF
     ENDDO
     call histbeg(nname,npstn,npplon,npplon_bounds, &
                            npplat,npplat_bounds, &
                            itau0,zjulian,dtime,nnhori,nnid_day,phys_domain_id)
    endif
!$OMP END MASTER

  end subroutine histbeg_phy_points
 
  subroutine histwrite2d_phy(nid,lpoint,name,itau,field)
  USE dimphy, only: klon
  USE mod_phys_lmdz_para
  USE ioipsl, only: histwrite
  USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat
  implicit none
    
    integer,intent(in) :: nid
    logical,intent(in) :: lpoint 
    character*(*), intent(IN) :: name
    integer, intent(in) :: itau
    real,dimension(:),intent(in) :: field
    REAL,dimension(klon_mpi) :: buffer_omp
    INTEGER, allocatable, dimension(:) :: index2d
    REAL :: Field2d(nbp_lon,jj_nb)

    integer :: ip
    real,allocatable,dimension(:) :: fieldok

    IF (size(field)/=klon) CALL abort_physic('iophy::histwrite2d','Field first dimension not equal to klon',1)
    
    CALL Gather_omp(field,buffer_omp)    
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,Field2d)
    if(.NOT.lpoint) THEN
     ALLOCATE(index2d(nbp_lon*jj_nb))
     ALLOCATE(fieldok(nbp_lon*jj_nb))
     CALL histwrite(nid,name,itau,Field2d,nbp_lon*jj_nb,index2d)
    else
     ALLOCATE(fieldok(npstn))
     ALLOCATE(index2d(npstn))

     if(is_sequential) then
!     klon_mpi_begin=1
!     klon_mpi_end=klon
      DO ip=1, npstn
       fieldok(ip)=buffer_omp(nptabij(ip))
      ENDDO
     else
      DO ip=1, npstn
!     print*,'histwrite2d is_sequential npstn ip name nptabij',npstn,ip,name,nptabij(ip)
       IF(nptabij(ip).GE.klon_mpi_begin.AND. &
          nptabij(ip).LE.klon_mpi_end) THEN
         fieldok(ip)=buffer_omp(nptabij(ip)-klon_mpi_begin+1)
       ENDIF
      ENDDO
     endif
     CALL histwrite(nid,name,itau,fieldok,npstn,index2d)
!
    endif
    deallocate(index2d)
    deallocate(fieldok)
!$OMP END MASTER    
  end subroutine histwrite2d_phy

  subroutine histwrite3d_phy(nid,lpoint,name,itau,field)
  USE dimphy, only: klon
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat
  use ioipsl, only: histwrite
  implicit none
    
    integer,intent(in) :: nid
    logical,intent(in) :: lpoint
    character*(*), intent(IN) :: name
    integer, intent(in) :: itau
    real,dimension(:,:),intent(in) :: field  ! --> field(klon,:)
    REAL,dimension(klon_mpi,size(field,2)) :: buffer_omp
    REAL :: Field3d(nbp_lon,jj_nb,size(field,2))
    INTEGER :: ip, n, nlev
    INTEGER, ALLOCATABLE, dimension(:) :: index3d
    real,allocatable, dimension(:,:) :: fieldok

    IF (size(field,1)/=klon) CALL abort_physic('iophy::histwrite3d','Field first dimension not equal to klon',1)
    nlev=size(field,2)

!   print*,'hist3d_phy mpi_rank npstn=',mpi_rank,npstn

!   DO ip=1, npstn
!    print*,'hist3d_phy mpi_rank nptabij',mpi_rank,nptabij(ip)
!   ENDDO

    CALL Gather_omp(field,buffer_omp)
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,field3d)
    if(.NOT.lpoint) THEN
     ALLOCATE(index3d(nbp_lon*jj_nb*nlev))
     ALLOCATE(fieldok(nbp_lon*jj_nb,nlev))
     CALL histwrite(nid,name,itau,Field3d,nbp_lon*jj_nb*nlev,index3d)
    else
      nlev=size(field,2)
      ALLOCATE(index3d(npstn*nlev))
      ALLOCATE(fieldok(npstn,nlev))

      if(is_sequential) then
!      klon_mpi_begin=1
!      klon_mpi_end=klon
       DO n=1, nlev
       DO ip=1, npstn
        fieldok(ip,n)=buffer_omp(nptabij(ip),n)
       ENDDO
       ENDDO
      else
       DO n=1, nlev
       DO ip=1, npstn
        IF(nptabij(ip).GE.klon_mpi_begin.AND. &
         nptabij(ip).LE.klon_mpi_end) THEN
         fieldok(ip,n)=buffer_omp(nptabij(ip)-klon_mpi_begin+1,n)
        ENDIF
       ENDDO
       ENDDO
      endif
      CALL histwrite(nid,name,itau,fieldok,npstn*nlev,index3d)
    endif 
  deallocate(index3d)
  deallocate(fieldok)
!$OMP END MASTER    
  end subroutine histwrite3d_phy
  
end module iophy
