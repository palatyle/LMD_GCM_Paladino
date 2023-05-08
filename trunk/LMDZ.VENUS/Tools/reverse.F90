subroutine reverse4d(nlon,nlat,np,nt,newf)
!==============================================================================
! Purpose:
! reverse lat, lon and vertical axis on a 3D+time field
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: nlon ! longitude size
integer,intent(in) :: nlat ! latitude size
integer,intent(in) :: np   ! vertical size
integer,intent(in) :: nt   ! time size
real,intent(inout),dimension(nlon,nlat,np,nt) :: newf ! field to be reversed

!==============================================================================
! Local variables:
!==============================================================================
integer :: i,j,k,l
real,dimension(nlon,nlat,np,nt) :: tmp,tmp2 ! reversed fields
include "planet.h"

! Vertical axis:
    do k=1,np
          tmp(:,:,k,:)=newf(:,:,np+1-k,:)
    enddo
! horizontal dimensions
if(planet.eq."Venus") then
  do l=1,nt
      do i=1,nlat
        do j=1,nlon
          tmp2(j,i,:,l)=tmp(nlon+1-j,nlat+1-i,:,l)
        enddo
      enddo
  enddo
else
  tmp2=tmp
endif

newf=tmp2

end subroutine reverse4d

!===============================================================================

subroutine reverse2d(nlon,nlat,newf)
!==============================================================================
! Purpose:
! reverse lat and lon on a 2D field
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: nlon ! longitude size
integer,intent(in) :: nlat ! latitude size
real,intent(inout),dimension(nlon,nlat) :: newf ! field to be reversed

!==============================================================================
! Local variables:
!==============================================================================
integer :: i,j
real,dimension(nlon,nlat) :: tmp ! reversed field

      do i=1,nlat
        do j=1,nlon
          tmp(j,i)=newf(nlon+1-j,nlat+1-i)
        enddo
      enddo
newf=tmp

end subroutine reverse2d

!===============================================================================

subroutine reverse3d(nlon,nlat,nt,newf)
!==============================================================================
! Purpose:
! reverse lat and lon on a 2D+time field
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: nlon ! longitude size
integer,intent(in) :: nlat ! latitude size
integer,intent(in) :: nt ! time size
real,intent(inout),dimension(nlon,nlat,nt) :: newf ! field to be reversed

!==============================================================================
! Local variables:
!==============================================================================
integer :: i,j,l
real,dimension(nlon,nlat,nt) :: tmp ! reversed field

  do l=1,nt
      do i=1,nlat
        do j=1,nlon
          tmp(j,i,l)=newf(nlon+1-j,nlat+1-i,l)
        enddo
      enddo
  enddo
newf=tmp

end subroutine reverse3d

!===============================================================================

subroutine reverselev(np,newf)
!==============================================================================
! Purpose:
! reverse vertical pressure axis 
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: np   ! vertical size
real,intent(inout),dimension(np) :: newf ! field to be reversed

!==============================================================================
! Local variables:
!==============================================================================
integer :: k
real,dimension(np) :: tmp ! reversed field

    do k=1,np
          tmp(k)=newf(np+1-k)
    enddo
newf=tmp

end subroutine reverselev
