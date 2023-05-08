subroutine getslopes(ngrid,geopot)
    
use geometry_mod, only: longitude, latitude ! in radians
use slope_mod, only: theta_sl, psi_sl
use comcstfi_h, only: g, rad, pi
use mod_phys_lmdz_para, only: is_parallel
use mod_grid_phy_lmdz, only: nbp_lon, nbp_lat
implicit none


! This routine computes slope inclination and orientation for the GCM (callslope=.true. in callphys.def)
! It works fine with a non-regular grid for zoomed simulations.
! slope inclination angle (deg)  0 == horizontal, 90 == vertical
! slope orientation angle (deg)  0 == Northward,  90 == Eastward, 180 == Southward, 270 == Westward
! TN 04/1013

integer,intent(in) :: ngrid ! nnumber of atmospheric columns
real,intent(in) :: geopot(ngrid)     ! geopotential on phy grid
real topogrid(nbp_lon,nbp_lat) ! topography on lat/lon grid with poles and only one -180/180 point
real latigrid(nbp_lon,nbp_lat),longgrid(nbp_lon,nbp_lat) ! meshgrid of latitude and longitude values (radians)
real theta_val ! slope inclination
real psi_val   ! slope orientation
real gradx(nbp_lon,nbp_lat) ! x: latitude-wise topography gradient,  increasing northward
real grady(nbp_lon,nbp_lat) ! y: longitude-wise topography gradient, increasing westward
integer i,j,ig0
integer id2,idm1 ! a trick to compile testphys1d with debug option

if (is_parallel) then
  ! This routine only works in serial mode so stop now.
  write(*,*) "getslopes Error: this routine is not designed to run in parallel"
  stop
endif

id2  = 2
idm1 = nbp_lon-1

! rearrange topography on a 2d array
do j=2,nbp_lat-1
   ig0= 1+(j-2)*nbp_lon
   do i=1,nbp_lon
      topogrid(i,j)=geopot(ig0+i)/g
      latigrid(i,j)=latitude(ig0+i)
      longgrid(i,j)=longitude(ig0+i)
   enddo
enddo
!poles :
topogrid(:,1) = geopot(1)/g
latigrid(:,1) = latitude(1)
longgrid(:,1) = longitude(1)
topogrid(:,nbp_lat) = geopot(ngrid)/g
latigrid(:,nbp_lat) = latitude(ngrid)
longgrid(:,nbp_lat) = longitude(ngrid)



! compute topography gradient
! topogrid and rad are both in meters
do j=2,nbp_lat-1
   do i=1,nbp_lon
     gradx(i,j) = (topogrid(i,j+1) - topogrid(i,j-1)) / (latigrid(i,j+1)-latigrid(i,j-1))
     gradx(i,j) = gradx(i,j) / rad
   enddo
   grady(1,j) = (topogrid(id2,j) - topogrid(nbp_lon,j)) / (2*pi+longgrid(id2,j)-longgrid(nbp_lon,j)) 
   grady(1,j) = grady(1,j) / rad
   grady(nbp_lon,j) = (topogrid(1,j) - topogrid(idm1,j)) / (2*pi+longgrid(1,j)-longgrid(idm1,j)) 
   grady(nbp_lon,j) = grady(nbp_lon,j) / rad
   do i=2,nbp_lon-1
     grady(i,j) = (topogrid(i+1,j) - topogrid(i-1,j)) / (longgrid(i+1,j)-longgrid(i-1,j)) 
     grady(i,j) = grady(i,j) / rad
   enddo
enddo
! poles :
gradx(:,1) = 0.
grady(:,1) = 0.
gradx(:,nbp_lat) = 0.
grady(:,nbp_lat) = 0.



! compute slope inclination and orientation :
theta_sl(:) = 0.
psi_sl(:)   = 0.
do j=2,nbp_lat-1
   do i=1,nbp_lon
   
     ig0= 1+(j-2)*nbp_lon
   
     theta_val=atan(sqrt( (gradx(i,j))**2 + (grady(i,j))**2 ))
     
     psi_val=0.
     if (gradx(i,j) .ne. 0.) psi_val= -pi/2. - atan(grady(i,j)/gradx(i,j))
     if (gradx(i,j) .ge. 0.) psi_val= psi_val - pi
     psi_val = 3*pi/2. - psi_val
     psi_val = psi_val*180./pi
     psi_val = MODULO(psi_val,360.)
     
     theta_sl(ig0+i) = theta_val
     psi_sl(ig0+i)   = psi_val
          
   enddo
enddo



end subroutine getslopes
