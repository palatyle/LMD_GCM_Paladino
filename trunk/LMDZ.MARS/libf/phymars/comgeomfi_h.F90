
       module comgeomfi_h

       implicit none

       REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: sinlon 
       REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: coslon
       REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: sinlat
       REAL,ALLOCATABLE,SAVE,DIMENSION(:) :: coslat

       contains

         subroutine ini_comgeomfi_h(ngrid)

         implicit none
         integer,intent(in) :: ngrid ! number of atmospheric columns

         allocate(sinlat(ngrid))
         allocate(coslat(ngrid))
         allocate(sinlon(ngrid))
         allocate(coslon(ngrid))

         end subroutine ini_comgeomfi_h


         subroutine end_comgeomfi_h

         implicit none

         if (allocated(sinlat)) deallocate(sinlat)
         if (allocated(coslat)) deallocate(coslat)
         if (allocated(sinlon)) deallocate(sinlon)
         if (allocated(coslon)) deallocate(coslon)

         end subroutine end_comgeomfi_h

         subroutine ini_fillgeom(ngrid,plat,plon,parea)

         implicit none
         INTEGER,INTENT(IN) :: ngrid ! number of atmospheric columns
         REAL,INTENT(IN) :: plat(ngrid),plon(ngrid),parea(ngrid)
         integer :: ig

         DO ig=1,ngrid
            sinlat(ig)=sin(plat(ig))
            coslat(ig)=cos(plat(ig))
            sinlon(ig)=sin(plon(ig))
            coslon(ig)=cos(plon(ig))
         ENDDO

         end subroutine ini_fillgeom

       end module comgeomfi_h

