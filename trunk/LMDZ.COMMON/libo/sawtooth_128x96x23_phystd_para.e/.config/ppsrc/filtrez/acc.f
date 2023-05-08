










!
! $Header$
!
        subroutine acc(vec,d,im)
        implicit none
        integer :: im
        real :: vec(im,im),d(im)
        integer :: i,j
        real ::sum
        real,external :: ssum
        do j=1,im
          do i=1,im
            d(i)=vec(i,j)*vec(i,j)
          enddo
          sum=ssum(im,d,1)
          sum=sqrt(sum)
          do i=1,im
            vec(i,j)=vec(i,j)/sum
          enddo
        enddo
        return
        end
