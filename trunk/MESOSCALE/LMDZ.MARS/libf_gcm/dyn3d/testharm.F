      PROGRAM testharm
      implicit none

#include "dimensions.h"
#include "paramet.h"

      integer i,j,ii,jj
      real pi
      real psi(jjp1,jjp1,iim) , fphi(iim,iim)
      character file*20, str2*2, str2b*2

#include "serre.h"
#include "fxyprim.h"

      pi=2.*asin(1.)

      call harmspher( fphi,psi )
      print*,'FIN de harmspher'

      do i=1,iim
         write(str2,'(i2.2)') i
         file='Phi'//str2
         open(10,file=file,form='formatted')
         print*,'i,file',i,file
         do ii=1,iim
            write(10,*) fx(float(ii))*180./pi,fphi(ii,i)
         enddo
         close(10)
      enddo

      do i=1,iim
         do j=1,jjp1
            write(str2,'(i2.2)') i
            write(str2b,'(i2.2)') j
            file='Psi.'//str2//'.'//str2b
            open(10,file=file,form='formatted')
            do jj=1,jjp1
               write(10,*) fy(float(jj))*180./pi,psi(jj,j,i)
            enddo
            close(10)
         enddo
      enddo

      end
