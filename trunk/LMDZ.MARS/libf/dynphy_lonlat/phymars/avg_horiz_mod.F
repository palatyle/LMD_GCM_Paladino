      MODULE avg_horiz_mod

      IMPLICIT NONE

      CONTAINS
      
      subroutine avg_horiz(imo,jmo,imn,jmn,rlonuo,rlatvo,rlonv,
     .             rlatu,zdata,outdata)

c=======================================================================
c Maximum/minimum value composite approach for interpolate the variables
c  from surface.nc (summit, base... ) 
c  Chao Wang
c=======================================================================

      implicit none
      !inputs/outputs
      integer, intent(in)  ::  imo,jmo
      integer, intent(in)  :: imn,jmn
      real, intent(in)  :: rlonv(imn+1)    
      real, intent(in)  :: rlatu(jmn+1)   
      real, intent(in)  :: zdata(imo*jmo)
      real, intent(out)  :: outdata(imn+1,jmn+1)

      !local variables:
      real pi
      integer i,j,l,m,n
      real tmplon(imn+2)
      real tmplat(jmn+2)
      real rlonuo(imo)    
      real rlatvo(jmo)   

      integer iik(imo,jmo),jjk(imo,jmo)
      integer tmp
      integer nmesh
      
c     initialize 
      iik(:,:)=0
      jjk(:,:)=0
      outdata(:,:)=0.
      pi=2.*ASIN(1.)

      tmplon(1)=rlonv(1)
      do i=2,imn+1
        tmplon(i)= (rlonv(i-1)+rlonv(i))/2.
      enddo
      tmplon(imn+2)=rlonv(imn+1)

      tmplat(1)=rlatu(1)
      do j=2,jmn+1
        tmplat(j)=(rlatu(j-1)+rlatu(j))/2.
      enddo
      tmplat(jmn+2)=rlatu(jmn+1)


      do i=1,imo 
        do j=1,jmo
          do m=1,imn+1
            do n=1,jmn+1
              
              if ( ((rlonuo(i) .ge. tmplon(m)) .and.
     &                       (rlonuo(i) .lt. tmplon(m+1))) .and.
     &           ((rlatvo(j) .lt. tmplat(n) ) .and.
     &            (rlatvo(j) .ge. tmplat(n+1) ))) then
                iik(i,j)=m
                jjk(i,j)=n
                goto 1000
              endif
              
            enddo 
          enddo
1000      continue
        enddo
      enddo

      do m=1,imn+1
        do n=1,jmn+1
          nmesh=0
          do i=1,imo 
            do j=1,jmo
               if ((iik(i,j) .eq. m) .and. (jjk(i,j) .eq. n)) then
                  outdata(m,n)=outdata(m,n)+zdata(i+imo*(j-1)) 
                  nmesh=nmesh+1
               endif       

            enddo 
          enddo
          outdata(m,n)=outdata(m,n)/nmesh
        enddo
      enddo
      
      END SUBROUTINE avg_horiz

      END MODULE avg_horiz_mod
