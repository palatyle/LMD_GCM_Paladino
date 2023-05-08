c***************************************************

      subroutine flujo(date)


c     fgg           nov 2002     first version
c***************************************************

      implicit none


c     common variables and constants

      include 'param.h'
      include 'param_v3.h'


c     Arguments

      real date


c     Local variable and constants

      integer i
      integer inter
      real    nada

c*************************************************

      if(date.lt.1985.) date=1985.
      if(date.gt.2001.) date=2001.
      
      do i=1,ninter
        if(i.le.24) then
          fluxtop(i)=(((ct1(i)+p1(i)*date)/2.)
     $                *sin(2.*3.1416/11.*(date-1985.-3.1416))
     $                +(ct2(i)+p2(i)*date)+1.)*fluxtophr(i)
        else
          fluxtop(i)=fluxtophr(i)
        end if
      end do
     
      return
      end
      


