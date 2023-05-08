c**********************************************************************

      subroutine phdisrate(zenit,i)

c     apr 2002       fgg           first version
c**********************************************************************


      implicit none

      include 'param.h'
      include 'param_v3.h'

c     arguments

      integer         i                !altitude
      real            zenit

c     local variables

      integer         inter,iz
      real            lambda
      real            jdis(nabs,ninter,nzmax)
      character*1     dn


c**********************************************************************

      if(zenit.gt.100.) then
         dn='n'
         else
         dn='d'
      end if

      if(dn.eq.'n') then
         jdistot(1,i) = 0.
         jdistot(2,i) = 0.
         jdistot(3,i) = 0.
         jdistot(4,i) = 0.
         jdistot(5,i) = 0.
         jdistot(6,i) = 0.
         jdistot_b(1,i) = 0.
         jdistot_b(2,i) = 0.
      else

c     photodissociation rate

         jdistot(1,i)=0.
         jdistot(2,i)=0.
         jdistot(4,i)=0.
         jdistot(6,i)=0.
         jdistot_b(1,i)=0.
         jdistot_b(2,i)=0.
         do inter=1,32
            jdis(1,inter,i) = jfotsout(inter,1,i) * fluxtop(inter) 
     $      * efdisco2(inter)
            if(inter.gt.29.and.inter.le.32) then
               jdistot(1,i) = jdistot(1,i) + jdis(1,inter,i)
            else if(inter.le.29) then
               jdistot_b(1,i) = jdistot_b(1,i) + jdis(1,inter,i)
            end if
         end do

         do inter=1,31
            jdis(2,inter,i) = jfotsout(inter,2,i) * fluxtop(inter)
     $      * efdiso2(inter)
            if(inter.ge.31) then
               jdistot(2,i) = jdistot(2,i) + jdis(2,inter,i)
            else if(inter.lt.31) then
               jdistot_b(2,i) = jdistot_b(2,i) + jdis(2,inter,i)
            end if
         end do

         jdistot(3,i) = 0.

         do inter=25,31
            jdis(4,inter,i) = jfotsout(inter,4,i) * fluxtop(inter)
     $      * efdish2o(inter)
            jdistot(4,i) = jdistot(4,i) + jdis(4,inter,i)
         end do

         do inter=1,15
            jdis(5,inter,i) = jfotsout(inter,5,i) * fluxtop(inter)
     $      * efdish2(inter)
            jdistot(5,i) = jdistot(5,i) + jdis(5,inter,i)
         end do

         do inter=25,33
            jdis(6,inter,i) = jfotsout(inter,6,i) * fluxtop(inter)
     $      * efdish2o2(inter) 
            jdistot(6,i) = jdistot(6,i) + jdis(6,inter,i)
         end do 
      end if


      return
      

      end










