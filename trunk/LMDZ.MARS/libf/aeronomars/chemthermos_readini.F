      subroutine chemthermos_readini

      use param_v4_h, only: rcoef
      implicit none

c     local variables
      real*8        rc1,rc2,rc3
      integer       ierr,i

c*************************+PROGRAM STARTS**************************

!Opens external file to read reaction coefficients
      open(unit=10,file='chemthermos_reactionrates.def',status='old',
     $     iostat=ierr)
      if(ierr.eq.0) then
         read(10,*)
         do i=1,61
            read(10,*)
            read(10,*)
            read(10,*) rc1,rc2,rc3
            rcoef(i,1)=rc1
            rcoef(i,2)=rc2
            rcoef(i,3)=rc3
         enddo
         close(10)
      else
         write(*,*)
         write(*,*) 'Cannot find file chemthermos_reactionrates.def'
         stop
      endif
      
      end
