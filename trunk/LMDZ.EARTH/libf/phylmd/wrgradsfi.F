!
! $Header$
!
      subroutine wrgradsfi(if,nl,fieldfi_p,name,titlevar)
      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      implicit none

c   Declarations

#include "dimensions.h"
cym#include "dimphy.h"

c   arguments
      integer if,nl
      real fieldfi_p(klon,nl)
      real fieldfi(klon_glo,nl)
      real fielddyn((iim+1)*(jjm+1),llm)
      character*10 name
      character*10 titlevar

c   local

      integer lm,l,lnblnk



c     print*,'Transformation pour ',name
      call Gather(fieldfi_p,fieldfi)

c$OMP MASTER      
      if (is_mpi_root) then
        call gr_fi_dyn(nl,klon,iim+1,jjm+1,fieldfi,fielddyn)
        
c     print*,'Transformation OK '
        call wrgrads(if,nl,fielddyn,name,titlevar)
c     print*,'Ecriture ok'
      endif
c$OMP END MASTER
      
      return
      end
