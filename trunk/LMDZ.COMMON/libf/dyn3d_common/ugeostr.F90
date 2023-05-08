!
! $Id: ugeostr.F90 1474 2011-01-14 11:04:45Z lguez $
!
subroutine ugeostr(phi,ucov)

  ! Calcul du vent covariant geostrophique a partir du champ de
  ! geopotentiel.
  ! We actually compute: (1 - cos^8 \phi) u_g
  ! to have a wind going smoothly to 0 at the equator.
  ! We assume that the surface pressure is uniform so that model
  ! levels are pressure levels.

  USE comconst_mod, ONLY: omeg,rad

  implicit none

  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"

  real ucov(iip1,jjp1,llm),phi(iip1,jjp1,llm)
  real um(jjm,llm),fact,u(iip1,jjm,llm)
  integer i,j,l

  real zlat

  um(:,:)=0 ! initialize um()

  DO j=1,jjm

     if (abs(sin(rlatv(j))).lt.1.e-4) then
        zlat=1.e-4
     else
        zlat=rlatv(j)
     endif
     fact=cos(zlat)
     fact=fact*fact
     fact=fact*fact
     fact=fact*fact
     fact=(1.-fact)/ &
          (2.*omeg*sin(zlat)*(rlatu(j+1)-rlatu(j)))
     fact=-fact/rad
     DO l=1,llm
        DO i=1,iim
           u(i,j,l)=fact*(phi(i,j+1,l)-phi(i,j,l))
           um(j,l)=um(j,l)+u(i,j,l)/REAL(iim)
        ENDDO
     ENDDO
  ENDDO

  !   calcul des champ de vent:

  DO l=1,llm
     DO i=1,iip1
        ucov(i,1,l)=0.
        ucov(i,jjp1,l)=0.
     end DO
     DO  j=2,jjm
        DO  i=1,iim
           ucov(i,j,l) = 0.5*(u(i,j,l)+u(i,j-1,l))*cu(i,j)
        end DO
        ucov(iip1,j,l)=ucov(1,j,l)
     end DO
  end DO

  print *, 301

end subroutine ugeostr
