      SUBROUTINE dustlift(ngrid,nlay,nq,rho,pcdh_true,pcdh,co2ice,
     $                  dqslift)
      IMPLICIT NONE

c=======================================================================
c
c  Dust lifting by surface winds
c  Computing flux to the middle of the first layer
c  (Called by vdifc)
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"

c
c   arguments:
c   ----------

c   INPUT
      integer ngrid, nlay, nq  
      real rho(ngrid)  ! density (kg.m-3) at surface
      real pcdh_true(ngrid) ! Cd 
      real pcdh(ngrid) ! Cd * |V| 
      real co2ice(ngrid)

c   OUTPUT
      real dqslift(ngrid,nq) !surface dust flux to mid-layer (<0 when lifing)
c     real pb(ngrid,nlay) ! diffusion to surface coeff.

c   local:
c   ------
      INTEGER ig,iq
      REAL fhoriz(ngridmx)  ! Horizontal dust flux 
      REAL ust,us
      REAL stress_seuil
      SAVE stress_seuil
c      DATA stress_seuil/0.0225/   ! stress seuil soulevement (N.m2)
!****WRF
!****WRF: additional ASCII file to define dust opacity
          REAL alpha  
          INTEGER ierr
          OPEN(99,file='stress.def',status='old',form='formatted'
     .     ,iostat=ierr)
          IF(ierr.NE.0) THEN
             stress_seuil = 0.0225
             alpha = 1.   
             write(*,*) 'No file stress.def - set ', stress_seuil, alpha
             !stop
          ELSE
             READ(99,*) stress_seuil
             READ(99,*) alpha
             write(*,*) 'definir seuil stress : ', stress_seuil, alpha
             CLOSE(99)
          ENDIF
          alpha_lift(1) = alpha          
!****WRF
!****WRF





c     ---------------------------------
c     Computing horizontal flux: fhoriz
c     ---------------------------------

      do ig=1,ngrid
          fhoriz(ig) = 0.      ! initialisation

c         Selection of points where surface dust is available
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         if (latid(ig).ge.80.) goto 99  ! N permanent  polar caps
c         if (latid(ig).le.-80.) goto 99 ! S polar deposits
c         if  ((longd(ig).ge.-141. .and. longd(ig).le.-127.)
c    &   .and.(latid(ig).ge.12.   .and. latid(ig).le.23.))goto 99 ! olympus
c         if  ((longd(ig).ge.-125. .and. longd(ig).le.-118.)
c    &   .and.(latid(ig).ge.-12.   .and. latid(ig).le.-6.))goto 99 ! Arsia
c         if  ((longd(ig).ge.-116. .and. longd(ig).le.-109.)
c    &   .and.(latid(ig).ge.-5.   .and. latid(ig).le. 5.))goto 99 ! pavonis
c         if  ((longd(ig).ge.-109. .and. longd(ig).le.-100.)
c    &   .and.(latid(ig).ge. 7.   .and. latid(ig).le. 16.))goto 99 ! ascraeus
c         if  ((longd(ig).ge.  61. .and. longd(ig).le.  63.)
c    &   .and.(latid(ig).ge. 63. .and. latid(ig).le. 64.))goto 99 !weird point
          if (co2ice(ig).gt.0.) goto 99


c         Is the wind strong enough ?
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ust = sqrt(stress_seuil/rho(ig))
          us = pcdh(ig) /  sqrt(pcdh_true(ig)) ! ustar=cd*v /sqrt(cd)

          if (us.gt.ust) then
c            If lifting ? 
c            Calcul du flux suivant Marticorena ( en fait white (1979))

             fhoriz(ig) = 2.61*(rho(ig)/g) *
     &      (us -ust) * (us + ust)**2
          end if
 99      continue
      end do


cc       write(*,*) 'caca', us, ust
cc       do iq=1,nq
cc       do ig=1,ngrid
cc          write(*,*) 'cacao', -alpha_lift(iq) * fhoriz(ig)
cc       enddo
cc       enddo



c     -------------------------------------
c     Computing vertical flux and diffusion
c     -------------------------------------
 
       do iq=1,nq
         do ig=1,ngrid
             dqslift(ig,iq)= -alpha_lift(iq)* fhoriz(ig)

cc  le  flux vertical remplace le terme de diffusion turb. qui est mis a zero
c            zb(ig,1) = 0.
cc           If surface deposition by turbulence diffusion (impaction...)
cc           if(fhoriz(ig).ne.0) then
cc           zb(ig,1) = zcdh(ig)*zb0(ig,1)
cc           AMount of Surface deposition !
cc           pdqs_dif(ig,iq)=pdqs_dif(ig,iq) +
cc    &      zb(ig,1)*zq(ig,1,iq)/ptimestep
cc          write(*,*) 'zb(1)  = ' ,  zb(ig,1),zcdh(ig),zb0(ig,1)
cc

         enddo
       enddo

      RETURN
      END

