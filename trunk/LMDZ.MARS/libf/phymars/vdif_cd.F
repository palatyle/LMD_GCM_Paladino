      SUBROUTINE vdif_cd(ngrid,nlay,pz0,
     &           pg,pz,pu,pv,wstar,pts,ph,pcdv,pcdh)
      USE comcstfi_h
      use turb_mod, only: turb_resolved
      IMPLICIT NONE
c=======================================================================
c
c   Subject: computation of the surface drag coefficient using the
c   -------  approch developed by Loui for ECMWF.
c
c   Author: Frederic Hourdin  15 /10 /93
c   Modified by : Arnaud Colaitis 03/08/11
c   -------
c
c   Arguments:
c   ----------
c
c   inputs:
c   ------
c     ngrid            size of the horizontal grid
c     pg               gravity (m s -2)
c     pz(ngrid,nlay)   height of layers
c     pu(ngrid,nlay)   u component of the wind
c     pv(ngrid,nlay)   v component of the wind
c     pts(ngrid)       surface temperature
c     ph(ngrid)        potential temperature T*(p/ps)^kappa
c
c   outputs:
c   --------
c     pcdv(ngrid)      Cd for the wind
c     pcdh(ngrid)      Cd for potential temperature
c
c=======================================================================
c
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "callkeys.h"

c   Arguments:
c   ----------

      INTEGER, INTENT(IN) :: ngrid,nlay
      REAL, INTENT(IN) :: pz0(ngrid)
      REAL, INTENT(IN) :: pg,pz(ngrid,nlay)
      REAL, INTENT(IN) :: pu(ngrid,nlay),pv(ngrid,nlay)
      REAL, INTENT(IN) :: pts(ngrid),ph(ngrid,nlay)
      REAL, INTENT(IN) :: wstar(ngrid)
      REAL, INTENT(OUT) :: pcdv(ngrid),pcdh(ngrid) ! momentum and heat drag coefficient

c   Local:
c   ------

      INTEGER ig

      REAL karman,nu    ! Von Karman constant and fluid kinematic viscosity
      LOGICAL firstcal
      DATA karman,nu/.41,0.001/
      DATA firstcal/.true./
      SAVE karman,nu

c    Local(2):
c    ---------
      REAL z1,zcd0

      REAL rib(ngrid)  ! Bulk Richardson number
      REAL rig(ngrid)  ! Gradient Richardson number
      REAL fm(ngrid) ! stability function for momentum
      REAL fh(ngrid) ! stability function for heat
      REAL z1z0,z1z0t ! ratios z1/z0 and z1/z0T

c phim = 1+betam*zeta   or   (1-bm*zeta)**am
c phih = alphah + betah*zeta    or   alphah(1.-bh*zeta)**ah
      REAL betam, betah, alphah, bm, bh, lambda
c ah and am are assumed to be -0.25 and -0.5 respectively

      REAL cdn(ngrid),chn(ngrid)  ! neutral momentum and heat drag coefficient
      REAL pz0t        ! initial thermal roughness length. (local)
      REAL ric         ! critical richardson number
      REAL reynolds(ngrid)    ! reynolds number for UBL
      REAL prandtl(ngrid)     ! prandtl number for UBL
      REAL pz0tcomp(ngrid)     ! computed z0t
      REAL ite
      REAL residual
      REAL zu2(ngrid)
c-----------------------------------------------------------------------
c   couche de surface:
c   ------------------

c Original formulation :

      if(.not.callrichsl) then

      DO ig=1,ngrid
         z1=1.E+0 + pz(ig,1)/pz0(ig)
         zcd0=karman/log(z1)
         zcd0=zcd0*zcd0
         pcdv(ig)=zcd0
         pcdh(ig)=zcd0
      ENDDO
     
!      print*,'old : cd,ch; ',pcdv,pcdh
      else

      reynolds(:)=0.

c New formulation (AC) :

c phim = 1+betam*zeta   or   (1-bm*zeta)**am
c phih = alphah + betah*zeta    or   alphah(1.-bh*zeta)**ah
c am=-0.25, ah=-0.5

      pz0t = 0.     ! for the sake of simplicity
      pz0tcomp(:) = 0.1*pz0(:)
      rib(:)=0.

      pcdv(:)=0.
      pcdh(:)=0.

c this formulation assumes alphah=1., implying betah=betam
c We use Dyer et al. parameters, as they cover a broad range of Richardson numbers :
      bm=16.            !UBL
      bh=16.            !UBL
      alphah=1.
      betam=5.         !SBL
      betah=5.         !SBL
      lambda=(sqrt(bh/bm))/alphah
      ric=betah/(betam**2)

      DO ig=1,ngrid

         ite=0.
         residual=abs(pz0tcomp(ig)-pz0t)

         do while((residual .gt. 0.01*pz0(ig)) .and.  (ite .lt. 10.))

         pz0t=pz0tcomp(ig)

         if ((pu(ig,1) .ne. 0.) .or. (pv(ig,1) .ne. 0.)) then

c Classical Richardson number formulation

c         rib(ig) = (pg/ph(ig,1))*pz(ig,1)*(ph(ig,1)-pts(ig))
c     &           /(pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1))

c Richardson number formulation proposed by D.E. England et al. (1995)

!         zu2=MAX(pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1),0.25*wstar(ig)**2)
!         zu2=pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1)
!         zu2(ig)=MAX(pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1),             &
!     &      (0.3*wstar(ig))**2)
          zu2(ig)=pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1) 
     &     + (log(1.+0.7*wstar(ig) + 2.3*wstar(ig)**2))**2 
          if(turb_resolved) then
             zu2(ig)=MAX(zu2(ig),1.)
          endif
!       zu2(ig)=pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1) + (0.5*wstar(ig))**2

               ! we add the wstar to simulate
               ! bulk Ri changes due to subgrid wind feeding the thermals

!          rig(ig) = (pg/ph(ig,1))*((pz(ig,1)-pz0(ig))**2
!     &         /(pz(ig,1)-pz0t))*(ph(ig,1)-pts(ig))
!     &         /zu2

          rib(ig) = (pg/pts(ig))
!     &      *pz(ig,1)*pz0(ig)/sqrt(pz(ig,1)*pz0t)
     &      *sqrt(pz(ig,1)*pz0(ig))
     &      *(((log(pz(ig,1)/pz0(ig)))**2)/(log(pz(ig,1)/pz0t)))
     &      *(ph(ig,1)-pts(ig))
     &  /zu2(ig)

         else
         print*,'warning, infinite Richardson at surface'
         print*,pu(ig,1),pv(ig,1)
         rib(ig) = ric          ! traiter ce cas !
         endif

         z1z0=pz(ig,1)/pz0(ig)
         z1z0t=pz(ig,1)/pz0t

         cdn(ig)=karman/log(z1z0)
         cdn(ig)=cdn(ig)*cdn(ig)
         chn(ig)=cdn(ig)*log(z1z0)/log(z1z0t) 

c Stable case :
      if (rib(ig) .gt. 0.) then
c From D.E. England et al. (95)
      prandtl(ig)=1.
         if(rib(ig) .lt. ric) then
c Assuming alphah=1. and bh=bm for stable conditions :
            fm(ig)=((ric-rib(ig))/ric)**2
            fh(ig)=((ric-rib(ig))/ric)**2
         else
c For Ri>Ric, we consider Ri->Infinity => no turbulent mixing at surface
!            fm(ig)=0.
!            fh(ig)=0.
            fm(ig)=1.
            fh(ig)=1.
         endif
c Unstable case :
      else
c From D.E. England et al. (95)
         fm(ig)=sqrt(1.-lambda*bm*rib(ig))
         fh(ig)=(1./alphah)*((1.-lambda*bh*rib(ig))**0.5)*
     &                     (1.-lambda*bm*rib(ig))**0.25
         prandtl(ig)=alphah*((1.-lambda*bm*rib(ig))**0.25)/
     &             ((1.-lambda*bh*rib(ig))**0.5)
      endif

       reynolds(ig)=karman*sqrt(fm(ig))
     &      *sqrt(zu2(ig))
c     &      *sqrt(pu(ig,1)**2 + pv(ig,1)**2)
     &       *pz0(ig)/(log(z1z0)*nu)
       pz0tcomp(ig)=pz0(ig)*exp(-karman*7.3*
     &              (reynolds(ig)**0.25)*(prandtl(ig)**0.5))

          
      residual = abs(pz0t-pz0tcomp(ig))
      ite = ite+1
!      print*, "iteration nnumber, residual",ite,residual

      enddo  ! of while

          pz0t=0.

c Drag computation :

         pcdv(ig)=cdn(ig)*fm(ig)
         pcdh(ig)=chn(ig)*fh(ig)
       
      ENDDO
!
!      print*,'new : cd,ch; ',pcdv,pcdh

! Some useful diagnostics :

!       call WRITEDIAGFI(ngrid,'RiB',
!     &              'Bulk Richardson nb','no units',
!     &                         2,rib)
!                call WRITEDIAGFI(ngrid,'RiG',
!     &              'Grad Richardson nb','no units',
!     &                         2,rig)
!        call WRITEDIAGFI(ngrid,'Pr',
!     &              'Prandtl nb','no units',
!     &                         0,prandtl)
!       call WRITEDIAGFI(ngrid,'Re',
!     &              'Reynolds nb','no units',
!     &                         0,reynolds)
!        call WRITEDIAGFI(ngrid,'z0tcomp',
!     &              'computed z0t','m', 
!     &                         2,pz0tcomp)


      endif !of if call richardson surface layer

      RETURN
      END
