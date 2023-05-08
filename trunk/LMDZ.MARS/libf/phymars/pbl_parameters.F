      SUBROUTINE pbl_parameters(ngrid,nlay,ps,pplay,pz0, 
     & pg,zzlay,zzlev,pu,pv,wstar_in,hfmax,zmax,pts,ph,z_out,n_out,
     & T_out,u_out,ustar,tstar,L_mo,vhf,vvv)
      USE comcstfi_h
      IMPLICIT NONE
!=======================================================================
!
!   Anlysis of the PBL from input temperature, wind field and thermals outputs.
!
!   -------  
!
!   Author: Arnaud Colaitis 09/01/12
!   -------
!
!   Arguments:
!   ----------
!
!   inputs:
!   ------
!     ngrid            size of the horizontal grid
!     nlay             size of the vertical grid
!     pz0(ngrid)       surface roughness length
!     pg               gravity (m s -2)
!     zzlay(ngrid,nlay)   height of mid-layers
!     zzlev(ngrid,nlay+1)   height of layers interfaces
!     pu(ngrid,nlay)   u component of the wind
!     pv(ngrid,nlay)   v component of the wind
!     wstar_in(ngrid)  free convection velocity in PBL
!     hfmax(ngrid)     maximum vertical turbulent heat flux in thermals
!     zmax(ngrid)      height reached by the thermals (pbl height)
!     pts(ngrid)       surface temperature
!     ph(ngrid,nlay)   potential temperature T*(p/ps)^kappa
!     z_out(n_out)     heights of interpolation
!     n_out            number of points for interpolation
!
!   outputs:
!   ------
!
!     Teta_out(ngrid,n_out)  interpolated teta
!     u_out(ngrid,n_out)     interpolated u
!     ustar(ngrid)     friction velocity
!     tstar(ngrid)     friction temperature
!     L_mo(ngrid)      monin_obukhov length
!
!
!=======================================================================
!
!-----------------------------------------------------------------------
!   Declarations:
!   -------------

#include "callkeys.h"

!   Arguments:
!   ----------

      INTEGER, INTENT(IN) :: ngrid,nlay,n_out
      REAL, INTENT(IN) :: pz0(ngrid),ps(ngrid),pplay(ngrid,nlay)
      REAL, INTENT(IN) :: pg,zzlay(ngrid,nlay),zzlev(ngrid,nlay)
      REAL, INTENT(IN) :: pu(ngrid,nlay),pv(ngrid,nlay)
      REAL, INTENT(IN) :: wstar_in(ngrid),hfmax(ngrid),zmax(ngrid)
      REAL, INTENT(IN) :: pts(ngrid),ph(ngrid,nlay)
      REAL, INTENT(IN) :: z_out(n_out)

!    Outputs:
!    --------

      REAL, INTENT(OUT) :: T_out(ngrid,n_out),u_out(ngrid,n_out)
      REAL Teta_out(ngrid,n_out)
      REAL, INTENT(OUT) :: ustar(ngrid), tstar(ngrid)
      REAL, INTENT(OUT) :: L_mo(ngrid)

!   Local:
!   ------

      INTEGER ig,k,n
      REAL karman,nu
      DATA karman,nu/.41,0.001/
      SAVE karman,nu

!    Local(2):
!    ---------
      
      REAL zout
      REAL rib(ngrid)  ! Bulk Richardson number
      REAL fm(ngrid) ! stability function for momentum
      REAL fh(ngrid) ! stability function for heat
      REAL z1z0,z1z0t ! ratios z1/z0 and z1/z0T
          ! phim = 1+betam*zeta   or   (1-bm*zeta)**am
          ! phih = alphah + betah*zeta    or   alphah(1.-bh*zeta)**ah
      REAL betam, betah, alphah, bm, bh, lambda
          ! ah and am are assumed to be -0.25 and -0.5 respectively
      REAL cdn(ngrid),chn(ngrid)  ! neutral momentum and heat drag coefficient
      REAL pz0t        ! initial thermal roughness length. (local)
      REAL ric         ! critical richardson number
      REAL reynolds(ngrid)    ! reynolds number for UBL
      REAL prandtl(ngrid)     ! prandtl number for UBL
      REAL pz0tcomp(ngrid)     ! computed z0t
      REAL ite
      REAL residual,zcd0,z1
      REAL pcdv(ngrid),pcdh(ngrid)
      REAL zu2(ngrid)                  ! Large-scale wind at first layer
      REAL pbl_teta(ngrid)             ! mixed-layer potential temperature
      INTEGER pbl_height_index(ngrid)  ! index of nearest vertical grid point for zmax
      REAL dteta(ngrid,nlay),x(ngrid)  ! potential temperature gradient and z/zi
      REAL dvhf(ngrid),dvvv(ngrid)     ! dimensionless vertical heat flux and 
                                       ! dimensionless vertical velocity variance
      REAL vhf(ngrid),vvv(ngrid)       ! vertical heat flux and vertical velocity variance
      INTEGER ii(1)
! temporary
      INTEGER dimout

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! PART I : RICHARDSON/REYNOLDS/THERMAL_ROUGHNESS/STABILITY_COEFFICIENTS
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      DO n=1,n_out

c Initialisation :

      L_mo(:)=0.
      ustar(:)=0.
      tstar(:)=0.
      zout=z_out(n)
      reynolds(:)=0.
      pz0t = 0.
      pz0tcomp(:) = 0.1*pz0(:)
      rib(:)=0.
      pcdv(:)=0.
      pcdh(:)=0.

      ! this formulation assumes alphah=1., implying betah=betam
      ! We use Dyer et al. parameters, as they cover a broad range of Richardson numbers :

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

       zu2(ig)=pu(ig,1)*pu(ig,1)+pv(ig,1)*pv(ig,1)
     &     + (log(1.+0.7*wstar_in(ig) + 2.3*wstar_in(ig)**2))**2

       DO WHILE((residual .gt. 0.01*pz0(ig)) .and.  (ite .lt. 10.))

         pz0t=pz0tcomp(ig)
         IF (zu2(ig) .ne. 0.) THEN
            ! Richardson number formulation proposed by D.E. England et al. (1995)
          rib(ig) = (pg/pts(ig))
     &        *sqrt(zzlay(ig,1)*pz0(ig))
     &        *(((log(zzlay(ig,1)/pz0(ig)))**2)/(log(zzlay(ig,1)/pz0t)))
     &        *(ph(ig,1)-pts(ig))/(pu(ig,1)*pu(ig,1)+pv(ig,1)*pv(ig,1))
         ELSE
            print*,'warning, infinite Richardson at surface'
            print*,pu(ig,1),pv(ig,1)
            rib(ig) = ric
         ENDIF

         z1z0=zzlay(ig,1)/pz0(ig)
         z1z0t=zzlay(ig,1)/pz0t

         cdn(ig)=karman/log(z1z0)
         cdn(ig)=cdn(ig)*cdn(ig)
         chn(ig)=cdn(ig)*log(z1z0)/log(z1z0t) 

         ! STABLE BOUNDARY LAYER :
         IF (rib(ig) .gt. 0.) THEN
            ! From D.E. England et al. (95)
            prandtl(ig)=1.
            if(rib(ig) .lt. ric) then
               ! Assuming alphah=1. and bh=bm for stable conditions :
               fm(ig)=((ric-rib(ig))/ric)**2
               fh(ig)=((ric-rib(ig))/ric)**2
            else
               ! For Ri>Ric, we consider Ri->Infinity => no turbulent mixing at surface
               fm(ig)=0.
               fh(ig)=0.
            endif
         ! UNSTABLE BOUNDARY LAYER :
         ELSE
            ! From D.E. England et al. (95)
            fm(ig)=sqrt(1.-lambda*bm*rib(ig))
            fh(ig)=(1./alphah)*((1.-lambda*bh*rib(ig))**0.5)*
     &                     (1.-lambda*bm*rib(ig))**0.25
            prandtl(ig)=alphah*((1.-lambda*bm*rib(ig))**0.25)/
     &             ((1.-lambda*bh*rib(ig))**0.5)
         ENDIF
 
        reynolds(ig)=karman*sqrt(fm(ig))
     &              *sqrt(zu2(ig))
     &              *pz0(ig)/(log(z1z0)*nu)
        pz0tcomp(ig)=pz0(ig)*exp(-karman*7.3*
     &              (reynolds(ig)**0.25)*(prandtl(ig)**0.5))
        residual = abs(pz0t-pz0tcomp(ig))
        ite = ite+1

       ENDDO  ! of while
       pz0t=0.

! Drag computation:

         pcdv(ig)=cdn(ig)*fm(ig)
         pcdh(ig)=chn(ig)*fh(ig)
       
      ENDDO ! of ngrid

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! PART II : USTAR/TSTAR/U_OUT/TETA_OUT COMPUTATIONS
!------------------------------------------------------------------------
!------------------------------------------------------------------------

! u* theta* computation

      DO ig=1,ngrid
         IF (rib(ig) .ge. ric) THEN
           ustar(ig)=0.
           tstar(ig)=0.
         ELSE
           ustar(ig)=sqrt(pcdv(ig))
     &        *sqrt(pu(ig,1)*pu(ig,1)+pv(ig,1)*pv(ig,1))
           tstar(ig)=-pcdh(ig)*(pts(ig)-ph(ig,1))
     &        /sqrt(pcdv(ig))
         ENDIF
      ENDDO

! Interpolation:

      DO ig=1,ngrid
        IF(zout .lt. pz0tcomp(ig)) THEN
           u_out(ig,n)=0.
           Teta_out(ig,n)=pts(ig)

        ELSE
          IF (rib(ig) .ge. ric) THEN ! ustar=tstar=0  (and fm=fh=0)
           u_out(ig,n)=0
           Teta_out(ig,n)=pts(ig)
          ELSE
           u_out(ig,n)= ustar(ig)*log(zout/pz0(ig))/
     &(karman*sqrt(fm(ig)))

           Teta_out(ig,n)=pts(ig)+(tstar(ig)*sqrt(fm(ig))*log(zout/
     & (pz0tcomp(ig)))/
     &(karman*fh(ig)))
          ENDIF
        ENDIF

        IF (zout .lt. pz0(ig)) THEN
           u_out(ig,n)=0.
        ENDIF 

      ENDDO

! when using convective adjustment without thermals, a vertical potential temperature
! profile is assumed up to the thermal roughness length. Hence, theoretically, theta
! interpolated at any height in the surface layer is theta at the first level.

      IF ((.not.calltherm).and.(calladj)) THEN
       Teta_out(:,n)=ph(:,1)
       u_out(:,n)=(sqrt(cdn(:))*sqrt(pu(:,1)*pu(:,1)+pv(:,1)*pv(:,1))
     &                                /karman)*log(zout/pz0(:))
      ENDIF
              T_out(:,n) = Teta_out(:,n)*(exp(
     &   (zout/zzlay(:,1))*(log(pplay(:,1)/ps))
     &                             )
     &                         )**rcp

      ENDDO   !of n=1,n_out


!------------------------------------------------------------------------
!------------------------------------------------------------------------
! PART III : WSTAR COMPUTATION
!------------------------------------------------------------------------
!------------------------------------------------------------------------

! Detection of the mixed-layer potential temperature
! ------------

! Nearest index for the pbl height

      IF (calltherm) THEN

      pbl_height_index(:)=1

      DO k=1,nlay-1
         DO ig=1, ngrid
            IF (abs(zmax(ig)-zzlay(ig,k)) .lt. 
     &              abs(zmax(ig)-zzlay(ig,pbl_height_index(ig)))) THEN
               pbl_height_index(ig)=k
            ENDIF
         ENDDO
      ENDDO

! Potential temperature gradient

      dteta(:,nlay)=0.
      DO k=1,nlay-1
         DO ig=1, ngrid
         dteta(ig,k) = (ph(ig,k+1)-ph(ig,k))/(zzlay(ig,k+1)-zzlay(ig,k))
         ENDDO
      ENDDO

! Computation of the pbl mixed layer temperature

      DO ig=1, ngrid
         ii=MINLOC(abs(dteta(ig,1:pbl_height_index(ig))))
         pbl_teta(ig) = ph(ig,ii(1))
      ENDDO


!------------------------------------------------------------------------
!------------------------------------------------------------------------
! PART IV : VERTICAL_VELOCITY_VARIANCE/VERTICAL_TURBULENT_FLUX PROFILES
!------------------------------------------------------------------------
!------------------------------------------------------------------------

! We follow Spiga et. al 2010 (QJRMS)
! ------------

      DO ig=1, ngrid
         IF (zmax(ig) .gt. 0.) THEN
            x(ig) = zout/zmax(ig)
         ELSE
            x(ig) = 999.
         ENDIF
      ENDDO

      DO ig=1, ngrid
         ! dimensionless vertical heat flux
         IF (x(ig) .le. 0.3) THEN
            dvhf(ig) = ((-3.85/log(x(ig)))+0.07*log(x(ig)))
     &                                       *exp(-4.61*x(ig))
         ELSEIF (x(ig) .le. 1.) THEN
            dvhf(ig) = -1.52*x(ig) + 1.24
         ELSE
            dvhf(ig) = 0.
         ENDIF
         ! dimensionless vertical velocity variance
         IF (x(ig) .le. 1.) THEN
            dvvv(ig) = 2.05*(x(ig)**(2./3.))*(1.-0.64*x(ig))**2
         ELSE
            dvvv(ig) = 0.
         ENDIF
      ENDDO

      vhf(:) = dvhf(:)*hfmax(:)
      vvv(:) = dvvv(:)*(wstar_in(:))**2

      ENDIF ! of if calltherm

#ifndef MESOSCALE
            call WRITEDIAGFI(ngrid,'rib_pbl',
     &   'Richardson in pbl parameter','m/s',2,rib)
            call WRITEDIAGFI(ngrid,'cdn_pbl',
     &   'neutral momentum coef','m/s',2,cdn)
            call WRITEDIAGFI(ngrid,'fm_pbl',
     &   'momentum stab function','m/s',2,fm)
            call WRITEDIAGFI(ngrid,'uv',
     &   'wind norm first layer','m/s',2,sqrt(zu2))
            call WRITEDIAGFI(ngrid,'uvtrue',
     &   'wind norm first layer','m/s',2,sqrt(pu(:,1)**2+pv(:,1)**2))
            call WRITEDIAGFI(ngrid,'chn_pbl',
     &   'neutral momentum coef','m/s',2,chn)
            call WRITEDIAGFI(ngrid,'fh_pbl',
     &   'momentum stab function','m/s',2,fh)
            call WRITEDIAGFI(ngrid,'B_pbl',
     &   'buoyancy','m/',2,(ph(:,1)-pts(:))/pts(:))
            call WRITEDIAGFI(ngrid,'zot_pbl',
     &   'buoyancy','ms',2,pz0tcomp)
       call WRITEDIAGFI(ngrid,'zz1','buoyancy','m',2,zzlay(:,1))
#endif

      RETURN
      END
