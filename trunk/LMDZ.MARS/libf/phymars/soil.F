      subroutine soil(ngrid,nsoil,firstcall,
     &          therm_i,
     &          timestep,tsurf,tsoil,
     &          capcal,fluxgrd)

      use comsoil_h, only: layer, mlayer, volcapa,
     &                     mthermdiff, thermdiff, coefq,
     &                     coefd, alph, beta, mu
      use surfdat_h, only: watercaptag, inert_h2o_ice

      implicit none

!-----------------------------------------------------------------------
!  Author: Ehouarn Millour
!
!  Purpose: Compute soil temperature using an implict 1st order scheme
!  
!  Note: depths of layers and mid-layers, soil thermal inertia and 
!        heat capacity are commons in comsoil_h
!-----------------------------------------------------------------------

#include "callkeys.h"

c-----------------------------------------------------------------------
!  arguments
!  ---------
!  inputs:
      integer ngrid	! number of (horizontal) grid-points 
      integer nsoil	! number of soil layers 
      logical firstcall ! identifier for initialization call 
      real therm_i(ngrid,nsoil) ! thermal inertia
      real timestep	    ! time step
      real tsurf(ngrid)   ! surface temperature
! outputs:
      real tsoil(ngrid,nsoil) ! soil (mid-layer) temperature
      real capcal(ngrid) ! surface specific heat
      real fluxgrd(ngrid) ! surface diffusive heat flux

! local variables:
      integer ig,ik

! 0. Initialisations and preprocessing step
      if (firstcall.or.tifeedback) then
      ! note: firstcall is set to .true. or .false. by the caller
      !       and not changed by soil.F 
! 0.1 Build mthermdiff(:), the mid-layer thermal diffusivities
      do ig=1,ngrid
        if (watercaptag(ig)) then
          do ik=0,nsoil-1
! If we have permanent ice, we use the water ice thermal inertia from ground to last layer.
              mthermdiff(ig,ik)=inert_h2o_ice*inert_h2o_ice/volcapa
          enddo
        else
          do ik=0,nsoil-1
           mthermdiff(ig,ik)=therm_i(ig,ik+1)*therm_i(ig,ik+1)/volcapa
          enddo
        endif
      enddo

#ifdef MESOSCALE
      do ig=1,ngrid
        if ( therm_i(ig,1) .ge. inert_h2o_ice ) then
          print *, "limit max TI ", therm_i(ig,1), inert_h2o_ice
          do ik=0,nsoil-1
               mthermdiff(ig,ik)=inert_h2o_ice*inert_h2o_ice/volcapa
          enddo
        endif
      enddo
#endif

! 0.2 Build thermdiff(:), the "interlayer" thermal diffusivities
      do ig=1,ngrid
        do ik=1,nsoil-1
      thermdiff(ig,ik)=((layer(ik)-mlayer(ik-1))*mthermdiff(ig,ik)
     &                +(mlayer(ik)-layer(ik))*mthermdiff(ig,ik-1))
     &                    /(mlayer(ik)-mlayer(ik-1))
!	write(*,*),'soil: ik: ',ik,' thermdiff:',thermdiff(ig,ik)
	enddo
      enddo

! 0.3 Build coefficients mu, q_{k+1/2}, d_k, alpha_k and capcal
      ! mu
      mu=mlayer(0)/(mlayer(1)-mlayer(0))

      ! q_{1/2}
      coefq(0)=volcapa*layer(1)/timestep
	! q_{k+1/2}
        do ik=1,nsoil-1
          coefq(ik)=volcapa*(layer(ik+1)-layer(ik))
     &                 /timestep
	enddo

      do ig=1,ngrid
	! d_k
	do ik=1,nsoil-1
	  coefd(ig,ik)=thermdiff(ig,ik)/(mlayer(ik)-mlayer(ik-1))
	enddo
	
	! alph_{N-1}
	alph(ig,nsoil-1)=coefd(ig,nsoil-1)/
     &                  (coefq(nsoil-1)+coefd(ig,nsoil-1))
        ! alph_k
        do ik=nsoil-2,1,-1
	  alph(ig,ik)=coefd(ig,ik)/(coefq(ik)+coefd(ig,ik+1)*
     &                              (1.-alph(ig,ik+1))+coefd(ig,ik))
	enddo

        ! capcal
! Cstar
        capcal(ig)=volcapa*layer(1)+
     &              (thermdiff(ig,1)/(mlayer(1)-mlayer(0)))*
     &              (timestep*(1.-alph(ig,1)))
! Cs
        capcal(ig)=capcal(ig)/(1.+mu*(1.0-alph(ig,1))*
     &                         thermdiff(ig,1)/mthermdiff(ig,0))
!      write(*,*)'soil: ig=',ig,' capcal(ig)=',capcal(ig)
      enddo ! of do ig=1,ngrid
            
      endif ! of if (firstcall.or.tifeedback)

!  1. Compute soil temperatures
      IF (.not.firstcall) THEN
! First layer:
      do ig=1,ngrid
        tsoil(ig,1)=(tsurf(ig)+mu*beta(ig,1)*
     &                         thermdiff(ig,1)/mthermdiff(ig,0))/
     &              (1.+mu*(1.0-alph(ig,1))*
     &               thermdiff(ig,1)/mthermdiff(ig,0))
      enddo
! Other layers:
      do ik=1,nsoil-1
        do ig=1,ngrid
	  tsoil(ig,ik+1)=alph(ig,ik)*tsoil(ig,ik)+beta(ig,ik)
	enddo
      enddo
      
      ENDIF! of if (.not.firstcall)

!  2. Compute beta coefficients (preprocessing for next time step)
! Bottom layer, beta_{N-1}
      do ig=1,ngrid
        beta(ig,nsoil-1)=coefq(nsoil-1)*tsoil(ig,nsoil)
     &                   /(coefq(nsoil-1)+coefd(ig,nsoil-1))
      enddo
! Other layers
      do ik=nsoil-2,1,-1
        do ig=1,ngrid
	  beta(ig,ik)=(coefq(ik)*tsoil(ig,ik+1)+
     &                 coefd(ig,ik+1)*beta(ig,ik+1))/
     &                 (coefq(ik)+coefd(ig,ik+1)*(1.0-alph(ig,ik+1))
     &                  +coefd(ig,ik))
	enddo
      enddo

!  3. Compute surface diffusive flux & calorific capacity
      do ig=1,ngrid
! Cstar
!        capcal(ig)=volcapa(ig,1)*layer(ig,1)+
!     &              (thermdiff(ig,1)/(mlayer(ig,1)-mlayer(ig,0)))*
!     &              (timestep*(1.-alph(ig,1)))
! Fstar
        fluxgrd(ig)=(thermdiff(ig,1)/(mlayer(1)-mlayer(0)))*
     &              (beta(ig,1)+(alph(ig,1)-1.0)*tsoil(ig,1))

!        mu=mlayer(ig,0)/(mlayer(ig,1)-mlayer(ig,0))
!        capcal(ig)=capcal(ig)/(1.+mu*(1.0-alph(ig,1))*
!     &                         thermdiff(ig,1)/mthermdiff(ig,0))
! Fs
        fluxgrd(ig)=fluxgrd(ig)+(capcal(ig)/timestep)*
     &              (tsoil(ig,1)*(1.+mu*(1.0-alph(ig,1))*
     &                         thermdiff(ig,1)/mthermdiff(ig,0))
     &               -tsurf(ig)-mu*beta(ig,1)*
     &                          thermdiff(ig,1)/mthermdiff(ig,0))
      enddo

      end

