      subroutine soil(ngrid,nsoil,firstcall,lastcall,
     &          therm_i,
     &          timestep,tsurf,tsoil,
     &          capcal,fluxgrd)

      use comsoil_h, only: layer, mlayer, volcapa, inertiedat
      use comcstfi_mod, only: pi
      use time_phylmdz_mod, only: daysec
      use planete_mod, only: year_day
      use geometry_mod, only: longitude, latitude ! in radians

      implicit none

!-----------------------------------------------------------------------
!  Author: Ehouarn Millour
!
!  Purpose: Compute soil temperature using an implict 1st order scheme
!  
!  Note: depths of layers and mid-layers, soil thermal inertia and 
!        heat capacity are commons in comsoil.h
!-----------------------------------------------------------------------

c-----------------------------------------------------------------------
!  arguments
!  ---------
!  inputs:
      integer,intent(in) :: ngrid	! number of (horizontal) grid-points 
      integer,intent(in) :: nsoil	! number of soil layers 
      logical,intent(in) :: firstcall ! identifier for initialization call 
      logical,intent(in) :: lastcall
      real,intent(in) :: therm_i(ngrid,nsoil) ! thermal inertia
      real,intent(in) :: timestep	    ! time step
      real,intent(in) :: tsurf(ngrid)   ! surface temperature
! outputs:
      real,intent(out) :: tsoil(ngrid,nsoil) ! soil (mid-layer) temperature
      real,intent(out) :: capcal(ngrid) ! surface specific heat
      real,intent(out) :: fluxgrd(ngrid) ! surface diffusive heat flux

! local saved variables:
      real,dimension(:,:),save,allocatable :: mthermdiff ! mid-layer thermal diffusivity
      real,dimension(:,:),save,allocatable :: thermdiff ! inter-layer thermal diffusivity
      real,dimension(:),save,allocatable :: coefq ! q_{k+1/2} coefficients
      real,dimension(:,:),save,allocatable :: coefd ! d_k coefficients
      real,dimension(:,:),save,allocatable :: alph ! alpha_k coefficients
      real,dimension(:,:),save,allocatable :: beta ! beta_k coefficients
      real,save :: mu
!$OMP THREADPRIVATE(mthermdiff,thermdiff,coefq,coefd,alph,beta,mu)
            
! local variables:
      integer ig,ik
      real :: inertia_min,inertia_max
      real :: diurnal_skin ! diurnal skin depth (m)
      real :: annual_skin ! anuual skin depth (m)

! 0. Initialisations and preprocessing step
      if (firstcall) then
      ! note: firstcall is set to .true. or .false. by the caller
      !       and not changed by soil.F 

      ALLOCATE(mthermdiff(ngrid,0:nsoil-1)) ! mid-layer thermal diffusivity
      ALLOCATE(thermdiff(ngrid,nsoil-1))    ! inter-layer thermal diffusivity
      ALLOCATE(coefq(0:nsoil-1))              ! q_{k+1/2} coefficients
      ALLOCATE(coefd(ngrid,nsoil-1))        ! d_k coefficients
      ALLOCATE(alph(ngrid,nsoil-1))         ! alpha_k coefficients
      ALLOCATE(beta(ngrid,nsoil-1))         ! beta_k coefficients

! 0.1 Build mthermdiff(:), the mid-layer thermal diffusivities
      do ig=1,ngrid
        do ik=0,nsoil-1
	  mthermdiff(ig,ik)=therm_i(ig,ik+1)*therm_i(ig,ik+1)/volcapa
!	  write(*,*),'soil: ik: ',ik,' mthermdiff:',mthermdiff(ig,ik)
	enddo
      enddo

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
      !write(*,*)'soil: ig=',ig,' capcal(ig)=',capcal(ig)
      enddo ! of do ig=1,ngrid
      
      ! Additional checks: is the vertical discretization sufficient
      ! to resolve diurnal and annual waves?
      do ig=1,ngrid
        ! extreme inertia for this column
        inertia_min=minval(inertiedat(ig,:))
        inertia_max=maxval(inertiedat(ig,:))
        ! diurnal and annual skin depth
        diurnal_skin=(inertia_min/volcapa)*sqrt(daysec/pi)
        annual_skin=(inertia_max/volcapa)*sqrt(year_day*daysec/pi)
        if (0.5*diurnal_skin<layer(1)) then
        ! one should have the fist layer be at least half of diurnal skin depth
          write(*,*) "soil Error: grid point ig=",ig
          write(*,*) "            longitude=",longitude(ig)*(180./pi)
          write(*,*) "             latitude=",latitude(ig)*(180./pi)
          write(*,*) "  first soil layer depth ",layer(1)
          write(*,*) "  not small enough for a diurnal skin depth of ",
     &                diurnal_skin
          write(*,*) " change soil layer distribution (comsoil_h.F90)"
          stop
        endif
        if (2.*annual_skin>layer(nsoil)) then
        ! one should have the full soil be at least twice the diurnal skin depth
          write(*,*) "soil Error: grid point ig=",ig
          write(*,*) "            longitude=",longitude(ig)*(180./pi)
          write(*,*) "             latitude=",latitude(ig)*(180./pi)
          write(*,*) "  total soil layer depth ",layer(nsoil)
          write(*,*) "  not large enough for an annual skin depth of ",
     &                annual_skin
          write(*,*) " change soil layer distribution (comsoil_h.F90)"
          stop
        endif
      enddo ! of do ig=1,ngrid
      
      else ! of if (firstcall)


!  1. Compute soil temperatures
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
      
      endif! of if (firstcall)

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

!         print*,'this far in soil 1'
!         print*,'thermdiff=',thermdiff(ig,1)
!         print*,'mlayer=',mlayer
!         print*,'beta=',beta(ig,1)
!         print*,'alph=',alph(ig,1)
!         print*,'tsoil=',tsoil(ig,1)

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

