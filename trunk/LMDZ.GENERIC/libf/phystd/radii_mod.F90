!==================================================================
module radii_mod
!==================================================================
!  module to centralize the radii calculations for aerosols
! OK for water but should be extended to other aerosols (CO2,...)
!==================================================================
      
!     water cloud optical properties

      use callkeys_mod, only: radfixed,Nmix_co2,                    &
                pres_bottom_tropo,pres_top_tropo,size_tropo,        &
                pres_bottom_strato,size_strato
      
      real, save ::  rad_h2o
      real, save ::  rad_h2o_ice
      real, save ::  Nmix_h2o
      real, save ::  Nmix_h2o_ice
!$OMP THREADPRIVATE(rad_h2o,rad_h2o_ice,Nmix_h2o,Nmix_h2o_ice)
      real, parameter ::  coef_chaud=0.13
      real, parameter ::  coef_froid=0.09


contains


!==================================================================
   subroutine su_aer_radii(ngrid,nlayer,reffrad,nueffrad)
!==================================================================
!     Purpose
!     -------
!     Compute the effective radii of liquid and icy water particles
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================
      use ioipsl_getin_p_mod, only: getin_p
      use radinc_h, only: naerkind
      use aerosol_mod, only: iaero_back2lay, iaero_co2, iaero_dust, &
                             iaero_h2o, iaero_h2so4,iaero_nh3,iaero_aurora
      use callkeys_mod, only: size_nh3_cloud

      Implicit none

      integer,intent(in) :: ngrid
      integer,intent(in) :: nlayer

      real, intent(out) :: reffrad(ngrid,nlayer,naerkind)      !aerosols radii (K)
      real, intent(out) :: nueffrad(ngrid,nlayer,naerkind)     !variance     

      logical, save :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)
      integer :: iaer   
      
      print*,'enter su_aer_radii'
          do iaer=1,naerkind
!     these values will change once the microphysics gets to work
!     UNLESS tracer=.false., in which case we should be working with
!     a fixed aerosol layer, and be able to define reffrad in a 
!     .def file. To be improved!

            if(iaer.eq.iaero_co2)then ! CO2 ice
               reffrad(1:ngrid,1:nlayer,iaer) = 1.e-4
               nueffrad(1:ngrid,1:nlayer,iaer) = 0.1 
            endif

            if(iaer.eq.iaero_h2o)then ! H2O ice
               reffrad(1:ngrid,1:nlayer,iaer) = 1.e-5
               nueffrad(1:ngrid,1:nlayer,iaer) = 0.1 
            endif

            if(iaer.eq.iaero_dust)then ! dust
               reffrad(1:ngrid,1:nlayer,iaer) = 1.e-5
               nueffrad(1:ngrid,1:nlayer,iaer) = 0.1 
            endif
 
            if(iaer.eq.iaero_h2so4)then ! H2O ice
               reffrad(1:ngrid,1:nlayer,iaer) = 1.e-6
               nueffrad(1:ngrid,1:nlayer,iaer) = 0.1 
            endif
            
            if(iaer.eq.iaero_back2lay)then ! Two-layer aerosols
               reffrad(1:ngrid,1:nlayer,iaer) = 2.e-6
               nueffrad(1:ngrid,1:nlayer,iaer) = 0.1 
            endif


	    if(iaer.eq.iaero_nh3)then ! Nh3 cloud
               reffrad(1:ngrid,1:nlayer,iaer) = size_nh3_cloud
               nueffrad(1:ngrid,1:nlayer,iaer) = 0.1 
            endif

	    if(iaer.eq.iaero_aurora)then ! Auroral aerosols
               reffrad(1:ngrid,1:nlayer,iaer) = 3.e-7
               nueffrad(1:ngrid,1:nlayer,iaer) = 0.1 
            endif


            if(iaer.gt.5)then
               print*,'Error in callcorrk, naerkind is too high (>5).'
               print*,'The code still needs generalisation to arbitrary'
               print*,'aerosol kinds and number.'
               call abort
            endif

         enddo


         if (radfixed) then

            write(*,*)"radius of H2O water particles:"
            rad_h2o=13. ! default value
            call getin_p("rad_h2o",rad_h2o)
            write(*,*)" rad_h2o = ",rad_h2o

            write(*,*)"radius of H2O ice particles:"
            rad_h2o_ice=35. ! default value
            call getin_p("rad_h2o_ice",rad_h2o_ice)
            write(*,*)" rad_h2o_ice = ",rad_h2o_ice

         else

            write(*,*)"Number mixing ratio of H2O water particles:"
            Nmix_h2o=1.e6 ! default value
            call getin_p("Nmix_h2o",Nmix_h2o)
            write(*,*)" Nmix_h2o = ",Nmix_h2o

            write(*,*)"Number mixing ratio of H2O ice particles:"
            Nmix_h2o_ice=Nmix_h2o ! default value
            call getin_p("Nmix_h2o_ice",Nmix_h2o_ice)
            write(*,*)" Nmix_h2o_ice = ",Nmix_h2o_ice
         endif

      print*,'exit su_aer_radii'

   end subroutine su_aer_radii
!==================================================================


!==================================================================
   subroutine h2o_reffrad(ngrid,nlayer,pq,pt,reffrad,nueffrad)
!==================================================================
!     Purpose
!     -------
!     Compute the effective radii of liquid and icy water particles
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================
      use watercommon_h, Only: T_h2O_ice_liq,T_h2O_ice_clouds,rhowater,rhowaterice
      use comcstfi_mod, only: pi
      Implicit none

      integer,intent(in) :: ngrid
      integer,intent(in) :: nlayer

      real, intent(in) :: pq(ngrid,nlayer) !water ice mixing ratios (kg/kg)
      real, intent(in) :: pt(ngrid,nlayer) !temperature (K)
      real, intent(out) :: reffrad(ngrid,nlayer)      !aerosol radii
      real, intent(out) :: nueffrad(ngrid,nlayer) ! dispersion      

      integer :: ig,l
      real zfice ,zrad,zrad_liq,zrad_ice
      real,external :: CBRT            
      

      if (radfixed) then
         do l=1,nlayer
            do ig=1,ngrid
               zfice = 1.0 - (pt(ig,l)-T_h2O_ice_clouds) / (T_h2O_ice_liq-T_h2O_ice_clouds)
               zfice = MIN(MAX(zfice,0.0),1.0)
               reffrad(ig,l)= rad_h2o * (1.-zfice) + rad_h2o_ice * zfice
               nueffrad(ig,l) = coef_chaud * (1.-zfice) + coef_froid * zfice
            enddo
         enddo
      else
         do l=1,nlayer
            do ig=1,ngrid
               zfice = 1.0 - (pt(ig,l)-T_h2O_ice_clouds) / (T_h2O_ice_liq-T_h2O_ice_clouds)
               zfice = MIN(MAX(zfice,0.0),1.0)
               zrad_liq  = CBRT( 3*pq(ig,l)/(4*Nmix_h2o*pi*rhowater) )
               zrad_ice  = CBRT( 3*pq(ig,l)/(4*Nmix_h2o_ice*pi*rhowaterice) )
               nueffrad(ig,l) = coef_chaud * (1.-zfice) + coef_froid * zfice
               zrad = zrad_liq * (1.-zfice) + zrad_ice * zfice

               reffrad(ig,l) = min(max(zrad,1.e-6),1000.e-6)
               enddo
            enddo      
      end if

   end subroutine h2o_reffrad
!==================================================================


!==================================================================
   subroutine h2o_cloudrad(ngrid,nlayer,pql,reffliq,reffice)
!==================================================================
!     Purpose
!     -------
!     Compute the effective radii of liquid and icy water particles
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================
      use watercommon_h, Only: rhowater,rhowaterice
      use comcstfi_mod, only: pi
      Implicit none

      integer,intent(in) :: ngrid
      integer,intent(in) :: nlayer

      real, intent(in) :: pql(ngrid,nlayer) !condensed water mixing ratios (kg/kg)
      real, intent(out) :: reffliq(ngrid,nlayer),reffice(ngrid,nlayer)     !liquid and ice water particle radii (m)

      real,external :: CBRT            
      integer :: i,k

      if (radfixed) then
         reffliq(1:ngrid,1:nlayer)= rad_h2o
         reffice(1:ngrid,1:nlayer)= rad_h2o_ice
      else
         do k=1,nlayer
           do i=1,ngrid
             reffliq(i,k) = CBRT(3*pql(i,k)/(4*Nmix_h2o*pi*rhowater))
             reffliq(i,k) = min(max(reffliq(i,k),1.e-6),1000.e-6)
           
             reffice(i,k) = CBRT(3*pql(i,k)/(4*Nmix_h2o_ice*pi*rhowaterice))
             reffice(i,k) = min(max(reffice(i,k),1.e-6),1000.e-6)
           enddo
         enddo
      endif

   end subroutine h2o_cloudrad
!==================================================================



!==================================================================
   subroutine co2_reffrad(ngrid,nlayer,nq,pq,reffrad)
!==================================================================
!     Purpose
!     -------
!     Compute the effective radii of co2 ice particles
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================
      USE tracer_h, only:igcm_co2_ice,rho_co2
      use comcstfi_mod, only: pi
      Implicit none

      integer,intent(in) :: ngrid,nlayer,nq

      real, intent(in) :: pq(ngrid,nlayer,nq) !tracer mixing ratios (kg/kg)
      real, intent(out) :: reffrad(ngrid,nlayer)      !co2 ice particles radii (m)

      integer :: ig,l
      real :: zrad   
      real,external :: CBRT            
            
      

      if (radfixed) then
         reffrad(1:ngrid,1:nlayer) = 5.e-5 ! CO2 ice
      else
         do l=1,nlayer
            do ig=1,ngrid
               zrad = CBRT( 3*pq(ig,l,igcm_co2_ice)/(4*Nmix_co2*pi*rho_co2) )
               reffrad(ig,l) = min(max(zrad,1.e-6),100.e-6)
            enddo
         enddo      
      end if

   end subroutine co2_reffrad
!==================================================================



!==================================================================
   subroutine dust_reffrad(ngrid,nlayer,reffrad)
!==================================================================
!     Purpose
!     -------
!     Compute the effective radii of dust particles
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================
      Implicit none

      integer,intent(in) :: ngrid
      integer,intent(in) :: nlayer

      real, intent(out) :: reffrad(ngrid,nlayer)      !dust particles radii (m)
            
      reffrad(1:ngrid,1:nlayer) = 2.e-6 ! dust

   end subroutine dust_reffrad
!==================================================================


!==================================================================
   subroutine h2so4_reffrad(ngrid,nlayer,reffrad)
!==================================================================
!     Purpose
!     -------
!     Compute the effective radii of h2so4 particles
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================
      Implicit none

      integer,intent(in) :: ngrid
      integer,intent(in) :: nlayer

      real, intent(out) :: reffrad(ngrid,nlayer)      !h2so4 particle radii (m)
                
      reffrad(1:ngrid,1:nlayer) = 1.e-6 ! h2so4

   end subroutine h2so4_reffrad
!==================================================================

!==================================================================
   subroutine back2lay_reffrad(ngrid,reffrad,nlayer,pplev)
!==================================================================
!     Purpose
!     -------
!     Compute the effective radii of particles in a 2-layer model
!
!     Authors
!     -------
!     Sandrine Guerlet (2013)
!
!==================================================================
 
      use aerosol_mod   !! Particle sizes and boundaries of aerosol layers defined there
     Implicit none

      integer,intent(in) :: ngrid

      real, intent(out) :: reffrad(ngrid,nlayer)      ! particle radii (m)
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
      INTEGER,INTENT(IN) :: nlayer ! number of atmospheric layers
      REAL :: expfactor
      INTEGER l,ig
            
      reffrad(:,:)=1e-6  !!initialization, not important
          DO ig=1,ngrid
            DO l=1,nlayer-1
              IF (pplev(ig,l) .le. pres_bottom_tropo .and. pplev(ig,l) .ge. pres_top_tropo) THEN
                reffrad(ig,l) = size_tropo
              ELSEIF (pplev(ig,l) .lt. pres_top_tropo .and. pplev(ig,l) .gt. pres_bottom_strato) THEN
                expfactor=log(size_strato/size_tropo) / log(pres_bottom_strato/pres_top_tropo)
                reffrad(ig,l)= size_tropo*((pplev(ig,l)/pres_top_tropo)**expfactor)
              ELSEIF (pplev(ig,l) .le. pres_bottom_strato) then
                reffrad(ig,l) = size_strato
              ENDIF
            ENDDO
          ENDDO

   end subroutine back2lay_reffrad
!==================================================================


end module radii_mod
!==================================================================
