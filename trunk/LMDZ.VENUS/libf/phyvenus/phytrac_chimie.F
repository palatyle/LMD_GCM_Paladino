      subroutine phytrac_chimie (
     $                    debutphy,
     $                    gmtime,
     $                    nqmax,
     $                    nlon,
     $                    lat,
     $                    lon,
     $                    nlev, 
     $                    pdtphys,
     $                    temp,
     $                    pplay,
     $                    trac,
     $                    d_tr_chem,
     $                    iter)

      use chemparam_mod
      use conc, only: mmean

      implicit none
      
#include "clesphys.h" 
#include "YOMCST.h"

!===================================================================
!     input
!===================================================================

      integer :: nlon, nlev                     ! number of gridpoints and levels
      integer :: nqmax                          ! number of tracers

      real :: gmtime                            ! day fraction
      real :: pdtphys                           ! phytrac_chimie timestep (s)
      real, dimension(nlon,nlev) :: temp        ! temperature (k)
      real, dimension(nlon,nlev) :: pplay       ! pressure (pa)
      real, dimension(nlon,nlev,nqmax) :: trac  ! tracer mass mixing ratio

      logical :: debutphy                       ! first call flag

!===================================================================
!     output
!===================================================================

      real, dimension(nlon,nlev,nqmax) :: d_tr_chem  ! chemical tendency for each tracer
      integer, dimension(nlon,nlev) :: iter          ! chemical iterations

!===================================================================
!     local
!===================================================================

      real :: sza_local   ! solar zenith angle (deg)
      real :: lon_sun

      integer :: i, iq
      integer :: ilon, ilev

      real  lat(nlon), lat_local(nlon)
      real  lon(nlon), lon_local(nlon)

      real, dimension(nlon,nlev) :: mrtwv, mrtsa ! total water and total sulfuric acid
      real, dimension(nlon,nlev) :: mrwv, mrsa   ! gas-phase water and gas-phase sulfuric acid
      real, dimension(nlon,nlev) :: trac_sum
      real, dimension(nlon,nlev,nqmax) :: ztrac  ! local tracer mixing ratio
     
!===================================================================
!     first call : initialisations
!===================================================================

      if (debutphy) then
      
!-------------------------------------------------------------------
!        case of tracers re-initialisation with chemistry
!-------------------------------------------------------------------
         if (reinit_trac .and. ok_chem) then

            print*, "Tracers are re-initialised"
            trac(:,:,:) = 1.0e-30

            if ((i_ocs /= 0) .and. (i_co /= 0) .and. (i_hcl /= 0) 
     $           .and. (i_so2 /= 0) .and. (i_h2o /= 0) .and. (i_n2/ = 0)
     $           .and. (i_co2 /= 0)) then

               trac(:,1:22,i_ocs) = 3.e-6
               trac(:,1:22,i_co)  = 25.e-6
               trac(:,:,i_hcl)    = 0.4e-6
               trac(:,1:22,i_so2) = 10.e-6
               trac(:,1:22,i_h2o) = 30.e-6
               trac(:,:,i_n2)     = 0.35e-1
   
!          ensure that sum of mixing ratios = 1

               trac_sum(:,:) = 0.

               do iq = 1,nqmax - nmicro
                  if (iq /= i_co2) then
                     trac_sum(:,:) = trac_sum(:,:) + trac(:,:,iq) 
                  end if
               end do

!          initialise co2

               trac(:,:,i_co2) = 1. - trac_sum(:,:)

            else
               write(*,*) "at least one tracer is missing : stop"
               stop
            end if
        
!           convert volume to mass mixing ratio

            do iq = 1,nqmax - nmicro
               trac(:,:,iq) = trac(:,:,iq)*m_tr(iq)/mmean(:,:)
            end do
    
         end if  ! reinit_trac

!-------------------------------------------------------------------
!        case of detailed microphysics without chemistry
!-------------------------------------------------------------------
         if (.not. ok_chem .and. ok_cloud .and. cl_scheme == 2) then

!           convert mass to volume mixing ratio

            do iq = 1,nqmax - nmicro
               ztrac(:,:,iq) = trac(:,:,iq)*mmean(:,:)/m_tr(iq)
            end do
          
!           initialise microphysics
 
            call vapors4muphy_ini(nlon,nlev,ztrac)

!           convert volume to mass mixing ratio

            do iq = 1,nqmax - nmicro
               trac(:,:,iq) = ztrac(:,:,iq)*m_tr(iq)/mmean(:,:)
            end do
    
         end if

      end if  ! debutphy 

!===================================================================
!     convert mass to volume mixing ratio : gas phase
!===================================================================

      do iq = 1,nqmax - nmicro
         ztrac(:,:,iq) = max(trac(:,:,iq)*mmean(:,:)/m_tr(iq), 1.e-30)
      end do

!===================================================================
!     microphysics: simplified scheme (phd aurelien stolzenbach)
!===================================================================

      if (ok_cloud .and. cl_scheme == 1) then

!     convert mass to volume mixing ratio : liquid phase

         ztrac(:,:,i_h2so4liq) = max(trac(:,:,i_h2so4liq)
     $                             *mmean(:,:)/m_tr(i_h2so4liq), 1.e-30)
         ztrac(:,:,i_h2oliq) = max(trac(:,:,i_h2oliq)
     $                             *mmean(:,:)/m_tr(i_h2oliq), 1.e-30)
              
!     total water and sulfuric acid (gas + liquid)

         mrtwv(:,:) = ztrac(:,:,i_h2o) + ztrac(:,:,i_h2oliq)
         mrtsa(:,:) = ztrac(:,:,i_h2so4) + ztrac(:,:,i_h2so4liq)

!     all water and sulfuric acid is put in the gas-phase

         mrwv(:,:) = mrtwv(:,:)
         mrsa(:,:) = mrtsa(:,:)

!     call microphysics

         call new_cloud_venus(nlev, nlon, temp, pplay, 
     $                        mrtwv, mrtsa, mrwv, mrsa)

!     update water vapour and sulfuric acid

         ztrac(:,:,i_h2o) = mrwv(:,:)
         ztrac(:,:,i_h2oliq) = mrtwv(:,:) - ztrac(:,:,i_h2o)
     
         ztrac(:,:,i_h2so4) = mrsa(:,:)
         ztrac(:,:,i_h2so4liq) = mrtsa(:,:) - ztrac(:,:,i_h2so4)

      end if  ! simplified scheme

!===================================================================
!     microphysics: detailed scheme (phd sabrina guilbon)
!     !!! to be confirmed whether mad_muphy expects mmr or vmr for h2o and h2so4
!===================================================================

      if (ok_cloud .and. cl_scheme == 2) then

         do iq = nqmax-nmicro+1,nqmax
            ztrac(:,:,iq) = trac(:,:,iq)
         end do

         do ilon = 1,nlon       
            do ilev = 1, nlev
               if (temp(ilon,ilev) < 500.) then
                  call mad_muphy(pdtphys,                               ! timestep
     $                           temp(ilon,ilev),pplay(ilon,ilev),      ! temperature and pressure
     $                           ztrac(ilon,ilev,i_h2o),
     $                           ztrac(ilon,ilev,i_h2so4),      
     $                           ztrac(ilon,ilev,i_m0_aer),
     $                           ztrac(ilon,ilev,i_m3_aer),   
     $                           ztrac(ilon,ilev,i_m0_mode1drop),
     $                           ztrac(ilon,ilev,i_m0_mode1ccn), 
     $                           ztrac(ilon,ilev,i_m3_mode1sa),
     $                           ztrac(ilon,ilev,i_m3_mode1w),    
     $                           ztrac(ilon,ilev,i_m3_mode1ccn),   
     $                           ztrac(ilon,ilev,i_m0_mode2drop),
     $                           ztrac(ilon,ilev,i_m0_mode2ccn),
     $                           ztrac(ilon,ilev,i_m3_mode2sa),
     $                           ztrac(ilon,ilev,i_m3_mode2w), 
     $                           ztrac(ilon,ilev,i_m3_mode2ccn))
               else
                  ztrac(ilon,ilev,i_m0_aer)       = 0. 
                  ztrac(ilon,ilev,i_m3_aer)       = 0.
                  ztrac(ilon,ilev,i_m0_mode1drop) = 0.
                  ztrac(ilon,ilev,i_m0_mode1ccn)  = 0.
                  ztrac(ilon,ilev,i_m3_mode1sa)   = 0.
                  ztrac(ilon,ilev,i_m3_mode1w)    = 0.
                  ztrac(ilon,ilev,i_m3_mode1ccn)  = 0.
                  ztrac(ilon,ilev,i_m0_mode2drop) = 0.
                  ztrac(ilon,ilev,i_m0_mode2ccn)  = 0.
                  ztrac(ilon,ilev,i_m3_mode2sa)   = 0.
                  ztrac(ilon,ilev,i_m3_mode2w)    = 0.
                  ztrac(ilon,ilev,i_m3_mode2ccn)  = 0.
               end if
            end do
         end do

      end if  ! detailed scheme
            
!===================================================================
!     photochemistry
!===================================================================

      if (ok_chem) then

         lon_sun = (0.5 - gmtime)*2.*rpi
         lon_local(:) = lon(:)*rpi/180.
         lat_local(:) = lat(:)*rpi/180.

         do ilon = 1,nlon

!           solar zenith angle

            sza_local = acos(cos(lat_local(ilon))*cos(lon_local(ilon))
     $                 *cos(lon_sun) + cos(lat_local(ilon))
     $                 *sin(lon_local(ilon))*sin(lon_sun))*180./rpi 
      
            call photochemistry_venus(nlev, nlon, pdtphys,
     $                                pplay(ilon,:)/100.,
     $                                temp(ilon,:),
     $                                ztrac(ilon,:,:),
     $                                mmean(ilon,:),
     $                                sza_local, nqmax, iter(ilon,:))

         end do

      end if  ! ok_chem

!===================================================================
!     compute tendencies
!===================================================================

!     gas phase

      do iq = 1,nqmax - nmicro
         ztrac(:,:,iq) = max(ztrac(:,:,iq)*m_tr(iq)/mmean(:,:), 
     $                       1.e-30)
         d_tr_chem(:,:,iq) = (ztrac(:,:,iq) - trac(:,:,iq))/pdtphys
      end do

!     liquid phase or moments

      if (ok_cloud .and. cl_scheme == 1) then
         ztrac(:,:,i_h2so4liq) = max(ztrac(:,:,i_h2so4liq)
     $                               *m_tr(i_h2so4liq)/mmean(:,:),
     $                               1.e-30)
         ztrac(:,:,i_h2oliq)   = max(ztrac(:,:,i_h2oliq)
     $                               *m_tr(i_h2oliq)/mmean(:,:),
     $                               1.e-30)
         d_tr_chem(:,:,i_h2so4liq) = (ztrac(:,:,i_h2so4liq) 
     $                              - trac(:,:,i_h2so4liq))/pdtphys
         d_tr_chem(:,:,i_h2oliq) = (ztrac(:,:,i_h2oliq) 
     $                            - trac(:,:,i_h2oliq))/pdtphys
      end if


      if (ok_cloud .and. cl_scheme == 2) then
         do iq = nqmax-nmicro+1,nqmax
            d_tr_chem(:,:,iq) = (ztrac(:,:,iq) - trac(:,:,iq))/pdtphys
         end do
      end if

      end subroutine phytrac_chimie
