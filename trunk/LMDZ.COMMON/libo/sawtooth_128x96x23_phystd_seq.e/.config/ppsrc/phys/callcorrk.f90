










MODULE callcorrk_mod

IMPLICIT NONE

CONTAINS

      subroutine callcorrk(ngrid,nlayer,pq,nq,qsurf,           &
          albedo,albedo_equivalent,emis,mu0,pplev,pplay,pt,    & 
          tsurf,fract,dist_star,aerosol,muvar,                 &
          dtlw,dtsw,fluxsurf_lw,                               &
          fluxsurf_sw,fluxsurfabs_sw,fluxtop_lw,               &
          fluxabs_sw,fluxtop_dn,                               &
          OLR_nu,OSR_nu,                                       &
          int_dtaui,int_dtauv,                                 &
          tau_col,cloudfrac,totcloudfrac,                      &
          clearsky,firstcall,lastcall)

      use mod_phys_lmdz_para, only : is_master
      use radinc_h, only: L_NSPECTV, L_NSPECTI, naerkind, banddir, corrkdir,&
                          L_LEVELS, L_NGAUSS, L_NLEVRAD, L_NLAYRAD, L_REFVAR
      use radcommon_h, only: wrefvar, Cmk, fzeroi, fzerov, gasi, gasv, &
                             glat_ig, gweight, pfgasref, pgasmax, pgasmin, &
                             pgasref, tgasmax, tgasmin, tgasref, scalep, &
                             ubari, wnoi, stellarf, glat, dwnv, dwni, tauray
      use watercommon_h, only: psat_water, epsi
      use datafile_mod, only: datadir
      use ioipsl_getin_p_mod, only: getin_p
      use gases_h, only: ngasmx
      use radii_mod, only : su_aer_radii,co2_reffrad,h2o_reffrad,dust_reffrad,h2so4_reffrad,back2lay_reffrad
      use aerosol_mod, only : iaero_co2,iaero_h2o,iaero_dust,iaero_h2so4, iaero_back2lay, iaero_nh3, iaero_aurora
      use tracer_h, only: igcm_h2o_vap, igcm_h2o_ice, igcm_co2_ice
      use comcstfi_mod, only: pi, mugaz, cpp
      use callkeys_mod, only: varactive,diurnal,tracer,water,varfixed,satval,diagdtau,    &
                              kastprof,strictboundcorrk,specOLR,CLFvarying,               &
                              tplanckmin,tplanckmax
      use optcv_mod, only: optcv
      use optci_mod, only: optci
      implicit none

!==================================================================
!
!     Purpose
!     -------
!     Solve the radiative transfer using the correlated-k method for
!     the gaseous absorption and the Toon et al. (1989) method for
!     scatttering due to aerosols.
!
!     Authors
!     ------- 
!     Emmanuel 01/2001, Forget 09/2001
!     Robin Wordsworth (2009)
!
!==================================================================

!-----------------------------------------------------------------------
!     Declaration of the arguments (INPUT - OUTPUT) on the LMD GCM grid
!     Layer #1 is the layer near the ground. 
!     Layer #nlayer is the layer at the top.
!-----------------------------------------------------------------------


      ! INPUT
      INTEGER,INTENT(IN) :: ngrid                  ! Number of atmospheric columns.
      INTEGER,INTENT(IN) :: nlayer                 ! Number of atmospheric layers.
      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq)       ! Tracers (kg/kg_of_air).
      INTEGER,INTENT(IN) :: nq                     ! Number of tracers.
      REAL,INTENT(IN) :: qsurf(ngrid,nq)           ! Tracers on surface (kg.m-2).
      REAL,INTENT(IN) :: albedo(ngrid,L_NSPECTV)   ! Spectral Short Wavelengths Albedo. By MT2015
      REAL,INTENT(IN) :: emis(ngrid)               ! Long Wave emissivity.
      REAL,INTENT(IN) :: mu0(ngrid)                ! Cosine of sun incident angle.
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1)     ! Inter-layer pressure (Pa).
      REAL,INTENT(IN) :: pplay(ngrid,nlayer)       ! Mid-layer pressure (Pa).
      REAL,INTENT(IN) :: pt(ngrid,nlayer)          ! Air temperature (K).
      REAL,INTENT(IN) :: tsurf(ngrid)              ! Surface temperature (K).
      REAL,INTENT(IN) :: fract(ngrid)              ! Fraction of day.
      REAL,INTENT(IN) :: dist_star                 ! Distance star-planet (AU).
      REAL,INTENT(IN) :: muvar(ngrid,nlayer+1)
      REAL,INTENT(IN) :: cloudfrac(ngrid,nlayer)   ! Fraction of clouds (%).
      logical,intent(in) :: clearsky
      logical,intent(in) :: firstcall              ! Signals first call to physics.
      logical,intent(in) :: lastcall               ! Signals last call to physics.
      
      ! OUTPUT
      REAL,INTENT(OUT) :: aerosol(ngrid,nlayer,naerkind) ! Aerosol tau.
      REAL,INTENT(OUT) :: dtlw(ngrid,nlayer)             ! Heating rate (K/s) due to LW radiation.
      REAL,INTENT(OUT) :: dtsw(ngrid,nlayer)             ! Heating rate (K/s) due to SW radiation.
      REAL,INTENT(OUT) :: fluxsurf_lw(ngrid)             ! Incident LW flux to surf (W/m2).
      REAL,INTENT(OUT) :: fluxsurf_sw(ngrid)             ! Incident SW flux to surf (W/m2)
      REAL,INTENT(OUT) :: fluxsurfabs_sw(ngrid)          ! Absorbed SW flux by the surface (W/m2). By MT2015.
      REAL,INTENT(OUT) :: fluxtop_lw(ngrid)              ! Outgoing LW flux to space (W/m2).
      REAL,INTENT(OUT) :: fluxabs_sw(ngrid)              ! SW flux absorbed by the planet (W/m2).
      REAL,INTENT(OUT) :: fluxtop_dn(ngrid)              ! Incident top of atmosphere SW flux (W/m2).
      REAL,INTENT(OUT) :: OLR_nu(ngrid,L_NSPECTI)        ! Outgoing LW radition in each band (Normalized to the band width (W/m2/cm-1).
      REAL,INTENT(OUT) :: OSR_nu(ngrid,L_NSPECTV)        ! Outgoing SW radition in each band (Normalized to the band width (W/m2/cm-1).
      REAL,INTENT(OUT) :: tau_col(ngrid)                 ! Diagnostic from aeropacity.
      REAL,INTENT(OUT) :: albedo_equivalent(ngrid)       ! Spectrally Integrated Albedo. For Diagnostic. By MT2015
      REAL,INTENT(OUT) :: totcloudfrac(ngrid)            ! Column Fraction of clouds (%).
      REAL,INTENT(OUT) :: int_dtaui(ngrid,nlayer,L_NSPECTI) ! VI optical thickness of layers within narrowbands for diags ().
      REAL,INTENT(OUT) :: int_dtauv(ngrid,nlayer,L_NSPECTV) ! IR optical thickness of layers within narrowbands for diags ().
      
      
      
      

      ! Globally varying aerosol optical properties on GCM grid ; not needed everywhere so not in radcommon_h.   
      REAL :: QVISsQREF3d(ngrid,nlayer,L_NSPECTV,naerkind)
      REAL :: omegaVIS3d(ngrid,nlayer,L_NSPECTV,naerkind)
      REAL :: gVIS3d(ngrid,nlayer,L_NSPECTV,naerkind)
      REAL :: QIRsQREF3d(ngrid,nlayer,L_NSPECTI,naerkind)
      REAL :: omegaIR3d(ngrid,nlayer,L_NSPECTI,naerkind)
      REAL :: gIR3d(ngrid,nlayer,L_NSPECTI,naerkind)

!      REAL :: omegaREFvis3d(ngrid,nlayer,naerkind)
!      REAL :: omegaREFir3d(ngrid,nlayer,naerkind) ! not sure of the point of these...

      REAL,ALLOCATABLE,SAVE :: reffrad(:,:,:)  ! aerosol effective radius (m)
      REAL,ALLOCATABLE,SAVE :: nueffrad(:,:,:) ! aerosol effective variance
!$OMP THREADPRIVATE(reffrad,nueffrad)

!-----------------------------------------------------------------------
!     Declaration of the variables required by correlated-k subroutines
!     Numbered from top to bottom (unlike in the GCM)
!-----------------------------------------------------------------------

      REAL*8 tmid(L_LEVELS),pmid(L_LEVELS)
      REAL*8 tlevrad(L_LEVELS),plevrad(L_LEVELS)

      ! Optical values for the optci/cv subroutines
      REAL*8 stel(L_NSPECTV),stel_fract(L_NSPECTV)
      ! NB: Arrays below are "save" to avoid reallocating them at every call
      ! not because their content needs be reused from call to the next
      REAL*8,allocatable,save :: dtaui(:,:,:)
      REAL*8,allocatable,save :: dtauv(:,:,:)
      REAL*8,allocatable,save :: cosbv(:,:,:)
      REAL*8,allocatable,save :: cosbi(:,:,:)
      REAL*8,allocatable,save :: wbari(:,:,:)
      REAL*8,allocatable,save :: wbarv(:,:,:)
!$OMP THREADPRIVATE(dtaui,dtauv,cosbv,cosbi,wbari,wbarv)
      REAL*8,allocatable,save :: tauv(:,:,:)
      REAL*8,allocatable,save :: taucumv(:,:,:)
      REAL*8,allocatable,save :: taucumi(:,:,:)
!$OMP THREADPRIVATE(tauv,taucumv,taucumi)
      REAL*8 tauaero(L_LEVELS,naerkind)
      REAL*8 nfluxtopv,nfluxtopi,nfluxtop,fluxtopvdn
      REAL*8 nfluxoutv_nu(L_NSPECTV)                 ! Outgoing band-resolved VI flux at TOA (W/m2).
      REAL*8 nfluxtopi_nu(L_NSPECTI)                 ! Net band-resolved IR flux at TOA (W/m2).
      REAL*8 fluxupi_nu(L_NLAYRAD,L_NSPECTI)         ! For 1D diagnostic.
      REAL*8 fmneti(L_NLAYRAD),fmnetv(L_NLAYRAD)
      REAL*8 fluxupv(L_NLAYRAD),fluxupi(L_NLAYRAD)
      REAL*8 fluxdnv(L_NLAYRAD),fluxdni(L_NLAYRAD)
      REAL*8 albi,acosz
      REAL*8 albv(L_NSPECTV)                         ! Spectral Visible Albedo.

      INTEGER ig,l,k,nw,iaer

      real,save :: szangle
      logical,save :: global1d
!$OMP THREADPRIVATE(szangle,global1d)

      real*8,allocatable,save :: taugsurf(:,:)
      real*8,allocatable,save :: taugsurfi(:,:)
!$OMP THREADPRIVATE(taugsurf,taugsurfi)
      real*8 qvar(L_LEVELS)   ! Mixing ratio of variable component (mol/mol).

      ! Local aerosol optical properties for each column on RADIATIVE grid.
      real*8,save,allocatable ::  QXVAER(:,:,:)
      real*8,save,allocatable ::  QSVAER(:,:,:)
      real*8,save,allocatable ::  GVAER(:,:,:)
      real*8,save,allocatable ::  QXIAER(:,:,:)
      real*8,save,allocatable ::  QSIAER(:,:,:)
      real*8,save,allocatable ::  GIAER(:,:,:)
!$OMP THREADPRIVATE(QXVAER,QSVAER,GVAER,QXIAER,QSIAER,GIAER)
      real, dimension(:,:,:), save, allocatable :: QREFvis3d
      real, dimension(:,:,:), save, allocatable :: QREFir3d
!$OMP THREADPRIVATE(QREFvis3d,QREFir3d)


      ! Miscellaneous :
      real*8  temp,temp1,temp2,pweight
      character(len=10) :: tmp1
      character(len=10) :: tmp2
      character(len=100) :: message
      character(len=10),parameter :: subname="callcorrk"

      ! For fixed water vapour profiles.
      integer i_var
      real RH
      real*8 pq_temp(nlayer)
! real(KIND=r8) :: pq_temp(nlayer) ! better F90 way.. DOESNT PORT TO F77!!!
      real psat,qsat

      logical OLRz
      real*8 NFLUXGNDV_nu(L_NSPECTV)

      ! Included by RW for runaway greenhouse 1D study.
      real vtmp(nlayer)
      REAL*8 muvarrad(L_LEVELS)
      
      ! Included by MT for albedo calculations.      
      REAL*8 albedo_temp(L_NSPECTV) ! For equivalent albedo calculation.
      REAL*8 surface_stellar_flux   ! Stellar flux reaching the surface. Useful for equivalent albedo calculation.


!===============================================================
!           I.a Initialization on first call
!===============================================================


      if(firstcall) then

        ! test on allocated necessary because of CLFvarying (two calls to callcorrk in physiq)
        if(.not.allocated(QXVAER)) allocate(QXVAER(L_LEVELS,L_NSPECTV,naerkind))
        if(.not.allocated(QSVAER)) allocate(QSVAER(L_LEVELS,L_NSPECTV,naerkind))
        if(.not.allocated(GVAER)) allocate(GVAER(L_LEVELS,L_NSPECTV,naerkind))
        if(.not.allocated(QXIAER)) allocate(QXIAER(L_LEVELS,L_NSPECTI,naerkind))
        if(.not.allocated(QSIAER)) allocate(QSIAER(L_LEVELS,L_NSPECTI,naerkind))
        if(.not.allocated(GIAER)) allocate(GIAER(L_LEVELS,L_NSPECTI,naerkind))

         !!! ALLOCATED instances are necessary because of CLFvarying (strategy to call callcorrk twice in physiq...)
         IF(.not.ALLOCATED(QREFvis3d)) ALLOCATE(QREFvis3d(ngrid,nlayer,naerkind))
         IF(.not.ALLOCATED(QREFir3d)) ALLOCATE(QREFir3d(ngrid,nlayer,naerkind))
         ! Effective radius and variance of the aerosols
         IF(.not.ALLOCATED(reffrad)) allocate(reffrad(ngrid,nlayer,naerkind))
         IF(.not.ALLOCATED(nueffrad)) allocate(nueffrad(ngrid,nlayer,naerkind))

         call system('rm -f surf_vals_long.out')

         if(naerkind.gt.4)then
            message='Code not general enough to deal with naerkind > 4 yet.'
            call abort_physic(subname,message,1)
         endif
         call su_aer_radii(ngrid,nlayer,reffrad,nueffrad)
         
         
!--------------------------------------------------
!             Set up correlated k
!--------------------------------------------------


         print*, "callcorrk: Correlated-k data base folder:",trim(datadir)
         call getin_p("corrkdir",corrkdir)
         print*, "corrkdir = ",corrkdir
         write( tmp1, '(i3)' ) L_NSPECTI
         write( tmp2, '(i3)' ) L_NSPECTV
         banddir=trim(adjustl(tmp1))//'x'//trim(adjustl(tmp2))
         banddir=trim(adjustl(corrkdir))//'/'//trim(adjustl(banddir))

         call setspi            ! Basic infrared properties.
         call setspv            ! Basic visible properties.
         call sugas_corrk       ! Set up gaseous absorption properties.
         call suaer_corrk       ! Set up aerosol optical properties.
        

         ! now that L_NGAUSS has been initialized (by sugas_corrk)
         ! allocate related arrays
         if(.not.allocated(dtaui)) ALLOCATE(dtaui(L_NLAYRAD,L_NSPECTI,L_NGAUSS))
         if(.not.allocated(dtauv)) ALLOCATE(dtauv(L_NLAYRAD,L_NSPECTV,L_NGAUSS))
         if(.not.allocated(cosbv)) ALLOCATE(cosbv(L_NLAYRAD,L_NSPECTV,L_NGAUSS))
         if(.not.allocated(cosbi)) ALLOCATE(cosbi(L_NLAYRAD,L_NSPECTI,L_NGAUSS))
         if(.not.allocated(wbari)) ALLOCATE(wbari(L_NLAYRAD,L_NSPECTI,L_NGAUSS))
         if(.not.allocated(wbarv)) ALLOCATE(wbarv(L_NLAYRAD,L_NSPECTV,L_NGAUSS))
         if(.not.allocated(tauv)) ALLOCATE(tauv(L_NLEVRAD,L_NSPECTV,L_NGAUSS))
         if(.not.allocated(taucumv)) ALLOCATE(taucumv(L_LEVELS,L_NSPECTV,L_NGAUSS))
         if(.not.allocated(taucumi)) ALLOCATE(taucumi(L_LEVELS,L_NSPECTI,L_NGAUSS))
         if(.not.allocated(taugsurf)) ALLOCATE(taugsurf(L_NSPECTV,L_NGAUSS-1))
         if(.not.allocated(taugsurfi)) ALLOCATE(taugsurfi(L_NSPECTI,L_NGAUSS-1))

         if((igcm_h2o_vap.eq.0) .and. varactive)then
            message='varactive in callcorrk but no h2o_vap tracer.'
            call abort_physic(subname,message,1)
         endif

         OLR_nu(:,:) = 0.
         OSR_nu(:,:) = 0.

         if (ngrid.eq.1) then
            PRINT*, 'Simulate global averaged conditions ?'
            global1d = .false. ! default value
            call getin_p("global1d",global1d)
            write(*,*) "global1d = ",global1d
            
            ! Test of incompatibility : if global1d is true, there should not be any diurnal cycle.
            if (global1d.and.diurnal) then
               message='if global1d is true, diurnal must be set to false'
               call abort_physic(subname,message,1)
            endif

            if (global1d) then
               PRINT *,'Solar Zenith angle (deg.) ?'
               PRINT *,'(assumed for averaged solar flux S/4)'
               szangle=60.0  ! default value
               call getin_p("szangle",szangle)
               write(*,*) "szangle = ",szangle
            endif
         endif

      end if ! of if (firstcall)

!=======================================================================
!          I.b  Initialization on every call   
!=======================================================================
 
      qxvaer(:,:,:)=0.0
      qsvaer(:,:,:)=0.0
      gvaer(:,:,:) =0.0

      qxiaer(:,:,:)=0.0
      qsiaer(:,:,:)=0.0
      giaer(:,:,:) =0.0

!--------------------------------------------------
!     Effective radius and variance of the aerosols
!--------------------------------------------------

      do iaer=1,naerkind

         if ((iaer.eq.iaero_co2).and.tracer.and.(igcm_co2_ice.gt.0)) then ! Treat condensed co2 particles.
            call co2_reffrad(ngrid,nlayer,nq,pq,reffrad(1,1,iaero_co2))
            if (is_master) then
	       print*,'Max. CO2 ice particle size = ',maxval(reffrad(1:ngrid,1:nlayer,iaer))/1.e-6,' um'
               print*,'Min. CO2 ice particle size = ',minval(reffrad(1:ngrid,1:nlayer,iaer))/1.e-6,' um'
            end if
	 end if
         
         if ((iaer.eq.iaero_h2o).and.water) then ! Treat condensed water particles. To be generalized for other aerosols ...
            call h2o_reffrad(ngrid,nlayer,pq(1,1,igcm_h2o_ice),pt, &
                             reffrad(1,1,iaero_h2o),nueffrad(1,1,iaero_h2o))
            if (is_master) then
               print*,'Max. H2O cloud particle size = ',maxval(reffrad(1:ngrid,1:nlayer,iaer))/1.e-6,' um'
               print*,'Min. H2O cloud particle size = ',minval(reffrad(1:ngrid,1:nlayer,iaer))/1.e-6,' um'
            end if
         endif
         
         if(iaer.eq.iaero_dust)then
            call dust_reffrad(ngrid,nlayer,reffrad(1,1,iaero_dust))
            if (is_master) then
               print*,'Dust particle size = ',reffrad(1,1,iaer)/1.e-6,' um'
            end if
         endif
         
         if(iaer.eq.iaero_h2so4)then
            call h2so4_reffrad(ngrid,nlayer,reffrad(1,1,iaero_h2so4))
            if (is_master) then
               print*,'H2SO4 particle size =',reffrad(1,1,iaer)/1.e-6,' um'
            end if
         endif
         
          if(iaer.eq.iaero_back2lay)then
            call back2lay_reffrad(ngrid,reffrad(1,1,iaero_back2lay),nlayer,pplev)
         endif
!         if(iaer.eq.iaero_nh3)then
!	    call nh3_reffrad(ngrid,nlayer,reffrad(1,1,iaero_nh3))
!         endif
!         if(iaer.eq.iaero_aurora)then
!	    call aurora_reffrad(ngrid,nlayer,reffrad(1,1,iaero_aurora))
!         endif
        
     end do !iaer=1,naerkind.


      ! How much light do we get ?
      do nw=1,L_NSPECTV
         stel(nw)=stellarf(nw)/(dist_star**2)
      end do

      ! Get 3D aerosol optical properties.
      call aeroptproperties(ngrid,nlayer,reffrad,nueffrad,         &
           QVISsQREF3d,omegaVIS3d,gVIS3d,                          &
           QIRsQREF3d,omegaIR3d,gIR3d,                             &
           QREFvis3d,QREFir3d)                                     

      ! Get aerosol optical depths.
      call aeropacity(ngrid,nlayer,nq,pplay,pplev,pq,aerosol,      &
           reffrad,QREFvis3d,QREFir3d,                             & 
           tau_col,cloudfrac,totcloudfrac,clearsky)                
          


!-----------------------------------------------------------------------    
      do ig=1,ngrid ! Starting Big Loop over every GCM column
!-----------------------------------------------------------------------


!=======================================================================
!              II.  Transformation of the GCM variables
!=======================================================================


!-----------------------------------------------------------------------
!    Aerosol optical properties Qext, Qscat and g.
!    The transformation in the vertical is the same as for temperature.
!-----------------------------------------------------------------------
           
           
            do iaer=1,naerkind
               ! Shortwave.
               do nw=1,L_NSPECTV 
               
                  do l=1,nlayer

                     temp1=QVISsQREF3d(ig,nlayer+1-l,nw,iaer)         &
                         *QREFvis3d(ig,nlayer+1-l,iaer)

                     temp2=QVISsQREF3d(ig,max(nlayer-l,1),nw,iaer)    &
                         *QREFvis3d(ig,max(nlayer-l,1),iaer)

                     qxvaer(2*l,nw,iaer)  = temp1
                     qxvaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                     temp1=temp1*omegavis3d(ig,nlayer+1-l,nw,iaer)
                     temp2=temp2*omegavis3d(ig,max(nlayer-l,1),nw,iaer)

                     qsvaer(2*l,nw,iaer)  = temp1
                     qsvaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                     temp1=gvis3d(ig,nlayer+1-l,nw,iaer)
                     temp2=gvis3d(ig,max(nlayer-l,1),nw,iaer)

                     gvaer(2*l,nw,iaer)  = temp1
                     gvaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                  end do ! nlayer

                  qxvaer(1,nw,iaer)=qxvaer(2,nw,iaer)
                  qxvaer(2*nlayer+1,nw,iaer)=0.

                  qsvaer(1,nw,iaer)=qsvaer(2,nw,iaer)
                  qsvaer(2*nlayer+1,nw,iaer)=0.

                  gvaer(1,nw,iaer)=gvaer(2,nw,iaer)
                  gvaer(2*nlayer+1,nw,iaer)=0.

               end do ! L_NSPECTV
             
               do nw=1,L_NSPECTI
                  ! Longwave
                  do l=1,nlayer

                     temp1=QIRsQREF3d(ig,nlayer+1-l,nw,iaer)         &
                          *QREFir3d(ig,nlayer+1-l,iaer)

                     temp2=QIRsQREF3d(ig,max(nlayer-l,1),nw,iaer)    &
                          *QREFir3d(ig,max(nlayer-l,1),iaer)

                     qxiaer(2*l,nw,iaer)  = temp1
                     qxiaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                     temp1=temp1*omegair3d(ig,nlayer+1-l,nw,iaer)
                     temp2=temp2*omegair3d(ig,max(nlayer-l,1),nw,iaer)

                     qsiaer(2*l,nw,iaer)  = temp1
                     qsiaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                     temp1=gir3d(ig,nlayer+1-l,nw,iaer)
                     temp2=gir3d(ig,max(nlayer-l,1),nw,iaer)

                     giaer(2*l,nw,iaer)  = temp1
                     giaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                  end do ! nlayer

                  qxiaer(1,nw,iaer)=qxiaer(2,nw,iaer)
                  qxiaer(2*nlayer+1,nw,iaer)=0.

                  qsiaer(1,nw,iaer)=qsiaer(2,nw,iaer)
                  qsiaer(2*nlayer+1,nw,iaer)=0.

                  giaer(1,nw,iaer)=giaer(2,nw,iaer)
                  giaer(2*nlayer+1,nw,iaer)=0.

               end do ! L_NSPECTI
               
            end do ! naerkind

            ! Test / Correct for freaky s. s. albedo values.
            do iaer=1,naerkind
               do k=1,L_LEVELS

                  do nw=1,L_NSPECTV
                     if(qsvaer(k,nw,iaer).gt.1.05*qxvaer(k,nw,iaer))then
                        message='Serious problems with qsvaer values' 
                        call abort_physic(subname,message,1)
                     endif
                     if(qsvaer(k,nw,iaer).gt.qxvaer(k,nw,iaer))then
                        qsvaer(k,nw,iaer)=qxvaer(k,nw,iaer)
                     endif
                  end do

                  do nw=1,L_NSPECTI 
                     if(qsiaer(k,nw,iaer).gt.1.05*qxiaer(k,nw,iaer))then
                        message='Serious problems with qsvaer values' 
                        call abort_physic(subname,message,1)
                     endif
                     if(qsiaer(k,nw,iaer).gt.qxiaer(k,nw,iaer))then
                        qsiaer(k,nw,iaer)=qxiaer(k,nw,iaer)
                     endif
                  end do

               end do ! L_LEVELS
            end do ! naerkind

!-----------------------------------------------------------------------
!     Aerosol optical depths
!-----------------------------------------------------------------------
            
         do iaer=1,naerkind     ! a bug was here           
            do k=0,nlayer-1
               
               pweight=(pplay(ig,L_NLAYRAD-k)-pplev(ig,L_NLAYRAD-k+1))/   &
                       (pplev(ig,L_NLAYRAD-k)-pplev(ig,L_NLAYRAD-k+1))
               temp=aerosol(ig,L_NLAYRAD-k,iaer)/QREFvis3d(ig,L_NLAYRAD-k,iaer)
               tauaero(2*k+2,iaer)=max(temp*pweight,0.d0)
               tauaero(2*k+3,iaer)=max(temp-tauaero(2*k+2,iaer),0.d0)

            end do
            ! boundary conditions
            tauaero(1,iaer)          = tauaero(2,iaer)
            !tauaero(1,iaer)          = 0.
            !JL18 at time of testing, the two above conditions gave the same results bit for bit. 
	    
         end do ! naerkind

         ! Albedo and Emissivity.
         albi=1-emis(ig)   ! Long Wave.
         DO nw=1,L_NSPECTV ! Short Wave loop.
            albv(nw)=albedo(ig,nw)
         ENDDO

      if ((ngrid.eq.1).and.(global1d)) then ! Fixed zenith angle 'szangle' in 1D simulations w/ globally-averaged sunlight.
         acosz = cos(pi*szangle/180.0)
         print*,'acosz=',acosz,', szangle=',szangle
      else
         acosz=mu0(ig) ! Cosine of sun incident angle : 3D simulations or local 1D simulations using latitude.
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Note by JL13 : In the following, some indices were changed in the interpolations,
!!!                so that the model results are less dependent on the number of layers !
!!!
!!!           ---  The older versions are commented with the comment !JL13index  ---
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!-----------------------------------------------------------------------
!     Water vapour (to be generalised for other gases eventually ...)
!-----------------------------------------------------------------------
      
      if(varactive)then

         i_var=igcm_h2o_vap
         do l=1,nlayer
            qvar(2*l)   = pq(ig,nlayer+1-l,i_var)
            qvar(2*l+1) = pq(ig,nlayer+1-l,i_var)    
!JL13index            qvar(2*l+1) = (pq(ig,nlayer+1-l,i_var)+pq(ig,max(nlayer-l,1),i_var))/2    
!JL13index            ! Average approximation as for temperature...
         end do
         qvar(1)=qvar(2)

      elseif(varfixed)then

         do l=1,nlayer ! Here we will assign fixed water vapour profiles globally.
            RH = satval * ((pplay(ig,l)/pplev(ig,1) - 0.02) / 0.98)
            if(RH.lt.0.0) RH=0.0
            
            call Psat_water(pt(ig,l),pplay(ig,l),psat,qsat)

            !pq_temp(l) = qsat      ! fully saturated everywhere
            pq_temp(l) = RH * qsat ! ~realistic profile (e.g. 80% saturation at ground)
         end do
         
         do l=1,nlayer
            qvar(2*l)   = pq_temp(nlayer+1-l)
            qvar(2*l+1) = (pq_temp(nlayer+1-l)+pq_temp(max(nlayer-l,1)))/2
         end do
         
         qvar(1)=qvar(2)

         ! Lowest layer of atmosphere
         RH = satval * (1 - 0.02) / 0.98
         if(RH.lt.0.0) RH=0.0

         qvar(2*nlayer+1)= RH * qsat ! ~realistic profile (e.g. 80% saturation at ground)
 
      else
         do k=1,L_LEVELS
            qvar(k) = 1.0D-7
         end do
      end if ! varactive/varfixed

      if(.not.kastprof)then
         ! IMPORTANT: Now convert from kg/kg to mol/mol.
         do k=1,L_LEVELS
            qvar(k) = qvar(k)/(epsi+qvar(k)*(1.-epsi))
         end do
      end if

!-----------------------------------------------------------------------
!     kcm mode only !
!-----------------------------------------------------------------------

      if(kastprof)then
      
         if(.not.global1d)then ! garde-fou/safeguard added by MT (to be removed in the future)
           message='You have to fix mu0, the cosinus of the solar angle'
           call abort_physic(subname,message,1)
         endif
         
         ! Initial values equivalent to mugaz.
         DO l=1,nlayer
            muvarrad(2*l)   = mugaz
            muvarrad(2*l+1) = mugaz
         END DO

         if(ngasmx.gt.1)then

            DO l=1,nlayer
               muvarrad(2*l)   =  muvar(ig,nlayer+2-l)
               muvarrad(2*l+1) = (muvar(ig,nlayer+2-l) + &
                                  muvar(ig,max(nlayer+1-l,1)))/2
            END DO
      
            muvarrad(1) = muvarrad(2)
            muvarrad(2*nlayer+1) = muvar(ig,1)

            print*,'Recalculating qvar with VARIABLE epsi for kastprof'
            print*,'Assumes that the variable gas is H2O!!!'
            print*,'Assumes that there is only one tracer'
            
            !i_var=igcm_h2o_vap
            i_var=1
            
            if(nq.gt.1)then
               message='Need 1 tracer only to run kcm1d.e' 
               call abort_physic(subname,message,1)
            endif
            
            do l=1,nlayer
               vtmp(l)=pq(ig,l,i_var)/(epsi+pq(ig,l,i_var)*(1.-epsi)) 
               !vtmp(l)=pq(ig,l,i_var)*muvar(ig,l+1)/mH2O !JL to be changed
            end do

            do l=1,nlayer
               qvar(2*l)   = vtmp(nlayer+1-l)
               qvar(2*l+1) = vtmp(nlayer+1-l)
!               qvar(2*l+1) = ( vtmp(nlayer+1-l) + vtmp(max(nlayer-l,1)) )/2
            end do
            qvar(1)=qvar(2)

            write(*,*)trim(subname),' :Warning: reducing qvar in callcorrk.'
            write(*,*)trim(subname),' :Temperature profile no longer consistent ', &
                   'with saturated H2O. qsat=',satval
                   
            do k=1,L_LEVELS
               qvar(k) = qvar(k)*satval
            end do

         endif
      else ! if kastprof
         DO l=1,nlayer
            muvarrad(2*l)   = muvar(ig,nlayer+2-l)
            muvarrad(2*l+1) = (muvar(ig,nlayer+2-l)+muvar(ig,max(nlayer+1-l,1)))/2
         END DO
      
         muvarrad(1) = muvarrad(2)
         muvarrad(2*nlayer+1)=muvar(ig,1)         
      endif ! if kastprof
      
      ! Keep values inside limits for which we have radiative transfer coefficients !!!
      if(L_REFVAR.gt.1)then ! (there was a bug here)
         do k=1,L_LEVELS
            if(qvar(k).lt.wrefvar(1))then
               qvar(k)=wrefvar(1)+1.0e-8
            elseif(qvar(k).gt.wrefvar(L_REFVAR))then
               qvar(k)=wrefvar(L_REFVAR)-1.0e-8
            endif
         end do
      endif

!-----------------------------------------------------------------------
!     Pressure and temperature
!-----------------------------------------------------------------------

      DO l=1,nlayer
         plevrad(2*l)   = pplay(ig,nlayer+1-l)/scalep
         plevrad(2*l+1) = pplev(ig,nlayer+1-l)/scalep
         tlevrad(2*l)   = pt(ig,nlayer+1-l)
         tlevrad(2*l+1) = (pt(ig,nlayer+1-l)+pt(ig,max(nlayer-l,1)))/2
      END DO
      
      plevrad(1) = 0.
!      plevrad(2) = 0.   !! JL18 enabling this line puts the radiative top at p=0 which was the idea before, but does not seem to perform best after all. 

      tlevrad(1) = tlevrad(2)
      tlevrad(2*nlayer+1)=tsurf(ig)
      
      pmid(1) = pplay(ig,nlayer)/scalep
      pmid(2) =  pmid(1)

      tmid(1) = tlevrad(2)
      tmid(2) = tmid(1)
    
      DO l=1,L_NLAYRAD-1
         tmid(2*l+1) = tlevrad(2*l+1)
         tmid(2*l+2) = tlevrad(2*l+1)
         pmid(2*l+1) = plevrad(2*l+1)
         pmid(2*l+2) = plevrad(2*l+1)
      END DO
      pmid(L_LEVELS) = plevrad(L_LEVELS)
      tmid(L_LEVELS) = tlevrad(L_LEVELS)

!!Alternative interpolation:
!         pmid(3) = pmid(1)
!         pmid(4) = pmid(1) 
!         tmid(3) = tmid(1)
!         tmid(4) = tmid(1)
!      DO l=2,L_NLAYRAD-1
!         tmid(2*l+1) = tlevrad(2*l)
!         tmid(2*l+2) = tlevrad(2*l)
!         pmid(2*l+1) = plevrad(2*l)
!         pmid(2*l+2) = plevrad(2*l)
!      END DO
!      pmid(L_LEVELS) = plevrad(L_LEVELS-1)
!      tmid(L_LEVELS) = tlevrad(L_LEVELS-1)

      ! Test for out-of-bounds pressure.
      if(plevrad(3).lt.pgasmin)then
         print*,'Minimum pressure is outside the radiative'
         print*,'transfer kmatrix bounds, exiting.'
         message="Minimum pressure outside of kmatrix bounds"
         call abort_physic(subname,message,1)
      elseif(plevrad(L_LEVELS).gt.pgasmax)then
         print*,'Maximum pressure is outside the radiative'
         print*,'transfer kmatrix bounds, exiting.'
         message="Minimum pressure outside of kmatrix bounds"
         call abort_physic(subname,message,1)
      endif

      ! Test for out-of-bounds temperature.
      ! -- JVO 20 : Also add a sanity test checking that tlevrad is
      !             within Planck function temperature boundaries,
      !             which would cause gfluxi/sfluxi to crash.
      do k=1,L_LEVELS

         if(tlevrad(k).lt.tgasmin)then
            print*,'Minimum temperature is outside the radiative'
            print*,'transfer kmatrix bounds'
            print*,"k=",k," tlevrad(k)=",tlevrad(k)
            print*,"tgasmin=",tgasmin
            if (strictboundcorrk) then
              message="Minimum temperature outside of kmatrix bounds"
              call abort_physic(subname,message,1)
            else
              print*,'***********************************************'
              print*,'we allow model to continue with tlevrad<tgasmin' 
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              !tlevrad(k)=tgasmin ! Used in the source function !
            endif
         elseif(tlevrad(k).gt.tgasmax)then
            print*,'Maximum temperature is outside the radiative'
            print*,'transfer kmatrix bounds, exiting.'
            print*,"k=",k," tlevrad(k)=",tlevrad(k)
            print*,"tgasmax=",tgasmax
            if (strictboundcorrk) then
              message="Maximum temperature outside of kmatrix bounds"
              call abort_physic(subname,message,1)
            else
              print*,'***********************************************'
              print*,'we allow model to continue with tlevrad>tgasmax'  
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              !tlevrad(k)=tgasmax ! Used in the source function !
            endif
         endif

         if (tlevrad(k).lt.tplanckmin) then
            print*,'Minimum temperature is outside the boundaries for'
            print*,'Planck function integration set in callphys.def, aborting.'
            print*,"k=",k," tlevrad(k)=",tlevrad(k)
            print*,"tplanckmin=",tplanckmin
            message="Minimum temperature outside Planck function bounds - Change tplanckmin in callphys.def"
            call abort_physic(subname,message,1)
          else if (tlevrad(k).gt.tplanckmax) then
            print*,'Maximum temperature is outside the boundaries for'
            print*,'Planck function integration set in callphys.def, aborting.'
            print*,"k=",k," tlevrad(k)=",tlevrad(k)
            print*,"tplanckmax=",tplanckmax
            message="Maximum temperature outside Planck function bounds - Change tplanckmax in callphys.def"
            call abort_physic(subname,message,1)
          endif

      enddo

      do k=1,L_NLAYRAD+1
         if(tmid(k).lt.tgasmin)then
            print*,'Minimum temperature is outside the radiative'
            print*,'transfer kmatrix bounds, exiting.'
            print*,"k=",k," tmid(k)=",tmid(k)
            print*,"tgasmin=",tgasmin
            if (strictboundcorrk) then
              message="Minimum temperature outside of kmatrix bounds"
              call abort_physic(subname,message,1)
            else
              print*,'***********************************************'
              print*,'we allow model to continue but with tmid=tgasmin'
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              tmid(k)=tgasmin
            endif
         elseif(tmid(k).gt.tgasmax)then
            print*,'Maximum temperature is outside the radiative'
            print*,'transfer kmatrix bounds, exiting.'
            print*,"k=",k," tmid(k)=",tmid(k)
            print*,"tgasmax=",tgasmax
            if (strictboundcorrk) then
              message="Maximum temperature outside of kmatrix bounds"
              call abort_physic(subname,message,1)
            else
              print*,'***********************************************'
              print*,'we allow model to continue but with tmid=tgasmax'
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              tmid(k)=tgasmax
            endif
         endif
      enddo

!=======================================================================
!          III. Calling the main radiative transfer subroutines
!=======================================================================


         Cmk= 0.01 * 1.0 / (glat(ig) * mugaz * 1.672621e-27) ! q_main=1.0 assumed.
         glat_ig=glat(ig)

!-----------------------------------------------------------------------
!        Short Wave Part
!-----------------------------------------------------------------------

         if(fract(ig) .ge. 1.0e-4) then ! Only during daylight.
            if((ngrid.eq.1).and.(global1d))then
               do nw=1,L_NSPECTV
                  stel_fract(nw)= stel(nw)* 0.25 / acosz ! globally averaged = divide by 4, and we correct for solar zenith angle
               end do
            else
               do nw=1,L_NSPECTV
                  stel_fract(nw)= stel(nw) * fract(ig)
               end do
            endif

            call optcv(dtauv,tauv,taucumv,plevrad,                 &
                 qxvaer,qsvaer,gvaer,wbarv,cosbv,tauray,tauaero,   &
                 tmid,pmid,taugsurf,qvar,muvarrad)

            call sfluxv(dtauv,tauv,taucumv,albv,dwnv,wbarv,cosbv,  &
                 acosz,stel_fract,                                 &
                 nfluxtopv,fluxtopvdn,nfluxoutv_nu,nfluxgndv_nu,   &
                 fmnetv,fluxupv,fluxdnv,fzerov,taugsurf)

         else ! During the night, fluxes = 0.
            nfluxtopv       = 0.0d0
            fluxtopvdn      = 0.0d0
            nfluxoutv_nu(:) = 0.0d0
            nfluxgndv_nu(:) = 0.0d0
            do l=1,L_NLAYRAD
               fmnetv(l)=0.0d0
               fluxupv(l)=0.0d0
               fluxdnv(l)=0.0d0
            end do
         end if


         ! Equivalent Albedo Calculation (for OUTPUT). MT2015
         if(fract(ig) .ge. 1.0e-4) then ! equivalent albedo makes sense only during daylight.       
            surface_stellar_flux=sum(nfluxgndv_nu(1:L_NSPECTV))      
            if(surface_stellar_flux .gt. 1.0e-3) then ! equivalent albedo makes sense only if the stellar flux received by the surface is positive.
               DO nw=1,L_NSPECTV                  
                  albedo_temp(nw)=albedo(ig,nw)*nfluxgndv_nu(nw)
               ENDDO
               albedo_temp(1:L_NSPECTV)=albedo_temp(1:L_NSPECTV)/surface_stellar_flux
               albedo_equivalent(ig)=sum(albedo_temp(1:L_NSPECTV))
            else
               albedo_equivalent(ig)=0.0 ! Spectrally Integrated Albedo not defined for non-irradiated grid points. So we arbitrary set the equivalent albedo to 0.
            endif
         else
            albedo_equivalent(ig)=0.0 ! Spectrally Integrated Albedo not defined for non-irradiated grid points. So we arbitrary set the equivalent albedo to 0.
         endif


!-----------------------------------------------------------------------
!        Long Wave Part
!-----------------------------------------------------------------------

         call optci(plevrad,tlevrad,dtaui,taucumi,                  &
              qxiaer,qsiaer,giaer,cosbi,wbari,tauaero,tmid,pmid,    &
              taugsurfi,qvar,muvarrad)

         call sfluxi(plevrad,tlevrad,dtaui,taucumi,ubari,albi,      &
              wnoi,dwni,cosbi,wbari,nfluxtopi,nfluxtopi_nu,         & 
              fmneti,fluxupi,fluxdni,fluxupi_nu,fzeroi,taugsurfi)

!-----------------------------------------------------------------------
!     Transformation of the correlated-k code outputs
!     (into dtlw, dtsw, fluxsurf_lw, fluxsurf_sw, fluxtop_lw, fluxtop_sw)

!     Flux incident at the top of the atmosphere
         fluxtop_dn(ig)=fluxtopvdn 

         fluxtop_lw(ig)  = real(nfluxtopi)
         fluxabs_sw(ig)  = real(-nfluxtopv)
         fluxsurf_lw(ig) = real(fluxdni(L_NLAYRAD))
         fluxsurf_sw(ig) = real(fluxdnv(L_NLAYRAD))
         
!        Flux absorbed by the surface. By MT2015.         
         fluxsurfabs_sw(ig) = fluxsurf_sw(ig)*(1.-albedo_equivalent(ig))

         if(fluxtop_dn(ig).lt.0.0)then
            print*,'Achtung! fluxtop_dn has lost the plot!'
            print*,'fluxtop_dn=',fluxtop_dn(ig)
            print*,'acosz=',acosz
            print*,'aerosol=',aerosol(ig,:,:)
            print*,'temp=   ',pt(ig,:)
            print*,'pplay=  ',pplay(ig,:)
            message="Achtung! fluxtop_dn has lost the plot!"
            call abort_physic(subname,message,1)
         endif

!     Spectral output, for exoplanet observational comparison
         if(specOLR)then
            do nw=1,L_NSPECTI 
               OLR_nu(ig,nw)=nfluxtopi_nu(nw)/DWNI(nw) !JL Normalize to the bandwidth
            end do
            do nw=1,L_NSPECTV 
               !GSR_nu(ig,nw)=nfluxgndv_nu(nw)
               OSR_nu(ig,nw)=nfluxoutv_nu(nw)/DWNV(nw) !JL Normalize to the bandwidth
            end do
         endif

!     Finally, the heating rates

         DO l=2,L_NLAYRAD
            dtsw(ig,L_NLAYRAD+1-l)=(fmnetv(l)-fmnetv(l-1))  &
                *glat(ig)/(cpp*scalep*(plevrad(2*l+1)-plevrad(2*l-1)))
            dtlw(ig,L_NLAYRAD+1-l)=(fmneti(l)-fmneti(l-1))  &
                *glat(ig)/(cpp*scalep*(plevrad(2*l+1)-plevrad(2*l-1)))
         END DO      

!     These are values at top of atmosphere
         dtsw(ig,L_NLAYRAD)=(fmnetv(1)-nfluxtopv)           &
             *glat(ig)/(cpp*scalep*(plevrad(3)-plevrad(2)))
         dtlw(ig,L_NLAYRAD)=(fmneti(1)-nfluxtopi)           &
             *glat(ig)/(cpp*scalep*(plevrad(3)-plevrad(2)))

      !  Optical thickness diagnostics (added by JVO)
      if (diagdtau) then
        do l=1,L_NLAYRAD
          do nw=1,L_NSPECTV
            int_dtauv(ig,l,nw) = 0.0d0
             DO k=1,L_NGAUSS
              ! Output exp(-tau) because gweight ponderates exp and not tau itself
              int_dtauv(ig,l,nw)= int_dtauv(ig,l,nw) + exp(-dtauv(l,nw,k))*gweight(k)
             ENDDO
          enddo
          do nw=1,L_NSPECTI
           int_dtaui(ig,l,nw) = 0.0d0
             DO k=1,L_NGAUSS
              ! Output exp(-tau) because gweight ponderates exp and not tau itself
              int_dtaui(ig,l,nw)= int_dtaui(ig,l,nw) + exp(-dtaui(l,nw,k))*gweight(k)
             ENDDO
          enddo
        enddo
      endif        


!-----------------------------------------------------------------------    
      end do ! End of big loop over every GCM column.
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!     Additional diagnostics
!-----------------------------------------------------------------------

      ! IR spectral output, for exoplanet observational comparison
      if(lastcall.and.(ngrid.eq.1))then  ! could disable the 1D output, they are in the diagfi and diagspec... JL12

         print*,'Saving scalar quantities in surf_vals.out...'
         print*,'psurf = ', pplev(1,1),' Pa'
         open(116,file='surf_vals.out')
         write(116,*) tsurf(1),pplev(1,1),fluxtop_dn(1),         &
                      real(-nfluxtopv),real(nfluxtopi) 
         close(116)


!          USEFUL COMMENT - Do Not Remove.
!
!           if(specOLR)then
!               open(117,file='OLRnu.out')
!               do nw=1,L_NSPECTI
!                  write(117,*) OLR_nu(1,nw)
!               enddo
!               close(117)
!
!               open(127,file='OSRnu.out')
!               do nw=1,L_NSPECTV
!                  write(127,*) OSR_nu(1,nw)
!               enddo
!               close(127)
!           endif

           ! OLR vs altitude: do it as a .txt file.
         OLRz=.false.
         if(OLRz)then
            print*,'saving IR vertical flux for OLRz...'
            open(118,file='OLRz_plevs.out')
            open(119,file='OLRz.out')
            do l=1,L_NLAYRAD
               write(118,*) plevrad(2*l)
               do nw=1,L_NSPECTI
                  write(119,*) fluxupi_nu(l,nw) 
               enddo
            enddo 
            close(118)
            close(119)
         endif

      endif

      ! See physiq.F for explanations about CLFvarying. This is temporary.
      if (lastcall .and. .not.CLFvarying) then
        IF( ALLOCATED( gasi ) ) DEALLOCATE( gasi )
        IF( ALLOCATED( gasv ) ) DEALLOCATE( gasv )
!$OMP BARRIER
!$OMP MASTER
        IF( ALLOCATED( pgasref ) ) DEALLOCATE( pgasref )
        IF( ALLOCATED( tgasref ) ) DEALLOCATE( tgasref )
        IF( ALLOCATED( wrefvar ) ) DEALLOCATE( wrefvar )
        IF( ALLOCATED( pfgasref ) ) DEALLOCATE( pfgasref )
        IF( ALLOCATED( gweight ) ) DEALLOCATE( gweight )
!$OMP END MASTER
!$OMP BARRIER
        IF ( ALLOCATED(reffrad)) DEALLOCATE(reffrad)
        IF ( ALLOCATED(nueffrad)) DEALLOCATE(nueffrad)
      endif


    end subroutine callcorrk

END MODULE callcorrk_mod
