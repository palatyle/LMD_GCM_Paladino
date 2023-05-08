      subroutine calchim_asis(ngrid,nlayer,nq,                              &
                         ptimestep,pplay,pplev,pt,pdt,dist_sol,mu0,fract,   &
                         zzlev,zzlay,zday,pq,pdq,dqchim,dqschim,            &
                         tauref,co2ice,                                     &
                         pu,pdu,pv,pdv,surfdust,surfice)

      use tracer_h, only:   igcm_co2, igcm_co, igcm_o, igcm_o1d, igcm_o2, &
                            igcm_o3, igcm_h, igcm_h2, igcm_oh, igcm_ho2,  &
                            igcm_h2o2, igcm_ch4, igcm_n2, igcm_h2o_vap,   &
                            igcm_n, igcm_no, igcm_no2, igcm_n2d,   &
                            mmol

      use conc_mod, only: mmean ! mean molecular mass of the atmosphere
      USE comcstfi_mod
      use callkeys_mod

      implicit none

!=======================================================================
!
!   subject:
!   --------
!
!  Prepare the call for the photochemical module, and send back the
!  tendencies from photochemistry in the chemical species mass mixing ratios
!
!   Author:   Sebastien Lebonnois (08/11/2002)
!   -------
!    update 12/06/2003 for water ice clouds and compatibility with dust
!    update 07/2003 for coupling with thermosphere (Monica Angelats-i-Coll)
!    update 03/05/2005 cosmetic changes (Franck Lefevre)
!    update sept. 2008 identify tracers by their names (Ehouarn Millour)
!    update 05/12/2011 synchronize with latest version of chemistry (Franck Lefevre)
!    update 16/03/2012 optimization (Franck Lefevre)
!    update 11/12/2013 optimization (Franck Lefevre)
!    update 20/10/2017 adapted to LMDZ GENERIC+cosmetic changes (Benjamin Charnay)
!
!   Arguments:
!   ----------
!
!  Input:
!
!    ptimestep                  timestep (s)
!    pplay(ngrid,nlayer)    Pressure at the middle of the layers (Pa)
!    pplev(ngrid,nlayer+1)  Intermediate pressure levels (Pa)
!    pt(ngrid,nlayer)       Temperature (K)
!    pdt(ngrid,nlayer)      Temperature tendency (K)
!    pu(ngrid,nlayer)       u component of the wind (ms-1)
!    pdu(ngrid,nlayer)      u component tendency (K)
!    pv(ngrid,nlayer)       v component of the wind (ms-1)
!    pdv(ngrid,nlayer)      v component tendency (K)
!    dist_sol                   distance of the sun (AU)
!    mu0(ngrid)               cos of solar zenith angle (=1 when sun at zenith)
!    fract(ngrid)           day fraction (for diurnal=false)
!    pq(ngrid,nlayer,nq)  Advected fields, ie chemical species here
!    pdq(ngrid,nlayer,nq) Previous tendencies on pq
!    tauref(ngrid)            Optical depth at 7 hPa
!    co2ice(ngrid)            co2 ice surface layer (kg.m-2)
!    surfdust(ngrid,nlayer) dust surface area (m2/m3)
!    surfice(ngrid,nlayer)  ice surface area (m2/m3)
!
!  Output:
!
!    dqchim(ngrid,nlayer,nq) ! tendencies on pq due to chemistry
!    dqschim(ngrid,nq)         ! tendencies on qsurf 
!
!=======================================================================

#include "chimiedata.h"

!     input:

      integer,intent(in) :: ngrid    ! number of atmospheric columns
      integer,intent(in) :: nlayer   ! number of atmospheric layers
      integer,intent(in) :: nq       ! number of tracers
      real :: ptimestep
      real :: pplay(ngrid,nlayer)    ! pressure at the middle of the layers
      real :: zzlay(ngrid,nlayer)    ! pressure at the middle of the layers
      real :: pplev(ngrid,nlayer+1)  ! intermediate pressure levels
      real :: zzlev(ngrid,nlayer+1)  ! altitude at layer boundaries
      real :: pt(ngrid,nlayer)       ! temperature
      real :: pdt(ngrid,nlayer)      ! temperature tendency
      real :: pu(ngrid,nlayer)       ! u component of the wind (m.s-1)
      real :: pdu(ngrid,nlayer)      ! u component tendency
      real :: pv(ngrid,nlayer)       ! v component of the wind (m.s-1)
      real :: pdv(ngrid,nlayer)      ! v component tendency
      real :: dist_sol               ! distance of the sun (AU)
      real :: mu0(ngrid)             ! cos of solar zenith angle (=1 when sun at zenith)
      real :: pq(ngrid,nlayer,nq)    ! tracers mass mixing ratio
      real :: pdq(ngrid,nlayer,nq)   ! previous tendencies
      real :: zday                   ! date (time since Ls=0, in martian days)
      real :: tauref(ngrid)          ! optical depth at 7 hPa
      real :: co2ice(ngrid)          ! co2 ice surface layer (kg.m-2)
      real :: surfdust(ngrid,nlayer) ! dust surface area (m2/m3)
      real :: surfice(ngrid,nlayer)  !  ice surface area (m2/m3)

!     output:

      real :: dqchim(ngrid,nlayer,nq)   ! tendencies on pq due to chemistry
      real :: dqschim(ngrid,nq)         ! tendencies on qsurf 

!     local variables:

      integer,save :: nbq                   ! number of tracers used in the chemistry
      integer,allocatable,save :: niq(:)    ! array storing the indexes of the tracers
      integer :: iloc(1)                    ! index of major species
      integer :: ig,l,i,iq,iqmax
      integer :: foundswitch, lswitch
      integer,save :: chemthermod

      integer,save :: i_co2  = 0
      integer,save :: i_co   = 0
      integer,save :: i_o    = 0
      integer,save :: i_o1d  = 0
      integer,save :: i_o2   = 0
      integer,save :: i_o3   = 0
      integer,save :: i_h    = 0
      integer,save :: i_h2   = 0
      integer,save :: i_oh   = 0
      integer,save :: i_ho2  = 0
      integer,save :: i_h2o2 = 0
      integer,save :: i_ch4  = 0
      integer,save :: i_n2   = 0
      integer,save :: i_h2o  = 0
      integer,save :: i_n    = 0
      integer,save :: i_no    = 0
      integer,save :: i_no2    = 0
      integer,save :: i_n2d    = 0

      integer :: ig_vl1

      real    :: latvl1, lonvl1
      real    :: zq(ngrid,nlayer,nq)   ! pq+pdq*ptimestep before chemistry
                                       ! new mole fraction after
      real    :: zt(ngrid,nlayer)      ! temperature
      real    :: zu(ngrid,nlayer)      ! u component of the wind
      real    :: zv(ngrid,nlayer)      ! v component of the wind
      real    :: taucol                ! optical depth at 7 hPa

      logical,save :: firstcall = .true.
      logical,save :: depos = .false.  ! switch for dry deposition

!     for each column of atmosphere:

      real :: zpress(nlayer)       ! Pressure (mbar)
      real :: zdens(nlayer)        ! Density  (cm-3)
      real :: ztemp(nlayer)        ! Temperature (K)
      real :: zlocal(nlayer)       ! Altitude (km)
      real :: zycol(nlayer,nq)     ! Composition (mole fractions)
      real :: zmmean(nlayer)       ! Mean molecular mass (g.mole-1)
      real :: szacol               ! Solar zenith angle
      real :: fract(ngrid)         ! day fraction
      real :: fractcol             ! day fraction
      real :: surfice1d(nlayer)    ! Ice surface area (cm2/cm3)
      real :: surfdust1d(nlayer)   ! Dust surface area (cm2/cm3)
      real :: jo3(nlayer)          ! Photodissociation rate O3->O1D (s-1)

      integer :: iter(nlayer)      !  Number of chemical iterations
                                   !  within one physical timestep 

!     for output:

      logical :: output             ! to issue calls to writediagfi and stats
      parameter (output = .true.)
      real :: jo3_3d(ngrid,nlayer)  ! Photodissociation rate O3->O1D (s-1)
      real :: iter_3d(ngrid,nlayer) ! Number of chemical iterations
                                    !  within one physical timestep

!=======================================================================
!     initialization of the chemistry (first call only)
!=======================================================================

      if (firstcall) then

         if (photochem) then
            print*,'calchim: Read photolysis lookup table'
            call read_phototable
         end if
         ! find index of chemical tracers to use
         allocate(niq(nq))
         ! Listed here are all tracers that can go into photochemistry
         nbq = 0 ! to count number of tracers
         ! Species ALWAYS present if photochem=.T. or thermochem=.T.
         i_co2 = igcm_co2
         if (i_co2 == 0) then
            write(*,*) "calchim: Error; no CO2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_co2
         end if
         i_co = igcm_co
         if (i_co == 0) then
            write(*,*) "calchim: Error; no CO tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_co
         end if
         i_o = igcm_o
         if (i_o == 0) then
            write(*,*) "calchim: Error; no O tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o
         end if
         i_o1d = igcm_o1d
         if (i_o1d == 0) then
            write(*,*) "calchim: Error; no O1D tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o1d
         end if
         i_o2 = igcm_o2
         if (i_o2 == 0) then
            write(*,*) "calchim: Error; no O2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o2
         end if
         i_o3 = igcm_o3
         if (i_o3 == 0) then
            write(*,*) "calchim: Error; no O3 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o3
         end if
         i_h = igcm_h
         if (i_h == 0) then
            write(*,*) "calchim: Error; no H tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h
         end if
         i_h2 = igcm_h2
         if (i_h2 == 0) then
            write(*,*) "calchim: Error; no H2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h2
         end if
         i_oh = igcm_oh
         if (i_oh == 0) then
            write(*,*) "calchim: Error; no OH tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_oh
         end if
         i_ho2 = igcm_ho2
         if (i_ho2 == 0) then
            write(*,*) "calchim: Error; no HO2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_ho2
         end if
         i_h2o2 = igcm_h2o2
         if (i_h2o2 == 0) then
            write(*,*) "calchim: Error; no H2O2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h2o2
         end if
         i_ch4 = igcm_ch4
         if (i_ch4 == 0) then
            write(*,*) "calchim: Error; no CH4 tracer !!!"
            write(*,*) "CH4 will be ignored in the chemistry"
         else
            nbq = nbq + 1
            niq(nbq) = i_ch4
         end if
         i_n2 = igcm_n2
         if (i_n2 == 0) then
            write(*,*) "calchim: Error; no N2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_n2
         end if
         i_n = igcm_n
         if (i_n == 0) then
            if (photochem) then
               write(*,*) "calchim: Error; no N tracer !!!"
               write(*,*) "N will be ignored in the chemistry"
            end if
         else
            nbq = nbq + 1
            niq(nbq) = i_n
         end if
         i_n2d = igcm_n2d
         if (i_n2d == 0) then
            if (photochem) then
               write(*,*) "calchim: Error; no N2D tracer !!!"
               write(*,*) "N2d will be ignored in the chemistry"
            end if
         else
            nbq = nbq + 1
            niq(nbq) = i_n2d
         end if
         i_no = igcm_no
         if (i_no == 0) then
            if (photochem) then
               write(*,*) "calchim: Error; no NO tracer !!!"
               write(*,*) "NO will be ignored in the chemistry"
            end if
         else
            nbq = nbq + 1
            niq(nbq) = i_no
         end if
         i_no2 = igcm_no2
         if (i_no2 == 0) then
            if (photochem) then
               write(*,*) "calchim: Error; no NO2 tracer !!!"
               write(*,*) "NO2 will be ignored in the chemistry"
            end if
         else
            nbq = nbq + 1
            niq(nbq) = i_no2
         end if
         i_h2o = igcm_h2o_vap
         if (i_h2o == 0) then
            write(*,*) "calchim: Error; no water vapor tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h2o
         end if
         !Check tracers needed for thermospheric chemistry
!         if(thermochem) then
!            chemthermod=0  !Default: C/O/H chemistry
!            !Nitrogen chemistry
!            !NO is used to determine if N chemistry is wanted
!            !chemthermod=2 -> N chemistry
!            if (i_no == 0) then
!               write(*,*) "calchim: no NO tracer"
!               write(*,*) "C/O/H themosp chemistry only "
!            else
!               chemthermod=2
!               write(*,*) "calchim: NO in traceur.def"
!               write(*,*) "Nitrogen chemistry included"
!            end if
!            ! N
!            if(chemthermod == 2) then
!               if (i_n == 0) then
!                  write(*,*) "calchim: Error; no N tracer !!!"
!                  write(*,*) "N is needed if NO is in traceur.def"
!                  stop
!               end if
!            ! NO2
!               if (i_no2 == 0) then
!                  write(*,*) "calchim: Error; no NO2 tracer !!!"
!                  write(*,*) "NO2 is needed if NO is in traceur.def"
!                  stop
!               end if
!            ! N(2D)
!               if (i_n2d == 0) then
!                  write(*,*) "calchim: Error; no N2D !!!"
!                  write(*,*) "N2D is needed if NO is in traceur.def"
!                  stop
!               end if
!            endif    !Of if(chemthermod == 2)
!         endif      !Of thermochem

         write(*,*) 'calchim: found nbq    = ',nbq,' tracers'
               
         firstcall = .false.
      end if ! if (firstcall)

! Initializations

      zycol(:,:)    = 0.
      dqchim(:,:,:) = 0.
      dqschim(:,:)  = 0.

!     latvl1= 22.27
!     lonvl1= -47.94
!     ig_vl1= 1+ int( (1.5-(latvl1-90.)*jjm/180.)  -2 )*iim +    &
!             int(1.5+(lonvl1+180)*iim/360.)

!=======================================================================
!     loop over grid
!=======================================================================

      do ig = 1,ngrid
         
         foundswitch = 0
         do l = 1,nlayer
            do i = 1,nbq
               iq = niq(i) ! get tracer index
               zq(ig,l,iq) = pq(ig,l,iq) + pdq(ig,l,iq)*ptimestep
               zycol(l,iq) = zq(ig,l,iq)*mmean(ig,l)/mmol(iq)
            end do
            zt(ig,l)  = pt(ig,l) + pdt(ig,l)*ptimestep
            zu(ig,l)  = pu(ig,l) + pdu(ig,l)*ptimestep
            zv(ig,l)  = pv(ig,l) + pdv(ig,l)*ptimestep
            zpress(l) = pplay(ig,l)/100.
            ztemp(l)  = zt(ig,l)
            zdens(l)  = zpress(l)/(kb*1.e4*ztemp(l))
            zlocal(l) = zzlay(ig,l)/1000.
            zmmean(l) = mmean(ig,l)

!           surfdust1d and surfice1d: conversion from m2/m3 to cm2/cm3

            surfdust1d(l) = surfdust(ig,l)*1.e-2
            surfice1d(l)  = surfice(ig,l)*1.e-2

!           search for switch index between regions

!            if (photochem .and. thermochem) then
!               if (foundswitch == 0 .and. pplay(ig,l) < 1.e-1) then
!                  lswitch = l
!                  foundswitch = 1
!               end if
!            end if
            if (.not. photochem) then
               lswitch = 22
            end if
!            if (.not. thermochem) then
               lswitch = min(50,nlayer+1)
!            end if

         end do ! of do l=1,nlayer

         szacol = acos(mu0(ig))*180./pi
         taucol = tauref(ig)*(700./610.)  ! provisoire en attente de nouveau jmars
         fractcol=fract(ig)

!=======================================================================
!     call chemical subroutines
!=======================================================================

!        chemistry in lower atmosphere

         if (photochem) then

            call photochemistry_asis(nlayer,nq,ngrid,                          &
                                ig,lswitch,zycol,szacol,fractcol,ptimestep,    &
                                zpress,ztemp,zdens,zmmean,dist_sol,            &
                                surfdust1d,surfice1d,jo3,taucol,iter)

!        ozone photolysis, for output

            do l = 1,nlayer
               jo3_3d(ig,l) = jo3(l)
               iter_3d(ig,l) = iter(l)
            end do

!        condensation of h2o2

!            call perosat(ngrid, nlayer, nq,                        &
!                         ig,ptimestep,pplev,pplay,                 &
!                         ztemp,zycol,dqcloud,dqscloud)
         end if

!        chemistry in upper atmosphere
        
!         if (thermochem) then
!            call chemthermos(ig,nlayer,lswitch,chemthermod,zycol,ztemp,&
!                             zdens,zpress,zlocal,szacol,ptimestep,zday)
!         end if

!        dry deposition

!         if (depos) then
!            call deposition(ngrid, nlayer, nq,                      &
!                            ig, ig_vl1, pplay, pplev, zzlay, zzlev, & 
!                            zu, zv, zt, zycol, ptimestep, co2ice)
!         end if

!=======================================================================
!     tendencies
!=======================================================================

!     index of the most abundant species at each level

!         major(:) = maxloc(zycol, dim = 2)

!     tendency for the most abundant species = - sum of others

         do l = 1,nlayer
            iloc=maxloc(zycol(l,:))
            iqmax=iloc(1)
            do i = 1,nbq
               iq = niq(i) ! get tracer index
               if (iq /= iqmax) then
                  dqchim(ig,l,iq) = (zycol(l,iq)*mmol(iq)/mmean(ig,l) - zq(ig,l,iq))/ptimestep
                  dqchim(ig,l,iq) = max(dqchim(ig,l,iq),-zq(ig,l,iq)/ptimestep)
                  dqchim(ig,l,iqmax) = dqchim(ig,l,iqmax) - dqchim(ig,l,iq) 
               end if
            end do
         end do ! of do l = 1,nlayer






!=======================================================================
!     end of loop over grid
!=======================================================================

      end do ! of do ig=1,ngrid

!=======================================================================
!     write outputs
!=======================================================================

! value of parameter 'output' to trigger writting of outputs
! is set above at the declaration of the variable.

      if (photochem .and. output) then
            call writediagfi(ngrid,'jo3','j o3->o1d',    &
                             's-1',3,jo3_3d(1,1))
            call writediagfi(ngrid,'iter','iterations',  &
                             ' ',3,iter_3d(1,1))
            if (callstats) then
               call wstats(ngrid,'jo3','j o3->o1d',       &
                           's-1',3,jo3_3d(1,1))
               call wstats(ngrid,'mmean','mean molecular mass',       &
                           'g.mole-1',3,mmean(1,1))
            endif
      end if ! of if (output)

      return
      end
