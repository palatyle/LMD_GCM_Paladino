      subroutine deposition_source_early_earth(ngrid, nlayer, nq,   &
                           ig, zzlay, zzlev,zdens,      &
                           zycol, ptimestep)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     dry deposition of chemical species
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      use gases_h
      use tracer_h, only:   igcm_co2, igcm_co, igcm_o, igcm_o1d, igcm_o2, &
                            igcm_o3, igcm_h, igcm_h2, igcm_oh, igcm_ho2,  &
                            igcm_h2o2, igcm_ch4, igcm_n2, igcm_h2o_vap,   &
                            igcm_no, igcm_n, igcm_no2, igcm_n2d,          &
                            igcm_ch3, igcm_ch, igcm_3ch2, igcm_1ch2,      &
                            igcm_cho, igcm_ch2o, igcm_ch3o,               &
                            igcm_c, igcm_c2, igcm_c2h, igcm_c2h2,         &
                            igcm_c2h3, igcm_c2h4, igcm_c2h6, igcm_ch2co,  &
                            igcm_ch3co, igcm_hcaer

      implicit none
!
!
!     input
! 
      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of tracers
      integer ig                     ! grid point index
      real    zzlay(ngrid,nlayer)    ! altitude at the middle of the layers (m)
      real    zzlev(ngrid,nlayer+1)  ! altitude at layer boundaries (m)
      real    zdens(nlayer)  ! density (cm^-3)
      real    zycol(nlayer,nq)       ! composition (volume mixing ratio)
      real    ptimestep              ! physical timestep (s)
!
!     local
!
      real    vd                         ! dry deposition velocity (cm.s-1)
      real    deltaz                     ! thickness of first layer (m)
      real    loss                       ! loss rate (s-1)                       
      real    prod                       ! production rate (s-1) 
      integer iq
      logical, save :: firstcall = .true.

      integer, allocatable, save ::  SF_mode(:)
      real, allocatable, save ::  SF_value(:)
      real, allocatable, save ::  prod_rate(:)

      if (.not.allocated(SF_mode)) allocate(SF_mode(nq))
      if (.not.allocated(SF_value)) allocate(SF_value(nq))
      if (.not.allocated(prod_rate)) allocate(prod_rate(nq))

if (firstcall) then
   print*,'photochemistry: initialize deposition/source'
!   print*,'SF_mode=1 (fixed mixing ratio) / 2 (fixed sedimentation velocity in cm/s and/or flux in molecules/m^3)'
   do iq=1,nq
     SF_mode(iq)=2
     SF_value(iq)=0.
     prod_rate(iq)=0.
   enddo 

! Cases SF_mode=1 (fixed mixing ratio)

   SF_mode(igcm_co2)=1
   SF_value(igcm_co2)=gfrac(igas_CO2)
   SF_mode(igcm_ch4)=1
   SF_value(igcm_ch4)=gfrac(igas_CH4)
!   SF_mode(igcm_o2)=1
!   SF_value(igcm_o2)=1.e-8
   SF_mode(igcm_h2)=1
   SF_value(igcm_h2)=1.e-3

! Cases SF_mode=2 with vsed /=0
   SF_value(igcm_o)=1.
   SF_value(igcm_o2)=1.
   SF_value(igcm_h)=1.
   SF_value(igcm_oh)=1.
   SF_value(igcm_ho2)=1.
   SF_value(igcm_h2o2)=0.2
!   SF_value(igcm_h2)=2.4e-4
   SF_value(igcm_co)=1.2e-4
   SF_value(igcm_cho)=1.
   SF_value(igcm_ch2o)=0.2
   SF_value(igcm_ch3)=1
   SF_value(igcm_o3)=0.07
   SF_value(igcm_hcaer)=0.01

! Cases SF_mode=2 (production flux in mol/m2/s)
!   prod_rate(igcm_h2)=3.5e13

firstcall=.false.
endif !first call


!        thickness of first layer (m)
deltaz = zzlev(ig,2) - zzlev(ig,1)


do iq=1,nq

  if(SF_mode(iq).eq.1) then
!    zycol(1,iq) = SF_value(iq) +(zycol(1,iq)-SF_value(iq))*exp(-ptimestep/8.64e4)
    zycol(1,iq) = SF_value(iq)

  elseif(SF_mode(iq).eq.2) then
!        loss/prod rate (s-1)
    loss = 0.01*SF_value(iq)/deltaz
    prod=prod_rate(iq)/(zdens(1)*1e6*deltaz)
    zycol(1,iq) = zycol(1,iq)*exp(-loss*ptimestep)+prod*ptimestep

  endif

enddo

      return
      end
