module tracer_mod

 implicit none
 
      ! number of tracers:
      integer,save :: nqmx ! initialized in conf_phys
   
      character*30,allocatable,save ::  noms(:)  ! name of the tracer
      real,allocatable,save :: mmol(:)           ! mole mass of tracer (g/mol-1) 
      real,allocatable,save :: radius(:)   ! dust and ice particle radius (m)
      real,allocatable,save :: rho_q(:)    ! tracer densities (kg.m-3)
      real,allocatable,save :: alpha_lift(:) ! saltation vertical flux/horiz flux ratio (m-1)
      real,allocatable,save :: alpha_devil(:) ! lifting coeeficient by dust devil

      real,save :: varian      ! Characteristic variance of log-normal distribution
      real,save :: r3n_q     ! used to compute r0 from number and mass mixing ratio
      real,save :: rho_dust     ! Mars dust density (kg.m-3)
      real,save :: rho_ice     ! Water ice density (kg.m-3)
      real,save :: nuice_ref   ! Effective variance of the water ice dist.
      real,save :: nuice_sed   ! Sedimentation effective variance of the water ice dist.
      real,save :: ref_r0        ! for computing reff=ref_r0*r0 (in log.n. distribution)
      real,save :: rho_ice_co2     ! co2 ice density (kg.m-3)
      real,save :: nuiceco2_sed   ! Sedimentation effective variance of the co2 ice dist.
      real,save :: nuiceco2_ref   ! Effective variance of the co2 ice dist.
      
      real,save :: ccn_factor  ! ratio of nuclei for water ice particles

      INTEGER,ALLOCATABLE,SAVE :: nqdust(:) ! to store the indexes of dust tracers (cf aeropacity)
      real,allocatable,save :: dryness(:)!"Dryness coefficient" for grnd water ice sublimation


! tracer indexes: these are initialized in initracer and should be 0 if the
!                 corresponding tracer does not exist
      ! dust
      integer,allocatable,save :: igcm_dustbin(:) ! for dustbin 'dust' tracers
      ! dust, special doubleq case
      integer,save :: igcm_dust_mass   ! dust mass mixing ratio
                                  !   (for transported dust)
      integer,save :: igcm_dust_number ! dust number mixing ratio
                                  !   (transported dust)
      integer,save :: igcm_ccn_mass   ! CCN mass mixing ratio
      integer,save :: igcm_ccn_number ! CCN number mixing ratio
      integer,save :: igcm_dust_submicron ! submicron dust mixing ratio
      integer,save :: igcm_stormdust_mass   !  storm dust mass mixing ratio
      integer,save :: igcm_stormdust_number !  storm dust number mixing ratio
      integer,save :: igcm_topdust_mass   !  topdust mass mixing ratio
      integer,save :: igcm_topdust_number !  topdust number mixing ratio

      integer,save :: igcm_ccnco2_mass   ! CCN (dust and/or water ice) for CO2 mass mixing ratio
      integer,save :: igcm_ccnco2_number ! CCN (dust and/or water ice) for CO2 number mixing ratio

      ! water
      integer,save :: igcm_h2o_vap ! water vapour
      integer,save :: igcm_h2o_ice ! water ice
      integer,save :: igcm_co2_ice ! co2 ice

      ! chemistry:
      integer,save :: igcm_co2
      integer,save :: igcm_co
      integer,save :: igcm_o
      integer,save :: igcm_o1d
      integer,save :: igcm_o2
      integer,save :: igcm_o3
      integer,save :: igcm_h
      integer,save :: igcm_h2
      integer,save :: igcm_oh
      integer,save :: igcm_ho2
      integer,save :: igcm_h2o2
      integer,save :: igcm_n2
      integer,save :: igcm_ar
      integer,save :: igcm_n
      integer,save :: igcm_no
      integer,save :: igcm_no2
      integer,save :: igcm_n2d
      integer,save :: igcm_he
      integer,save :: igcm_ch4
      ! Ions
      integer,save :: igcm_co2plus
      integer,save :: igcm_oplus
      integer,save :: igcm_o2plus
      integer,save :: igcm_coplus
      integer,save :: igcm_cplus
      integer,save :: igcm_nplus
      integer,save :: igcm_noplus
      integer,save :: igcm_n2plus 
      integer,save :: igcm_hplus
      integer,save :: igcm_hco2plus
      integer,save :: igcm_hcoplus
      integer,save :: igcm_elec
      ! other tracers
      integer,save :: igcm_ar_n2 ! for simulations using co2 +neutral gas


!-----------------------------------------------------------------------

  contains
  
    subroutine ini_tracer_mod(nq,tname)
      implicit none
      
      integer,intent(in) :: nq ! number of tracers
      character(len=*),intent(in) :: tname(nq) ! tracer names
      
      integer :: iq, count
      character(len=20) :: txt ! to store some text
      
      ! set dimension and tracer names
      nqmx=nq
      allocate(noms(nq))
      do iq=1,nq
        noms(iq)=tname(iq)
        write(*,*) "tracer_mod names : ", trim(noms(iq))
      enddo
     
#ifndef MESOSCALE 
      ! check if tracers have 'old' names
      count=0
      do iq=1,nq
        txt=" "
        write(txt,'(a1,i2.2)') 'q',iq
        if (txt.eq.tname(iq)) then
          count=count+1
        endif
      enddo ! of do iq=1,nq
      
      if ((count.eq.nq).and.(nq.ne.0)) then
        write(*,*) "ini_tracer_mod: tracers seem to follow old naming ", &
                   "convention (q01,q02,...)"
        write(*,*) "you should run newstart to rename them"
        call abort_physic("ini_tracer_mod","tracer name issue",1)
      endif
#endif
            
      ! allocate module arrays:
      ! -- not domain-dependent
      allocate(mmol(nq))
      allocate(radius(nq))
      allocate(rho_q(nq))
      allocate(alpha_lift(nq))
      allocate(alpha_devil(nq))
      allocate(igcm_dustbin(nq))
      allocate(nqdust(nq))
      
    end subroutine ini_tracer_mod

    subroutine end_tracer_mod

    implicit none

      if (allocated(noms)) deallocate(noms)
      if (allocated(mmol)) deallocate(mmol)
      if (allocated(radius)) deallocate(radius)
      if (allocated(rho_q)) deallocate(rho_q)
      if (allocated(alpha_lift)) deallocate(alpha_lift)
      if (allocated(alpha_devil)) deallocate(alpha_devil)
      if (allocated(igcm_dustbin)) deallocate(igcm_dustbin)
      if (allocated(nqdust)) deallocate(nqdust)

    end subroutine end_tracer_mod

end module tracer_mod
