      SUBROUTINE initracer(ngrid,nq,nametrac)

      use surfdat_h, ONLY: dryness, watercaptag
      USE tracer_h
      USE callkeys_mod, only: water, rho_volc
      use ioipsl_getin_p_mod, only: getin_p
      IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Initialization related to tracer 
c   (transported dust, water, chemical species, ice...)
c
c   Name of the tracer
c
c   Test of dimension :
c   Initialize COMMON tracer in tracer.h, using tracer names provided
c   by the argument nametrac
c
c   author: F.Forget
c   ------
c            Ehouarn Millour (oct. 2008) identify tracers by their names
c=======================================================================

      integer,intent(in) :: ngrid,nq
      character(len=30),intent(in) :: nametrac(nq) ! name of the tracer from dynamics

      character(len=30) :: txt ! to store some text
      integer iq,ig,count
      real r0_lift , reff_lift



c-----------------------------------------------------------------------
c  radius(nq)      ! aerosol particle radius (m)
c  rho_q(nq)       ! tracer densities (kg.m-3)
c  qext(nq)        ! Single Scat. Extinction coeff at 0.67 um
c  alpha_lift(nq)  ! saltation vertical flux/horiz flux ratio (m-1)
c  alpha_devil(nq) ! lifting coeeficient by dust devil
c  rho_dust          ! Mars dust density
c  rho_ice           ! Water ice density
c  rho_volc          ! Volcanic ash density
c  doubleq           ! if method with mass (iq=1) and number(iq=2) mixing ratio
c  varian            ! Characteristic variance of log-normal distribution
c-----------------------------------------------------------------------

       nqtot=nq
       !! we allocate once for all arrays in common in tracer_h.F90
       !! (supposedly those are not used before call to initracer)
       IF (.NOT.ALLOCATED(noms))         ALLOCATE(noms(nq))
       IF (.NOT.ALLOCATED(mmol))         ALLOCATE(mmol(nq))
       IF (.NOT.ALLOCATED(radius))       ALLOCATE(radius(nq))
       IF (.NOT.ALLOCATED(rho_q))        ALLOCATE(rho_q(nq))
       IF (.NOT.ALLOCATED(qext))         ALLOCATE(qext(nq))
       IF (.NOT.ALLOCATED(alpha_lift))   ALLOCATE(alpha_lift(nq))
       IF (.NOT.ALLOCATED(alpha_devil))  ALLOCATE(alpha_devil(nq))
       IF (.NOT.ALLOCATED(qextrhor))     ALLOCATE(qextrhor(nq))
       IF (.NOT.ALLOCATED(igcm_dustbin)) ALLOCATE(igcm_dustbin(nq))

!might not need this
!       ALLOCATE(igcm_volc_1(nq))
!       ALLOCATE(igcm_volc_2(nq))
!       ALLOCATE(igcm_volc_3(nq))
!       ALLOCATE(igcm_volc_4(nq))
!       ALLOCATE(igcm_volc_5(nq))
!       ALLOCATE(igcm_volc_6(nq))
!       ALLOCATE(igcm_h2so4(nq))
       !! initialization
       alpha_lift(:)=0.
       alpha_devil(:)=0.
       
       ! Added by JVO 2017 : these arrays are handled later
       ! -> initialization is the least we can do, please !!!
       radius(:)=0.
       qext(:)=0.


! Initialization: copy tracer names from dynamics
        do iq=1,nq
          noms(iq)=nametrac(iq)
          write(*,*)"initracer: iq=",iq,"noms(iq)=",trim(noms(iq))
        enddo


! Identify tracers by their names: (and set corresponding values of mmol)
      ! 0. initialize tracer indexes to zero:
      ! NB: igcm_* indexes are commons in 'tracer.h'
      do iq=1,nq
        igcm_dustbin(iq)=0
      enddo
      igcm_dust_mass=0
      igcm_dust_number=0
      igcm_h2o_vap=0
      igcm_h2o_ice=0
      igcm_co2=0
      igcm_co=0
      igcm_o=0
      igcm_o1d=0
      igcm_o2=0
      igcm_o3=0
      igcm_h=0
      igcm_h2=0
      igcm_oh=0
      igcm_ho2=0
      igcm_h2o2=0
      igcm_n2=0
      igcm_n=0
      igcm_n2d=0
      igcm_no=0
      igcm_no2=0
      igcm_ar=0
      igcm_ar_n2=0
      igcm_co2_ice=0
      igcm_volc_1=0
      igcm_volc_2=0
      igcm_volc_3=0
      igcm_volc_4=0
      igcm_volc_5=0
      igcm_volc_6=0
      igcm_h2so4=0

      igcm_ch4=0
      igcm_ch3=0
      igcm_ch=0
      igcm_3ch2=0
      igcm_1ch2=0
      igcm_cho=0
      igcm_ch2o=0
      igcm_ch3o=0
      igcm_c=0
      igcm_c2=0
      igcm_c2h=0
      igcm_c2h2=0
      igcm_c2h3=0
      igcm_c2h4=0
      igcm_c2h6=0
      igcm_ch2co=0
      igcm_ch3co=0
      igcm_hcaer=0



      ! 1. find dust tracers
      count=0

      ! 2. find chemistry and water tracers
      do iq=1,nq
        if (noms(iq).eq."co2") then
          igcm_co2=iq
          mmol(igcm_co2)=44.
          count=count+1
!          write(*,*) 'co2: count=',count
        endif
        if (noms(iq).eq."co2_ice") then
          igcm_co2_ice=iq
          mmol(igcm_co2_ice)=44.
          count=count+1
!          write(*,*) 'co2_ice: count=',count
        endif
        if (noms(iq).eq."h2o_vap") then
          igcm_h2o_vap=iq
          mmol(igcm_h2o_vap)=18.
          count=count+1
!          write(*,*) 'h2o_vap: count=',count
        endif
        if (noms(iq).eq."h2o_ice") then
          igcm_h2o_ice=iq
          mmol(igcm_h2o_ice)=18.
          count=count+1
!          write(*,*) 'h2o_ice: count=',count
        endif
        if (noms(iq).eq."co") then
          igcm_co=iq
          mmol(igcm_co)=28.
          count=count+1
        endif
        if (noms(iq).eq."o") then
          igcm_o=iq
          mmol(igcm_o)=16.
          count=count+1
        endif
        if (noms(iq).eq."o1d") then
          igcm_o1d=iq
          mmol(igcm_o1d)=16.
          count=count+1
        endif
        if (noms(iq).eq."o2") then
          igcm_o2=iq
          mmol(igcm_o2)=32.
          count=count+1
        endif
        if (noms(iq).eq."o3") then
          igcm_o3=iq
          mmol(igcm_o3)=48.
          count=count+1
        endif
        if (noms(iq).eq."h") then
          igcm_h=iq
          mmol(igcm_h)=1.
          count=count+1
        endif
        if (noms(iq).eq."h2") then
          igcm_h2=iq
          mmol(igcm_h2)=2.
          count=count+1
        endif
        if (noms(iq).eq."oh") then
          igcm_oh=iq
          mmol(igcm_oh)=17.
          count=count+1
        endif
        if (noms(iq).eq."ho2") then
          igcm_ho2=iq
          mmol(igcm_ho2)=33.
          count=count+1
        endif
        if (noms(iq).eq."h2o2") then
          igcm_h2o2=iq
          mmol(igcm_h2o2)=34.
          count=count+1
        endif
        if (noms(iq).eq."n2") then
          igcm_n2=iq
          mmol(igcm_n2)=28.
          count=count+1
        endif
        if (noms(iq).eq."ch4") then
          igcm_ch4=iq
          mmol(igcm_ch4)=16.
          count=count+1
        endif
        if (noms(iq).eq."ar") then
          igcm_ar=iq
          mmol(igcm_ar)=40.
          count=count+1
        endif
        if (noms(iq).eq."n") then
          igcm_n=iq
          mmol(igcm_n)=14.
          count=count+1
        endif
        if (noms(iq).eq."no") then
          igcm_no=iq
          mmol(igcm_no)=30.
          count=count+1
        endif
        if (noms(iq).eq."no2") then
          igcm_no2=iq
          mmol(igcm_no2)=46.
          count=count+1
        endif
        if (noms(iq).eq."n2d") then
          igcm_n2d=iq
          mmol(igcm_n2d)=28.
          count=count+1
        endif
        if (noms(iq).eq."ch3") then
          igcm_ch3=iq
          mmol(igcm_ch3)=15.
          count=count+1
        endif
        if (noms(iq).eq."ch") then
          igcm_ch=iq
          mmol(igcm_ch)=13.
          count=count+1
        endif
        if (noms(iq).eq."3ch2") then
          igcm_3ch2=iq
          mmol(igcm_3ch2)=14.
          count=count+1
        endif
        if (noms(iq).eq."1ch2") then
          igcm_1ch2=iq
          mmol(igcm_1ch2)=14.
          count=count+1
        endif
        if (noms(iq).eq."cho") then
          igcm_cho=iq
          mmol(igcm_cho)=29.
          count=count+1
        endif
        if (noms(iq).eq."ch2o") then
          igcm_ch2o=iq
          mmol(igcm_ch2o)=30.
          count=count+1
        endif
        if (noms(iq).eq."ch3o") then
          igcm_ch3o=iq
          mmol(igcm_ch3o)=31.
          count=count+1
        endif
        if (noms(iq).eq."c") then
          igcm_c=iq
          mmol(igcm_c)=12.
          count=count+1
        endif
        if (noms(iq).eq."c2") then
          igcm_c2=iq
          mmol(igcm_c2)=24.
          count=count+1
        endif
        if (noms(iq).eq."c2h") then
          igcm_c2h=iq
          mmol(igcm_c2h)=25.
          count=count+1
        endif
        if (noms(iq).eq."c2h2") then
          igcm_c2h2=iq
          mmol(igcm_c2h2)=26.
          count=count+1
        endif
        if (noms(iq).eq."c2h3") then
          igcm_c2h3=iq
          mmol(igcm_c2h3)=27.
          count=count+1
        endif
        if (noms(iq).eq."c2h4") then
          igcm_c2h4=iq
          mmol(igcm_c2h4)=28.
          count=count+1
        endif
        if (noms(iq).eq."c2h6") then
          igcm_c2h6=iq
          mmol(igcm_c2h6)=30.
          count=count+1
        endif
        if (noms(iq).eq."ch2co") then
          igcm_ch2co=iq
          mmol(igcm_ch2co)=42.
          count=count+1
        endif
        if (noms(iq).eq."ch3co") then
          igcm_ch3co=iq
          mmol(igcm_ch3co)=43.
          count=count+1
        endif
        if (noms(iq).eq."hcaer") then
          igcm_hcaer=iq
          mmol(igcm_hcaer)=50.
          count=count+1
        endif
      enddo ! of do iq=1,nq
      
      ! check that we identified all tracers:
      if (count.ne.nq) then
        write(*,*) "initracer: found only ",count," tracers"
        write(*,*) "               expected ",nq
        do iq=1,count
          write(*,*)'      ',iq,' ',trim(noms(iq))
       enddo
!        stop
      else
        write(*,*) "initracer: found all expected tracers, namely:"
        do iq=1,nq
          write(*,*)'      ',iq,' ',trim(noms(iq))
        enddo
      endif


c------------------------------------------------------------
c     Initialisation tracers ....
c------------------------------------------------------------
      rho_q(1:nq)=0

!      call zerophys(nq,rho_q)
      ! Paladino
      write(*,*) "Volcanish Ash Density?"
      rho_volc=1000
      call getin_p("rho_volc",rho_volc)
      write(*,*) rho_volc ! Volcanic ash density (kg.m-3)

      rho_dust=2500.  ! Mars dust density (kg.m-3)
      rho_ice=920.    ! Water ice density (kg.m-3)
      rho_co2=1620.   ! CO2 ice density (kg.m-3)
      ! rho_volc=700.   ! Volcanic ash density (kg.m-3)
      rho_h2so4=1840. ! H2SO4 density (kg.m-3)
        
        write(*,*) 'rho_q before allocation', rho_q

       do iq=1,nq
        if (noms(iq).eq."volc_1") then
          igcm_volc_1=iq
          radius(igcm_volc_1) = 7.8e-6 ! 7.8 um
          rho_q(igcm_volc_1) = rho_volc
          count=count+1
        write(*,*) 'VOLCANISM : Initialization'
        endif
        
        write(*,*), 'rho_q after volc1', rho_q

        if (noms(iq).eq."volc_2") then
          igcm_volc_2=iq
         radius(igcm_volc_2) = 15.6e-6 ! 15.6 um
         rho_q(igcm_volc_2) = rho_volc
          count=count+1
        endif

        if (noms(iq).eq."volc_3") then
          igcm_volc_3=iq
          rho_q(igcm_volc_3) = rho_volc
          radius(igcm_volc_3) = 125.e-6 ! 125 um
          count=count+1
        endif

        if (noms(iq).eq."volc_4") then
          igcm_volc_4=iq
          rho_q(igcm_volc_4) = rho_volc
          radius(igcm_volc_4) = 1.e-3 ! 1 mm 
          count=count+1
        endif

        if (noms(iq).eq."volc_5") then
          igcm_volc_5=iq
          rho_q(igcm_volc_5) = rho_volc
          radius(igcm_volc_5) = 760.e-6 ! 760 um
          count=count+1
        endif

        if (noms(iq).eq."volc_6") then
          igcm_volc_6=iq
          rho_q(igcm_volc_6) = rho_volc
          radius(igcm_volc_6) = 3045.e-6 ! 3045 um
          count=count+1
        endif
        if (noms(iq).eq."h2so4") then
          igcm_h2so4=iq
          rho_q(igcm_h2so4) = rho_h2so4
          radius(igcm_h2so4) = 1.e-6 ! 1 um
          mmol(igcm_h2so4) = 98.079 ! g/mol
         count=count+1
        endif

       enddo ! of do iq=1,nq



c     Initialization for water vapor
c     ------------------------------
      if(water) then
         radius(igcm_h2o_vap)=0.
         Qext(igcm_h2o_vap)=0.
         alpha_lift(igcm_h2o_vap) =0.
         alpha_devil(igcm_h2o_vap)=0.
	 qextrhor(igcm_h2o_vap)= 0.

         !! this is defined in surfdat_h.F90
         IF (.not.ALLOCATED(dryness)) ALLOCATE(dryness(ngrid))
         IF (.not.ALLOCATED(watercaptag)) ALLOCATE(watercaptag(ngrid))

         do ig=1,ngrid
           if (ngrid.ne.1) watercaptag(ig)=.false.
           dryness(ig) = 1.
         enddo


           radius(igcm_h2o_ice)=3.e-6
           rho_q(igcm_h2o_ice)=rho_ice
           Qext(igcm_h2o_ice)=0.
!           alpha_lift(igcm_h2o_ice) =0.
!           alpha_devil(igcm_h2o_ice)=0.
           qextrhor(igcm_h2o_ice)= (3./4.)*Qext(igcm_h2o_ice) 
     $       / (rho_ice*radius(igcm_h2o_ice))


      end if  ! (water)


!
!     some extra (possibly redundant) sanity checks for tracers:
!     ---------------------------------------------------------
       if (water) then
       ! verify that we indeed have h2o_vap and h2o_ice tracers
         if (igcm_h2o_vap.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use water option without ",
     &                "an h2o_vap tracer !"
           stop
         endif
         if (igcm_h2o_ice.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use water option without ",
     &                "an h2o_ice tracer !"
           stop
         endif
       endif


c     Output for records:
c     ~~~~~~~~~~~~~~~~~~
      write(*,*)
      Write(*,*) '******** initracer : dust transport parameters :'
      write(*,*) 'alpha_lift = ', alpha_lift
      write(*,*) 'alpha_devil = ', alpha_devil
      write(*,*) 'radius  = ', radius
      write(*,*) 'Qext  = ', qext 
      write(*,*)

      end
