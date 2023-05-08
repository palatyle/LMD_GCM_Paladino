      SUBROUTINE initracer(ngrid,nq,qsurf)

       use tracer_mod
       USE comcstfi_h
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
c   Initialize tracer related data in tracer_mod, using tracer names provided
c   by the dynamics in "infotrac"
c
c
c   author: F.Forget
c   ------
c    Modifs: Franck Montmessin, Sebastien Lebonnois (june 2003)
c            Ehouarn Millour (oct. 2008) identify tracers by their names
c=======================================================================


      include "callkeys.h"

      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nq ! number of tracers
      real,intent(out) :: qsurf(ngrid,nq) ! tracer on surface (e.g.  kg.m-2)

      integer iq,ig,count
      real r0_lift , reff_lift, nueff_lift
      real r0_storm,reff_storm
c     Ratio of small over large dust particles (used when both 
c       doubleq and the submicron mode are active); In Montmessin
c       et al. (2002), a value of 25 has been deduced;
      real, parameter :: popratio = 25.
      character(len=30) :: txt ! to store some text

c-----------------------------------------------------------------------
c  radius(nq)      ! aerosol particle radius (m)
c  rho_q(nq)       ! tracer densities (kg.m-3)
c  alpha_lift(nq)  ! saltation vertical flux/horiz flux ratio (m-1)
c  alpha_devil(nq) ! lifting coeeficient by dust devil
c  rho_dust          ! Mars dust density
c  rho_ice           ! Water ice density
c  nuice_ref         ! Effective variance nueff of the
c                    !   water-ice size distributions
c  doubleq           ! if method with mass (iq=1) and number(iq=2) mixing ratio
c  varian            ! Characteristic variance of log-normal distribution
c-----------------------------------------------------------------------


c------------------------------------------------------------
c         NAME and molar mass of the tracer 
c------------------------------------------------------------
   
! Identify tracers by their names: (and set corresponding values of mmol)
      ! 0. initialize tracer indexes to zero:
      igcm_dustbin(1:nq)=0
      igcm_co2_ice=0
      igcm_ccnco2_mass=0
      igcm_ccnco2_number=0
      igcm_dust_mass=0
      igcm_dust_number=0
      igcm_ccn_mass=0
      igcm_ccn_number=0
      igcm_dust_submicron=0
      igcm_h2o_vap=0
      igcm_h2o_ice=0
      igcm_stormdust_mass=0
      igcm_stormdust_number=0
      igcm_topdust_mass=0
      igcm_topdust_number=0
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
      igcm_ch4=0
      igcm_n2=0
      igcm_ar=0
      igcm_ar_n2=0
      igcm_n=0
      igcm_no=0
      igcm_no2=0
      igcm_n2d=0
      igcm_he=0
      igcm_co2plus=0
      igcm_oplus=0
      igcm_o2plus=0
      igcm_coplus=0
      igcm_cplus=0
      igcm_nplus=0
      igcm_noplus=0
      igcm_n2plus=0
      igcm_hplus=0
      igcm_hco2plus=0
      igcm_hcoplus=0
      igcm_elec=0

      ! 1. find dust tracers
      count=0
      if (dustbin.gt.0) then
        do iq=1,nq
          txt=" "
          write(txt,'(a4,i2.2)')'dust',count+1
          if (noms(iq).eq.txt) then
            count=count+1
            igcm_dustbin(count)=iq
            mmol(iq)=100.
          endif
        enddo !do iq=1,nq
      endif ! of if (dustbin.gt.0)
      if (doubleq) then
        do iq=1,nq
          if (noms(iq).eq."dust_mass") then
            igcm_dust_mass=iq
            count=count+1
          endif
          if (noms(iq).eq."dust_number") then
            igcm_dust_number=iq
            count=count+1
          endif
        enddo
      endif ! of if (doubleq)
      if (microphys) then
        do iq=1,nq
          if (noms(iq).eq."ccn_mass") then
            igcm_ccn_mass=iq
            count=count+1
          endif
          if (noms(iq).eq."ccn_number") then
            igcm_ccn_number=iq
            count=count+1
          endif
        enddo
      endif ! of if (microphys)
      if (submicron) then
        do iq=1,nq
          if (noms(iq).eq."dust_submicron") then
            igcm_dust_submicron=iq
            mmol(iq)=100.
            count=count+1
          endif
        enddo
      endif ! of if (submicron)
       if (rdstorm) then
        do iq=1,nq
          if (noms(iq).eq."stormdust_mass") then
            igcm_stormdust_mass=iq
            count=count+1
          endif
          if (noms(iq).eq."stormdust_number") then
            igcm_stormdust_number=iq
            count=count+1
          endif
        enddo
      endif ! of if (rdstorm)
       if (slpwind) then
        do iq=1,nq
          if (noms(iq).eq."topdust_mass") then
            igcm_topdust_mass=iq
            count=count+1
          endif
          if (noms(iq).eq."topdust_number") then
            igcm_topdust_number=iq
            count=count+1
          endif
        enddo
      endif ! of if (slpwind)    
      ! 2. find chemistry and water tracers
      do iq=1,nq
        if (noms(iq).eq."co2") then
          igcm_co2=iq
          mmol(igcm_co2)=44.
          count=count+1
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
        if (noms(iq).eq."he") then
          igcm_he=iq
          mmol(igcm_he)=4.
          count=count+1
        endif
        if (noms(iq).eq."co2plus") then
          igcm_co2plus=iq
          mmol(igcm_co2plus)=44.
          count=count+1
        endif
        if (noms(iq).eq."oplus") then
          igcm_oplus=iq
          mmol(igcm_oplus)=16.
          count=count+1
        endif
        if (noms(iq).eq."o2plus") then
          igcm_o2plus=iq
          mmol(igcm_o2plus)=32.
          count=count+1
        endif
        if (noms(iq).eq."coplus") then
          igcm_coplus=iq
          mmol(igcm_coplus)=28.
          count=count+1
        endif
        if (noms(iq).eq."cplus") then
          igcm_cplus=iq
          mmol(igcm_cplus)=12.
          count=count+1
        endif
        if (noms(iq).eq."nplus") then
          igcm_nplus=iq
          mmol(igcm_nplus)=14.
          count=count+1
        endif
        if (noms(iq).eq."noplus") then
          igcm_noplus=iq
          mmol(igcm_noplus)=30.
          count=count+1
        endif
        if (noms(iq).eq."n2plus") then
          igcm_n2plus=iq
          mmol(igcm_n2plus)=28.
          count=count+1
        endif
        if (noms(iq).eq."hplus") then
          igcm_hplus=iq
          mmol(igcm_hplus)=1.
          count=count+1
        endif
        if (noms(iq).eq."hco2plus") then
          igcm_hco2plus=iq
          mmol(igcm_hco2plus)=45.
          count=count+1
        endif
        if (noms(iq).eq."hcoplus") then
          igcm_hcoplus=iq
          mmol(igcm_hcoplus)=29.
          count=count+1
        endif
        if (noms(iq).eq."elec") then
          igcm_elec=iq
          mmol(igcm_elec)=1./1822.89
          count=count+1
        endif
        if (noms(iq).eq."h2o_vap") then
          igcm_h2o_vap=iq
          mmol(igcm_h2o_vap)=18.
          count=count+1
        endif
        if (noms(iq).eq."co2_ice") then
          igcm_co2_ice=iq
          mmol(igcm_co2_ice)=44.
          count=count+1
        endif
        if (noms(iq).eq."h2o_ice") then
          igcm_h2o_ice=iq
          mmol(igcm_h2o_ice)=18.
          count=count+1
        endif
        ! Other stuff: e.g. for simulations using co2 + neutral gaz
        if (noms(iq).eq."Ar_N2") then
          igcm_ar_n2=iq
          mmol(igcm_ar_n2)=30.
          count=count+1
        endif
        if (co2clouds) then
           if (noms(iq).eq."ccnco2_mass") then
              igcm_ccnco2_mass=iq
              count=count+1
           endif
           if (noms(iq).eq."ccnco2_number") then
              igcm_ccnco2_number=iq
              count=count+1
           endif
        endif
      enddo                     ! of do iq=1,nq
     
      ! check that we identified all tracers:
      if (count.ne.nq) then
        write(*,*) "initracer: found only ",count," tracers"
        write(*,*) "               expected ",nq
        do iq=1,count
          write(*,*)'      ',iq,' ',trim(noms(iq))
        enddo
        stop
      else
        write(*,*) "initracer: found all expected tracers, namely:"
        do iq=1,nq
          write(*,*)'      ',iq,' ',trim(noms(iq))
        enddo
      endif

      ! if water cycle but iceparty=.false., there will nevertheless be
      ! water ice at the surface (iceparty is not used anymore, but this
      ! part is still relevant, as we want to stay compatible with the
      ! older versions).
      if (water.and.(igcm_h2o_ice.eq.0)) then
        igcm_h2o_ice=igcm_h2o_vap ! so that qsurf(i_h2o_ice) is identified
                                  ! even though there is no q(i_h2o_ice)
      else
       ! surface ice qsurf(i_h2o_ice) was loaded twice by phyetat0,
       ! as qsurf(i_h2o_vap) & as qsurf(i_h2o_ice), so to be clean:
       if (igcm_h2o_vap.ne.0) then
         qsurf(1:ngrid,igcm_h2o_vap)=0
       endif
      endif

c------------------------------------------------------------
c     Initialize tracers .... (in tracer_mod)
c------------------------------------------------------------
      ! start by setting everything to (default) zero
      rho_q(1:nq)=0     ! tracer density (kg.m-3)
      radius(1:nq)=0.   ! tracer particle radius (m)
      alpha_lift(1:nq) =0.  ! tracer saltation vertical flux/horiz flux ratio (m-1)
      alpha_devil(1:nq)=0.  ! tracer lifting coefficient by dust devils


      ! some reference values
      rho_dust=2500.  ! Mars dust density (kg.m-3)
      rho_ice=920.    ! Water ice density (kg.m-3)
      rho_ice_co2=1650. 
      !Mangan et al., Icarus 2017 :CO2 density = 1.72391-2.53×10−4T – 2.87×10−6T^2
      nuice_ref=0.1   ! Effective variance nueff of the
                      ! water-ice size distribution
      !!!nuice_sed=0.45   ! Sedimentation effective variance
                      ! of the water-ice size distribution

      if (doubleq) then
c       "doubleq" technique 
c       -------------------
c      (transport of mass and number mixing ratio)
c       iq=1: Q mass mixing ratio, iq=2: N number mixing ratio

        if( (nq.lt.2).or.(water.and.(nq.lt.4)) ) then
          write(*,*)'initracer: nq is too low : nq=', nq
          write(*,*)'water= ',water,' doubleq= ',doubleq   
        end if

        nueff_lift = 0.5
        varian=sqrt(log(1.+nueff_lift))

        rho_q(igcm_dust_mass)=rho_dust
        rho_q(igcm_dust_number)=rho_dust

c       Intermediate calcul for computing geometric mean radius r0
c       as a function of mass and number mixing ratio Q and N
c       (r0 = (r3n_q * Q/ N)^(1/3))
        r3n_q = exp(-4.5*varian**2)*(3./4.)/(pi*rho_dust)

c       Intermediate calcul for computing effective radius reff
c       from geometric mean radius r0
c       (reff = ref_r0 * r0)
        ref_r0 = exp(2.5*varian**2)
        
c       lifted dust :
c       '''''''''''
        reff_lift = 3.0e-6 !3.e-6 !Effective radius of lifted dust (m)
        alpha_devil(igcm_dust_mass)=9.e-9   !  dust devil lift mass coeff
c       alpha_lift(igcm_dust_mass)=3.0e-15  !  Lifted mass coeff

!! default lifting settings
!! -- GCM: alpha_lift not zero because large-scale lifting by default
!! -- MESOSCALE: alpha_lift zero because no lifting at all in mesoscale by default
#ifdef MESOSCALE
        alpha_lift(igcm_dust_mass)=0.0 
#else
        alpha_lift(igcm_dust_mass)=1.e-6 !1.e-6 !Lifted mass coeff
        IF (dustinjection.ge.1) THEN
        	reff_lift = 3.0e-6 ! Effective radius of lifted dust (m)
        	alpha_lift(igcm_dust_mass)=(4/3.)*reff_lift*rho_dust
     &           				/2.4
     	ENDIF
#endif

        r0_lift = reff_lift/ref_r0
        alpha_devil(igcm_dust_number)=r3n_q*
     &                        alpha_devil(igcm_dust_mass)/r0_lift**3
        alpha_lift(igcm_dust_number)=r3n_q*
     &                        alpha_lift(igcm_dust_mass)/r0_lift**3

        radius(igcm_dust_mass) = reff_lift
        radius(igcm_dust_number) = reff_lift

        write(*,*) "initracer: doubleq_param reff_lift:", reff_lift
        write(*,*) "initracer: doubleq_param nueff_lift:", nueff_lift
        write(*,*) "initracer: doubleq_param alpha_lift:",
     &    alpha_lift(igcm_dust_mass)
!c ----------------------------------------------------------------------
!c rocket dust storm scheme
!c lifting tracer stormdust using same distribution than
!c normal dust 
        if (rdstorm) then
          reff_storm=3.e-6 ! reff_lift !3.e-6
          r0_storm=reff_storm/ref_r0
          rho_q(igcm_stormdust_mass)=rho_dust
          rho_q(igcm_stormdust_number)=rho_dust

          alpha_devil(igcm_stormdust_mass)=9.e-9   ! dust devil lift mass coeff
          alpha_lift(igcm_stormdust_mass)=4./3./2.4*reff_storm*rho_dust

          write(*,*) 'alpha_lift(rds):',alpha_lift(igcm_stormdust_mass)
  
          alpha_devil(igcm_stormdust_number)=r3n_q*
     &                      alpha_devil(igcm_stormdust_mass)/r0_storm**3
          alpha_lift(igcm_stormdust_number)=r3n_q*
     &                       alpha_lift(igcm_stormdust_mass)/r0_storm**3
  
          radius(igcm_stormdust_mass) = reff_storm
          radius(igcm_stormdust_number) = reff_storm
        end if !(rdstorm)
!c ----------------------------------------------------------------------
!c slope wind scheme
!c you need a radius value for topdust to active its sedimentation
!c we take the same value as for the normal dust
        if (slpwind) then
          rho_q(igcm_topdust_mass)=rho_dust
          rho_q(igcm_topdust_number)=rho_dust
          radius(igcm_topdust_mass) = 3.e-6
          radius(igcm_topdust_number) = 3.e-6
        end if !(slpwind)
!c ----------------------------------------------------------------------
      
      else

       ! initialize varian, which may be used (e.g. by surfacearea)
       ! even with conrath dust
       nueff_lift = 0.5
       varian=sqrt(log(1.+nueff_lift))

       if (dustbin.gt.1) then
        print*,'initracer: STOP!',
     $   ' properties of dust need to be set in initracer !!!'
        stop

       else if (dustbin.eq.1) then

c       This will be used for 1 dust particle size:
c       ------------------------------------------
        radius(igcm_dustbin(1))=3.e-6
        alpha_lift(igcm_dustbin(1))=0.0e-6
        alpha_devil(igcm_dustbin(1))=7.65e-9
        rho_q(igcm_dustbin(1))=rho_dust

       endif
      end if    ! (doubleq)


c     Scavenging of dust particles by H2O clouds:
c     ------------------------------------------
c     Initialize the two tracers used for the CCNs
      if (water.AND.doubleq.AND.scavenging) then
        radius(igcm_ccn_mass) = radius(igcm_dust_mass)
        alpha_lift(igcm_ccn_mass) = 1e-30
        alpha_devil(igcm_ccn_mass) = 1e-30
        rho_q(igcm_ccn_mass) = rho_dust

        radius(igcm_ccn_number) = radius(igcm_ccn_mass)
        alpha_lift(igcm_ccn_number) = alpha_lift(igcm_ccn_mass)
        alpha_devil(igcm_ccn_number) = alpha_devil(igcm_ccn_mass)
        rho_q(igcm_ccn_number) = rho_q(igcm_ccn_mass)
      endif ! of if (water.AND.doubleq.AND.scavenging)

c     Submicron dust mode:
c     --------------------

      if (submicron) then
        radius(igcm_dust_submicron)=0.1e-6
        rho_q(igcm_dust_submicron)=rho_dust
        if (doubleq) then
c         If doubleq is also active, we use the population ratio:
          alpha_lift(igcm_dust_submicron) = 
     &      alpha_lift(igcm_dust_number)*popratio*
     &      rho_q(igcm_dust_submicron)*4./3.*pi*
     &      radius(igcm_dust_submicron)**3.
          alpha_devil(igcm_dust_submicron)=1.e-30
        else
          alpha_lift(igcm_dust_submicron)=1e-6
          alpha_devil(igcm_dust_submicron)=1.e-30
        endif ! (doubleq)
      end if  ! (submicron)

c     Initialization for water vapor
c     ------------------------------
      if(water) then
         radius(igcm_h2o_vap)=0.
         alpha_lift(igcm_h2o_vap) =0.
         alpha_devil(igcm_h2o_vap)=0.
         if(water.and.(nq.ge.2)) then
           radius(igcm_h2o_ice)=3.e-6
           rho_q(igcm_h2o_ice)=rho_ice
           alpha_lift(igcm_h2o_ice) =0.
           alpha_devil(igcm_h2o_ice)=0.
         elseif(water.and.(nq.lt.2)) then
            write(*,*) 'nq is too low : nq=', nq
            write(*,*) 'water= ',water
         endif

      end if  ! (water)
      
! Initialisation for CO2 clouds
      if (co2clouds ) then 
        radius(igcm_ccnco2_mass) = radius(igcm_dust_mass)
        alpha_lift(igcm_ccnco2_mass) = 1e-30
        alpha_devil(igcm_ccnco2_mass) = 1e-30
        rho_q(igcm_ccnco2_mass) = rho_dust
        radius(igcm_ccnco2_number) = radius(igcm_ccnco2_mass)
        alpha_lift(igcm_ccnco2_number) = alpha_lift(igcm_ccnco2_mass)
        alpha_devil(igcm_ccnco2_number) = alpha_devil(igcm_ccnco2_mass)
        rho_q(igcm_ccnco2_number) = rho_q(igcm_ccnco2_mass)
     
        radius(igcm_co2)=0.
        alpha_lift(igcm_co2) =0.
        alpha_devil(igcm_co2)=0.
        radius(igcm_co2_ice)=1.e-8
        rho_q(igcm_co2_ice)=rho_ice_co2
        alpha_lift(igcm_co2_ice) =0.
        alpha_devil(igcm_co2_ice)=0.

      endif 
      
c     Output for records:
c     ~~~~~~~~~~~~~~~~~~
      write(*,*)
      Write(*,*) '******** initracer : dust transport parameters :'
      write(*,*) 'alpha_lift = ', alpha_lift
      write(*,*) 'alpha_devil = ', alpha_devil
      write(*,*) 'radius  = ', radius
      if(doubleq) then
        write(*,*) 'reff_lift (um) =  ', reff_lift
        write(*,*) 'size distribution variance  = ', varian
        write(*,*) 'r3n_q , ref_r0 : ', r3n_q , ref_r0
      end if

!
!     some extra (possibly redundant) sanity checks for tracers: 
!     ---------------------------------------------------------

       if (doubleq) then 
       ! verify that we indeed have dust_mass and dust_number tracers 
         if (igcm_dust_mass.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use doubleq option without ",
     &                "a dust_mass tracer !"
           stop
         endif
         if (igcm_dust_number.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use doubleq option without ",
     &                "a dust_number tracer !"
           stop
         endif
       endif

       if ((.not.doubleq).and.(dustbin.gt.0)) then
       ! verify that we indeed have 'dustbin' dust tracers
         count=0
         do iq=1,dustbin
           if (igcm_dustbin(iq).ne.0) then
             count=count+1
           endif
         enddo
         if (count.ne.dustbin) then
           write(*,*) "initracer: error !!"
           write(*,*) "  dusbin is set to ",dustbin,
     &                " but we only have the following dust tracers:"
           do iq=1,count
             write(*,*)"   ",trim(noms(igcm_dustbin(iq)))
           enddo
           stop
         endif
       endif

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

       if (co2clouds) then
          !verify that we have co2_ice and co2 tracers
          if (igcm_co2 .eq. 0) then 
             write(*,*) "initracer: error !!"
             write(*,*) "  cannot use co2 clouds option without ",
     &            "a co2 tracer !"
          stop
          endif
          if (igcm_co2_ice .eq. 0) then 
             write(*,*) "initracer: error !!"
             write(*,*) "  cannot use co2 clouds option without ",
     &            "a co2_ice tracer !"
             stop
          endif 
       endif
 
       if (rdstorm) then 
       ! verify that we indeed have stormdust_mass and stormdust_number tracers 
         if (igcm_stormdust_mass.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use rdstorm option without ",
     &                "a stormdust_mass tracer !"
           stop
         endif
         if (igcm_stormdust_number.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use rdstorm option without ",
     &                "a stormdust_number tracer !"
           stop
         endif
       endif

       if (slpwind) then 
       ! verify that we indeed have topdust_mass and topdust_number tracers 
         if (igcm_topdust_mass.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use slpwind option without ",
     &                "a topdust_mass tracer !"
           stop
         endif
         if (igcm_topdust_number.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use slpwind option without ",
     &                "a topdust_number tracer !"
           stop
         endif
       endif
      
       if (callnlte) then ! NLTE requirements
         if (nltemodel.ge.1) then
           ! check that co2, co, o and n2 tracers are available
           if (igcm_co2.eq.0) then
             write(*,*) "initracer: error !!"
             write(*,*) "  with nltemodel>0, we need the co2 tracer!"
             stop
           endif
           if (igcm_co.eq.0) then
             write(*,*) "initracer: error !!"
             write(*,*) "  with nltemodel>0, we need the co tracer!"
             stop
           endif
           if (igcm_o.eq.0) then
             write(*,*) "initracer: error !!"
             write(*,*) "  with nltemodel>0, we need the o tracer!"
             stop
           endif
           if (igcm_n2.eq.0) then
             write(*,*) "initracer: error !!"
             write(*,*) "  with nltemodel>0, we need the n2 tracer!"
             stop
           endif
         endif 
       endif

       if (scavenging) then
       ! verify that we indeed have ccn_mass and ccn_number tracers
         if (igcm_ccn_mass.eq.0 .and. igcm_ccnco2_mass.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use scavenging option without ",
     &                "a ccn_mass or ccnco2_mass tracer !"
           stop
         endif
         if (igcm_ccn_number.eq.0 .and. igcm_ccnco2_number.eq.0 ) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use scavenging option without ",
     &                "a ccn_number or ccnco2_number tracer !"
           stop
         endif
       endif ! of if (scavenging)

       if (photochem .or. callthermos) then
       ! verify that we indeed have the chemistry tracers
         if (igcm_co2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a co2 tracer !"
         stop
         endif
         if (igcm_co.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a co tracer !"
         stop
         endif
         if (igcm_o.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a o tracer !"
         stop
         endif
         if (igcm_o1d.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a o1d tracer !"
         stop
         endif
         if (igcm_o2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "an o2 tracer !"
         stop
         endif
         if (igcm_o3.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "an o3 tracer !"
         stop
         endif
         if (igcm_h.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a h tracer !"
         stop
         endif
         if (igcm_h2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a h2 tracer !"
         stop
         endif
         if (igcm_oh.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "an oh tracer !"
         stop
         endif
         if (igcm_ho2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a ho2 tracer !"
         stop
      endif
         if (igcm_h2o2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a h2o2 tracer !"
         stop
         endif
         if (igcm_n2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a n2 tracer !"
         stop
         endif
         if (igcm_ar.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "an ar tracer !"
         stop
         endif
       endif ! of if (photochem .or. callthermos)

      end
