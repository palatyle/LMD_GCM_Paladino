      SUBROUTINE initracer(qsurf,co2ice)

       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Initialisation related to tracer 
c (transported dust, water, chemical species, ice...)
c
c   Name of the tracer
c
c   Test of dimension :
c   Initialise COMMON tracer in tracer.h
c
c   If water=T : q(iq=nqmx) is the water mass mixing ratio
c   If water=T and iceparty=T : q(iq=nqmx-1) is the ice mass mixing ratio

c   If there is transported dust, it uses iq=1 to iq=dustbin
c   If there is no transported dust : dustbin=0
c   If doubleq=T : q(iq=1) is the dust mass mixing ratio
c                  q(iq=2) is the dust number mixing ratio 

c   If (photochem.or.thermochem) there is "ncomp" chemical species (ncomp
c   is set in aeronomars/chimiedata.h) using the ncomp iq values starting at
c      iq=nqchem_min = dustbin+1   (nqchem_min is defined in inifis.F)
c   
c
c   author: F.Forget
c   ------
c    Modifs: Franck Montmessin, Sebastien Lebonnois (june 2003)
c
c=======================================================================


#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "tracer.h"

#include "comgeomfi.h"
#include "watercap.h"
#include "aerice.h"
#include "fisice.h"
#include "chimiedata.h"


      real qsurf(ngridmx,nqmx)       ! tracer on surface (e.g.  kg.m-2)
      real co2ice(ngridmx)           ! co2 ice mass on surface (e.g.  kg.m-2)
      integer iq,ig
      real r0_lift , reff_lift

c-----------------------------------------------------------------------
c  radius(nqmx)      ! aerosol particle radius (m)
c  qext(nqmx)        ! Single Scat. Extinction coeff at 0.67 um
c  alpha_lift(nqmx)  ! saltation vertical flux/horiz flux ratio (m-1)
c  alpha_devil(nqmx) ! lifting coeeficient by dust devil
c  rho_dust          ! Mars dust density
c  rho_ice           ! Water ice density
c  doubleq           ! if method with mass (iq=1) and number(iq=2) mixing ratio
c  varian            ! Characteristic variance of log-normal distribution
c-----------------------------------------------------------------------


c------------------------------------------------------------
c     Test Dimensions tracers
c------------------------------------------------------------

      if(photochem.or.thermochem) then
          if (iceparty) then
              if ((nqchem_min+ncomp+1).ne.nqmx) then
                 print*,'********* Dimension problem! ********'
                 print*,"nqchem_min+ncomp+1).ne.nqmx"
                 print*,"ncomp: ",ncomp
                 print*,"nqchem_min: ",nqchem_min
                 print*,"nqmx: ",nqmx
                 print*,'Change ncomp in chimiedata.h'
               endif
          else
               if ((nqchem_min+ncomp).ne.nqmx) then
                 print*,'********* Dimension problem! ********'
                 print*,"nqchem_min+ncomp).ne.nqmx"
                 print*,"ncomp: ",ncomp
                 print*,"nqchem_min: ",nqchem_min
                 print*,"nqmx: ",nqmx
                 print*,'Change ncomp in chimiedata.h'
                 STOP
               endif
            endif
          endif

c------------------------------------------------------------
c         NAME and molar mass of the tracer 
c------------------------------------------------------------

c noms and mmol vectors:
      if (water) then 
         mmol(nqmx) = 18.
         noms(nqmx)   = 'h2o'
      end if
      if (iceparty) then
         noms(nqmx-1) = 'ice'
         mmol(nqmx-1) = 18.
      end if
      if(photochem.or.thermochem) then
         do iq=nqchem_min, nqchem_min+ncomp-1
             noms(iq) = nomchem(iq-nqchem_min+1)
             mmol(iq) = mmolchem(iq-nqchem_min+1)
         enddo
      end if
      if (dustbin.ge.1) then
          do iq=1,dustbin
            noms(iq) = 'dust'
            mmol(iq) = 100.
          enddo
          if (doubleq) then
           noms(1) = 'dust mass mix. ratio'
           noms(dustbin) = 'dust number mix. ratio'
          end if
      end if

c     Simulation of CO2 + neutral gaz
      if ((dustbin.eq.0).and.(.not.water)) then
         noms(1) = 'co2'
         mmol(1) = 44
         if (nqmx.eq.2)then
            noms(nqmx) = 'Ar_N2'
            mmol(nqmx) = 30
         end if 
      end if 
   

c------------------------------------------------------------
c     Initialisation tracers ....
c------------------------------------------------------------

      rho_dust=2500.  ! Mars dust density (kg.m-3)
      rho_ice=920.    ! Water ice density (kg.m-3)

      if (doubleq) then
c       "doubleq" technique 
c       -------------------
c      (transport of mass and number mixing ratio)
c       iq=1: Q mass mixing ratio, iq=2: N number mixing ratio

        if( (nqmx.lt.2).or.(water.and.(nqmx.lt.3)) ) then
            write(*,*) 'nqmx is too low : nqmx=', nqmx
            write(*,*) 'water= ',water,' doubleq= ',doubleq   
        end if

        varian=0.637    ! Characteristic variance    
        qext(1)=3.04    ! reference extinction at 0.67 um for ref dust
        qext(2)=3.04    ! reference extinction at 0.67 um for ref dust

c       Intermediate calcul for computing geometric mean radius r0
c       as a function of mass and number mixing ratio Q and N
c       (r0 = (r3n_q * Q/ N)
        r3n_q = exp(-4.5*varian**2)*(3./4.)/(pi*rho_dust)

c       Intermediate calcul for computing effective radius reff
c       from geometric mean radius r0
c       (reff = ref_r0 * r0)
        ref_r0 = exp(2.5*varian**2)
        
c       lifted dust :
c       '''''''''''
        reff_lift = 3.e-6      !  Effective radius of lifted dust (m)
        alpha_devil(1)=9.e-9   !  dust devil lift mass coeff
        alpha_lift(1)=3.0e-15  !  Lifted mass coeff

        r0_lift = reff_lift/ref_r0
        alpha_devil(2)= r3n_q * alpha_devil(1)/r0_lift**3
        alpha_lift(2)= r3n_q * alpha_lift(1)/r0_lift**3

c       Not used:
        radius(1) = 0.
	radius(2) = 0.

      else

       if (dustbin.gt.1) then
        print*,'ATTENTION:',
     $   ' properties of dust need input in initracer !!!'
	stop

       else if (dustbin.eq.1) then

c       This will be used for 1 dust particle size:
c       ------------------------------------------
cc        radius(1)=3.e-6
cc        Qext(1)=3.04
cc        alpha_lift(1)=0.0e-6
cc        alpha_devil(1)=7.65e-9
        radius(1)=1.e-6
        Qext(1)=3.04        
        alpha_lift(1)=1.
        alpha_devil(1)=0.0e-9        
	qextrhor(1)= (3./4.)*Qext(1) / (rho_dust*radius(1))

       endif
      end if    ! (doubleq)

c     Initialization for photochemistry:
c     ---------------------------------
      if (photochem) then
        do iq=nqchem_min,nqmx
	 radius(iq)=0.
         Qext(iq)=0.
         alpha_lift(iq) =0.
         alpha_devil(iq)=0.
	 qextrhor(iq)= 0.
        enddo
      endif

c     Initialization for water vapor
c     ------------------------------
      if(water) then
         radius(nqmx)=0.
         Qext(nqmx)=0.
         alpha_lift(nqmx) =0.
         alpha_devil(nqmx)=0.
	 qextrhor(nqmx)= 0.

c       "Dryness coefficient" controlling the evaporation and
c        sublimation from the ground water ice (close to 1)
c        HERE, the goal is to correct for the fact
c        that the simulated permanent water ice polar caps
c        is larger than the actual cap and the atmospheric
c        opacity not always realistic.

         do ig=1,ngridmx
           if (ngridmx.ne.1) watercaptag(ig)=.false.
           dryness(ig) = 1.
           if (activice) pclc(ig)=1.
         enddo

         IF (caps) THEN
c Perennial H20 north cap defined by watercaptag=true (allows surface to be
c hollowed by sublimation in vdifc).
c Cloud area fraction (pclc) is defined here.
         do ig=1,ngridmx
           if (lati(ig)*180./pi.gt.84) then
             if (ngridmx.ne.1) watercaptag(ig)=.true.
             dryness(ig) = 1.
             if (activice)then
               pclc(ig)=1.
               print*,'Cloud area ratio : ',pclc(ig),' at lat '
     $        ,lati(ig)*180./pi
             endif
c Use the following cap definition for high spatial resolution (latitudinal bin <= 5 deg)
c	      if (lati(ig)*180./pi.lt.85.and.long(ig).ge.0) then
c	        if (ngridmx.ne.1) watercaptag(ig)=.true.
c	        pclc(ig)=.3
c               dryness(ig) = 1.
c               print*,'Cloud area ratio : ',pclc(ig),' at lat '
c	      endif
c             if (lati(ig)*180./pi.ge.85) then
c               if (ngridmx.ne.1) watercaptag(ig)=.true.
c               dryness(ig) = 1.
c	        pclc(ig)=.3
c               print*,'Cloud area ratio : ',pclc(ig),' at lat '
c	      endif
           endif  ! (lati>80 deg)
         end do ! (ngridmx)
        ENDIF ! (caps)

         if(iceparty.and.nqmx.ge.2) then
           radius(nqmx-1)=3.e-6
           Qext(nqmx-1)=0.
           alpha_lift(nqmx-1) =0.
           alpha_devil(nqmx-1)=0.
           if (activice) then
             radius(nqmx-1)=rcrystal
             Qext(nqmx-1)=qrefice
           endif
           qextrhor(nqmx-1)= (3./4.)*Qext(nqmx-1) 
     $       / (rho_ice*radius(nqmx-1))
         elseif(iceparty.and.nqmx.lt.2) then
            write(*,*) 'nqmx is too low : nqmx=', nqmx
            write(*,*) 'water= ',water,' iceparty= ',iceparty
         endif

      end if  ! (water)

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
      write(*,*) 'Qext  = ', qext 
      write(*,*)

      end
