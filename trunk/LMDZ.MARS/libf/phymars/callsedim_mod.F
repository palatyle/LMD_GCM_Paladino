      MODULE callsedim_mod

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE callsedim(ngrid,nlay,ptimestep,
     &                pplev,zlev,zlay,pt,pdt,
     &                rdust,rstormdust,rtopdust,
     &                rice,rsedcloud,rhocloud,
     &                pq,pdqfi,pdqsed,pdqs_sed,nq, 
     &                tau,tauscaling)

      USE ioipsl_getincom, only: getin
      USE updaterad, only: updaterdust,updaterice_micro,updaterice_typ
      USE tracer_mod, only: noms, igcm_dust_mass, igcm_dust_number,
     &                      rho_dust, rho_q, radius, varian,
     &                      igcm_ccn_mass, igcm_ccn_number,
     &                      igcm_h2o_ice, nuice_sed, nuice_ref,
     &                      igcm_ccnco2_mass,igcm_ccnco2_number,
     &                      igcm_co2_ice, igcm_stormdust_mass, 
     &                      igcm_stormdust_number,igcm_topdust_mass, 
     &                      igcm_topdust_number
      USE newsedim_mod, ONLY: newsedim
      USE comcstfi_h, ONLY: g
      USE dimradmars_mod, only: naerkind
      IMPLICIT NONE

c=======================================================================
c      Sedimentation of the  Martian aerosols
c      depending on their density and radius
c
c      F.Forget 1999
c
c      Modified by J.-B. Madeleine 2010: Now includes the doubleq
c        technique in order to have only one call to callsedim in
c        physiq.F.
c
c      Modified by J. Audouard 09/16: Now includes the co2clouds case
c        If the co2 microphysics is on, then co2 theice & ccn tracers 
c        are being sedimented in the microtimestep (co2cloud.F), not 
c        in this routine.
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------
      
      include "callkeys.h"

c
c   arguments:
c   ----------

      integer,intent(in) :: ngrid  ! number of horizontal grid points
      integer,intent(in) :: nlay   ! number of atmospheric layers
      real,intent(in) :: ptimestep ! physics time step (s)
      real,intent(in) :: pplev(ngrid,nlay+1) ! pressure at inter-layers (Pa)
      real,intent(in) :: zlev(ngrid,nlay+1) ! altitude at layer boundaries
      real,intent(in) :: zlay(ngrid,nlay)   ! altitude at the middle of the layers
      real,intent(in) :: pt(ngrid,nlay) ! temperature at mid-layer (K)
      real,intent(in) :: pdt(ngrid,nlay) ! tendency on temperature, from
                                         ! previous processes (K/s)
c    Aerosol radius provided by the water ice microphysical scheme:
      real,intent(out) :: rdust(ngrid,nlay) ! Dust geometric mean radius (m)
      real,intent(out) :: rstormdust(ngrid,nlay) ! Stormdust geometric mean radius (m) 
      real,intent(out) :: rtopdust(ngrid,nlay) ! topdust geometric mean radius (m)
      real,intent(out) :: rice(ngrid,nlay)  ! H2O Ice geometric mean radius (m)
c     Sedimentation radius of water ice
      real,intent(in) :: rsedcloud(ngrid,nlay)
c     Cloud density (kg.m-3)
      real,intent(inout) :: rhocloud(ngrid,nlay)
c    Traceurs :
      real,intent(in) :: pq(ngrid,nlay,nq)  ! tracers (kg/kg)
      real,intent(in) :: pdqfi(ngrid,nlay,nq)  ! tendency before sedimentation (kg/kg.s-1)
      real,intent(out) :: pdqsed(ngrid,nlay,nq) ! tendency due to sedimentation (kg/kg.s-1)
      real,intent(out) :: pdqs_sed(ngrid,nq)    ! flux at surface (kg.m-2.s-1)
      integer,intent(in) :: nq  ! number of tracers
      real,intent(in) :: tau(ngrid,naerkind) ! dust opacity
      real,intent(in) :: tauscaling(ngrid)
      
c   local:
c   ------

      INTEGER l,ig, iq

      real zqi(ngrid,nlay,nq) ! to locally store tracers
      real zt(ngrid,nlay) ! to locally store temperature
      real masse (ngrid,nlay) ! Layer mass (kg.m-2)
      real epaisseur (ngrid,nlay) ! Layer thickness (m)
      real wq(ngrid,nlay+1) ! displaced tracer mass (kg.m-2)
      real r0(ngrid,nlay) ! geometric mean radius used for
                                !   sedimentation (m)
      real r0dust(ngrid,nlay) ! geometric mean radius used for
                                    !   dust (m)
      real r0stormdust(ngrid,nlay) ! Geometric mean radius used for stormdust (m)
!                                    !   CCNs (m)
      real r0topdust(ngrid,nlay) ! Geometric mean radius used for topdust (m)
!                                    !   CCNs (m)
      real,save :: beta ! correction for the shape of the ice particles (cf. newsedim)
c     for ice radius computation
      REAL Mo,No
      REAl ccntyp



c     Discrete size distributions (doubleq)
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       1) Parameters used to represent the changes in fall
c          velocity as a function of particle size;
      integer ir
      integer,parameter :: nr=12 !(nr=7) ! number of bins
      real,save :: rd(nr)
      real qr(ngrid,nlay,nr)
      real,save :: rdi(nr+1)    ! extreme and intermediate radii
      real Sq(ngrid,nlay)
      real,parameter :: rdmin=1.e-8
      real,parameter :: rdmax=30.e-6
      real,parameter :: rdimin=1.e-8 ! 1.e-7
      real,parameter :: rdimax=1.e-4

c       2) Second size distribution for the log-normal integration
c          (the mass mixing ratio is computed for each radius)

      integer iint
      integer,parameter :: ninter=4 ! number of points between each rdi radii
      real,save :: rr(ninter,nr)
      integer radpower
      real sigma0

c       3) Other local variables used in doubleq

      INTEGER,SAVE :: idust_mass  ! index of tracer containing dust mass
                                  !   mix. ratio
      INTEGER,SAVE :: idust_number ! index of tracer containing dust number
                                   !   mix. ratio
      INTEGER,SAVE :: iccn_mass  ! index of tracer containing CCN mass
                                 !   mix. ratio
      INTEGER,SAVE :: iccn_number ! index of tracer containing CCN number
                                  !   mix. ratio
      INTEGER,SAVE :: istormdust_mass  !  index of tracer containing
                                       !stormdust mass mix. ratio
      INTEGER,SAVE :: istormdust_number !  index of tracer containing
                                        !stormdust number mix. ratio 
      INTEGER,SAVE :: itopdust_mass  !  index of tracer containing
                                       !topdust mass mix. ratio
      INTEGER,SAVE :: itopdust_number !  index of tracer containing
                                        !topdust number mix. ratio                        
      INTEGER,SAVE :: iccnco2_number ! index of tracer containing CCN number
      INTEGER,SAVE :: iccnco2_mass ! index of tracer containing CCN number
      INTEGER,SAVE :: ico2_ice ! index of tracer containing CCN number


      LOGICAL,SAVE :: firstcall=.true.



c    ** un petit test de coherence
c       --------------------------
      ! AS: firstcall OK absolute
      IF (firstcall) THEN
         
c       Doubleq: initialization
        IF (doubleq) THEN
         do ir=1,nr
             rd(ir)= rdmin*(rdmax/rdmin)**(float(ir-1)/float(nr-1))
         end do
         rdi(1)=rdimin
         do ir=2,nr
           rdi(ir)= sqrt(rd(ir-1)*rd(ir))
         end do
         rdi(nr+1)=rdimax

         do ir=1,nr
           do iint=1,ninter
             rr(iint,ir)=
     &        rdi(ir)*
     &        (rdi(ir+1)/rdi(ir))**(float(iint-1)/float(ninter-1))
c             write(*,*) rr(iint,ir)
           end do
         end do

      ! identify tracers corresponding to mass mixing ratio and
      ! number mixing ratio
        idust_mass=0      ! dummy initialization
        idust_number=0    ! dummy initialization

        do iq=1,nq
          if (noms(iq).eq."dust_mass") then
            idust_mass=iq
            write(*,*)"callsedim: idust_mass=",idust_mass
          endif
          if (noms(iq).eq."dust_number") then
            idust_number=iq
            write(*,*)"callsedim: idust_number=",idust_number
          endif
        enddo

        ! check that we did find the tracers
        if ((idust_mass.eq.0).or.(idust_number.eq.0)) then
          write(*,*) 'callsedim: error! could not identify'
          write(*,*) ' tracers for dust mass and number mixing'
          write(*,*) ' ratio and doubleq is activated!'
          stop
        endif
        ENDIF !of if (doubleq)

        IF (microphys) THEN
          iccn_mass=0
          iccn_number=0
          do iq=1,nq
            if (noms(iq).eq."ccn_mass") then
              iccn_mass=iq
              write(*,*)"callsedim: iccn_mass=",iccn_mass
            endif
            if (noms(iq).eq."ccn_number") then
              iccn_number=iq
              write(*,*)"callsedim: iccn_number=",iccn_number
            endif
          enddo
          ! check that we did find the tracers
          if ((iccn_mass.eq.0).or.(iccn_number.eq.0)) then
            write(*,*) 'callsedim: error! could not identify'
            write(*,*) ' tracers for ccn mass and number mixing'
            write(*,*) ' ratio and microphys is activated!'
            stop
          endif
        ENDIF !of if (microphys)

        IF (co2clouds) THEN
          iccnco2_mass=0
          iccnco2_number=0
          ico2_ice=0
          do iq=1,nq
            if (noms(iq).eq."ccnco2_mass") then
              iccnco2_mass=iq
              write(*,*)"callsedim: iccnco2_mass=",iccnco2_mass
            endif
            if (noms(iq).eq."co2_ice") then
              ico2_ice=iq
              write(*,*)"callsedim: ico2_ice=",ico2_ice
            endif
            if (noms(iq).eq."ccnco2_number") then
              iccnco2_number=iq
              write(*,*)"callsedim: iccnco2_number=",iccnco2_number
            endif
          enddo
          ! check that we did find the tracers
          if ((iccnco2_mass.eq.0).or.(iccnco2_number.eq.0)) then
            write(*,*) 'callsedim: error! could not identify'
            write(*,*) ' tracers for ccn co2 mass and number mixing'
            write(*,*) ' ratio and co2clouds are activated!'
            stop
          endif
       ENDIF                    !of if (co2clouds)

       IF (water) THEN
         write(*,*) "correction for the shape of the ice particles ?"
         beta=0.75 ! default value
         call getin("ice_shape",beta)
         write(*,*) " ice_shape = ",beta

          write(*,*) "water_param nueff Sedimentation:", nuice_sed
          IF (activice) THEN
            write(*,*) "water_param nueff Radiative:", nuice_ref
          ENDIF
       ENDIF

       IF (rdstorm) THEN ! identifying stormdust tracers for sedimentation
           istormdust_mass=0      ! dummy initialization
           istormdust_number=0    ! dummy initialization

           do iq=1,nq
             if (noms(iq).eq."stormdust_mass") then
               istormdust_mass=iq
               write(*,*)"callsedim: istormdust_mass=",istormdust_mass
             endif
             if (noms(iq).eq."stormdust_number") then
               istormdust_number=iq
               write(*,*)"callsedim: istormdust_number=",
     &                                           istormdust_number
             endif
           enddo

           ! check that we did find the tracers
           if ((istormdust_mass.eq.0).or.(istormdust_number.eq.0)) then
             write(*,*) 'callsedim: error! could not identify'
             write(*,*) ' tracers for stormdust mass and number mixing'
             write(*,*) ' ratio and rdstorm is activated!'
             stop
           endif
       ENDIF !of if (rdstorm)

       IF (slpwind) THEN ! identifying topdust tracers for sedimentation
           itopdust_mass=0      ! dummy initialization
           itopdust_number=0    ! dummy initialization

           do iq=1,nq
             if (noms(iq).eq."topdust_mass") then
               itopdust_mass=iq
               write(*,*)"callsedim: itopdust_mass=",itopdust_mass
             endif
             if (noms(iq).eq."topdust_number") then
               itopdust_number=iq
               write(*,*)"callsedim: itopdust_number=",
     &                                           itopdust_number
             endif
           enddo

           ! check that we did find the tracers
           if ((itopdust_mass.eq.0).or.(itopdust_number.eq.0)) then
             write(*,*) 'callsedim: error! could not identify'
             write(*,*) ' tracers for topdust mass and number mixing'
             write(*,*) ' ratio and slpwind is activated!'
             stop
           endif
       ENDIF !of if (slpwind)

        firstcall=.false.
      ENDIF ! of IF (firstcall)

c-----------------------------------------------------------------------
c    1. Initialization
c    -----------------

!      zqi(1:ngrid,1:nlay,1:nqmx) = 0.
c     Update the mass mixing ratio and temperature with the tendencies coming
c       from other parameterizations:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      zqi(1:ngrid,1:nlay,1:nq)=pq(1:ngrid,1:nlay,1:nq)
     &                         +pdqfi(1:ngrid,1:nlay,1:nq)*ptimestep
      zt(1:ngrid,1:nlay)=pt(1:ngrid,1:nlay)
     &                         +pdt(1:ngrid,1:nlay)*ptimestep

c    Computing the different layer properties
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c    Mass (kg.m-2), thickness(m), crossing time (s)  etc.

      do  l=1,nlay
        do ig=1, ngrid
          masse(ig,l)=(pplev(ig,l) - pplev(ig,l+1)) /g 
          epaisseur(ig,l)= zlev(ig,l+1) - zlev(ig,l)
        end do
      end do

c =================================================================
c     Compute the geometric mean radius used for sedimentation

      if (doubleq) then
        do l=1,nlay
          do ig=1, ngrid
     
         call updaterdust(zqi(ig,l,igcm_dust_mass),
     &                    zqi(ig,l,igcm_dust_number),r0dust(ig,l),
     &                    tauscaling(ig))
          
          end do
        end do
      endif
c     rocket dust storm
      if (rdstorm) then
        do l=1,nlay
          do ig=1, ngrid
     
         call updaterdust(zqi(ig,l,igcm_stormdust_mass),
     &               zqi(ig,l,igcm_stormdust_number),r0stormdust(ig,l),
     &               tauscaling(ig))
          
          end do
        end do
      endif
c     entrainment by slope wind
      if (slpwind) then
        do l=1,nlay
          do ig=1, ngrid
     
         call updaterdust(zqi(ig,l,igcm_topdust_mass),
     &               zqi(ig,l,igcm_topdust_number),r0topdust(ig,l),
     &               tauscaling(ig))
          
          end do
        end do
      endif
c =================================================================
      do iq=1,nq
        if(radius(iq).gt.1.e-9 .and.(iq.ne.ico2_ice) .and.
     &        (iq .ne. iccnco2_mass) .and. (iq .ne. 
     &        iccnco2_number)) then   ! no sedim for gaz or CO2 clouds  (done in microtimestep)

c -----------------------------------------------------------------
c         DOUBLEQ CASE
c -----------------------------------------------------------------
          if ( doubleq.and.
     &     ((iq.eq.idust_mass).or.(iq.eq.idust_number).or.
     &     (iq.eq.istormdust_mass).or.(iq.eq.istormdust_number).or.
     &     (iq.eq.itopdust_mass).or.(iq.eq.itopdust_number)) ) then
     
c           Computing size distribution:
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if ((iq.eq.idust_mass).or.(iq.eq.idust_number)) then
              do  l=1,nlay
                do ig=1, ngrid
                  r0(ig,l)=r0dust(ig,l)
                end do
              end do
            else if ((iq.eq.istormdust_mass).or.
     &                                (iq.eq.istormdust_number)) then
              do  l=1,nlay
                do ig=1, ngrid
                  r0(ig,l)=r0stormdust(ig,l)
                end do
              end do
            else if ((iq.eq.itopdust_mass).or.
     &                                (iq.eq.itopdust_number)) then
              do  l=1,nlay
                do ig=1, ngrid
                  r0(ig,l)=r0topdust(ig,l)
                end do
              end do
            endif
            sigma0 = varian

c        Computing mass mixing ratio for each particle size
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF ((iq.EQ.idust_mass).or.(iq.EQ.istormdust_mass)
     &                          .or.(iq.EQ.itopdust_mass)) then
            radpower = 2
          ELSE  ! number
            radpower = -1
          ENDIF
          Sq(1:ngrid,1:nlay) = 0.
          do ir=1,nr
            do l=1,nlay
              do ig=1,ngrid
c                ****************
c                Size distribution integration
c                (Trapezoid Integration Method)
                 qr(ig,l,ir)=0.5*(rr(2,ir)-rr(1,ir))*
     &             (rr(1,ir)**radpower)*
     &             exp(-(log(rr(1,ir)/r0(ig,l)))**2/(2*sigma0**2))
                 do iint=2,ninter-1
                   qr(ig,l,ir)=qr(ig,l,ir) +
     &             0.5*(rr(iint+1,ir)-rr(iint-1,ir))*
     &             (rr(iint,ir)**radpower)*
     &             exp(-(log(rr(iint,ir)/r0(ig,l)))**2/
     &             (2*sigma0**2))
                 end do
                 qr(ig,l,ir)=qr(ig,l,ir) +
     &             0.5*(rr(ninter,ir)-rr(ninter-1,ir))*
     &             (rr(ninter,ir)**radpower)*
     &             exp(-(log(rr(ninter,ir)/r0(ig,l)))**2/
     &             (2*sigma0**2))

c                **************** old method (not recommended!)
c                qr(ig,l,ir)=(rd(ir)**(5-3*iq))*
c    &           exp( -(log(rd(ir)/r0(ig,l)))**2 / (2*sigma0**2) )
c                ******************************

                 Sq(ig,l)=Sq(ig,l)+qr(ig,l,ir)
              enddo
            enddo
          enddo

          do ir=1,nr
            do l=1,nlay
              do ig=1,ngrid
                 qr(ig,l,ir) = zqi(ig,l,iq)*qr(ig,l,ir)/Sq(ig,l)
              enddo
            enddo
          enddo

c         Computing sedimentation for each tracer
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          zqi(1:ngrid,1:nlay,iq) = 0.
          pdqs_sed(1:ngrid,iq) = 0.

          do ir=1,nr
               call newsedim(ngrid,nlay,1,1,ptimestep,
     &         pplev,masse,epaisseur,zt,rd(ir),(/rho_dust/),qr(1,1,ir),
     &         wq,0.5)

c            Tendencies
c            ~~~~~~~~~~
             do ig=1,ngrid
               pdqs_sed(ig,iq) = pdqs_sed(ig,iq)
     &                                + wq(ig,1)/ptimestep
             end do
             DO l = 1, nlay
               DO ig=1,ngrid
                 zqi(ig,l,iq)=zqi(ig,l,iq)+qr(ig,l,ir)
               ENDDO
             ENDDO            
          enddo ! of do ir=1,nr
c -----------------------------------------------------------------
c         WATER CYCLE CASE
c -----------------------------------------------------------------
           else if ((iq .eq. iccn_mass) .or. (iq .eq. iccn_number)
     &       .or. (iq .eq. igcm_h2o_ice)) then
            if (microphys) then
              ! water ice sedimentation
              call newsedim(ngrid,nlay,ngrid*nlay,ngrid*nlay,
     &        ptimestep,pplev,masse,epaisseur,zt,rsedcloud,rhocloud,
     &        zqi(1,1,iq),wq,beta)
            else
              ! water ice sedimentation
              call newsedim(ngrid,nlay,ngrid*nlay,1,
     &        ptimestep,pplev,masse,epaisseur,zt,rsedcloud,rho_q(iq),
     &        zqi(1,1,iq),wq,beta)
            endif ! of if (microphys)
c           Tendencies
c           ~~~~~~~~~~
            do ig=1,ngrid 
              pdqs_sed(ig,iq)=wq(ig,1)/ptimestep
            end do
c -----------------------------------------------------------------
c         GENERAL CASE
c -----------------------------------------------------------------
          else
            call newsedim(ngrid,nlay,1,1,ptimestep,
     &      pplev,masse,epaisseur,zt,radius(iq),rho_q(iq),
     &      zqi(1,1,iq),wq,1.0)
c           Tendencies
c           ~~~~~~~~~~
            do ig=1,ngrid 
              pdqs_sed(ig,iq)=wq(ig,1)/ptimestep
            end do
          endif ! of if doubleq and if water
c -----------------------------------------------------------------

c         Compute the final tendency:
c         ---------------------------
          DO l = 1, nlay
            DO ig=1,ngrid
              pdqsed(ig,l,iq)=(zqi(ig,l,iq)-
     $        (pq(ig,l,iq) + pdqfi(ig,l,iq)*ptimestep))/ptimestep
            ENDDO
          ENDDO

        endif ! of if(radius(iq).gt.1.e-9)
c =================================================================
      enddo ! of do iq=1,nq

c     Update the dust particle size "rdust"
c     -------------------------------------
      if (doubleq) then
       DO l = 1, nlay
        DO ig=1,ngrid
        
     
         call updaterdust(zqi(ig,l,igcm_dust_mass),
     &                    zqi(ig,l,igcm_dust_number),rdust(ig,l),
     &                    tauscaling(ig))     

          
        ENDDO
       ENDDO
      endif ! of if (doubleq)

      if (rdstorm) then
       DO l = 1, nlay
        DO ig=1,ngrid
         call updaterdust(zqi(ig,l,igcm_stormdust_mass),
     &                zqi(ig,l,igcm_stormdust_number),rstormdust(ig,l),
     &                tauscaling(ig))    
        ENDDO
       ENDDO
      endif ! of if (rdstorm)

      if (slpwind) then
       DO l = 1, nlay
        DO ig=1,ngrid
         call updaterdust(zqi(ig,l,igcm_topdust_mass),
     &                zqi(ig,l,igcm_topdust_number),rtopdust(ig,l),
     &                tauscaling(ig))    
        ENDDO
       ENDDO
      endif ! of if (slpwind)
 
c     Update the ice particle size "rice"
c     -------------------------------------
      if (water) then
       IF(microphys) THEN 
       
       
        DO l = 1, nlay
          DO ig=1,ngrid

         call updaterice_micro(zqi(ig,l,igcm_h2o_ice),
     &    zqi(ig,l,igcm_ccn_mass),zqi(ig,l,igcm_ccn_number),
     &    tauscaling(ig),rice(ig,l),rhocloud(ig,l))
           
          ENDDO
        ENDDO
        
       ELSE
       
        DO l = 1, nlay
          DO ig=1,ngrid
          
            call updaterice_typ(zqi(ig,l,igcm_h2o_ice),
     &                      tau(ig,1),zlay(ig,l),rice(ig,l))

          ENDDO
        ENDDO
       ENDIF ! of IF(microphys)
      endif ! of if (water)

      END SUBROUTINE callsedim
      
      END MODULE callsedim_mod

