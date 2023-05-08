      SUBROUTINE callsedim(ngrid,nlay, ptimestep,
     &                pplev,zlev, pt, pdt,
     &                pq, pdqfi, pdqsed,pdqs_sed,nq)

      use radinc_h, only : naerkind
      use radii_mod, only: h2o_reffrad
      use aerosol_mod, only : iaero_h2o
      USE tracer_h, only : igcm_co2_ice,igcm_h2o_ice,radius,rho_q
      use comcstfi_mod, only: g
      use callkeys_mod, only : water

      IMPLICIT NONE

!==================================================================
!     
!     Purpose
!     -------
!     Calculates sedimentation of aerosols depending on their
!     density and radius.
!     
!     Authors
!     -------
!     F. Forget (1999)
!     Tracer generalisation by E. Millour (2009)
!     
!==================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

c   arguments:
c   ----------

      integer,intent(in):: ngrid ! number of horizontal grid points
      integer,intent(in):: nlay  ! number of atmospheric layers
      real,intent(in):: ptimestep  ! physics time step (s)
      real,intent(in):: pplev(ngrid,nlay+1) ! pressure at inter-layers (Pa)
      real,intent(in):: pt(ngrid,nlay)      ! temperature at mid-layer (K)
      real,intent(in):: pdt(ngrid,nlay) ! tendency on temperature
      real,intent(in):: zlev(ngrid,nlay+1)  ! altitude at layer boundaries
      integer,intent(in) :: nq ! number of tracers
      real,intent(in) :: pq(ngrid,nlay,nq)  ! tracers (kg/kg)
      real,intent(in) :: pdqfi(ngrid,nlay,nq)  ! tendency on tracers before
                                               ! sedimentation (kg/kg.s-1)
      
      real,intent(out) :: pdqsed(ngrid,nlay,nq) ! tendency due to sedimentation (kg/kg.s-1)
      real,intent(out) :: pdqs_sed(ngrid,nq)    ! flux at surface (kg.m-2.s-1)
      
c   local:
c   ------

      INTEGER l,ig, iq

      ! for particles with varying radii:
      real reffrad(ngrid,nlay,naerkind) ! particle radius (m) 
      real nueffrad(ngrid,nlay,naerkind) ! aerosol effective radius variance

      real zqi(ngrid,nlay,nq) ! to locally store tracers
      real zt(ngrid,nlay) ! to locally store temperature (K)
      real masse (ngrid,nlay) ! Layer mass (kg.m-2)
      real epaisseur (ngrid,nlay) ! Layer thickness (m)
      real wq(ngrid,nlay+1) ! displaced tracer mass (kg.m-2)
c      real dens(ngrid,nlay) ! Mean density of the ice part. accounting for dust core


      LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)

c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
        firstcall=.false.
        ! add some tests on presence of required tracers/aerosols:
        if (water) then
          if (igcm_h2o_ice.eq.0) then
            write(*,*) "callsedim error: water=.true.",
     &                 " but igcm_h2o_ice=0"
          stop
          endif
          if (iaero_h2o.eq.0) then
            write(*,*) "callsedim error: water=.true.",
     &                 " but iaero_ho2=0"
          stop
          endif
        endif
      ENDIF ! of IF (firstcall)
      
!=======================================================================
!     Preliminary calculation of the layer characteristics
!     (mass (kg.m-2), thickness (m), etc.

      do  l=1,nlay
        do ig=1, ngrid
          masse(ig,l)=(pplev(ig,l) - pplev(ig,l+1)) /g 
          epaisseur(ig,l)= zlev(ig,l+1) - zlev(ig,l)
          zt(ig,l)=pt(ig,l)+pdt(ig,l)*ptimestep
        end do
      end do

!======================================================================
! Calculate the transport due to sedimentation for each tracer
! [This has been rearranged by L. Kerber to allow the sedimentation
!  of general tracers.]
 
      do iq=1,nq
       if((radius(iq).gt.1.e-9).and.(iq.ne.igcm_co2_ice)) then ! JVO 08/2017 : be careful radius was tested uninitialized (fixed) ... 
       
!         (no sedim for gases, and co2_ice sedim is done in condense_co2)      

! store locally updated tracers

          do l=1,nlay 
            do ig=1, ngrid
              zqi(ig,l,iq)=pq(ig,l,iq)+pdqfi(ig,l,iq)*ptimestep
            enddo
          enddo ! of do l=1,nlay
        
!======================================================================
! Sedimentation 
!======================================================================
! Water          
          if (water.and.(iq.eq.igcm_h2o_ice)) then
            ! compute radii for h2o_ice 
             call h2o_reffrad(ngrid,nlay,zqi(1,1,igcm_h2o_ice),zt,
     &                reffrad(1,1,iaero_h2o),nueffrad(1,1,iaero_h2o))
            ! call sedimentation for h2o_ice
             call newsedim(ngrid,nlay,ngrid*nlay,ptimestep,
     &            pplev,masse,epaisseur,zt,reffrad(1,1,iaero_h2o),
     &            rho_q(iq),zqi(1,1,igcm_h2o_ice),wq)

! General Case
          else 
             call newsedim(ngrid,nlay,1,ptimestep,
     &            pplev,masse,epaisseur,zt,radius(iq),rho_q(iq),
     &            zqi(1,1,iq),wq)
          endif

!=======================================================================
!     Calculate the tendencies
!======================================================================

          do ig=1,ngrid
!     Ehouarn: with new way of tracking tracers by name, this is simply
            pdqs_sed(ig,iq) = wq(ig,1)/ptimestep
          end do

          DO l = 1, nlay
            DO ig=1,ngrid
              pdqsed(ig,l,iq)=(zqi(ig,l,iq)-
     &        (pq(ig,l,iq) + pdqfi(ig,l,iq)*ptimestep))/ptimestep
            ENDDO
          ENDDO
       endif ! of no gases no co2_ice
      enddo ! of do iq=1,nq
      return
      end
